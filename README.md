# F.L.U.I.D — Framework for Ligand-target Unified In-silico Discovery

**Pipeline agnóstico de descoberta de alvos terapêuticos e docking molecular in silico.**

F.L.U.I.D é uma plataforma computacional para triagem virtual de fármacos que integra análise evolutiva (dN/dS), busca semântica na literatura científica (RAG) e simulações de docking molecular (AutoDock Vina) em um pipeline assíncrono de alto desempenho. A arquitetura é agnóstica ao patógeno — capaz de processar genomas e estruturas de organismos que vão de vírus de RNA a bactérias e protozoários.

---

## Arquitetura do Sistema

```
┌──────────────┐     ┌──────────────┐     ┌──────────────┐     ┌──────────────┐
│ DataIngestor │────>│  RAGEngine   │────>│ TargetScorer │────>│DockingManager│
│  O(N × L)    │     │  O(log n)    │     │ T = Σw·f×Φ×Ψ│     │  Vina + GC   │
└──────────────┘     └──────────────┘     └──────────────┘     └──────────────┘
       │                    │                    │                     │
       └────────────────────┴────────────────────┴─────────────────────┘
                              Celery Workers (Redis)
                                     │
                              ┌──────┴──────┐
                              │  FastAPI    │
                              │  REST API   │
                              └──────┬──────┘
                                     │
                         ┌───────────┴───────────┐
                         │  Desktop GUI          │
                         │  CustomTkinter + JS   │
                         └───────────────────────┘
```

---

## Módulos do Backend

### 1. DataIngestor — Pré-clustering de Sequências

O módulo de ingestão resolve o gargalo computacional do cálculo pairwise de dN/dS, cuja complexidade nativa é **O(N² × L)** (N sequências, L códons por sequência).

**Otimização implementada:**
- **Greedy Centroid Clustering** (inspirado no CD-HIT): agrupa sequências por similaridade, mantendo apenas os centroides representativos.
- **Pré-filtro k-mer Jaccard** (O(L/k)): descarta pares de sequências distantes sem realizar alinhamento, evitando comparações O(L²) desnecessárias.

| Etapa | Complexidade | Descrição |
|-------|-------------|-----------|
| Extração de k-mers | O(N × L) | Hash sets por sequência |
| Pré-filtro Jaccard | O(N × C × L/k) | Descarta pares distintos |
| Alinhamento exato | O(N × C × L) | Apenas para candidatos |
| **Total amortizado** | **O(N × L)** | C << N em dados biológicos |

### 2. RAGEngine — Busca Semântica Estruturada

Motor de Retrieval-Augmented Generation para orientação da Grid Box do docking molecular.

- **Ingestão**: artigos Open Access via Europe PMC REST API → chunking recursivo → embedding (MiniLM-L6-v2) → indexação em **ChromaDB** (busca vetorial HNSW com complexidade amortizada O(log n)).
- **Output estruturado**: a LLM (GPT-4o-mini) é forçada a retornar JSON válido via `PydanticOutputParser` contendo `enzyme_class`, `grid_center` (x, y, z em Å), `grid_size` e `confidence`. Texto livre não é permitido.
- **Fallback**: se a confiança da extração for inferior a 0.3, o sistema calcula o centroide geométrico dos resíduos catalíticos citados na literatura, evitando blind docking de complexidade O(n³).

### 3. TargetScorer — Ranqueamento Evolutivo

Implementa a fórmula de *targetability score* para ranqueamento de alvos terapêuticos:

```
T_raw = [ w_c · C + w_d · D + w_e · E ] × Φ_metal × Ψ_disorder
T_final = normalização min-max no batch
```

| Componente | Símbolo | Descrição |
|-----------|---------|-----------|
| Conservação | C | Sigmoide invertida de ω (dN/dS normalizado pela taxa basal μ) |
| Druggabilidade | D | Score de fpocket/DoGSiteScorer |
| Essencialidade | E | Derivado de dados de mutagênese/KO |
| Multiplicador metálico | Φ | 2.0 para Zn²⁺/Mn²⁺, 1.5 para outros íons, 1.0 sem metais |
| Penalização por desordem | Ψ | Função por partes (threshold 0.3, floor 0.1) |

A **normalização adaptativa** utiliza a taxa mutacional basal (μ) do patógeno para padronizar o ω entre classes taxonômicas distintas (vírus RNA vs bactérias).

### 4. DockingManager — AutoDock Vina Wrapper

- **Preparação de ligantes**: SMILES → RDKit (AddHs + ETKDG v3 + MMFF94) → Meeko → PDBQT.
- **Preparação de receptor**: PDB → `mk_prepare_receptor` (Meeko) com fallback OpenBabel.
- **Execução**: Vina via biblioteca Python ou subprocess, com Grid Box orientada pelo RAGEngine.
- **Garbage Collection rigoroso**: todo I/O em `tempfile.mkdtemp()`, rastreio de diretórios temporários, `shutil.rmtree()` pós-batch, safety net no destrutor. Nenhum arquivo `.pdbqt`, `.log` ou `.txt` persiste após a extração de ΔG.

O processamento é delegado a **Celery Workers** com Redis como message broker, garantindo que simulações CPU-bound não travem a API.

---

## Arquitetura do Frontend

### GUI Desktop (CustomTkinter)

A interface é construída com `customtkinter` em modo dark e utiliza **threading rigoroso** para evitar bloqueio do event loop:

- **Regra absoluta**: chamadas `requests.get/post` **jamais** executam na thread principal do `mainloop()`.
- **Submissão de jobs**: `threading.Thread(daemon=True)` para chamadas de rede.
- **Polling de status**: `self.after(2000, callback)` — timer nativo do Tk, thread-safe.
- **Atualizações de UI**: sempre via `self.after(0, callback)` para garantir execução no mainloop.

### Simulação de Fluidos 2D (matter.js + p5.js)

O motor de física 2D simula a interação fármaco-proteína em tempo real:

- **Broadphase Grid (Spatial Hashing)**: detecção de colisão em **O(n)** amortizado — jamais O(n²).
- **Movimento browniano**: forças aleatórias O(n) simulam difusão térmica no citosol.
- **Atração termodinâmica**: quando ΔG < -7.0 kcal/mol, uma `Constraint` (mola invisível) do matter.js atrai o ligante para a fenda catalítica da proteína, simulando o encaixe molecular.

---

## Instalação e Deploy

### Pré-requisitos

- Python 3.11+
- Docker e Docker Compose (para deploy em produção)

### Opção 1: Desenvolvimento Local

```bash
# Clonar o repositório
git clone https://github.com/albert-wb/F.L.U.I.D.git
cd F.L.U.I.D

# Instalar dependências do backend
pip install -r requirements.txt

# Instalar dependências do desktop
pip install -r requirements_desktop.txt

# Iniciar API + GUI Desktop
python run_local.py
```

O launcher sobe automaticamente a API FastAPI (porta 8000) e a GUI desktop.

| Comando | Descrição |
|---------|-----------|
| `python run_local.py` | API + GUI Desktop |
| `python run_local.py --api` | Apenas a API |
| `python run_local.py --gui` | Apenas a GUI |
| `python run_local.py --dash` | Dashboard web (Dash/Plotly) |

### Opção 2: Docker Compose (Produção)

```bash
# Configurar variáveis de ambiente
cp .env.example .env
# Editar .env com OPENAI_API_KEY, NCBI_API_KEY, etc.

# Buildar e subir os 3 serviços
docker compose up --build -d

# Verificar saúde
curl http://localhost:8000/api/v1/health
```

**Serviços Docker:**

| Serviço | Porta | Descrição |
|---------|-------|-----------|
| `fluid-api` | 8000 | FastAPI + Uvicorn (2 workers) |
| `fluid-redis` | 6379 | Redis 7 (broker + backend) |
| `fluid-worker` | — | Celery (concurrency=4, Vina + RDKit) |

### Endpoints da API

| Método | Rota | Descrição |
|--------|------|-----------|
| `POST` | `/api/v1/discover` | Submete job de descoberta |
| `POST` | `/api/v1/upload-receptor` | Upload de arquivo PDB |
| `POST` | `/api/v1/fetch-receptor` | Baixa PDB do RCSB por ID |
| `GET` | `/api/v1/jobs/{id}` | Consulta status do job |
| `GET` | `/api/v1/health` | Health check |
| `GET` | `/docs` | Documentação Swagger |

---

## Testes

```bash
# Executar suíte de testes
pytest tests/ -v --tb=short

# Com cobertura
pytest tests/ --cov=src --cov-report=term-missing
```

---

## Estrutura do Projeto

```
F.L.U.I.D/
├── Dockerfile
├── docker-compose.yml
├── requirements.txt
├── requirements_desktop.txt
├── run_local.py
│
├── src/
│   ├── api/
│   │   └── main.py              # FastAPI + Celery chain
│   ├── ingestion/
│   │   └── data_ingestor.py     # Clustering + IDR filter
│   ├── rag/
│   │   └── rag_engine.py        # LangChain + ChromaDB
│   ├── scoring/
│   │   └── target_scorer.py     # dN/dS + T = [w·f] × Φ × Ψ
│   ├── docking/
│   │   └── docking_manager.py   # Vina + RDKit + GC
│   ├── dashboard/
│   │   ├── visualizer.py        # Dash (web)
│   │   └── assets/
│   │       └── fluid_engine.js  # matter.js + p5.js
│   ├── desktop/
│   │   ├── api_client.py        # Cliente HTTP thread-safe
│   │   └── main_gui.py          # CustomTkinter GUI
│   └── workers/
│       ├── celery_app.py        # Config + stub inline
│       ├── docking_tasks.py
│       └── scoring_tasks.py
│
└── tests/
    ├── test_data_ingestor.py
    ├── test_target_scorer.py
    └── test_integration.py
```

---

## Licença

Este projeto é estritamente acadêmico e computacional. Consulte o arquivo `LICENSE` para detalhes.
