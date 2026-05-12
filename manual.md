# F.L.U.I.D — Manual do Usuário

**Guia de uso para pesquisadores e bioinformatas.**

---

## 1. O que é o F.L.U.I.D?

O F.L.U.I.D (Framework for Ligand-target Unified In-silico Discovery) é uma plataforma computacional para **descoberta de alvos terapêuticos** e **triagem virtual de fármacos**. O software automatiza o pipeline de bioinformática estrutural — desde a análise evolutiva de proteínas-alvo até a simulação de docking molecular com o AutoDock Vina.

**O que o F.L.U.I.D faz por você:**

1. Recebe a estrutura 3D de uma proteína-alvo (arquivo PDB ou ID do RCSB).
2. Analisa a conservação evolutiva e a vulnerabilidade farmacológica do alvo.
3. Consulta a literatura científica (via Inteligência Artificial) para identificar sítios catalíticos e orientar a simulação.
4. Executa o docking molecular com os ligantes que você fornecer.
5. Ranqueia os resultados e exibe uma **animação termodinâmica** dos encaixes mais favoráveis.

O sistema é **agnóstico ao patógeno**: funciona com vírus de RNA, vírus de DNA, bactérias, protozoários e fungos. Basta selecionar a classe correta na interface.

---

## 2. Instalação Rápida

```bash
# 1. Instalar dependências
pip install -r requirements.txt
pip install -r requirements_desktop.txt

# 2. Iniciar a plataforma
python run_local.py
```

A plataforma inicia duas aplicações:
- **API Backend** em `http://localhost:8000` (processamento dos dados)
- **Interface Gráfica Desktop** (janela do CustomTkinter)

---

## 3. Navegando pela Interface

A interface é dividida em duas regiões principais:

```
┌─────────────────┬────────────────────────────────────────┐
│                 │  [Barra de Status + Progresso]          │
│   SIDEBAR       ├────────────────────────────────────────┤
│   (Controles)   │  [Tabela de Resultados]                │
│                 ├────────────────────────────────────────┤
│                 │  [Simulação de Fluidos 2D]             │
└─────────────────┴────────────────────────────────────────┘
```

### 3.1 Sidebar (Painel de Controle — Esquerda)

A sidebar concentra todos os campos de entrada para configurar sua análise.

#### Receptor PDB

Há duas formas de fornecer a estrutura 3D da proteína-alvo:

| Método | Como usar |
|--------|-----------|
| **PDB ID** | Digite o identificador de 4 caracteres do RCSB (ex: `6LU7` para a protease principal do SARS-CoV-2). O sistema baixa automaticamente a estrutura do banco de dados RCSB. |
| **Upload local** | Clique em "📂 Upload PDB Local" para selecionar um arquivo `.pdb` do seu computador. Útil para estruturas obtidas por homologia, AlphaFold ou modelagem experimental. |

#### Classe do Patógeno

Selecione a família taxonômica do organismo-alvo no dropdown. Esta informação é **essencial** para a calibração da métrica evolutiva.

| Classe | Taxa basal de mutação (μ) | Exemplo de organismos |
|--------|--------------------------|----------------------|
| Vírus RNA | Alta (~10⁻⁴ sub/sítio/ciclo) | SARS-CoV-2, HIV, Influenza |
| Vírus DNA | Moderada (~10⁻⁶) | Herpesvírus, Adenovírus |
| Bactéria | Baixa (~10⁻⁹) | *M. tuberculosis*, *S. aureus* |
| Protozoário | Variável (~10⁻⁹) | *Plasmodium*, *Leishmania* |
| Fungo | Baixa (~10⁻⁹) | *Candida*, *Aspergillus* |

**Por que isso importa?** A razão dN/dS (taxa de mutações não-sinônimas vs. sinônimas) é o principal indicador de conservação evolutiva. No entanto, um valor de ω = 0.1 em um vírus RNA tem significado biológico distinto de ω = 0.1 em uma bactéria, pois as taxas basais de mutação diferem em ordens de magnitude. O F.L.U.I.D normaliza automaticamente essa métrica com base na classe selecionada.

#### Nome da Proteína (Opcional)

Campo livre para identificar a proteína nos resultados (ex: "3CLpro", "RdRp", "InhA").

#### Ligantes (SMILES)

Insira as representações SMILES dos compostos candidatos, **um por linha**:

```
CC(=O)Oc1ccccc1C(=O)O
c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34
CC(C)Cc1ccc(cc1)C(C)C(=O)O
```

Se nenhum ligante for fornecido, o pipeline executará apenas a análise de *targetability* sem docking.

#### Botão "Run F.L.U.I.D Pipeline"

Clique para submeter a análise. O botão ficará desabilitado durante o processamento para evitar submissões duplicadas.

---

### 3.2 Barra de Status (Topo — Direita)

Após a submissão, a barra de status exibe o progresso do pipeline em tempo real:

| Mensagem | Significado |
|----------|-------------|
| ⏳ Aguardando submissão... | Nenhum job foi submetido ainda |
| 📡 Submetendo job... | A requisição está sendo enviada para o servidor |
| ✅ Job submetido: abc123... | O servidor aceitou o job e começou a processar |
| 📡 Scoring in progress (25%) | O ranqueamento evolutivo está em andamento |
| 📡 Docking in progress (60%) | O AutoDock Vina está executando as simulações |
| ✅ Pipeline concluído! | Todos os resultados estão prontos |
| ❌ [Mensagem de erro] | Ocorreu uma falha no processamento |

A barra de progresso preenche gradativamente conforme as etapas são concluídas. O sistema faz polling automático a cada 2 segundos — não é necessária nenhuma ação do usuário.

> **Nota técnica:** O processamento pesado (docking molecular) ocorre em workers assíncronos (Celery). A interface gráfica permanece responsiva durante todo o processo.

---

### 3.3 Tabela de Resultados (Meio — Direita)

Quando o pipeline conclui, a tabela é preenchida com os ligantes ranqueados:

| Coluna | Descrição | Como interpretar |
|--------|-----------|------------------|
| **Ligand ID** | Identificador do composto | Corresponde à ordem dos SMILES inseridos |
| **Targetability Score** | Score de vulnerabilidade farmacológica (0 a 1) | **Quanto mais próximo de 0**: mais conservado evolucionariamente → mais difícil de adquirir resistência → melhor alvo terapêutico |
| **ΔG (kcal/mol)** | Energia Livre de Gibbs do encaixe | **Quanto mais negativo**: maior afinidade de ligação |
| **Status** | Indicador visual de match | ✅ Match = ΔG < -7.0 kcal/mol |

**Interpretação do Targetability Score:**

O score é baseado na razão dN/dS (ω), que mede a pressão de seleção evolutiva sobre a proteína:

- **ω < 1.0** → Seleção purificadora (conservação): a proteína é essencial e mutações são letais para o organismo. Alvos ideais.
- **ω ≈ 1.0** → Evolução neutra: sem pressão seletiva significativa.
- **ω > 1.0** → Seleção positiva (diversificação): a proteína muda rapidamente. Risco de resistência.

**Interpretação do ΔG:**

| Faixa de ΔG | Interpretação |
|-------------|---------------|
| < -9.0 kcal/mol | Afinidade muito alta (candidato prioritário) |
| -9.0 a -7.0 kcal/mol | Afinidade moderada a alta |
| -7.0 a -5.0 kcal/mol | Afinidade fraca |
| > -5.0 kcal/mol | Sem afinidade significativa |

Ligantes com ΔG < -7.0 kcal/mol são destacados em **verde** na tabela.

---

## 4. Animação Termodinâmica (Simulação de Fluidos 2D)

O quadro inferior da interface contém um botão para abrir a **simulação visual** do docking. Esta visualização utiliza um motor de física 2D baseado em `matter.js` e `p5.js`.

### Como abrir

Clique em **"🚀 Abrir Simulação de Fluidos"**. A simulação abrirá no seu navegador padrão.

### O que você verá

A simulação representa, de forma análoga, o ambiente molecular dentro da fenda catalítica:

1. **A proteína-alvo** aparece como uma esfera central maior, pulsante.
2. **Os ligantes** aparecem como esferas menores coloridas, flutuando ao redor.

### Comportamento das partículas

#### Movimento Browniano (estado basal)

Todos os ligantes exibem **movimento browniano** — um deslocamento aleatório contínuo que simula a difusão molecular no citosol. Na simulação, forças aleatórias são aplicadas a cada frame, reproduzindo a cinética de colisão térmica em solução aquosa.

Este é o estado padrão de todas as partículas antes que os resultados de docking sejam recebidos.

#### Atração Termodinâmica (encaixe favorável)

Quando o AutoDock Vina calcula uma **Energia Livre de Gibbs (ΔG) inferior a -7.0 kcal/mol**, o sistema interpreta isso como uma interação termodinamicamente favorável e cria uma **Constraint** (mola invisível) entre o ligante e a proteína.

Na prática visual:

- O ligante passa a ser **atraído magneticamente** em direção à fenda da proteína.
- A velocidade de atração é proporcional à magnitude do ΔG: afinidades mais negativas geram atração mais rápida.
- Quando o ligante atinge a superfície da proteína (distância < 15 pixels), o sistema o trava em posição fixa com uma **animação de glow**, indicando o encaixe molecular.

#### Ligantes sem afinidade

Ligantes cujo ΔG é superior a -7.0 kcal/mol permanecem em movimento browniano livre, sem interação com a proteína. Isso representa a ausência de encaixe termodinamicamente viável.

### Resumo visual

| Comportamento | Significado | Critério |
|--------------|-------------|----------|
| Flutuação aleatória | Difusão no citosol | Estado inicial |
| Atração para a proteína | Afinidade de ligação | ΔG < -7.0 kcal/mol |
| Travamento + glow | Encaixe molecular | Distância < 15 px do alvo |
| Sem interação | Ligante ineficaz | ΔG ≥ -7.0 kcal/mol |

---

## 5. Fluxo de Trabalho Recomendado

1. **Preparação**: obtenha o PDB ID da proteína-alvo no [RCSB Protein Data Bank](https://www.rcsb.org).
2. **Ligantes**: prepare a lista de SMILES dos compostos candidatos (ferramentas como PubChem ou ChEMBL podem ajudar).
3. **Submissão**: insira o PDB ID, selecione a classe do patógeno, cole os SMILES e clique em "Run".
4. **Aguarde**: o pipeline leva de 30 segundos (scoring only) a vários minutos (com docking), dependendo do número de ligantes.
5. **Análise**: examine a tabela de resultados. Priorize candidatos com baixo Targetability Score e ΔG < -7.0 kcal/mol.
6. **Visualização**: abra a simulação de fluidos para uma verificação visual do encaixe.

---

## 6. Solução de Problemas

| Problema | Causa provável | Solução |
|----------|---------------|---------|
| "API: offline" na sidebar | O backend não está rodando | Execute `python run_local.py --api` |
| Porta 8000 já em uso | Execução anterior não foi encerrada | O launcher libera a porta automaticamente. Se persistir, feche manualmente pelo Gerenciador de Tarefas |
| Timeout no docking | Muitos ligantes ou grid grande | Reduza o número de SMILES por submissão |
| ΔG = 0 para todos os ligantes | Receptor PDB sem estrutura adequada | Verifique se o PDB contém a cadeia correta |
| Simulação de fluidos não abre | Navegador bloqueou pop-up | Permita pop-ups para `file://` no navegador |

---

## 7. Referências Técnicas

- **dN/dS (ω)**: Nei, M. & Gojobori, T. (1986). Simple methods for estimating the numbers of synonymous and nonsynonymous nucleotide substitutions. *Molecular Biology and Evolution*, 3(5), 418-426.
- **AutoDock Vina**: Eberhardt, J. et al. (2021). AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. *J. Chem. Inf. Model.*, 61(8), 3891-3898.
- **Energia Livre de Gibbs**: valores de ΔG reportados em kcal/mol conforme convenção padrão do campo de docking molecular.
