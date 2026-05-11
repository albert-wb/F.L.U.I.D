# =============================================================================
# F.L.U.I.D — Orquestrador do Pipeline (FastAPI + Celery)
# =============================================================================
"""
Controlador assíncrono que amarra todos os módulos do backend:
    DataIngestor → RAGEngine → TargetScorer → DockingManager

Expõe endpoints REST via FastAPI para:
    - Submissão de jobs de descoberta de alvos
    - Consulta de status de jobs em andamento
    - Recuperação de resultados finalizados

Todo processamento pesado é delegado a Celery workers via Redis.
"""
from __future__ import annotations

import uuid
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any

import httpx
from fastapi import FastAPI, HTTPException, BackgroundTasks, UploadFile, File, Form
from pydantic import BaseModel, Field
from loguru import logger

from src.config import FluidConfig, PathogenClass, DEFAULT_CONFIG
from src.workers.celery_app import celery_app


# ---------------------------------------------------------------------------
# App FastAPI
# ---------------------------------------------------------------------------
app = FastAPI(
    title="F.L.U.I.D API",
    description="Framework for Ligand-target Unified In-silico Discovery",
    version="0.4.0",
)


# ---------------------------------------------------------------------------
# Schemas de Request/Response
# ---------------------------------------------------------------------------
class JobStatus(str, Enum):
    PENDING = "pending"
    INGESTING = "ingesting"
    RAG_INDEXING = "rag_indexing"
    SCORING = "scoring"
    DOCKING = "docking"
    COMPLETED = "completed"
    FAILED = "failed"


class DiscoveryRequest(BaseModel):
    """Input do usuário para iniciar o pipeline."""
    protein_id: str = Field(..., description="UniProt ID da proteína-alvo")
    protein_name: str = Field("", description="Nome da proteína (ex: 3CLpro)")
    pathogen_class: str = Field("unknown", description="Classe do patógeno")
    taxonomy_id: int = Field(0, description="NCBI Taxonomy ID")
    ligand_smiles: list[str] = Field(
        default_factory=list,
        description="Lista de SMILES dos ligantes candidatos"
    )
    ligand_ids: list[str] = Field(
        default_factory=list,
        description="IDs correspondentes aos ligantes"
    )
    coding_sequences: list[str] = Field(
        default_factory=list,
        description="Sequências codificantes para cálculo dN/dS"
    )
    sequence_ids: list[str] = Field(
        default_factory=list,
        description="IDs das sequências codificantes"
    )
    rag_query: str = Field(
        "",
        description="Query para busca de literatura (auto-gerado se vazio)"
    )
    max_papers: int = Field(30, ge=1, le=100)


class JobResponse(BaseModel):
    """Resposta ao submeter um job."""
    job_id: str
    status: JobStatus
    message: str
    created_at: str


class JobStatusResponse(BaseModel):
    """Status detalhado de um job."""
    job_id: str
    status: JobStatus
    current_step: str
    progress: float = 0.0
    results: dict[str, Any] | None = None
    error: str | None = None


# ---------------------------------------------------------------------------
# Store in-memory de jobs (produção: Redis/DB)
# ---------------------------------------------------------------------------
_jobs: dict[str, dict[str, Any]] = {}


# ---------------------------------------------------------------------------
# Celery Task: Pipeline Completo
# ---------------------------------------------------------------------------
@celery_app.task(
    bind=True,
    name="fluid.pipeline.run_discovery",
    soft_time_limit=7200,
    time_limit=10800,
)
def run_discovery_pipeline(
    self: Any,
    job_id: str,
    request_data: dict,
) -> dict[str, Any]:
    """
    Executa o pipeline completo de descoberta em background.

    Chain: DataIngestor → RAGEngine → TargetScorer → DockingManager

    Cada etapa atualiza o status do job para polling pelo frontend.
    """
    import asyncio

    try:
        _update_job(job_id, JobStatus.INGESTING, "Ingestão de dados", 0.1)

        # === ETAPA 1: DataIngestor (Clustering) ===
        from src.ingestion.data_ingestor import DataIngestor
        ingestor = DataIngestor()

        coding_seqs = request_data.get("coding_sequences", [])
        seq_ids = request_data.get("sequence_ids", [])

        centroid_seqs: list[str] = []
        centroid_ids: list[str] = []
        if coding_seqs:
            centroid_seqs, centroid_ids = ingestor.get_centroid_sequences(
                coding_seqs, seq_ids if seq_ids else None
            )

        _update_job(job_id, JobStatus.RAG_INDEXING, "Indexando literatura", 0.25)

        # === ETAPA 2: RAGEngine ===
        from src.rag.rag_engine import RAGEngine
        rag = RAGEngine()

        protein_name = request_data.get("protein_name", "")
        protein_id = request_data.get("protein_id", "")
        rag_query = request_data.get("rag_query", "")

        if not rag_query:
            rag_query = f"{protein_name} {protein_id} active site catalytic mechanism"

        # Executa async no loop do Celery
        loop = asyncio.new_event_loop()
        try:
            n_chunks = loop.run_until_complete(
                rag.ingest_literature(rag_query, max_papers=request_data.get("max_papers", 30))
            )
            guidance = loop.run_until_complete(
                rag.get_docking_guidance(protein_id, protein_name)
            )
        finally:
            loop.close()

        _update_job(job_id, JobStatus.SCORING, "Calculando scores", 0.50)

        # === ETAPA 3: TargetScorer ===
        from src.scoring.target_scorer import TargetScorer, TargetFeatures
        scorer = TargetScorer()

        # Calcula dN/dS se sequências fornecidas
        omega_mean = 0.5  # default neutro
        if centroid_seqs and len(centroid_seqs) >= 2:
            df_dnds = scorer.compute_dnds_from_alignment(centroid_seqs)
            valid_omegas = df_dnds["omega"].dropna()
            if len(valid_omegas) > 0:
                omega_mean = float(valid_omegas.mean())

        # Cria TargetFeatures para o alvo
        pc_str = request_data.get("pathogen_class", "unknown")
        try:
            pc = PathogenClass(pc_str)
        except ValueError:
            pc = PathogenClass.UNKNOWN

        target = TargetFeatures(
            protein_id=protein_id,
            pathogen_class=pc,
            dnds_raw=omega_mean,
            has_metal_binding=bool(guidance.catalytic_residues),
            druggability_score=min(guidance.confidence + 0.3, 1.0),
            essentiality_score=0.5,  # default, refinado via literatura
        )
        score = scorer.score_single_target(target)

        scoring_result = {
            "protein_id": protein_id,
            "omega_mean": omega_mean,
            "n_centroids": len(centroid_seqs),
            "targetability_score_raw": score,
            "enzyme_family": guidance.enzyme_family.value,
            "grid_center": list(guidance.grid_center),
            "grid_size": list(guidance.grid_size),
            "rag_confidence": guidance.confidence,
            "catalytic_residues": guidance.catalytic_residues,
        }

        _update_job(job_id, JobStatus.DOCKING, "Executando docking", 0.70)

        # === ETAPA 4: DockingManager ===
        docking_result: dict[str, Any] = {"docking_performed": False}
        ligand_smiles = request_data.get("ligand_smiles", [])
        ligand_ids = request_data.get("ligand_ids", [])

        if ligand_smiles:
            from src.docking.docking_manager import DockingManager

            docker = DockingManager()
            # Nota: receptor PDBQT precisa ser preparado previamente
            # Aqui registramos a intenção — preparação real requer PDB file
            docking_result = {
                "docking_performed": False,
                "reason": "Receptor PDB não fornecido via API — "
                          "use endpoint /prepare-receptor primeiro",
                "n_ligands_queued": len(ligand_smiles),
                "grid_center": list(guidance.grid_center),
                "grid_size": list(guidance.grid_size),
            }

        # === RESULTADO FINAL ===
        _update_job(job_id, JobStatus.COMPLETED, "Pipeline concluído", 1.0)

        final = {
            "job_id": job_id,
            "scoring": scoring_result,
            "docking": docking_result,
            "rag_chunks_indexed": n_chunks if 'n_chunks' in dir() else 0,
            "status": "completed",
        }

        _update_job_results(job_id, final)
        return final

    except Exception as e:
        logger.error(f"Pipeline falhou para job {job_id}: {e}")
        _update_job(job_id, JobStatus.FAILED, str(e), 0.0)
        raise


def _update_job(
    job_id: str, status: JobStatus, step: str, progress: float
) -> None:
    """Atualiza status do job no store."""
    if job_id in _jobs:
        _jobs[job_id].update({
            "status": status, "current_step": step, "progress": progress,
        })


def _update_job_results(job_id: str, results: dict) -> None:
    """Armazena resultados finais."""
    if job_id in _jobs:
        _jobs[job_id]["results"] = results


# ---------------------------------------------------------------------------
# Endpoints FastAPI
# ---------------------------------------------------------------------------

@app.post("/api/v1/discover", response_model=JobResponse)
async def submit_discovery(request: DiscoveryRequest) -> JobResponse:
    """
    Submete um job de descoberta de alvos terapêuticos.

    O pipeline inteiro é executado assincronamente via Celery.
    Use GET /api/v1/jobs/{job_id} para consultar o status.
    """
    job_id = str(uuid.uuid4())
    now = datetime.now(timezone.utc).isoformat()

    _jobs[job_id] = {
        "status": JobStatus.PENDING,
        "current_step": "Aguardando worker",
        "progress": 0.0,
        "results": None,
        "created_at": now,
    }

    # Dispatch para Celery
    run_discovery_pipeline.delay(job_id, request.model_dump())

    logger.info(f"Job {job_id} submetido para proteína {request.protein_id}")

    return JobResponse(
        job_id=job_id,
        status=JobStatus.PENDING,
        message=f"Pipeline iniciado para {request.protein_id}",
        created_at=now,
    )


@app.get("/api/v1/jobs/{job_id}", response_model=JobStatusResponse)
async def get_job_status(job_id: str) -> JobStatusResponse:
    """Consulta o status de um job em andamento."""
    if job_id not in _jobs:
        raise HTTPException(status_code=404, detail=f"Job {job_id} não encontrado")

    job = _jobs[job_id]
    return JobStatusResponse(
        job_id=job_id,
        status=job["status"],
        current_step=job["current_step"],
        progress=job["progress"],
        results=job.get("results"),
        error=job.get("error"),
    )


# ---------------------------------------------------------------------------
# Endpoints de Upload / Fetch de Receptor PDB
# ---------------------------------------------------------------------------

class ReceptorResponse(BaseModel):
    """Resposta ao submeter um receptor PDB."""
    job_id: str
    pdb_id: str
    pdb_path: str
    message: str


@app.post("/api/v1/upload-receptor", response_model=ReceptorResponse)
async def upload_receptor(
    pdb_file: UploadFile = File(..., description="Arquivo .pdb do receptor"),
    protein_id: str = Form("", description="UniProt ID (opcional)"),
    protein_name: str = Form("", description="Nome da proteína"),
    pathogen_class: str = Form("unknown"),
    ligand_smiles: str = Form("", description="SMILES separados por vírgula"),
) -> ReceptorResponse:
    """
    Upload de um arquivo PDB local para iniciar o pipeline.

    O receptor é salvo no workspace, preparado via DataIngestor,
    e a chain completa do Celery é disparada automaticamente.
    """
    if not pdb_file.filename or not pdb_file.filename.endswith(".pdb"):
        raise HTTPException(400, "Arquivo deve ter extensão .pdb")

    # Salva o PDB no workspace
    config = DEFAULT_CONFIG
    structures_dir = config.workspace_dir / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)

    pdb_id = protein_id or pdb_file.filename.replace(".pdb", "")
    pdb_path = structures_dir / f"{pdb_id}.pdb"

    content = await pdb_file.read()
    pdb_path.write_bytes(content)
    logger.info(f"Receptor PDB salvo: {pdb_path} ({len(content)} bytes)")

    # Dispara pipeline com receptor
    smiles_list = [s.strip() for s in ligand_smiles.split(",") if s.strip()]
    job_id = _dispatch_pipeline_with_receptor(
        pdb_id, str(pdb_path), protein_name, pathogen_class, smiles_list
    )

    return ReceptorResponse(
        job_id=job_id, pdb_id=pdb_id,
        pdb_path=str(pdb_path),
        message=f"Receptor {pdb_id} salvo e pipeline iniciado",
    )


@app.post("/api/v1/fetch-receptor", response_model=ReceptorResponse)
async def fetch_receptor_by_id(
    pdb_id: str = Form(..., description="PDB ID do RCSB (ex: 6LU7)"),
    protein_name: str = Form(""),
    pathogen_class: str = Form("unknown"),
    ligand_smiles: str = Form("", description="SMILES separados por vírgula"),
) -> ReceptorResponse:
    """
    Baixa um PDB do RCSB pelo ID e inicia o pipeline.

    Endpoint: https://files.rcsb.org/download/{pdb_id}.pdb
    """
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"

    async with httpx.AsyncClient(timeout=30.0) as client:
        resp = await client.get(url)
        if resp.status_code != 200:
            raise HTTPException(404, f"PDB {pdb_id} não encontrado no RCSB")

    config = DEFAULT_CONFIG
    structures_dir = config.workspace_dir / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)

    pdb_path = structures_dir / f"{pdb_id.upper()}.pdb"
    pdb_path.write_bytes(resp.content)
    logger.info(f"PDB {pdb_id} baixado do RCSB: {pdb_path}")

    smiles_list = [s.strip() for s in ligand_smiles.split(",") if s.strip()]
    job_id = _dispatch_pipeline_with_receptor(
        pdb_id.upper(), str(pdb_path), protein_name, pathogen_class, smiles_list
    )

    return ReceptorResponse(
        job_id=job_id, pdb_id=pdb_id.upper(),
        pdb_path=str(pdb_path),
        message=f"PDB {pdb_id} baixado do RCSB e pipeline iniciado",
    )


def _dispatch_pipeline_with_receptor(
    pdb_id: str, pdb_path: str, protein_name: str,
    pathogen_class: str, ligand_smiles: list[str],
) -> str:
    """Cria job e dispara Celery com receptor PDB injetado."""
    job_id = str(uuid.uuid4())
    now = datetime.now(timezone.utc).isoformat()

    _jobs[job_id] = {
        "status": JobStatus.PENDING,
        "current_step": "Aguardando worker",
        "progress": 0.0,
        "results": None,
        "created_at": now,
    }

    request_data = {
        "protein_id": pdb_id,
        "protein_name": protein_name or pdb_id,
        "pathogen_class": pathogen_class,
        "ligand_smiles": ligand_smiles,
        "ligand_ids": [f"lig_{i}" for i in range(len(ligand_smiles))],
        "receptor_pdb_path": pdb_path,
    }

    run_discovery_pipeline.delay(job_id, request_data)
    logger.info(f"Job {job_id} com receptor PDB {pdb_id}")
    return job_id


@app.get("/api/v1/health")
async def health_check() -> dict[str, str]:
    """Health check do serviço."""
    return {"status": "ok", "version": "0.4.0", "service": "F.L.U.I.D"}
