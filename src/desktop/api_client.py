# =============================================================================
# F.L.U.I.D — API Client (Desacoplado da UI)
# =============================================================================
"""
Cliente HTTP para comunicação com a API REST do F.L.U.I.D.

Toda lógica de rede está encapsulada aqui — a GUI nunca faz requests
diretamente. Isso garante separação de responsabilidades e testabilidade.

Thread-safety: todos os métodos são stateless (exceto base_url),
portanto seguros para chamada a partir de threads secundárias.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import requests
from requests.exceptions import ConnectionError, Timeout


# ---------------------------------------------------------------------------
# DTOs de Resposta
# ---------------------------------------------------------------------------
@dataclass
class JobSubmission:
    """Resultado de uma submissão de job."""
    success: bool
    job_id: str = ""
    message: str = ""
    error: str = ""


@dataclass
class JobStatus:
    """Status de um job em andamento."""
    job_id: str
    status: str = "unknown"
    current_step: str = ""
    progress: float = 0.0
    results: dict[str, Any] | None = None
    error: str | None = None


# ---------------------------------------------------------------------------
# Cliente
# ---------------------------------------------------------------------------
class FluidAPIClient:
    """
    Cliente HTTP para a API F.L.U.I.D.

    Uso:
        client = FluidAPIClient("http://localhost:8000")
        job = client.submit_discovery("P0DTD1", "3CLpro", "rna_virus", [...])
        status = client.poll_job(job.job_id)
    """

    def __init__(self, base_url: str = "http://localhost:8000") -> None:
        self.base_url = base_url.rstrip("/")
        self._timeout = 15  # segundos

    def health_check(self) -> bool:
        """Verifica se a API está online."""
        try:
            r = requests.get(
                f"{self.base_url}/api/v1/health", timeout=5
            )
            return r.status_code == 200
        except (ConnectionError, Timeout):
            return False

    def submit_discovery(
        self,
        protein_id: str,
        protein_name: str = "",
        pathogen_class: str = "unknown",
        ligand_smiles: list[str] | None = None,
        coding_sequences: list[str] | None = None,
    ) -> JobSubmission:
        """Submete job via POST /api/v1/discover."""
        payload = {
            "protein_id": protein_id,
            "protein_name": protein_name or protein_id,
            "pathogen_class": pathogen_class,
            "ligand_smiles": ligand_smiles or [],
            "ligand_ids": [f"lig_{i}" for i in range(len(ligand_smiles or []))],
            "coding_sequences": coding_sequences or [],
        }
        try:
            r = requests.post(
                f"{self.base_url}/api/v1/discover",
                json=payload, timeout=self._timeout,
            )
            if r.status_code == 200:
                d = r.json()
                return JobSubmission(True, d["job_id"], d.get("message", ""))
            return JobSubmission(False, error=f"HTTP {r.status_code}: {r.text[:200]}")
        except (ConnectionError, Timeout) as e:
            return JobSubmission(False, error=f"Conexão falhou: {e}")

    def upload_receptor(
        self,
        pdb_path: Path,
        protein_id: str = "",
        protein_name: str = "",
        pathogen_class: str = "unknown",
        ligand_smiles: str = "",
    ) -> JobSubmission:
        """Upload de PDB via POST /api/v1/upload-receptor."""
        try:
            with open(pdb_path, "rb") as f:
                r = requests.post(
                    f"{self.base_url}/api/v1/upload-receptor",
                    files={"pdb_file": (pdb_path.name, f, "chemical/x-pdb")},
                    data={
                        "protein_id": protein_id or pdb_path.stem,
                        "protein_name": protein_name,
                        "pathogen_class": pathogen_class,
                        "ligand_smiles": ligand_smiles,
                    },
                    timeout=self._timeout,
                )
            if r.status_code == 200:
                d = r.json()
                return JobSubmission(True, d["job_id"], d.get("message", ""))
            return JobSubmission(False, error=f"HTTP {r.status_code}")
        except (ConnectionError, Timeout, OSError) as e:
            return JobSubmission(False, error=str(e))

    def fetch_receptor(
        self,
        pdb_id: str,
        protein_name: str = "",
        pathogen_class: str = "unknown",
        ligand_smiles: str = "",
    ) -> JobSubmission:
        """Fetch PDB do RCSB via POST /api/v1/fetch-receptor."""
        try:
            r = requests.post(
                f"{self.base_url}/api/v1/fetch-receptor",
                data={
                    "pdb_id": pdb_id,
                    "protein_name": protein_name,
                    "pathogen_class": pathogen_class,
                    "ligand_smiles": ligand_smiles,
                },
                timeout=30,
            )
            if r.status_code == 200:
                d = r.json()
                return JobSubmission(True, d["job_id"], d.get("message", ""))
            return JobSubmission(False, error=f"HTTP {r.status_code}: {r.text[:200]}")
        except (ConnectionError, Timeout) as e:
            return JobSubmission(False, error=str(e))

    def poll_job(self, job_id: str) -> JobStatus:
        """Consulta status via GET /api/v1/jobs/{id}."""
        try:
            r = requests.get(
                f"{self.base_url}/api/v1/jobs/{job_id}", timeout=self._timeout,
            )
            if r.status_code == 200:
                d = r.json()
                return JobStatus(
                    job_id=job_id,
                    status=d.get("status", "unknown"),
                    current_step=d.get("current_step", ""),
                    progress=d.get("progress", 0.0),
                    results=d.get("results"),
                    error=d.get("error"),
                )
            return JobStatus(job_id, status="error", error=f"HTTP {r.status_code}")
        except (ConnectionError, Timeout) as e:
            return JobStatus(job_id, status="error", error=str(e))
