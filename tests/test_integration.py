# =============================================================================
# F.L.U.I.D — Testes de Integração End-to-End
# =============================================================================
"""
Testa o pipeline completo simulando requisições HTTP ao FastAPI.

Valida:
    1. Submissão de job via POST /discover retorna job_id
    2. Celery task é registrada corretamente
    3. Resultado final contém targetability_score e ΔG
    4. Upload de receptor PDB funciona
    5. Formato JSON compatível com frontend Dash
"""
from __future__ import annotations

from unittest.mock import patch, MagicMock, AsyncMock
from pathlib import Path
import json

import pytest
from fastapi.testclient import TestClient

from src.api.main import app, _jobs, JobStatus


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture
def client() -> TestClient:
    """TestClient do FastAPI (síncrono, sem server real)."""
    return TestClient(app)


@pytest.fixture(autouse=True)
def clean_jobs():
    """Limpa o store de jobs entre testes."""
    _jobs.clear()
    yield
    _jobs.clear()


# ---------------------------------------------------------------------------
# Testes de Endpoint: Health Check
# ---------------------------------------------------------------------------
class TestHealthCheck:
    def test_health_returns_ok(self, client: TestClient):
        resp = client.get("/api/v1/health")
        assert resp.status_code == 200
        data = resp.json()
        assert data["status"] == "ok"
        assert data["service"] == "F.L.U.I.D"


# ---------------------------------------------------------------------------
# Testes de Endpoint: Submissão de Job
# ---------------------------------------------------------------------------
class TestDiscoverySubmission:
    @patch("src.api.main.run_discovery_pipeline")
    def test_submit_returns_job_id(self, mock_task, client: TestClient):
        """POST /discover deve retornar job_id e status pending."""
        mock_task.delay = MagicMock()

        payload = {
            "protein_id": "P0DTD1",
            "protein_name": "3CLpro",
            "pathogen_class": "rna_virus",
            "ligand_smiles": ["CC(=O)Oc1ccccc1C(=O)O"],
            "ligand_ids": ["aspirin"],
        }

        resp = client.post("/api/v1/discover", json=payload)
        assert resp.status_code == 200

        data = resp.json()
        assert "job_id" in data
        assert data["status"] == "pending"
        assert len(data["job_id"]) == 36  # UUID format

        # Verifica que Celery foi chamado
        mock_task.delay.assert_called_once()

    @patch("src.api.main.run_discovery_pipeline")
    def test_submit_without_protein_id_fails(self, mock_task, client: TestClient):
        """Requisição sem protein_id deve retornar 422."""
        resp = client.post("/api/v1/discover", json={})
        assert resp.status_code == 422

    @patch("src.api.main.run_discovery_pipeline")
    def test_submit_creates_job_in_store(self, mock_task, client: TestClient):
        """Job deve existir no store após submissão."""
        mock_task.delay = MagicMock()

        resp = client.post("/api/v1/discover", json={"protein_id": "TEST1"})
        job_id = resp.json()["job_id"]

        assert job_id in _jobs
        assert _jobs[job_id]["status"] == JobStatus.PENDING


# ---------------------------------------------------------------------------
# Testes de Endpoint: Consulta de Status
# ---------------------------------------------------------------------------
class TestJobStatus:
    def test_nonexistent_job_returns_404(self, client: TestClient):
        resp = client.get("/api/v1/jobs/fake-uuid-1234")
        assert resp.status_code == 404

    def test_completed_job_returns_results(self, client: TestClient):
        """Job concluído deve retornar scoring + docking no JSON."""
        job_id = "test-job-123"
        _jobs[job_id] = {
            "status": JobStatus.COMPLETED,
            "current_step": "Pipeline concluído",
            "progress": 1.0,
            "results": {
                "scoring": {
                    "protein_id": "P0DTD1",
                    "targetability_score_raw": 0.847,
                    "omega_mean": 0.12,
                    "enzyme_family": "protease",
                },
                "docking": {
                    "docking_performed": True,
                    "results": [
                        {"ligand_id": "lig_0", "affinity_kcal": -8.3},
                        {"ligand_id": "lig_1", "affinity_kcal": -6.1},
                    ],
                },
            },
        }

        resp = client.get(f"/api/v1/jobs/{job_id}")
        assert resp.status_code == 200

        data = resp.json()
        assert data["status"] == "completed"
        assert data["progress"] == 1.0

        # Valida formato esperado pelo frontend Dash
        results = data["results"]
        assert "scoring" in results
        assert "docking" in results
        assert isinstance(results["scoring"]["targetability_score_raw"], float)
        assert results["scoring"]["targetability_score_raw"] == 0.847

        # Verifica ΔG mesclado
        docking_results = results["docking"]["results"]
        assert len(docking_results) == 2
        assert docking_results[0]["affinity_kcal"] == -8.3


# ---------------------------------------------------------------------------
# Testes de Endpoint: Upload de Receptor PDB
# ---------------------------------------------------------------------------
class TestReceptorUpload:
    @patch("src.api.main.run_discovery_pipeline")
    def test_upload_pdb_file(self, mock_task, client: TestClient, tmp_path: Path):
        """Upload de arquivo .pdb deve salvar e disparar pipeline."""
        mock_task.delay = MagicMock()

        pdb_content = b"ATOM      1  N   ALA A   1       1.000   2.000   3.000\nEND\n"

        resp = client.post(
            "/api/v1/upload-receptor",
            files={"pdb_file": ("test_receptor.pdb", pdb_content, "chemical/x-pdb")},
            data={
                "protein_id": "TEST_REC",
                "protein_name": "Test Receptor",
                "pathogen_class": "bacterium",
                "ligand_smiles": "CCO,CC(=O)O",
            },
        )

        assert resp.status_code == 200
        data = resp.json()
        assert data["pdb_id"] == "TEST_REC"
        assert "job_id" in data
        mock_task.delay.assert_called_once()

    @patch("src.api.main.run_discovery_pipeline")
    def test_upload_non_pdb_fails(self, mock_task, client: TestClient):
        """Arquivo sem extensão .pdb deve retornar 400."""
        resp = client.post(
            "/api/v1/upload-receptor",
            files={"pdb_file": ("file.txt", b"not a pdb", "text/plain")},
            data={"protein_id": "X"},
        )
        assert resp.status_code == 400


# ---------------------------------------------------------------------------
# Testes de Formato: Compatibilidade com Frontend Dash
# ---------------------------------------------------------------------------
class TestDashFrontendFormat:
    def test_result_json_serializable(self):
        """Resultado final deve ser JSON-serializável para o Dash."""
        result = {
            "job_id": "abc-123",
            "scoring": {
                "protein_id": "P0DTD1",
                "omega_mean": 0.15,
                "targetability_score_raw": 0.92,
                "enzyme_family": "protease",
                "grid_center": [10.5, 20.3, 15.7],
                "grid_size": [22.0, 22.0, 22.0],
                "rag_confidence": 0.75,
                "catalytic_residues": ["His41", "Cys145"],
            },
            "docking": {
                "docking_performed": True,
                "results": [
                    {
                        "ligand_id": "lig_0",
                        "protein_id": "P0DTD1",
                        "affinity_kcal": -9.1,
                        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                        "rmsd_lb": 0.0,
                        "rmsd_ub": 0.0,
                        "mode_rank": 1,
                        "success": True,
                    },
                ],
            },
            "status": "completed",
        }

        # Deve serializar sem erro
        json_str = json.dumps(result)
        parsed = json.loads(json_str)

        # Score e ΔG devem estar presentes
        assert parsed["scoring"]["targetability_score_raw"] == 0.92
        assert parsed["docking"]["results"][0]["affinity_kcal"] == -9.1

        # Formato deve ter as chaves que o Dash espera
        assert "scoring" in parsed
        assert "docking" in parsed
        assert "targetability_score_raw" in parsed["scoring"]

    def test_affinity_threshold_logic(self):
        """ΔG < -7.0 deve ativar atração no componente de fluido."""
        threshold = -7.0
        results = [
            {"ligand_id": "lig_0", "affinity_kcal": -9.1},  # Match
            {"ligand_id": "lig_1", "affinity_kcal": -5.2},  # No match
            {"ligand_id": "lig_2", "affinity_kcal": -7.5},  # Match
        ]
        matches = [r for r in results if r["affinity_kcal"] < threshold]
        assert len(matches) == 2
        assert all(r["affinity_kcal"] < threshold for r in matches)


# ---------------------------------------------------------------------------
# Testes do Pipeline (Unidade do Celery Task)
# ---------------------------------------------------------------------------
class TestPipelineUnit:
    def test_scoring_produces_valid_output(self):
        """TargetScorer deve produzir score entre 0 e 2 (antes de min-max)."""
        from src.scoring.target_scorer import TargetScorer, TargetFeatures
        from src.config import PathogenClass

        scorer = TargetScorer()
        target = TargetFeatures(
            protein_id="P0DTD1",
            pathogen_class=PathogenClass.RNA_VIRUS,
            dnds_raw=0.1,
            has_metal_binding=True,
            metal_ions=["Zn2+"],
            druggability_score=0.8,
            essentiality_score=0.9,
        )
        score = scorer.score_single_target(target)
        assert 0.0 <= score <= 2.0  # Multiplicador metálico pode levar > 1

    def test_clustering_reduces_input(self):
        """DataIngestor deve reduzir sequências redundantes."""
        from src.ingestion.data_ingestor import DataIngestor

        ingestor = DataIngestor()
        seqs = ["ATGATGATGATGATGATG"] * 20 + ["GCGCGCGCGCGCGCGCGC"] * 10
        centroids, ids = ingestor.get_centroid_sequences(seqs)
        assert len(centroids) == 2
        assert len(centroids) < len(seqs)
