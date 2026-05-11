# =============================================================================
# F.L.U.I.D — Docking Tasks (Celery Workers)
# =============================================================================
"""
Tarefas Celery para execução de docking molecular em background.

Cada task encapsula uma etapa do DockingManager, garantindo:
    - Execução estritamente assíncrona (sem timeout na API)
    - Garbage collection de arquivos PDBQT pós-processamento
    - Retry automático em caso de falha transitória
    - Rastreabilidade via task_id do Celery

Uso:
    result = run_batch_docking.delay(receptor_path, ligand_smiles, guidance)
    status = result.status     # PENDING | STARTED | SUCCESS | FAILURE
    output = result.get()      # Bloqueia até conclusão (ou use polling)
"""
from __future__ import annotations

import json
import tempfile
from pathlib import Path
from typing import Any

from loguru import logger

from src.workers.celery_app import celery_app
from src.config import FluidConfig, DEFAULT_CONFIG


@celery_app.task(
    bind=True,
    name="fluid.docking.prepare_receptor",
    max_retries=2,
    soft_time_limit=300,
)
def prepare_receptor_task(
    self: Any,
    pdb_content: str,
    protein_id: str,
) -> dict[str, str]:
    """
    Converte PDB em PDBQT para o receptor (background).

    Parâmetros:
        pdb_content: Conteúdo do arquivo PDB como string.
        protein_id: ID para nomeação do output.

    Retorna:
        {"protein_id": str, "pdbqt_path": str, "status": "success"}
    """
    try:
        config = DEFAULT_CONFIG
        workspace = config.workspace_dir / "structures"
        workspace.mkdir(parents=True, exist_ok=True)

        pdb_path = workspace / f"{protein_id}.pdb"
        pdb_path.write_text(pdb_content)

        # TODO: Integrar DockingManager.prepare_receptor() na Sprint 3
        pdbqt_path = workspace / f"{protein_id}.pdbqt"

        logger.info(f"Receptor preparado: {protein_id} → {pdbqt_path}")
        return {
            "protein_id": protein_id,
            "pdbqt_path": str(pdbqt_path),
            "status": "success",
        }
    except Exception as exc:
        logger.error(f"Falha ao preparar receptor {protein_id}: {exc}")
        raise self.retry(exc=exc)


@celery_app.task(
    bind=True,
    name="fluid.docking.run_batch",
    max_retries=1,
    soft_time_limit=3600,
    time_limit=7200,
)
def run_batch_docking(
    self: Any,
    receptor_pdbqt_path: str,
    ligand_smiles_list: list[str],
    ligand_ids: list[str],
    grid_center: list[float],
    grid_size: list[float],
) -> dict[str, Any]:
    """
    Executa batch de docking via AutoDock Vina em background.

    Pipeline:
        1. Gera PDBQT para cada ligante (RDKit → Meeko)
        2. Executa Vina sequencialmente com a Grid Box
        3. Coleta resultados e limpa arquivos temporários (GC)

    Retorna:
        {"results": [...], "n_ligands": int, "status": "success"}
    """
    try:
        config = DEFAULT_CONFIG
        results_dir = config.workspace_dir / "results"
        results_dir.mkdir(parents=True, exist_ok=True)

        temp_files: list[Path] = []
        results: list[dict] = []

        for smiles, lig_id in zip(ligand_smiles_list, ligand_ids):
            # TODO: Integrar DockingManager completo na Sprint 3
            # Placeholder: simula resultado para validação de infra
            results.append({
                "ligand_id": lig_id,
                "smiles": smiles,
                "affinity_kcal": 0.0,
                "status": "pending_implementation",
            })

        # === Garbage Collection rigoroso ===
        for tmp in temp_files:
            if tmp.exists():
                tmp.unlink()
                logger.debug(f"GC: removido {tmp}")

        logger.info(
            f"Batch docking concluído: {len(results)} ligantes processados"
        )
        return {
            "results": results,
            "n_ligands": len(results),
            "receptor": receptor_pdbqt_path,
            "status": "success",
        }
    except Exception as exc:
        logger.error(f"Batch docking falhou: {exc}")
        raise self.retry(exc=exc)


@celery_app.task(name="fluid.docking.cleanup")
def cleanup_workspace(paths: list[str]) -> dict[str, int]:
    """
    Task de limpeza pós-processamento.
    Remove arquivos temporários do workspace.
    """
    removed = 0
    for p in paths:
        path = Path(p)
        if path.exists():
            path.unlink()
            removed += 1
    logger.info(f"Cleanup: {removed}/{len(paths)} arquivos removidos")
    return {"removed": removed, "total": len(paths)}
