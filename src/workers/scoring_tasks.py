# =============================================================================
# F.L.U.I.D — Scoring Tasks (Celery Workers)
# =============================================================================
"""
Tarefas Celery para scoring em background.
Útil para datasets grandes onde o clustering + dN/dS
pairwise pode levar minutos.
"""
from __future__ import annotations

from typing import Any

from loguru import logger

from src.workers.celery_app import celery_app


@celery_app.task(
    bind=True,
    name="fluid.scoring.compute_dnds_batch",
    soft_time_limit=1800,
)
def compute_dnds_batch(
    self: Any,
    sequences: list[str],
    sequence_ids: list[str],
    pathogen_class: str,
    identity_threshold: float = 0.95,
) -> dict[str, Any]:
    """
    Pipeline completo: Clustering → dN/dS pairwise → Omega médio.

    Executa em background para não travar a API.
    """
    try:
        from src.ingestion.data_ingestor import DataIngestor
        from src.scoring.target_scorer import TargetScorer

        # 1. Clustering para reduzir N
        ingestor = DataIngestor()
        centroid_seqs, centroid_ids = ingestor.get_centroid_sequences(
            sequences, sequence_ids,
            identity_threshold=identity_threshold,
        )

        # 2. dN/dS pairwise nos centroids
        scorer = TargetScorer()
        df = scorer.compute_dnds_from_alignment(centroid_seqs)

        # 3. Estatísticas
        valid = df["omega"].dropna()
        result = {
            "n_input": len(sequences),
            "n_centroids": len(centroid_seqs),
            "centroid_ids": centroid_ids,
            "omega_mean": float(valid.mean()) if len(valid) > 0 else None,
            "omega_median": float(valid.median()) if len(valid) > 0 else None,
            "omega_std": float(valid.std()) if len(valid) > 0 else None,
            "n_pairs": len(df),
            "pathogen_class": pathogen_class,
            "pairwise_data": df.to_dict(orient="records"),
            "status": "success",
        }

        logger.info(
            f"dN/dS batch: {len(sequences)} seqs → "
            f"{len(centroid_seqs)} centroids → "
            f"{len(df)} pares, ω̄={result['omega_mean']:.4f}"
        )
        return result

    except Exception as exc:
        logger.error(f"Scoring batch falhou: {exc}")
        raise self.retry(exc=exc)


@celery_app.task(
    bind=True,
    name="fluid.scoring.rank_targets",
    soft_time_limit=600,
)
def rank_targets_task(
    self: Any,
    targets_data: list[dict],
) -> dict[str, Any]:
    """
    Ranqueia alvos a partir de dicts serializados (JSON-safe).
    """
    try:
        from src.config import PathogenClass
        from src.scoring.target_scorer import TargetScorer, TargetFeatures

        targets = [
            TargetFeatures(
                protein_id=d["protein_id"],
                pathogen_class=PathogenClass(d.get("pathogen_class", "unknown")),
                dn=d.get("dn", 0.0),
                ds=d.get("ds", 0.0),
                dnds_raw=d.get("dnds_raw", 0.0),
                has_metal_binding=d.get("has_metal_binding", False),
                metal_ions=d.get("metal_ions", []),
                druggability_score=d.get("druggability_score", 0.0),
                essentiality_score=d.get("essentiality_score", 0.0),
                disorder_fraction=d.get("disorder_fraction", 0.0),
            )
            for d in targets_data
        ]

        scorer = TargetScorer()
        df = scorer.score_targets(targets)

        return {
            "ranking": df.to_dict(orient="records"),
            "n_targets": len(targets),
            "top_target": df.iloc[0]["protein_id"] if not df.empty else None,
            "status": "success",
        }

    except Exception as exc:
        logger.error(f"Ranking falhou: {exc}")
        raise self.retry(exc=exc)
