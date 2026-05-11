# =============================================================================
# F.L.U.I.D — Módulo 1: DataIngestor (Parser Agnóstico + Clustering)
# =============================================================================
"""
Ingestão de dados biológicos com pré-clustering de sequências para
mitigar o custo O(N²×L) do cálculo pairwise dN/dS downstream.

Clustering: greedy centroid com pré-filtro k-mer Jaccard.
Complexidade amortizada: O(N × L).
"""
from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from loguru import logger

from src.config import FluidConfig, PathogenClass, DEFAULT_CONFIG

# Tags do UniProt para limpeza agnóstica
UNIPROT_FEATURE_TAGS: set[str] = {
    "Active site", "Metal binding", "Binding site",
    "Mutagenesis", "Glycosylation", "Disulfide bond",
    "Modified residue", "Lipidation",
}

DEFAULT_IDENTITY_THRESHOLD: float = 0.95
DEFAULT_KMER_SIZE: int = 6


# ---------------------------------------------------------------------------
# DTOs
# ---------------------------------------------------------------------------
@dataclass
class ProteinRecord:
    """Registro normalizado de uma proteína-alvo."""
    uniprot_id: str
    pdb_ids: list[str] = field(default_factory=list)
    organism: str = ""
    pathogen_class: PathogenClass = PathogenClass.UNKNOWN
    taxonomy_id: int = 0
    sequence: str = ""
    active_sites: list[dict] = field(default_factory=list)
    metal_binding_sites: list[dict] = field(default_factory=list)
    glycosylation_sites: list[dict] = field(default_factory=list)
    mutagenesis_sites: list[dict] = field(default_factory=list)
    disorder_regions: list[tuple[int, int]] = field(default_factory=list)
    length: int = 0


@dataclass
class LigandRecord:
    """Registro normalizado de um ligante candidato."""
    compound_id: str
    source: str = ""
    smiles: str = ""
    inchi_key: str = ""
    molecular_weight: float = 0.0
    logp: float = 0.0


@dataclass
class SequenceCluster:
    """Resultado de um cluster de sequências."""
    centroid_index: int
    centroid_id: str
    centroid_sequence: str
    member_indices: list[int] = field(default_factory=list)
    member_ids: list[str] = field(default_factory=list)
    size: int = 1


# ---------------------------------------------------------------------------
# Classe Principal
# ---------------------------------------------------------------------------
class DataIngestor:
    """
    Parser agnóstico com pré-clustering greedy centroid (estilo CD-HIT).

    Pipeline de clustering:
        1. Ordena sequências por comprimento (desc)
        2. Primeira seq → centroid do cluster 0
        3. Para cada seq seguinte:
           a. Pré-filtro Jaccard de k-mers — O(L/k)
           b. Se Jaccard ≥ 0.60 → alinhamento exato O(L)
           c. Se identidade ≥ threshold → agrupa; senão → novo cluster
        4. Retorna centroids para pipeline downstream

    Complexidade: O(N × C × L) onde C = clusters. Como C << N para
    datasets biológicos redundantes → O(N × L) amortizado.
    """

    def __init__(self, config: FluidConfig = DEFAULT_CONFIG) -> None:
        self._config = config

    # === Clustering ===

    @staticmethod
    def _extract_kmers(seq: str, k: int = DEFAULT_KMER_SIZE) -> set[str]:
        """Extrai k-mers únicos. O(L)."""
        if len(seq) < k:
            return {seq}
        return {seq[i:i + k] for i in range(len(seq) - k + 1)}

    @staticmethod
    def _jaccard(a: set[str], b: set[str]) -> float:
        """Similaridade de Jaccard. O(min(|A|,|B|))."""
        if not a or not b:
            return 0.0
        return len(a & b) / len(a | b)

    @staticmethod
    def _pairwise_identity(a: str, b: str) -> float:
        """Identidade percentual (alinhamento truncado). O(L)."""
        max_len = max(len(a), len(b))
        if max_len == 0:
            return 1.0
        matches = sum(1 for i in range(min(len(a), len(b))) if a[i] == b[i])
        return matches / max_len

    def cluster_sequences(
        self,
        sequences: list[str],
        sequence_ids: list[str] | None = None,
        *,
        identity_threshold: float = DEFAULT_IDENTITY_THRESHOLD,
        kmer_size: int = DEFAULT_KMER_SIZE,
        jaccard_prefilter: float = 0.60,
    ) -> list[SequenceCluster]:
        """
        Clusteriza sequências removendo redundância ≥ identity_threshold.
        Retorna lista de SequenceCluster com centroids e membros.
        """
        n = len(sequences)
        if n == 0:
            return []

        if sequence_ids is None:
            sequence_ids = [f"seq_{i}" for i in range(n)]

        # Ordena por comprimento desc (centroid = mais longa)
        order = sorted(range(n), key=lambda i: len(sequences[i]), reverse=True)

        # Pré-computa k-mers — O(N × L)
        kmer_cache: list[set[str]] = [
            self._extract_kmers(sequences[i], kmer_size) for i in range(n)
        ]

        clusters: list[SequenceCluster] = []
        assigned: set[int] = set()

        for idx in order:
            if idx in assigned:
                continue

            seq = sequences[idx]
            kmers = kmer_cache[idx]
            best_cluster: int | None = None
            best_id: float = 0.0

            for c_i, cl in enumerate(clusters):
                # Pré-filtro rápido via Jaccard
                jac = self._jaccard(kmers, kmer_cache[cl.centroid_index])
                if jac < jaccard_prefilter:
                    continue

                # Alinhamento exato
                ident = self._pairwise_identity(seq, cl.centroid_sequence)
                if ident >= identity_threshold and ident > best_id:
                    best_id = ident
                    best_cluster = c_i

            if best_cluster is not None:
                clusters[best_cluster].member_indices.append(idx)
                clusters[best_cluster].member_ids.append(sequence_ids[idx])
                clusters[best_cluster].size += 1
            else:
                clusters.append(SequenceCluster(
                    centroid_index=idx,
                    centroid_id=sequence_ids[idx],
                    centroid_sequence=seq,
                    member_indices=[idx],
                    member_ids=[sequence_ids[idx]],
                    size=1,
                ))
            assigned.add(idx)

        logger.info(
            f"Clustering: {n} seqs → {len(clusters)} clusters "
            f"(threshold={identity_threshold:.0%}, "
            f"redução={1 - len(clusters)/max(n,1):.1%})"
        )
        return clusters

    def get_centroid_sequences(
        self,
        sequences: list[str],
        sequence_ids: list[str] | None = None,
        *,
        identity_threshold: float = DEFAULT_IDENTITY_THRESHOLD,
    ) -> tuple[list[str], list[str]]:
        """Retorna (centroid_sequences, centroid_ids) para uso direto."""
        clusters = self.cluster_sequences(
            sequences, sequence_ids, identity_threshold=identity_threshold
        )
        if not clusters:
            return [], []
        seqs, ids = zip(*[(c.centroid_sequence, c.centroid_id) for c in clusters])
        return list(seqs), list(ids)

    # === Filtro de IDRs ===

    @staticmethod
    def compute_disorder_fraction(
        sequence: str, disorder_regions: list[tuple[int, int]]
    ) -> float:
        """Fração de resíduos em regiões de desordem intrínseca [0,1]."""
        if not sequence or not disorder_regions:
            return 0.0
        total = len(sequence)
        residues: set[int] = set()
        for s, e in disorder_regions:
            for p in range(max(1, s), min(e + 1, total + 1)):
                residues.add(p)
        return len(residues) / total

    @staticmethod
    def filter_disordered_targets(
        records: list[ProteinRecord], *, max_disorder: float = 0.70
    ) -> list[ProteinRecord]:
        """Remove proteínas com desordem acima do threshold."""
        return [
            r for r in records
            if DataIngestor.compute_disorder_fraction(
                r.sequence, r.disorder_regions
            ) <= max_disorder
        ]

    # === Parser UniProt ===

    @staticmethod
    def parse_uniprot_features(raw: list[dict]) -> dict[str, list[dict]]:
        """Agrupa features do UniProt por tag padronizada."""
        grouped: dict[str, list[dict]] = defaultdict(list)
        for feat in raw:
            ft = feat.get("type", "")
            if ft in UNIPROT_FEATURE_TAGS:
                loc = feat.get("location", {})
                grouped[ft].append({
                    "type": ft,
                    "description": feat.get("description", ""),
                    "start": loc.get("start", {}).get("value", 0),
                    "end": loc.get("end", {}).get("value", 0),
                })
        return dict(grouped)

    @staticmethod
    def extract_metal_ions(metal_features: list[dict]) -> list[str]:
        """Extrai símbolos de íons metálicos das features Metal binding."""
        pattern = re.compile(
            r"(Zinc|Zn|Manganese|Mn|Magnesium|Mg|Iron|Fe|"
            r"Cobalt|Co|Copper|Cu|Calcium|Ca|Nickel|Ni)", re.IGNORECASE
        )
        names = {
            "zinc": "Zn", "zn": "Zn", "manganese": "Mn", "mn": "Mn",
            "magnesium": "Mg", "mg": "Mg", "iron": "Fe", "fe": "Fe",
            "cobalt": "Co", "co": "Co", "copper": "Cu", "cu": "Cu",
            "calcium": "Ca", "ca": "Ca", "nickel": "Ni", "ni": "Ni",
        }
        ions: set[str] = set()
        for f in metal_features:
            for m in pattern.findall(f.get("description", "")):
                ions.add(names.get(m.lower(), m))
        return sorted(ions)

    # === API Stubs ===

    async def classify_organism(self, taxonomy_id: int) -> PathogenClass:
        """Consulta NCBI Taxonomy → PathogenClass."""
        raise NotImplementedError("Sprint 3")

    async def fetch_protein(self, uniprot_id: str) -> ProteinRecord:
        """Busca registro UniProt com filtro de IDRs."""
        raise NotImplementedError("Sprint 3")

    async def fetch_ligands(
        self, target_chembl_id: str, *, limit: int = 500
    ) -> pd.DataFrame:
        """Busca compostos ativos no ChEMBL."""
        raise NotImplementedError("Sprint 3")

    async def fetch_pdb_structure(self, pdb_id: str) -> bytes:
        """Baixa PDB/mmCIF do RCSB."""
        raise NotImplementedError("Sprint 3")
