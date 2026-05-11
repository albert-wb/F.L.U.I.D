# =============================================================================
# Testes — DataIngestor (Clustering)
# =============================================================================
import pytest
from src.ingestion.data_ingestor import DataIngestor, SequenceCluster


@pytest.fixture
def ingestor() -> DataIngestor:
    return DataIngestor()


class TestClustering:
    def test_identical_sequences_cluster_together(self, ingestor: DataIngestor):
        seqs = ["ATGATGATG"] * 5
        ids = [f"s{i}" for i in range(5)]
        clusters = ingestor.cluster_sequences(seqs, ids)
        assert len(clusters) == 1
        assert clusters[0].size == 5

    def test_distinct_sequences_separate(self, ingestor: DataIngestor):
        s1 = "ATGATGATGATGATGATG"
        s2 = "GCGCGCGCGCGCGCGCGC"
        clusters = ingestor.cluster_sequences([s1, s2], ["a", "b"])
        assert len(clusters) == 2

    def test_high_identity_clusters(self, ingestor: DataIngestor):
        base = "ATGATGATGATGATGATGATGATGATG"
        # 1 mutação → >95% identidade
        mut = base[:5] + "C" + base[6:]
        clusters = ingestor.cluster_sequences(
            [base, mut], ["ref", "mut"], identity_threshold=0.90
        )
        assert len(clusters) == 1

    def test_empty_input(self, ingestor: DataIngestor):
        assert ingestor.cluster_sequences([]) == []

    def test_get_centroids_returns_subset(self, ingestor: DataIngestor):
        seqs = ["ATGATGATGATGATGATG"] * 10 + ["GCGCGCGCGCGCGCGCGC"]
        c_seqs, c_ids = ingestor.get_centroid_sequences(seqs)
        assert len(c_seqs) == 2
        assert len(c_ids) == 2


class TestDisorderFilter:
    def test_fraction_calculation(self):
        seq = "A" * 100
        regions = [(1, 30)]
        assert abs(DataIngestor.compute_disorder_fraction(seq, regions) - 0.30) < 0.01

    def test_empty_regions(self):
        assert DataIngestor.compute_disorder_fraction("ACGT", []) == 0.0


class TestMetalIonExtraction:
    def test_zinc_extraction(self):
        features = [{"description": "Zinc; catalytic"}]
        ions = DataIngestor.extract_metal_ions(features)
        assert "Zn" in ions

    def test_multiple_ions(self):
        features = [
            {"description": "Zinc ion"},
            {"description": "Manganese cofactor"},
        ]
        ions = DataIngestor.extract_metal_ions(features)
        assert "Zn" in ions and "Mn" in ions
