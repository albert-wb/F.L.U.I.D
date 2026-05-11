# =============================================================================
# Testes — TargetScorer (Sprint 2 — Fórmula Exata)
# =============================================================================
import math
import numpy as np
import pytest
from src.config import PathogenClass
from src.scoring.target_scorer import TargetScorer, TargetFeatures


@pytest.fixture
def scorer() -> TargetScorer:
    return TargetScorer()


class TestDnDsPairwise:
    def test_identical(self, scorer: TargetScorer):
        seq = "ATGATGATGATG"
        dn, ds, w = scorer.compute_dnds_pairwise(seq, seq)
        assert dn == 0.0 and ds == 0.0 and math.isnan(w)

    def test_synonymous(self, scorer: TargetScorer):
        dn, ds, w = scorer.compute_dnds_pairwise("TTT", "TTC")
        assert dn == 0.0 and ds > 0.0 and w == 0.0

    def test_nonsynonymous(self, scorer: TargetScorer):
        dn, ds, w = scorer.compute_dnds_pairwise("AAA", "GAA")
        assert dn > 0.0


class TestNormalization:
    def test_rna_amplifies(self, scorer: TargetScorer):
        norm = scorer.normalize_omega(0.3, PathogenClass.RNA_VIRUS)
        assert norm > 0.3

    def test_bacterium_identity(self, scorer: TargetScorer):
        norm = scorer.normalize_omega(
            0.5, PathogenClass.BACTERIUM,
            reference_class=PathogenClass.BACTERIUM,
        )
        assert abs(norm - 0.5) < 1e-10

    def test_nan_propagates(self, scorer: TargetScorer):
        assert math.isnan(scorer.normalize_omega(float("nan"), PathogenClass.RNA_VIRUS))


class TestConservation:
    def test_low_omega_high_conservation(self, scorer: TargetScorer):
        assert scorer.compute_conservation_score(0.01) > 0.95

    def test_high_omega_low_conservation(self, scorer: TargetScorer):
        assert scorer.compute_conservation_score(2.0) < 0.05


class TestMetalMultiplier:
    def test_zinc_max(self, scorer: TargetScorer):
        t = TargetFeatures("x", has_metal_binding=True, metal_ions=["Zn2+"])
        assert scorer.compute_metal_multiplier(t) == 2.0

    def test_other_metal(self, scorer: TargetScorer):
        t = TargetFeatures("x", has_metal_binding=True, metal_ions=["Fe"])
        assert scorer.compute_metal_multiplier(t) == 1.5

    def test_no_metal(self, scorer: TargetScorer):
        t = TargetFeatures("x", has_metal_binding=False)
        assert scorer.compute_metal_multiplier(t) == 1.0


class TestDisorderPenalty:
    def test_low_disorder_no_penalty(self):
        assert TargetScorer.compute_disorder_penalty(0.2) == 1.0

    def test_high_disorder_penalized(self):
        psi = TargetScorer.compute_disorder_penalty(0.6)
        assert 0.0 < psi < 1.0

    def test_extreme_disorder_floor(self):
        assert TargetScorer.compute_disorder_penalty(0.95) == 0.1


class TestScoreFormula:
    def test_ranking_order(self, scorer: TargetScorer):
        targets = [
            TargetFeatures(
                "weak", pathogen_class=PathogenClass.BACTERIUM,
                dnds_raw=1.5, druggability_score=0.2, essentiality_score=0.1,
            ),
            TargetFeatures(
                "strong", pathogen_class=PathogenClass.BACTERIUM,
                dnds_raw=0.05, has_metal_binding=True, metal_ions=["Zn2+"],
                druggability_score=0.9, essentiality_score=0.95,
            ),
        ]
        df = scorer.score_targets(targets)
        assert df.iloc[0]["protein_id"] == "strong"

    def test_metal_multiplier_boosts(self, scorer: TargetScorer):
        base = TargetFeatures(
            "no_metal", pathogen_class=PathogenClass.BACTERIUM,
            dnds_raw=0.1, druggability_score=0.8, essentiality_score=0.8,
        )
        metal = TargetFeatures(
            "with_zn", pathogen_class=PathogenClass.BACTERIUM,
            dnds_raw=0.1, has_metal_binding=True, metal_ions=["Zn"],
            druggability_score=0.8, essentiality_score=0.8,
        )
        assert scorer.score_single_target(metal) > scorer.score_single_target(base)

    def test_min_max_normalization(self, scorer: TargetScorer):
        targets = [
            TargetFeatures("a", dnds_raw=0.1, druggability_score=0.5, essentiality_score=0.5),
            TargetFeatures("b", dnds_raw=0.9, druggability_score=0.1, essentiality_score=0.1),
        ]
        df = scorer.score_targets(targets)
        assert df["targetability_score"].max() == 1.0
        assert df["targetability_score"].min() == 0.0
