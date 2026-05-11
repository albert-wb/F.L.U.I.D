# =============================================================================
# Testes unitários — TargetScorer
# =============================================================================
import math
import numpy as np
import pandas as pd
import pytest

from src.config import PathogenClass, FluidConfig
from src.scoring.target_scorer import (
    TargetScorer,
    TargetFeatures,
    BASAL_MUTATION_RATES,
)


@pytest.fixture
def scorer() -> TargetScorer:
    return TargetScorer()


class TestDnDsPairwise:
    """Valida o cálculo de dN/dS entre pares de sequências."""

    def test_identical_sequences(self, scorer: TargetScorer):
        seq = "ATGATGATGATG"
        dn, ds, omega = scorer.compute_dnds_pairwise(seq, seq)
        assert dn == 0.0
        assert ds == 0.0
        assert math.isnan(omega)

    def test_synonymous_only(self, scorer: TargetScorer):
        # TTT (Phe) → TTC (Phe) — sinônima
        seq_a = "TTT"
        seq_b = "TTC"
        dn, ds, omega = scorer.compute_dnds_pairwise(seq_a, seq_b)
        assert dn == 0.0
        assert ds > 0.0
        assert omega == 0.0

    def test_nonsynonymous_mutation(self, scorer: TargetScorer):
        # AAA (Lys) → GAA (Glu) — não-sinônima
        seq_a = "AAA"
        seq_b = "GAA"
        dn, ds, omega = scorer.compute_dnds_pairwise(seq_a, seq_b)
        assert dn > 0.0

    def test_empty_sequences(self, scorer: TargetScorer):
        dn, ds, omega = scorer.compute_dnds_pairwise("", "")
        assert math.isnan(omega)


class TestNormalization:
    """Valida a normalização de omega por taxa mutacional."""

    def test_rna_virus_amplifies_omega(self, scorer: TargetScorer):
        omega_raw = 0.3
        norm = scorer.normalize_omega(omega_raw, PathogenClass.RNA_VIRUS)
        # RNA vírus tem taxa ~10⁴× maior que bactéria (ref)
        assert norm > omega_raw

    def test_bacterium_is_identity(self, scorer: TargetScorer):
        omega_raw = 0.5
        norm = scorer.normalize_omega(
            omega_raw, PathogenClass.BACTERIUM,
            reference_class=PathogenClass.BACTERIUM,
        )
        assert abs(norm - omega_raw) < 1e-10

    def test_nan_propagates(self, scorer: TargetScorer):
        assert math.isnan(
            scorer.normalize_omega(float("nan"), PathogenClass.RNA_VIRUS)
        )


class TestConservationScore:
    """Valida a sigmoide invertida de conservação."""

    def test_low_omega_high_conservation(self, scorer: TargetScorer):
        assert scorer.compute_conservation_score(0.01) > 0.95

    def test_high_omega_low_conservation(self, scorer: TargetScorer):
        assert scorer.compute_conservation_score(2.0) < 0.05

    def test_inflection_point(self, scorer: TargetScorer):
        assert abs(scorer.compute_conservation_score(0.5) - 0.5) < 0.01


class TestMetalBonus:
    """Valida a bonificação por coordenação metálica."""

    def test_zinc_max_bonus(self, scorer: TargetScorer):
        t = TargetFeatures(
            protein_id="test", has_metal_binding=True, metal_ions=["Zn2+"]
        )
        assert scorer.compute_metal_bonus(t) == 1.0

    def test_other_metal_partial(self, scorer: TargetScorer):
        t = TargetFeatures(
            protein_id="test", has_metal_binding=True, metal_ions=["Fe"]
        )
        assert scorer.compute_metal_bonus(t) == 0.7

    def test_no_metal_zero(self, scorer: TargetScorer):
        t = TargetFeatures(protein_id="test", has_metal_binding=False)
        assert scorer.compute_metal_bonus(t) == 0.0


class TestScoreTargets:
    """Valida o ranqueamento integrado."""

    def test_ranking_order(self, scorer: TargetScorer):
        targets = [
            TargetFeatures(
                protein_id="weak",
                pathogen_class=PathogenClass.BACTERIUM,
                dnds_raw=1.5,
                druggability_score=0.2,
                essentiality_score=0.1,
            ),
            TargetFeatures(
                protein_id="strong",
                pathogen_class=PathogenClass.BACTERIUM,
                dnds_raw=0.05,
                has_metal_binding=True,
                metal_ions=["Zn2+"],
                druggability_score=0.9,
                essentiality_score=0.95,
            ),
        ]
        df = scorer.score_targets(targets)
        assert df.iloc[0]["protein_id"] == "strong"
        assert df.iloc[0]["targetability_score"] > df.iloc[1]["targetability_score"]

    def test_disorder_penalty(self, scorer: TargetScorer):
        base = TargetFeatures(
            protein_id="ordered",
            pathogen_class=PathogenClass.BACTERIUM,
            dnds_raw=0.1,
            druggability_score=0.8,
            essentiality_score=0.8,
            disorder_fraction=0.1,
        )
        disordered = TargetFeatures(
            protein_id="disordered",
            pathogen_class=PathogenClass.BACTERIUM,
            dnds_raw=0.1,
            druggability_score=0.8,
            essentiality_score=0.8,
            disorder_fraction=0.8,
        )
        s1 = scorer.score_single_target(base)
        s2 = scorer.score_single_target(disordered)
        assert s1 > s2
