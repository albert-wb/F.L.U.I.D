# =============================================================================
# F.L.U.I.D — Módulo 3: TargetScorer
# Métrica dN/dS + Normalização Adaptativa + Score Final Especificado
# =============================================================================
"""
Avalia viabilidade de alvos terapêuticos com fórmula matemática exata:

    T = [ Σᵢ(wᵢ · fᵢ) ] × Φ_metal × Ψ_disorder

Onde:
    fᵢ ∈ {C, D, E}  — features normalizadas [0,1]
    Φ_metal          — multiplicador metálico [1.0, 1.5, 2.0]
    Ψ_disorder       — fator de penalização por IDRs [0.1, 1.0]

O resultado final é normalizado para [0, 1] via min-max no batch.

Complexidade:
    - compute_dnds_from_alignment: O(N² × L) (mitigado por clustering)
    - score_targets: O(T) linear
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from src.config import FluidConfig, PathogenClass, DEFAULT_CONFIG


# ---------------------------------------------------------------------------
# Constantes
# ---------------------------------------------------------------------------
BASAL_MUTATION_RATES: dict[PathogenClass, float] = {
    PathogenClass.RNA_VIRUS:  1e-3,
    PathogenClass.DNA_VIRUS:  1e-5,
    PathogenClass.BACTERIUM:  1e-7,
    PathogenClass.PROTOZOAN:  5e-8,
    PathogenClass.FUNGUS:     1e-8,
    PathogenClass.UNKNOWN:    1e-5,
}

_CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


# ---------------------------------------------------------------------------
# DTOs
# ---------------------------------------------------------------------------
@dataclass
class TargetFeatures:
    """Vetor de features de um alvo terapêutico candidato."""
    protein_id: str
    pathogen_class: PathogenClass = PathogenClass.UNKNOWN

    # Pressão seletiva
    dn: float = 0.0
    ds: float = 0.0
    dnds_raw: float = 0.0

    # Estrutural
    has_metal_binding: bool = False
    metal_ions: list[str] = field(default_factory=list)
    binding_pocket_volume: float = 0.0
    druggability_score: float = 0.0   # [0,1]

    # Funcional
    essentiality_score: float = 0.0   # [0,1]
    disorder_fraction: float = 0.0    # [0,1]

    # Resultado (preenchido pelo scorer)
    targetability_score: float = 0.0


# ---------------------------------------------------------------------------
# Classe Principal
# ---------------------------------------------------------------------------
class TargetScorer:
    """
    Calcula o Targetability Score com fórmula matemática exata.

    Score Final — Especificação Formal:
    ─────────────────────────────────────────────────────────────
    Seja o vetor de features f = (C, D, E) e pesos w = (w_c, w_d, w_e):

        S_raw = Σᵢ (wᵢ · fᵢ)          ... soma ponderada base

    Multiplicador metálico (Φ):
        Φ = 2.0  se {Zn²⁺, Mn²⁺} ∈ metal_ions
        Φ = 1.5  se outros metais ∈ metal_ions
        Φ = 1.0  se não há coordenação metálica

    Fator de penalização por desordem (Ψ):
        Ψ = 1.0                  se disorder_frac ≤ 0.3
        Ψ = 1.0 - (d - 0.3)     se 0.3 < disorder_frac ≤ 0.9
        Ψ = 0.1                  se disorder_frac > 0.9  (floor)

    Score final antes de normalização:
        T_raw = S_raw × Φ × Ψ

    Normalização min-max no batch:
        T = (T_raw - T_min) / (T_max - T_min)
    ─────────────────────────────────────────────────────────────
    """

    def __init__(self, config: FluidConfig = DEFAULT_CONFIG) -> None:
        self._config = config

    # =================================================================
    # SEÇÃO 1: dN/dS (Nei-Gojobori simplificado)
    # =================================================================

    @staticmethod
    def _classify_site(codon: str) -> tuple[float, float]:
        """Classifica sítios de um códon em (N_sites, S_sites)."""
        if codon not in _CODON_TABLE or _CODON_TABLE[codon] == "*":
            return 0.0, 0.0

        aa = _CODON_TABLE[codon]
        n_sites, s_sites = 0.0, 0.0

        for pos in range(3):
            ns, sy = 0, 0
            for base in "ACGT":
                if base == codon[pos]:
                    continue
                mut = codon[:pos] + base + codon[pos + 1:]
                if mut not in _CODON_TABLE:
                    continue
                m_aa = _CODON_TABLE[mut]
                if m_aa == "*":
                    continue
                if m_aa != aa:
                    ns += 1
                else:
                    sy += 1
            total = ns + sy
            if total > 0:
                n_sites += ns / total
                s_sites += sy / total

        return n_sites, s_sites

    @staticmethod
    def _count_differences(ca: str, cb: str) -> tuple[float, float]:
        """Conta diferenças (nd, sd) entre dois códons."""
        if ca == cb:
            return 0.0, 0.0
        aa_a = _CODON_TABLE.get(ca, "*")
        aa_b = _CODON_TABLE.get(cb, "*")
        if aa_a == "*" or aa_b == "*":
            return 0.0, 0.0

        diffs = sum(1 for i in range(3) if ca[i] != cb[i])
        if diffs == 1:
            return (1.0, 0.0) if aa_a != aa_b else (0.0, 1.0)
        # Múltiplas diferenças — aproximação proporcional
        if aa_a != aa_b:
            return float(diffs) * 0.7, float(diffs) * 0.3
        return float(diffs) * 0.3, float(diffs) * 0.7

    @staticmethod
    def _jukes_cantor(p: float) -> float:
        """Correção JC69: d = -3/4 × ln(1 - 4p/3)."""
        if p >= 0.75:
            return float("inf")
        if p <= 0.0:
            return 0.0
        return -0.75 * np.log(1.0 - (4.0 * p / 3.0))

    def compute_dnds_pairwise(
        self, seq_a: str, seq_b: str
    ) -> tuple[float, float, float]:
        """
        dN/dS entre duas sequências codificantes alinhadas.
        Retorna (dN, dS, omega). Complexidade: O(L).
        """
        seq_a = seq_a.upper().replace("-", "")
        seq_b = seq_b.upper().replace("-", "")
        num_codons = min(len(seq_a), len(seq_b)) // 3

        if num_codons == 0:
            return 0.0, 0.0, float("nan")

        t_n, t_s, t_nd, t_sd = 0.0, 0.0, 0.0, 0.0

        for i in range(num_codons):
            ca = seq_a[i*3:i*3+3]
            cb = seq_b[i*3:i*3+3]
            if len(ca) != 3 or len(cb) != 3:
                continue

            na, sa = self._classify_site(ca)
            nb, sb = self._classify_site(cb)
            t_n += (na + nb) / 2.0
            t_s += (sa + sb) / 2.0

            nd, sd = self._count_differences(ca, cb)
            t_nd += nd
            t_sd += sd

        pn = t_nd / t_n if t_n > 0 else 0.0
        ps = t_sd / t_s if t_s > 0 else 0.0

        dn = self._jukes_cantor(pn)
        ds = self._jukes_cantor(ps)

        if ds == 0.0 or np.isinf(ds):
            omega = float("nan")
        else:
            omega = dn / ds

        return float(dn), float(ds), float(omega)

    def compute_dnds_from_alignment(
        self, sequences: list[str]
    ) -> pd.DataFrame:
        """
        dN/dS pairwise para N sequências.
        Complexidade: O(N² × L). Use clustering antes para reduzir N.
        """
        records = []
        n = len(sequences)
        for i in range(n):
            for j in range(i + 1, n):
                dn, ds, omega = self.compute_dnds_pairwise(
                    sequences[i], sequences[j]
                )
                records.append({
                    "seq_i": i, "seq_j": j,
                    "dN": dn, "dS": ds, "omega": omega,
                })
        return pd.DataFrame(records)

    # =================================================================
    # SEÇÃO 2: Normalização Adaptativa
    # =================================================================

    @staticmethod
    def get_basal_rate(pc: PathogenClass) -> float:
        """Taxa mutacional basal μ para a classe do patógeno."""
        return BASAL_MUTATION_RATES.get(
            pc, BASAL_MUTATION_RATES[PathogenClass.UNKNOWN]
        )

    def normalize_omega(
        self,
        omega_raw: float,
        pathogen_class: PathogenClass,
        *,
        reference_class: PathogenClass = PathogenClass.BACTERIUM,
    ) -> float:
        """
        ω_norm = ω_raw × (μ_patógeno / μ_referência)

        RNA vírus → ω amplificado (se ainda baixo = conservação forte).
        Bactéria → ω ≈ ω_raw (referência).
        """
        if np.isnan(omega_raw) or np.isinf(omega_raw):
            return float("nan")
        mu_t = self.get_basal_rate(pathogen_class)
        mu_r = self.get_basal_rate(reference_class)
        return omega_raw * (mu_t / mu_r)

    def compute_conservation_score(self, omega_norm: float) -> float:
        """
        Sigmoide invertida: C = 1 / (1 + e^(k×(ω-ω₀)))
        k=10, ω₀=0.5. Mapeia ω → [0,1] (alto ω = baixa conservação).
        """
        if np.isnan(omega_norm):
            return 0.0
        return float(1.0 / (1.0 + np.exp(10.0 * (omega_norm - 0.5))))

    # =================================================================
    # SEÇÃO 3: Score Final (Fórmula Matemática Exata)
    # =================================================================

    @staticmethod
    def compute_metal_multiplier(target: TargetFeatures) -> float:
        """
        Multiplicador metálico Φ:
            Φ = 2.0 se {Zn²⁺, Mn²⁺} ∈ metal_ions (prioridade catalítica)
            Φ = 1.5 se outros metais coordenadores presentes
            Φ = 1.0 se sem coordenação metálica
        """
        if not target.has_metal_binding or not target.metal_ions:
            return 1.0

        normalized = {
            ion.upper().replace("+", "").replace("2", "").strip()
            for ion in target.metal_ions
        }
        if normalized & {"ZN", "MN"}:
            return 2.0
        return 1.5

    @staticmethod
    def compute_disorder_penalty(disorder_fraction: float) -> float:
        """
        Fator de penalização por desordem intrínseca Ψ:
            Ψ = 1.0                    se d ≤ 0.3
            Ψ = 1.0 - (d - 0.3)       se 0.3 < d ≤ 0.9
            Ψ = 0.1                    se d > 0.9  (floor 10%)
        """
        if disorder_fraction <= 0.3:
            return 1.0
        if disorder_fraction > 0.9:
            return 0.1
        return 1.0 - (disorder_fraction - 0.3)

    def score_single_target(self, target: TargetFeatures) -> float:
        """
        Score Final (antes de normalização min-max):

            T_raw = [ w_c·C + w_d·D + w_e·E ] × Φ × Ψ

        Onde:
            C = conservation_score(normalize_omega(ω_raw, pathogen_class))
            D = druggability_score              [0,1]
            E = essentiality_score              [0,1]
            Φ = metal_multiplier                [1.0, 1.5, 2.0]
            Ψ = disorder_penalty                [0.1, 1.0]
        """
        cfg = self._config

        # Feature vector
        omega_norm = self.normalize_omega(
            target.dnds_raw, target.pathogen_class
        )
        C = self.compute_conservation_score(omega_norm)
        D = target.druggability_score
        E = target.essentiality_score

        # Soma ponderada base
        s_raw = (
            cfg.w_conservation * C
            + cfg.w_druggability * D
            + cfg.w_essentiality * E
        )

        # Multiplicadores
        phi = self.compute_metal_multiplier(target)
        psi = self.compute_disorder_penalty(target.disorder_fraction)

        return float(s_raw * phi * psi)

    def score_targets(
        self, targets: list[TargetFeatures]
    ) -> pd.DataFrame:
        """
        Ranqueia alvos com normalização min-max no batch.

            T = (T_raw - T_min) / (T_max - T_min)

        Retorna DataFrame ordenado por targetability_score desc.
        Complexidade: O(T) linear.
        """
        records = []
        raw_scores: list[float] = []

        for t in targets:
            omega_norm = self.normalize_omega(t.dnds_raw, t.pathogen_class)
            C = self.compute_conservation_score(omega_norm)
            phi = self.compute_metal_multiplier(t)
            psi = self.compute_disorder_penalty(t.disorder_fraction)
            t_raw = self.score_single_target(t)
            raw_scores.append(t_raw)

            records.append({
                "protein_id": t.protein_id,
                "pathogen_class": t.pathogen_class.value,
                "dN": t.dn,
                "dS": t.ds,
                "omega_raw": t.dnds_raw,
                "omega_normalized": omega_norm,
                "conservation_C": C,
                "druggability_D": t.druggability_score,
                "essentiality_E": t.essentiality_score,
                "metal_multiplier_Phi": phi,
                "metal_ions": ", ".join(t.metal_ions),
                "disorder_fraction": t.disorder_fraction,
                "disorder_penalty_Psi": psi,
                "score_raw": t_raw,
            })

        df = pd.DataFrame(records)

        if not df.empty:
            # Normalização min-max no batch
            s_min = df["score_raw"].min()
            s_max = df["score_raw"].max()
            denom = s_max - s_min

            if denom > 1e-12:
                df["targetability_score"] = (df["score_raw"] - s_min) / denom
            else:
                # Todos os scores iguais → 1.0
                df["targetability_score"] = 1.0

            df = df.sort_values(
                "targetability_score", ascending=False
            ).reset_index(drop=True)
            df.index.name = "rank"
            df.index += 1

        return df

    # =================================================================
    # SEÇÃO 4: Utilitários Estatísticos
    # =================================================================

    @staticmethod
    def bootstrap_omega_ci(
        sequences: list[str],
        n_bootstrap: int = 1000,
        ci: float = 0.95,
    ) -> dict[str, float]:
        """
        Intervalo de confiança bootstrap para omega médio.
        Retorna {"mean", "ci_lower", "ci_upper"}.
        """
        scorer = TargetScorer()
        n = len(sequences)
        if n < 2:
            return {"mean": float("nan"), "ci_lower": 0.0, "ci_upper": 0.0}

        omegas: list[float] = []
        rng = np.random.default_rng(seed=42)

        for _ in range(n_bootstrap):
            idx = rng.choice(n, size=n, replace=True)
            batch = []
            for i in range(len(idx)):
                for j in range(i + 1, len(idx)):
                    if idx[i] == idx[j]:
                        continue
                    _, _, w = scorer.compute_dnds_pairwise(
                        sequences[idx[i]], sequences[idx[j]]
                    )
                    if not np.isnan(w) and not np.isinf(w):
                        batch.append(w)
            if batch:
                omegas.append(float(np.mean(batch)))

        if not omegas:
            return {"mean": float("nan"), "ci_lower": 0.0, "ci_upper": 0.0}

        arr = np.array(omegas)
        alpha = (1.0 - ci) / 2.0
        return {
            "mean": float(np.mean(arr)),
            "ci_lower": float(np.percentile(arr, alpha * 100)),
            "ci_upper": float(np.percentile(arr, (1 - alpha) * 100)),
        }
