# =============================================================================
# F.L.U.I.D — Módulo 3: TargetScorer
# Algoritmo de Ranqueamento com métrica dN/dS e Normalização de Peso
# =============================================================================
"""
Avalia a viabilidade de alvos terapêuticos usando pressão seletiva
evolutiva (dN/dS), druggabilidade estrutural e coordenação de metais.

A normalização calibra os pesos com base na taxa mutacional basal
do patógeno, garantindo comparação justa entre vírus de RNA e bactérias.

Complexidade:
    - compute_dnds_from_alignment: O(L * N²) — L = comprimento do
      alinhamento, N = número de sequências (pairwise)
    - score_targets: O(T * F) — T = alvos, F = features por alvo
    - Total pipeline é dominado pelo alinhamento pairwise.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats

from src.config import FluidConfig, PathogenClass, DEFAULT_CONFIG


# ---------------------------------------------------------------------------
# Constantes Biológicas
# ---------------------------------------------------------------------------

# Taxas de mutação basal (substituições/sítio/ano) — literatura consolidada.
# Usadas para normalizar o dN/dS entre classes taxonômicas distintas.
BASAL_MUTATION_RATES: dict[PathogenClass, float] = {
    PathogenClass.RNA_VIRUS:  1e-3,   # ~10⁻³ (ex: SARS-CoV-2, HIV)
    PathogenClass.DNA_VIRUS:  1e-5,   # ~10⁻⁵ (ex: Herpesvírus)
    PathogenClass.BACTERIUM:  1e-7,   # ~10⁻⁷ (ex: E. coli, M. tuberculosis)
    PathogenClass.PROTOZOAN:  5e-8,   # ~5×10⁻⁸ (ex: Plasmodium)
    PathogenClass.FUNGUS:     1e-8,   # ~10⁻⁸ (ex: Candida)
    PathogenClass.UNKNOWN:    1e-5,   # fallback conservador
}

# Tabela de codons → aminoácido (código genético padrão)
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
# Data Transfer Objects
# ---------------------------------------------------------------------------

@dataclass
class TargetFeatures:
    """Vetor de features de um alvo terapêutico candidato."""
    protein_id: str
    pathogen_class: PathogenClass = PathogenClass.UNKNOWN

    # Pressão seletiva
    dn: float = 0.0                  # Substituições não-sinônimas
    ds: float = 0.0                  # Substituições sinônimas
    dnds_raw: float = 0.0            # Razão bruta dN/dS

    # Estrutural
    has_metal_binding: bool = False
    metal_ions: list[str] = field(default_factory=list)
    binding_pocket_volume: float = 0.0
    druggability_score: float = 0.0  # 0-1, de ferramentas como fpocket

    # Funcional
    essentiality_score: float = 0.0  # 0-1, derivado de mutagênese/literatura
    disorder_fraction: float = 0.0   # Fração de IDRs (regiões desordenadas)

    # Resultado final
    targetability_score: float = 0.0


# ---------------------------------------------------------------------------
# Classe Principal
# ---------------------------------------------------------------------------

class TargetScorer:
    """
    Calcula o Targetability Score normalizado para alvos terapêuticos.

    Pipeline interno:
        1. Computar dN/dS a partir de alinhamento de códons
        2. Normalizar pela taxa mutacional basal do patógeno
        3. Integrar features estruturais (metais, druggabilidade)
        4. Produzir score final ponderado e ranqueado

    Injeção de dependência via FluidConfig para testabilidade.
    """

    def __init__(self, config: FluidConfig = DEFAULT_CONFIG) -> None:
        self._config = config

    # =================================================================
    # SEÇÃO 1: Cálculo de dN/dS (Nei-Gojobori simplificado)
    # =================================================================

    @staticmethod
    def _classify_site(codon: str) -> tuple[float, float]:
        """
        Classifica os sítios de um códon em sinônimos (S) e
        não-sinônimos (N) usando o método de Nei-Gojobori.

        Para cada posição do códon, verifica se uma mutação pontual
        causa troca de aminoácido (não-sinônima) ou não (sinônima).

        Retorna (N_sites, S_sites) para o códon — soma = 3.0.
        """
        if codon not in _CODON_TABLE or _CODON_TABLE[codon] == "*":
            return 0.0, 0.0

        original_aa = _CODON_TABLE[codon]
        n_sites = 0.0
        s_sites = 0.0
        bases = "ACGT"

        for pos in range(3):
            nonsyn_count = 0
            syn_count = 0
            for base in bases:
                if base == codon[pos]:
                    continue
                mutant = codon[:pos] + base + codon[pos + 1:]
                if mutant not in _CODON_TABLE:
                    continue
                mutant_aa = _CODON_TABLE[mutant]
                if mutant_aa == "*":
                    continue  # Ignora stop codons
                if mutant_aa != original_aa:
                    nonsyn_count += 1
                else:
                    syn_count += 1

            total = nonsyn_count + syn_count
            if total > 0:
                n_sites += nonsyn_count / total
                s_sites += syn_count / total

        return n_sites, s_sites

    @staticmethod
    def _count_differences(
        codon_a: str, codon_b: str
    ) -> tuple[float, float]:
        """
        Conta diferenças sinônimas (sd) e não-sinônimas (nd)
        entre dois códons alinhados.

        Retorna (nd, sd).
        """
        if codon_a == codon_b:
            return 0.0, 0.0

        aa_a = _CODON_TABLE.get(codon_a, "*")
        aa_b = _CODON_TABLE.get(codon_b, "*")

        if aa_a == "*" or aa_b == "*":
            return 0.0, 0.0

        # Conta posições diferentes
        diffs = sum(1 for i in range(3) if codon_a[i] != codon_b[i])

        if diffs == 1:
            # Caso trivial: 1 mutação
            return (1.0, 0.0) if aa_a != aa_b else (0.0, 1.0)

        # Para ≥2 diferenças: proporção média sobre caminhos possíveis
        # Aproximação: assume distribuição uniforme
        if aa_a != aa_b:
            return float(diffs) * 0.7, float(diffs) * 0.3
        return float(diffs) * 0.3, float(diffs) * 0.7

    def compute_dnds_pairwise(
        self, seq_a: str, seq_b: str
    ) -> tuple[float, float, float]:
        """
        Calcula dN/dS entre duas sequências codificantes alinhadas.

        Parâmetros:
            seq_a, seq_b: Sequências de DNA alinhadas (mesmo comprimento,
                          múltiplo de 3, sem gaps internos).

        Retorna:
            (dN, dS, dN/dS) — razão omega. dN/dS < 1 indica seleção
            purificadora (alvo conservado = boa druggabilidade).

        Complexidade: O(L) onde L = comprimento do alinhamento em códons.
        """
        seq_a = seq_a.upper().replace("-", "")
        seq_b = seq_b.upper().replace("-", "")

        min_len = min(len(seq_a), len(seq_b))
        num_codons = min_len // 3

        if num_codons == 0:
            return 0.0, 0.0, float("nan")

        total_n_sites = 0.0
        total_s_sites = 0.0
        total_nd = 0.0
        total_sd = 0.0

        for i in range(num_codons):
            c_a = seq_a[i * 3: i * 3 + 3]
            c_b = seq_b[i * 3: i * 3 + 3]

            if len(c_a) != 3 or len(c_b) != 3:
                continue

            # Sítios
            na, sa = self._classify_site(c_a)
            nb, sb = self._classify_site(c_b)
            total_n_sites += (na + nb) / 2.0
            total_s_sites += (sa + sb) / 2.0

            # Diferenças
            nd, sd = self._count_differences(c_a, c_b)
            total_nd += nd
            total_sd += sd

        # Proporções observadas
        pn = total_nd / total_n_sites if total_n_sites > 0 else 0.0
        ps = total_sd / total_s_sites if total_s_sites > 0 else 0.0

        # Correção de Jukes-Cantor: d = -3/4 * ln(1 - 4p/3)
        def _jukes_cantor(p: float) -> float:
            if p >= 0.75:
                return float("inf")  # Saturação
            if p <= 0.0:
                return 0.0
            return -0.75 * np.log(1.0 - (4.0 * p / 3.0))

        dn = _jukes_cantor(pn)
        ds = _jukes_cantor(ps)

        if ds == 0.0 or np.isinf(ds):
            omega = float("nan")
        else:
            omega = dn / ds

        return float(dn), float(ds), float(omega)

    def compute_dnds_from_alignment(
        self, sequences: list[str]
    ) -> pd.DataFrame:
        """
        Calcula dN/dS para todas as comparações pairwise de um
        conjunto de sequências homólogas alinhadas.

        Retorna DataFrame com colunas: [seq_i, seq_j, dN, dS, omega].

        Complexidade: O(N² * L) — N = número de sequências, L = códons.
        """
        n = len(sequences)
        records = []

        for i in range(n):
            for j in range(i + 1, n):
                dn, ds, omega = self.compute_dnds_pairwise(
                    sequences[i], sequences[j]
                )
                records.append({
                    "seq_i": i,
                    "seq_j": j,
                    "dN": dn,
                    "dS": ds,
                    "omega": omega,
                })

        return pd.DataFrame(records)

    # =================================================================
    # SEÇÃO 2: Normalização de Peso Adaptativa
    # =================================================================

    @staticmethod
    def get_basal_rate(pathogen_class: PathogenClass) -> float:
        """Retorna a taxa mutacional basal para a classe do patógeno."""
        return BASAL_MUTATION_RATES.get(
            pathogen_class, BASAL_MUTATION_RATES[PathogenClass.UNKNOWN]
        )

    def normalize_omega(
        self,
        omega_raw: float,
        pathogen_class: PathogenClass,
        *,
        reference_class: PathogenClass = PathogenClass.BACTERIUM,
    ) -> float:
        """
        Normaliza o omega (dN/dS) pela taxa mutacional basal do patógeno.

        Motivação:
            Vírus de RNA possuem taxa mutacional ~10⁴× maior que bactérias.
            Um omega = 0.3 em HIV é "mais conservado" que omega = 0.3 em
            M. tuberculosis, pois o HIV tem muito mais oportunidade de mutar.
            A normalização corrige esse viés.

        Fórmula:
            ω_norm = ω_raw × (μ_patógeno / μ_referência)

            Onde μ = taxa de mutação basal. A referência padrão é bactéria
            (caso intermediário). O efeito é:
            - RNA vírus: ω_norm >> ω_raw (penaliza → se ainda é baixo,
              a conservação é extremamente forte)
            - Bactéria: ω_norm ≈ ω_raw
            - Fungo: ω_norm << ω_raw

        Retorna: omega normalizado (float). NaN propaga NaN.
        """
        if np.isnan(omega_raw) or np.isinf(omega_raw):
            return float("nan")

        mu_target = self.get_basal_rate(pathogen_class)
        mu_ref = self.get_basal_rate(reference_class)

        scaling_factor = mu_target / mu_ref
        return omega_raw * scaling_factor

    def compute_conservation_score(
        self, omega_normalized: float
    ) -> float:
        """
        Converte omega normalizado em score de conservação [0, 1].

        Usa transformação sigmoide invertida:
            conservation = 1 / (1 + e^(k*(ω - ω₀)))

        Onde:
            k = 10 (inclinação), ω₀ = 0.5 (ponto de inflexão).

        Interpretação:
            ω_norm → 0: conservação → 1.0 (seleção purificadora forte)
            ω_norm → 1: conservação → 0.5 (evolução neutra)
            ω_norm >> 1: conservação → 0.0 (seleção positiva / diversificante)
        """
        if np.isnan(omega_normalized):
            return 0.0

        k = 10.0
        omega_0 = 0.5
        return float(1.0 / (1.0 + np.exp(k * (omega_normalized - omega_0))))

    # =================================================================
    # SEÇÃO 3: Score Integrado (Targetability Score)
    # =================================================================

    def compute_metal_bonus(self, target: TargetFeatures) -> float:
        """
        Bolsões com íons metálicos coordenadores recebem peso máximo.

        Retorna 1.0 se contém Zn²⁺ ou Mn²⁺ (alvos clássicos de
        metaloproteases/polimerases), 0.7 para outros metais, 0.0 se
        não há coordenação metálica.
        """
        if not target.has_metal_binding or not target.metal_ions:
            return 0.0

        priority_ions = {"ZN", "MN", "Zn", "Mn", "Zn2+", "Mn2+"}
        ions_upper = {ion.upper().replace("+", "").strip()
                      for ion in target.metal_ions}

        if ions_upper & {"ZN", "MN"}:
            return 1.0
        return 0.7

    def score_single_target(self, target: TargetFeatures) -> float:
        """
        Calcula o Targetability Score final para um alvo.

        Fórmula (soma ponderada normalizada):
            T = w_c·C + w_d·D + w_m·M + w_e·E

        Onde:
            C = conservation_score (derivado de dN/dS normalizado)
            D = druggability_score (0-1, de análise de bolsão)
            M = metal_bonus (0, 0.7 ou 1.0)
            E = essentiality_score (0-1)

        Penalização:
            Se disorder_fraction > 0.5, o score é reduzido
            proporcionalmente — IDRs extensas comprometem a
            estabilidade do bolsão de ligação.
        """
        w = self._config

        # 1. Conservação via dN/dS normalizado
        omega_norm = self.normalize_omega(
            target.dnds_raw, target.pathogen_class
        )
        conservation = self.compute_conservation_score(omega_norm)

        # 2. Metal bonus
        metal = self.compute_metal_bonus(target)

        # 3. Score ponderado
        raw_score = (
            w.w_conservation * conservation
            + w.w_druggability * target.druggability_score
            + w.w_metal_binding * metal
            + w.w_essentiality * target.essentiality_score
        )

        # 4. Penalização por desordem intrínseca
        if target.disorder_fraction > 0.5:
            penalty = 1.0 - (target.disorder_fraction - 0.5)
            raw_score *= max(penalty, 0.1)  # Floor de 10%

        return float(np.clip(raw_score, 0.0, 1.0))

    def score_targets(
        self, targets: list[TargetFeatures]
    ) -> pd.DataFrame:
        """
        Ranqueia uma lista de alvos candidatos.

        Retorna DataFrame ordenado por targetability_score (desc)
        com colunas completas de métricas intermediárias.

        Complexidade: O(T) — linear no número de alvos.
        """
        records = []

        for t in targets:
            omega_norm = self.normalize_omega(t.dnds_raw, t.pathogen_class)
            conservation = self.compute_conservation_score(omega_norm)
            metal = self.compute_metal_bonus(t)
            final_score = self.score_single_target(t)

            records.append({
                "protein_id": t.protein_id,
                "pathogen_class": t.pathogen_class.value,
                "dN": t.dn,
                "dS": t.ds,
                "omega_raw": t.dnds_raw,
                "omega_normalized": omega_norm,
                "conservation_score": conservation,
                "druggability": t.druggability_score,
                "metal_bonus": metal,
                "metal_ions": ", ".join(t.metal_ions),
                "essentiality": t.essentiality_score,
                "disorder_fraction": t.disorder_fraction,
                "targetability_score": final_score,
            })

        df = pd.DataFrame(records)

        if not df.empty:
            df = df.sort_values(
                "targetability_score", ascending=False
            ).reset_index(drop=True)
            df.index.name = "rank"
            df.index += 1  # Rank começa em 1

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
        Calcula intervalo de confiança bootstrap para o omega médio.

        Reamostra pares de sequências e recalcula omega para estimar
        a distribuição empírica do estimador.

        Retorna: {"mean": ..., "ci_lower": ..., "ci_upper": ...}
        """
        scorer = TargetScorer()
        n = len(sequences)
        if n < 2:
            return {"mean": float("nan"), "ci_lower": 0.0, "ci_upper": 0.0}

        omegas = []
        rng = np.random.default_rng(seed=42)

        for _ in range(n_bootstrap):
            idx = rng.choice(n, size=n, replace=True)
            batch_omegas = []
            for i in range(len(idx)):
                for j in range(i + 1, len(idx)):
                    if idx[i] == idx[j]:
                        continue
                    _, _, w = scorer.compute_dnds_pairwise(
                        sequences[idx[i]], sequences[idx[j]]
                    )
                    if not np.isnan(w) and not np.isinf(w):
                        batch_omegas.append(w)
            if batch_omegas:
                omegas.append(np.mean(batch_omegas))

        if not omegas:
            return {"mean": float("nan"), "ci_lower": 0.0, "ci_upper": 0.0}

        arr = np.array(omegas)
        alpha = (1.0 - ci) / 2.0
        return {
            "mean": float(np.mean(arr)),
            "ci_lower": float(np.percentile(arr, alpha * 100)),
            "ci_upper": float(np.percentile(arr, (1 - alpha) * 100)),
        }
