# =============================================================================
# F.L.U.I.D — Configurações Globais
# =============================================================================
from __future__ import annotations

import os
from pathlib import Path
from enum import Enum
from dataclasses import dataclass, field

from dotenv import load_dotenv

load_dotenv()


# ---------------------------------------------------------------------------
# Enums de Domínio
# ---------------------------------------------------------------------------
class PathogenClass(Enum):
    """Classificação taxonômica de alto nível do patógeno."""
    RNA_VIRUS = "rna_virus"
    DNA_VIRUS = "dna_virus"
    BACTERIUM = "bacterium"
    PROTOZOAN = "protozoan"
    FUNGUS = "fungus"
    UNKNOWN = "unknown"


class EnzymeFamilyType(Enum):
    """Famílias enzimáticas suportadas pelo pipeline."""
    PROTEASE = "protease"
    POLYMERASE = "polymerase"
    HELICASE = "helicase"
    INTEGRASE = "integrase"
    KINASE = "kinase"
    TRANSFERASE = "transferase"
    LIGASE = "ligase"
    HYDROLASE = "hydrolase"
    OTHER = "other"


# ---------------------------------------------------------------------------
# Configuração Central
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class FluidConfig:
    """Objeto imutável de configuração injetável em todos os módulos."""

    # Paths
    workspace_dir: Path = field(
        default_factory=lambda: Path(os.getenv("WORKSPACE_DIR", "./workspace"))
    )
    vina_binary: Path = field(
        default_factory=lambda: Path(
            os.getenv("VINA_BINARY_PATH", "vina")
        )
    )

    # API Keys
    openai_api_key: str = field(
        default_factory=lambda: os.getenv("OPENAI_API_KEY", "")
    )
    ncbi_api_key: str = field(
        default_factory=lambda: os.getenv("NCBI_API_KEY", "")
    )

    # Redis / Celery
    redis_url: str = field(
        default_factory=lambda: os.getenv("REDIS_URL", "redis://localhost:6379/0")
    )

    # Logging
    log_level: str = field(
        default_factory=lambda: os.getenv("LOG_LEVEL", "INFO")
    )

    # Docking defaults
    exhaustiveness: int = 32
    num_modes: int = 9
    energy_range: float = 3.0

    # Scoring weights (defaults — override via normalização)
    w_conservation: float = 0.35
    w_druggability: float = 0.25
    w_metal_binding: float = 0.20
    w_essentiality: float = 0.20

    def __post_init__(self) -> None:
        """Garante que o workspace exista."""
        self.workspace_dir.mkdir(parents=True, exist_ok=True)


# Singleton de conveniência
DEFAULT_CONFIG = FluidConfig()
