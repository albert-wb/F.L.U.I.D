# =============================================================================
# F.L.U.I.D — Dashboard: Visualizer
# =============================================================================
"""
Componente de visualização científica e animação de fluidos 2D.
Stub para implementação futura.
"""
from __future__ import annotations

from src.config import FluidConfig, DEFAULT_CONFIG


class Visualizer:
    """
    Gerencia o dashboard Dash/Plotly e injeta o componente JavaScript
    de simulação de fluido 2D (matter.js / p5.js).

    Responsabilidades:
        - Layout do dashboard com Dash Bootstrap Components
        - Componente customizado de animação de docking em fluido
        - Gráficos interativos (Plotly) de ranqueamento e resultados
    """

    def __init__(self, config: FluidConfig = DEFAULT_CONFIG) -> None:
        self._config = config

    def create_app(self):
        """Inicializa a aplicação Dash."""
        raise NotImplementedError

    def register_callbacks(self):
        """Registra callbacks reativos do Dash."""
        raise NotImplementedError
