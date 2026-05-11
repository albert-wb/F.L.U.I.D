# =============================================================================
# F.L.U.I.D — Dashboard Dash (Visualizer)
# =============================================================================
"""
Dashboard científico com:
    - Simulação de fluido 2D (matter.js + p5.js) via Iframe component
    - Gráficos interativos de ranqueamento (Plotly)
    - Formulário de submissão de jobs
    - Polling de status via API

Injeção do componente JavaScript:
    O motor de física é carregado via html.Iframe apontando para um
    HTML estático que importa matter.js, p5.js e fluid_engine.js.
    A comunicação Dash ↔ JS é feita via window.postMessage.
"""
from __future__ import annotations

from pathlib import Path

import dash
from dash import html, dcc, Input, Output, State, callback, no_update
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import requests
from loguru import logger

from src.config import FluidConfig, DEFAULT_CONFIG


# ---------------------------------------------------------------------------
# HTML do componente de fluido (carregado via Iframe)
# ---------------------------------------------------------------------------
FLUID_HTML = """<!DOCTYPE html>
<html lang="pt-BR">
<head>
<meta charset="UTF-8">
<title>F.L.U.I.D — Simulação de Docking</title>
<style>
  * { margin: 0; padding: 0; box-sizing: border-box; }
  body { background: #0a0e1a; overflow: hidden; }
  canvas { display: block; }
</style>
<script src="https://cdnjs.cloudflare.com/ajax/libs/matter-js/0.19.0/matter.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/p5.js/1.9.4/p5.min.js"></script>
</head>
<body>
<div id="fluid-container"></div>
<script>
%s
</script>
<script>
// Listener para comunicação com o Dash via postMessage
window.addEventListener('message', function(event) {
    try {
        const msg = event.data;
        if (!msg || !msg.action) return;

        switch (msg.action) {
            case 'addProtein':
                window.fluidAPI.addProtein(
                    msg.id, msg.label,
                    msg.x || 450, msg.y || 300, msg.size || 50
                );
                break;
            case 'addLigand':
                window.fluidAPI.addLigand(
                    msg.id, msg.label, msg.svg,
                    msg.x, msg.y, msg.radius
                );
                break;
            case 'submitResults':
                window.fluidAPI.submitResults(msg.results);
                break;
            case 'reset':
                window.fluidAPI.reset();
                break;
            case 'getState':
                window.parent.postMessage({
                    action: 'stateUpdate',
                    state: window.fluidAPI.getState()
                }, '*');
                break;
        }
    } catch (e) { console.error('Message handler error:', e); }
});
</script>
</body>
</html>"""


class Visualizer:
    """
    Dashboard Dash com visualização científica e animação de fluidos 2D.

    Componentes:
        1. Painel de submissão (protein ID, pathogen class, SMILES)
        2. Simulação de fluido (matter.js + p5.js via Iframe)
        3. Tabela de ranqueamento (Plotly)
        4. Gráfico de afinidades (scatter plot ΔG vs Score)
    """

    def __init__(self, config: FluidConfig = DEFAULT_CONFIG) -> None:
        self._config = config
        self._api_base = "http://localhost:8000"
        self._app: dash.Dash | None = None

    def create_app(self) -> dash.Dash:
        """Inicializa a aplicação Dash com layout completo."""
        self._app = dash.Dash(
            __name__,
            external_stylesheets=[dbc.themes.CYBORG],
            title="F.L.U.I.D — Drug Discovery Platform",
            suppress_callback_exceptions=True,
        )

        self._app.layout = self._build_layout()
        self._register_callbacks()
        return self._app

    def _get_fluid_html(self) -> str:
        """Lê o JS da engine e injeta no HTML do Iframe."""
        js_path = Path(__file__).parent / "assets" / "fluid_engine.js"
        try:
            js_code = js_path.read_text(encoding="utf-8")
        except FileNotFoundError:
            js_code = "// fluid_engine.js não encontrado"
            logger.error(f"fluid_engine.js não encontrado em {js_path}")
        return FLUID_HTML % js_code

    def _build_layout(self) -> dbc.Container:
        """Constrói o layout do dashboard."""
        return dbc.Container([
            # --- Header ---
            dbc.Row([
                dbc.Col([
                    html.H1(
                        "F.L.U.I.D",
                        className="text-center mt-3",
                        style={"color": "#4fc3f7", "fontWeight": "300",
                               "letterSpacing": "8px", "fontSize": "2.5rem"},
                    ),
                    html.P(
                        "Framework for Ligand-target Unified In-silico Discovery",
                        className="text-center text-muted mb-4",
                        style={"letterSpacing": "2px", "fontSize": "0.85rem"},
                    ),
                ])
            ]),

            # --- Main Content ---
            dbc.Row([
                # Coluna Esquerda: Controles + Status
                dbc.Col([
                    # Formulário de submissão
                    dbc.Card([
                        dbc.CardHeader("🧬 Submissão de Alvo", style={"color": "#4fc3f7"}),
                        dbc.CardBody([
                            dbc.Label("UniProt ID", html_for="input-protein"),
                            dbc.Input(id="input-protein", placeholder="ex: P0DTD1",
                                      type="text", className="mb-2"),
                            dbc.Label("Nome da Proteína", html_for="input-name"),
                            dbc.Input(id="input-name", placeholder="ex: 3CLpro",
                                      type="text", className="mb-2"),
                            dbc.Label("Classe do Patógeno", html_for="input-pathogen"),
                            dbc.Select(id="input-pathogen", options=[
                                {"label": "Vírus RNA", "value": "rna_virus"},
                                {"label": "Vírus DNA", "value": "dna_virus"},
                                {"label": "Bactéria", "value": "bacterium"},
                                {"label": "Protozoário", "value": "protozoan"},
                                {"label": "Fungo", "value": "fungus"},
                            ], value="rna_virus", className="mb-2"),
                            dbc.Label("SMILES (um por linha)", html_for="input-smiles"),
                            dbc.Textarea(id="input-smiles",
                                         placeholder="CC(=O)Oc1ccccc1C(=O)O\nC1=CC=CC=C1",
                                         style={"height": "80px"}, className="mb-3"),
                            dbc.Button("🚀 Iniciar Pipeline", id="btn-submit",
                                       color="info", className="w-100"),
                        ])
                    ], className="mb-3", style={"backgroundColor": "#111827"}),

                    # Status do Job
                    dbc.Card([
                        dbc.CardHeader("📊 Status", style={"color": "#66bb6a"}),
                        dbc.CardBody([
                            html.Div(id="job-status", children="Aguardando submissão..."),
                            dbc.Progress(id="job-progress", value=0,
                                         className="mt-2", color="info"),
                        ])
                    ], style={"backgroundColor": "#111827"}),

                ], width=3),

                # Coluna Central: Simulação de Fluido
                dbc.Col([
                    dbc.Card([
                        dbc.CardHeader("🌊 Simulação de Docking em Fluido",
                                       style={"color": "#4fc3f7"}),
                        dbc.CardBody([
                            html.Iframe(
                                id="fluid-iframe",
                                srcDoc=self._get_fluid_html(),
                                style={
                                    "width": "100%", "height": "620px",
                                    "border": "1px solid #1e293b",
                                    "borderRadius": "8px",
                                },
                            ),
                        ], style={"padding": "8px"}),
                    ], style={"backgroundColor": "#111827"}),
                ], width=6),

                # Coluna Direita: Resultados
                dbc.Col([
                    dbc.Card([
                        dbc.CardHeader("🎯 Ranqueamento", style={"color": "#ffa726"}),
                        dbc.CardBody([
                            dcc.Graph(id="graph-ranking",
                                      style={"height": "280px"}),
                        ])
                    ], className="mb-3", style={"backgroundColor": "#111827"}),

                    dbc.Card([
                        dbc.CardHeader("⚡ Afinidades (ΔG)", style={"color": "#ef5350"}),
                        dbc.CardBody([
                            dcc.Graph(id="graph-affinity",
                                      style={"height": "280px"}),
                        ])
                    ], style={"backgroundColor": "#111827"}),
                ], width=3),
            ]),

            # Polling interval
            dcc.Interval(id="interval-poll", interval=3000, disabled=True),
            dcc.Store(id="store-job-id"),

        ], fluid=True, style={
            "backgroundColor": "#0a0e1a", "minHeight": "100vh",
            "fontFamily": "'Inter', sans-serif",
        })

    def _register_callbacks(self) -> None:
        """Registra callbacks reativos do Dash."""
        app = self._app

        @app.callback(
            [Output("store-job-id", "data"),
             Output("interval-poll", "disabled"),
             Output("job-status", "children")],
            Input("btn-submit", "n_clicks"),
            [State("input-protein", "value"),
             State("input-name", "value"),
             State("input-pathogen", "value"),
             State("input-smiles", "value")],
            prevent_initial_call=True,
        )
        def submit_job(n_clicks, protein_id, name, pathogen, smiles_text):
            if not protein_id:
                return no_update, True, "⚠️ Insira um UniProt ID"

            smiles_list = [s.strip() for s in (smiles_text or "").split("\n") if s.strip()]
            ligand_ids = [f"lig_{i}" for i in range(len(smiles_list))]

            payload = {
                "protein_id": protein_id,
                "protein_name": name or protein_id,
                "pathogen_class": pathogen or "unknown",
                "ligand_smiles": smiles_list,
                "ligand_ids": ligand_ids,
            }

            try:
                resp = requests.post(
                    f"{self._api_base}/api/v1/discover",
                    json=payload, timeout=10,
                )
                data = resp.json()
                return data["job_id"], False, f"✅ Job {data['job_id'][:8]}... submetido"
            except Exception as e:
                return no_update, True, f"❌ Erro: {e}"

        @app.callback(
            [Output("job-status", "children", allow_duplicate=True),
             Output("job-progress", "value"),
             Output("graph-ranking", "figure"),
             Output("graph-affinity", "figure"),
             Output("interval-poll", "disabled", allow_duplicate=True)],
            Input("interval-poll", "n_intervals"),
            State("store-job-id", "data"),
            prevent_initial_call=True,
        )
        def poll_status(n, job_id):
            if not job_id:
                return no_update, 0, no_update, no_update, True

            try:
                resp = requests.get(
                    f"{self._api_base}/api/v1/jobs/{job_id}", timeout=5
                )
                data = resp.json()
            except Exception:
                return "⏳ Polling...", 0, no_update, no_update, False

            progress = int(data.get("progress", 0) * 100)
            status_text = f"📡 {data.get('current_step', '...')} ({progress}%)"

            ranking_fig = _empty_fig("Aguardando dados...")
            affinity_fig = _empty_fig("Aguardando docking...")
            stop_poll = False

            if data.get("status") == "completed":
                status_text = "✅ Pipeline concluído!"
                stop_poll = True

                results = data.get("results", {})
                scoring = results.get("scoring", {})

                if scoring:
                    ranking_fig = _build_ranking_chart(scoring)

                docking = results.get("docking", {})
                if docking and docking.get("docking_performed"):
                    affinity_fig = _build_affinity_chart(docking)

            elif data.get("status") == "failed":
                status_text = f"❌ Falhou: {data.get('error', 'desconhecido')}"
                stop_poll = True

            return status_text, progress, ranking_fig, affinity_fig, stop_poll

    def run(self, host: str = "0.0.0.0", port: int = 8050, debug: bool = True) -> None:
        """Inicia o dashboard."""
        app = self.create_app()
        app.run(host=host, port=port, debug=debug)


# ---------------------------------------------------------------------------
# Helpers de gráficos
# ---------------------------------------------------------------------------
def _empty_fig(msg: str) -> go.Figure:
    """Retorna figura vazia com mensagem centralizada."""
    fig = go.Figure()
    fig.add_annotation(text=msg, showarrow=False, font=dict(size=14, color="#666"))
    fig.update_layout(
        paper_bgcolor="#111827", plot_bgcolor="#111827",
        xaxis=dict(visible=False), yaxis=dict(visible=False),
        margin=dict(l=10, r=10, t=10, b=10),
    )
    return fig


def _build_ranking_chart(scoring: dict) -> go.Figure:
    """Constrói gráfico de barras do targetability score."""
    fig = go.Figure(go.Bar(
        x=[scoring.get("protein_id", "N/A")],
        y=[scoring.get("targetability_score_raw", 0)],
        marker_color="#4fc3f7",
        text=[f"{scoring.get('targetability_score_raw', 0):.3f}"],
        textposition="outside",
    ))
    fig.update_layout(
        paper_bgcolor="#111827", plot_bgcolor="#111827",
        font_color="#e0e0e0", yaxis_title="Targetability Score",
        margin=dict(l=40, r=10, t=10, b=30),
        yaxis=dict(range=[0, 1.5], gridcolor="#1e293b"),
    )
    return fig


def _build_affinity_chart(docking: dict) -> go.Figure:
    """Constrói scatter plot de afinidades ΔG."""
    results = docking.get("results", [])
    if not results:
        return _empty_fig("Sem resultados de docking")

    df = pd.DataFrame(results)
    fig = px.scatter(
        df, x="ligand_id", y="affinity_kcal",
        color="affinity_kcal", color_continuous_scale="RdYlGn_r",
    )
    fig.add_hline(y=-7.0, line_dash="dash", line_color="#76ff03",
                  annotation_text="Threshold: -7.0 kcal/mol")
    fig.update_layout(
        paper_bgcolor="#111827", plot_bgcolor="#111827",
        font_color="#e0e0e0", yaxis_title="ΔG (kcal/mol)",
        margin=dict(l=40, r=10, t=10, b=30),
    )
    return fig


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    viz = Visualizer()
    viz.run()
