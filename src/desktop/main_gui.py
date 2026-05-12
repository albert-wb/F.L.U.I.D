# =============================================================================
# F.L.U.I.D — Cliente Desktop (CustomTkinter)
# =============================================================================
"""
GUI Desktop com:
  - Sidebar: inputs PDB/patógeno/SMILES + botão de execução
  - Workspace: barra de progresso, tabela de resultados, viewer de fluidos
  - Threading: chamadas de rede NUNCA bloqueiam o mainloop()
  - Polling: .after(2000ms) consulta status até completed/failed
  - Webview: motor matter.js+p5.js embutido via pywebview
"""
from __future__ import annotations

import json
import threading
import webbrowser
from pathlib import Path
from tkinter import filedialog, ttk
from typing import Any

import customtkinter as ctk

from src.desktop.api_client import FluidAPIClient, JobStatus
from src.desktop.terminal_gui import LogTerminalFrame

# --- Tema ---
ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("dark-blue")

# Cores do design system
C = {
    "bg": "#0a0e1a", "surface": "#111827", "border": "#1e293b",
    "accent": "#4fc3f7", "green": "#66bb6a", "red": "#ef5350",
    "orange": "#ffa726", "match": "#76ff03", "text": "#e0e0e0",
    "muted": "#9e9e9e",
}

PATHOGEN_OPTIONS = [
    "rna_virus", "dna_virus", "bacterium", "protozoan", "fungus",
]
PATHOGEN_LABELS = {
    "rna_virus": "Vírus RNA", "dna_virus": "Vírus DNA",
    "bacterium": "Bactéria", "protozoan": "Protozoário", "fungus": "Fungo",
}


# ============================================================================
# SIDEBAR — Painel de Controle (25%)
# ============================================================================
class SidebarFrame(ctk.CTkFrame):
    """Painel de entrada com campos PDB, patógeno e SMILES."""

    def __init__(self, master: "App", **kw: Any) -> None:
        super().__init__(master, corner_radius=0, fg_color=C["surface"], **kw)
        self.app = master
        self._pdb_path: Path | None = None
        self._build()

    def _build(self) -> None:
        # Logo
        ctk.CTkLabel(
            self, text="F.L.U.I.D", font=("Consolas", 28, "bold"),
            text_color=C["accent"],
        ).pack(pady=(20, 0))
        ctk.CTkLabel(
            self, text="Drug Discovery Platform",
            font=("Segoe UI", 10), text_color=C["muted"],
        ).pack(pady=(0, 20))

        sep = ctk.CTkFrame(self, height=1, fg_color=C["border"])
        sep.pack(fill="x", padx=15, pady=5)

        # --- PDB Input ---
        ctk.CTkLabel(self, text="🧬  Receptor PDB", anchor="w",
                      font=("Segoe UI", 12, "bold")).pack(padx=20, anchor="w", pady=(15, 2))

        self.pdb_entry = ctk.CTkEntry(self, placeholder_text="PDB ID (ex: 6LU7)")
        self.pdb_entry.pack(padx=20, fill="x", pady=2)

        ctk.CTkButton(
            self, text="📂 Upload PDB Local", height=28,
            fg_color="transparent", border_width=1, border_color=C["border"],
            command=self._browse_pdb,
        ).pack(padx=20, fill="x", pady=(4, 2))

        self.file_label = ctk.CTkLabel(
            self, text="", font=("Segoe UI", 9), text_color=C["muted"],
        )
        self.file_label.pack(padx=20, anchor="w")

        # --- Patógeno ---
        ctk.CTkLabel(self, text="🦠  Classe do Patógeno", anchor="w",
                      font=("Segoe UI", 12, "bold")).pack(padx=20, anchor="w", pady=(15, 2))

        self.pathogen_var = ctk.StringVar(value="rna_virus")
        self.pathogen_menu = ctk.CTkOptionMenu(
            self, values=list(PATHOGEN_LABELS.values()),
            variable=self.pathogen_var, command=self._on_pathogen_change,
        )
        self.pathogen_menu.pack(padx=20, fill="x", pady=2)
        self.pathogen_var.set("Vírus RNA")

        # --- Proteína ---
        ctk.CTkLabel(self, text="🏷️  Nome da Proteína", anchor="w",
                      font=("Segoe UI", 12, "bold")).pack(padx=20, anchor="w", pady=(15, 2))
        self.name_entry = ctk.CTkEntry(self, placeholder_text="ex: 3CLpro (opcional)")
        self.name_entry.pack(padx=20, fill="x", pady=2)

        # --- SMILES ---
        ctk.CTkLabel(self, text="💊  Ligantes (SMILES)", anchor="w",
                      font=("Segoe UI", 12, "bold")).pack(padx=20, anchor="w", pady=(15, 2))

        self.smiles_text = ctk.CTkTextbox(self, height=80, font=("Consolas", 11))
        self.smiles_text.pack(padx=20, fill="x", pady=2)
        self.smiles_text.insert("0.0", "CC(=O)Oc1ccccc1C(=O)O")

        # --- Botão de Execução ---
        self.run_btn = ctk.CTkButton(
            self, text="▶  Run F.L.U.I.D Pipeline", height=42,
            font=("Segoe UI", 14, "bold"), fg_color="#0288d1",
            hover_color="#0277bd", command=self.app.on_run_pipeline,
        )
        self.run_btn.pack(padx=20, fill="x", pady=(25, 5))

        # --- Status da API ---
        self.api_status = ctk.CTkLabel(
            self, text="● API: verificando...",
            font=("Segoe UI", 10), text_color=C["muted"],
        )
        self.api_status.pack(padx=20, anchor="w", pady=(15, 5), side="bottom")

    def _browse_pdb(self) -> None:
        path = filedialog.askopenfilename(
            title="Selecionar arquivo PDB",
            filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")],
        )
        if path:
            self._pdb_path = Path(path)
            self.file_label.configure(text=f"📄 {self._pdb_path.name}")
            self.pdb_entry.delete(0, "end")
            self.pdb_entry.insert(0, self._pdb_path.stem)

    def _on_pathogen_change(self, value: str) -> None:
        pass  # Variável já atualizada via StringVar

    def get_inputs(self) -> dict[str, Any]:
        """Coleta todos os inputs do formulário."""
        label = self.pathogen_var.get()
        pc = next((k for k, v in PATHOGEN_LABELS.items() if v == label), "unknown")
        smiles_raw = self.smiles_text.get("0.0", "end").strip()
        smiles_list = [s.strip() for s in smiles_raw.split("\n") if s.strip()]
        return {
            "pdb_id": self.pdb_entry.get().strip(),
            "pdb_path": self._pdb_path,
            "protein_name": self.name_entry.get().strip(),
            "pathogen_class": pc,
            "ligand_smiles": smiles_list,
        }

    def set_api_status(self, online: bool) -> None:
        if online:
            self.api_status.configure(text="● API: online", text_color=C["green"])
        else:
            self.api_status.configure(text="● API: offline", text_color=C["red"])


# ============================================================================
# RESULT TABLE — Tabela de Ligantes Ranqueados
# ============================================================================
class ResultTableFrame(ctk.CTkFrame):
    """Tabela Treeview estilizada com resultados de docking."""

    COLUMNS = ("ligand_id", "score", "affinity", "status")
    HEADERS = ("Ligand ID", "Targetability Score", "ΔG (kcal/mol)", "Status")
    WIDTHS = (140, 150, 130, 100)

    def __init__(self, master: ctk.CTkFrame, **kw: Any) -> None:
        super().__init__(master, fg_color=C["surface"], corner_radius=8, **kw)
        self._build()

    def _build(self) -> None:
        ctk.CTkLabel(
            self, text="📊  Resultados de Docking", anchor="w",
            font=("Segoe UI", 13, "bold"), text_color=C["accent"],
        ).pack(padx=12, pady=(10, 5), anchor="w")

        # Estilo dark para Treeview
        style = ttk.Style()
        style.theme_use("default")
        style.configure("Dark.Treeview", background="#1a1f2e",
                         foreground="#e0e0e0", fieldbackground="#1a1f2e",
                         rowheight=28, font=("Consolas", 10))
        style.configure("Dark.Treeview.Heading", background="#111827",
                         foreground="#4fc3f7", font=("Segoe UI", 10, "bold"))
        style.map("Dark.Treeview", background=[("selected", "#0288d1")])

        container = ctk.CTkFrame(self, fg_color="transparent")
        container.pack(fill="both", expand=True, padx=10, pady=(0, 10))

        self.tree = ttk.Treeview(
            container, columns=self.COLUMNS, show="headings",
            style="Dark.Treeview", height=8,
        )
        for col, hdr, w in zip(self.COLUMNS, self.HEADERS, self.WIDTHS):
            self.tree.heading(col, text=hdr)
            self.tree.column(col, width=w, anchor="center")

        scrollbar = ttk.Scrollbar(container, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=scrollbar.set)
        self.tree.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

    def clear(self) -> None:
        for item in self.tree.get_children():
            self.tree.delete(item)

    def populate(self, results: dict[str, Any]) -> None:
        """Preenche tabela com resultados do pipeline."""
        self.clear()
        scoring = results.get("scoring", {})
        docking = results.get("docking", {})
        score_raw = scoring.get("targetability_score_raw", 0.0)

        docking_list = docking.get("results", [])
        if docking_list:
            for r in docking_list:
                aff = r.get("affinity_kcal", 0.0)
                tag = "match" if aff < -7.0 else ""
                self.tree.insert("", "end", values=(
                    r.get("ligand_id", "?"),
                    f"{score_raw:.4f}",
                    f"{aff:.1f}",
                    "✅ Match" if aff < -7.0 else "—",
                ), tags=(tag,))
            self.tree.tag_configure("match", foreground=C["match"])
        else:
            self.tree.insert("", "end", values=(
                scoring.get("protein_id", "?"),
                f"{score_raw:.4f}", "N/A",
                "Scoring only",
            ))


# ============================================================================
# FLUID VIEWER — Motor de Física 2D (pywebview bridge)
# ============================================================================
class FluidViewerFrame(ctk.CTkFrame):
    """
    Wrapper para o motor matter.js+p5.js.

    Estratégia de lançamento:
        pywebview exige rodar na main thread no Windows, o que é
        incompatível com o mainloop do CustomTkinter. Portanto:

        1. Salva o HTML com a engine de fluidos em arquivo local
        2. Abre no navegador padrão do sistema
        3. Injeta resultados via arquivo JSON que o JS lê via polling

        Isso garante que a GUI NUNCA congela.
    """

    def __init__(self, master: ctk.CTkFrame, **kw: Any) -> None:
        super().__init__(master, fg_color=C["surface"], corner_radius=8, **kw)
        self._html_path: Path | None = None
        self._data_path: Path | None = None
        self._build()

    def _build(self) -> None:
        ctk.CTkLabel(
            self, text="🌊  Simulação de Docking em Fluido", anchor="w",
            font=("Segoe UI", 13, "bold"), text_color=C["accent"],
        ).pack(padx=12, pady=(10, 5), anchor="w")

        self.launch_btn = ctk.CTkButton(
            self, text="🚀 Abrir Simulação de Fluidos", height=36,
            command=self._launch_viewer,
        )
        self.launch_btn.pack(padx=12, pady=10, fill="x")

        self.info_label = ctk.CTkLabel(
            self, text="Abre no navegador padrão com matter.js + p5.js",
            font=("Segoe UI", 10), text_color=C["muted"],
        )
        self.info_label.pack(padx=12, pady=(0, 10))

    def _get_fluid_html(self) -> str:
        """Lê o JS e gera HTML com polling de dados via localStorage."""
        js_path = Path(__file__).parent.parent / "dashboard" / "assets" / "fluid_engine.js"
        try:
            js_code = js_path.read_text(encoding="utf-8")
        except FileNotFoundError:
            js_code = "// fluid_engine.js não encontrado"

        # JS extra: polling de localStorage para receber dados do Python
        bridge_js = """
// Bridge: Python salva JSON em localStorage, JS consome
setInterval(function() {
    try {
        var data = localStorage.getItem('fluid_docking_data');
        if (data) {
            var parsed = JSON.parse(data);
            if (parsed._consumed) return;
            // Adiciona proteína
            if (parsed.protein_id) {
                window.fluidAPI.addProtein(
                    parsed.protein_id, parsed.protein_id, 450, 300, 55
                );
            }
            // Adiciona ligantes
            if (parsed.ligands) {
                parsed.ligands.forEach(function(l) {
                    window.fluidAPI.addLigand(l.id, l.id, null);
                });
            }
            // Submete resultados de docking
            if (parsed.docking_results) {
                window.fluidAPI.submitResults(parsed.docking_results);
            }
            parsed._consumed = true;
            localStorage.setItem('fluid_docking_data', JSON.stringify(parsed));
        }
    } catch(e) { console.error('Bridge error:', e); }
}, 1500);
"""

        return f"""<!DOCTYPE html>
<html><head><meta charset="UTF-8"><title>F.L.U.I.D — Simulação</title>
<style>*{{margin:0;padding:0}}body{{background:#0a0e1a;overflow:hidden}}</style>
<script src="https://cdnjs.cloudflare.com/ajax/libs/matter-js/0.19.0/matter.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/p5.js/1.9.4/p5.min.js"></script>
</head><body><div id="fluid-container"></div>
<script>{js_code}</script>
<script>{bridge_js}</script>
</body></html>"""

    def _launch_viewer(self) -> None:
        """Salva HTML e abre no navegador. Nunca bloqueia o mainloop."""
        viewer_dir = Path("workspace")
        viewer_dir.mkdir(exist_ok=True)

        self._html_path = viewer_dir / "fluid_viewer.html"
        self._data_path = viewer_dir / "fluid_data.json"

        self._html_path.write_text(self._get_fluid_html(), encoding="utf-8")
        webbrowser.open(self._html_path.resolve().as_uri())

        self.info_label.configure(
            text="✅ Simulação aberta no navegador",
            text_color=C["green"],
        )

    def inject_results(self, results: dict[str, Any]) -> None:
        """Injeta resultados para a simulação via arquivo JSON."""
        scoring = results.get("scoring", {})
        docking = results.get("docking", {})
        dock_list = docking.get("results", [])

        data = {
            "protein_id": scoring.get("protein_id", "target"),
            "ligands": [{"id": r.get("ligand_id", f"lig_{i}")}
                        for i, r in enumerate(dock_list)],
            "docking_results": dock_list,
            "_consumed": False,
        }

        # Salva JSON para o bridge JS
        data_path = Path("workspace/fluid_data.json")
        data_path.write_text(json.dumps(data), encoding="utf-8")

        # Também abre automaticamente se viewer não estava aberto
        if not self._html_path or not self._html_path.exists():
            self._launch_viewer()


# ============================================================================
# APP PRINCIPAL
# ============================================================================
class App(ctk.CTk):
    """
    Janela principal do F.L.U.I.D Desktop.

    Threading:
        - Chamadas de rede: threading.Thread (nunca no mainloop)
        - Polling: self.after(2000) — timer do Tk, thread-safe
        - Atualizações de UI: sempre no mainloop via after()
    """

    def __init__(self) -> None:
        super().__init__()
        self.title("F.L.U.I.D — Drug Discovery Platform")
        self.geometry("1280x780")
        self.minsize(1000, 600)
        self.configure(fg_color=C["bg"])

        self._client = FluidAPIClient("http://localhost:8000")
        self._current_job_id: str | None = None
        self._polling = False

        self._build_layout()
        self._check_api_health()

    def _build_layout(self) -> None:
        """Grid: sidebar 25% | workspace 75%."""
        self.grid_columnconfigure(0, weight=1)   # sidebar
        self.grid_columnconfigure(1, weight=3)   # workspace
        self.grid_rowconfigure(0, weight=1)

        # --- Sidebar ---
        self.sidebar = SidebarFrame(self)
        self.sidebar.grid(row=0, column=0, sticky="nsew")

        # --- Workspace ---
        workspace = ctk.CTkFrame(self, fg_color=C["bg"], corner_radius=0)
        workspace.grid(row=0, column=1, sticky="nsew", padx=(2, 0))
        workspace.grid_rowconfigure(1, weight=1)
        workspace.grid_columnconfigure(0, weight=1)

        # Top Bar: Status + Progress
        top = ctk.CTkFrame(workspace, fg_color=C["surface"], height=70, corner_radius=8)
        top.grid(row=0, column=0, sticky="ew", padx=10, pady=(10, 5))
        top.grid_columnconfigure(0, weight=1)

        self.status_label = ctk.CTkLabel(
            top, text="⏳ Aguardando submissão...",
            font=("Segoe UI", 13), text_color=C["text"], anchor="w",
        )
        self.status_label.grid(row=0, column=0, padx=15, pady=(10, 2), sticky="w")

        self.progress_bar = ctk.CTkProgressBar(top, progress_color=C["accent"])
        self.progress_bar.grid(row=1, column=0, padx=15, pady=(2, 10), sticky="ew")
        self.progress_bar.set(0)

        # --- Tabview: Pipeline | Terminal ---
        self.tabview = ctk.CTkTabview(
            workspace, fg_color=C["bg"],
            segmented_button_fg_color=C["surface"],
            segmented_button_selected_color="#0288d1",
            segmented_button_unselected_color=C["surface"],
        )
        self.tabview.grid(row=1, column=0, sticky="nsew", padx=10, pady=(5, 10))

        # Aba 1: Pipeline (tabela + fluidos)
        tab_pipeline = self.tabview.add("📊 Pipeline")
        tab_pipeline.grid_rowconfigure(0, weight=1)
        tab_pipeline.grid_rowconfigure(1, weight=1)
        tab_pipeline.grid_columnconfigure(0, weight=1)

        self.result_table = ResultTableFrame(tab_pipeline)
        self.result_table.grid(row=0, column=0, sticky="nsew", pady=(0, 5))

        self.fluid_viewer = FluidViewerFrame(tab_pipeline)
        self.fluid_viewer.grid(row=1, column=0, sticky="nsew", pady=(5, 0))

        # Aba 2: Terminal de observabilidade
        tab_terminal = self.tabview.add("📟 Terminal")
        tab_terminal.grid_rowconfigure(0, weight=1)
        tab_terminal.grid_columnconfigure(0, weight=1)

        self.log_terminal = LogTerminalFrame(
            tab_terminal, show_start_button=False,
        )
        self.log_terminal.grid(row=0, column=0, sticky="nsew")

    # ----------------------------------------------------------------
    # API Health Check (thread)
    # ----------------------------------------------------------------
    def _check_api_health(self) -> None:
        def _check() -> None:
            ok = self._client.health_check()
            self.after(0, lambda: self.sidebar.set_api_status(ok))
        threading.Thread(target=_check, daemon=True).start()
        self.after(15000, self._check_api_health)  # re-check every 15s

    # ----------------------------------------------------------------
    # Pipeline Execution (thread + polling)
    # ----------------------------------------------------------------
    def on_run_pipeline(self) -> None:
        """Chamado pelo botão Run — dispara em thread secundária."""
        inputs = self.sidebar.get_inputs()
        pdb_id = inputs["pdb_id"]

        if not pdb_id:
            self.status_label.configure(text="⚠️ Insira um PDB ID", text_color=C["orange"])
            return

        self.sidebar.run_btn.configure(state="disabled", text="⏳ Processando...")
        self.status_label.configure(text="📡 Submetendo job...", text_color=C["text"])
        self.progress_bar.set(0)
        self.result_table.clear()

        # Log no terminal
        tlog = self.log_terminal.get_logger()
        tlog.info("=" * 48)
        tlog.info(f"JOB SUBMETIDO: PDB={inputs['pdb_id']} classe={inputs['pathogen_class']}")
        tlog.info(f"  Ligantes: {len(inputs.get('ligand_smiles', []))}")
        tlog.info("=" * 48)

        threading.Thread(target=self._submit_job, args=(inputs,), daemon=True).start()

    def _submit_job(self, inputs: dict[str, Any]) -> None:
        """Executa na thread secundária — NUNCA no mainloop."""
        pdb_path = inputs.get("pdb_path")
        smiles_csv = ",".join(inputs.get("ligand_smiles", []))

        if pdb_path and pdb_path.exists():
            result = self._client.upload_receptor(
                pdb_path, inputs["pdb_id"], inputs["protein_name"],
                inputs["pathogen_class"], smiles_csv,
            )
        else:
            result = self._client.fetch_receptor(
                inputs["pdb_id"], inputs["protein_name"],
                inputs["pathogen_class"], smiles_csv,
            )

        # Atualiza UI via after() — thread-safe
        if result.success:
            self._current_job_id = result.job_id
            self.after(0, lambda: self._on_submitted(result.job_id))
            self.after(2000, self._poll_job)
        else:
            self.after(0, lambda: self._on_error(result.error))

    def _on_submitted(self, job_id: str) -> None:
        self.status_label.configure(
            text=f"✅ Job submetido: {job_id[:12]}...", text_color=C["green"],
        )
        self.progress_bar.set(0.05)
        self._polling = True
        self.log_terminal.get_logger().info(f"Job aceito pelo servidor: {job_id}")

    def _on_error(self, error: str) -> None:
        self.status_label.configure(text=f"❌ {error}", text_color=C["red"])
        self.sidebar.run_btn.configure(state="normal", text="▶  Run F.L.U.I.D Pipeline")
        self.log_terminal.get_logger().error(f"ERRO: {error}")

    # ----------------------------------------------------------------
    # Polling — .after(2000) — thread-safe, roda no mainloop
    # ----------------------------------------------------------------
    def _poll_job(self) -> None:
        """Consulta status a cada 2s via thread + after()."""
        if not self._current_job_id or not self._polling:
            return

        def _fetch() -> None:
            status = self._client.poll_job(self._current_job_id)
            self.after(0, lambda: self._handle_poll_result(status))

        threading.Thread(target=_fetch, daemon=True).start()

    def _handle_poll_result(self, status: JobStatus) -> None:
        """Processa resultado do polling — roda no mainloop via after()."""
        self.progress_bar.set(status.progress)

        step = status.current_step or status.status
        self.status_label.configure(
            text=f"📡 {step} ({int(status.progress * 100)}%)",
            text_color=C["text"],
        )
        self.log_terminal.get_logger().info(
            f"Polling: {step} | progress={status.progress:.0%}"
        )

        if status.status == "completed":
            self._polling = False
            self.status_label.configure(
                text="✅ Pipeline concluído!", text_color=C["green"],
            )
            self.progress_bar.set(1.0)
            self.sidebar.run_btn.configure(state="normal", text="▶  Run F.L.U.I.D Pipeline")

            if status.results:
                self.result_table.populate(status.results)
                self.fluid_viewer.inject_results(status.results)

        elif status.status == "failed":
            self._polling = False
            self._on_error(status.error or "Pipeline falhou")

        else:
            # Continua polling
            self.after(2000, self._poll_job)


# ============================================================================
# ENTRY POINT
# ============================================================================
def main() -> None:
    app = App()
    app.mainloop()


if __name__ == "__main__":
    main()
