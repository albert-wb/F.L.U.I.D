# =============================================================================
# F.L.U.I.D — Módulo de Observabilidade (Log Terminal)
# =============================================================================
"""
Terminal de log embutido no CustomTkinter para monitoramento do pipeline.

Arquitetura:
    ┌─────────────────────────────────────────────────┐
    │  Top Bar   [ ▶ Start Pipeline ] [ 📄 Export ]   │
    ├─────────────────────────────────────────────────┤
    │                                                 │
    │   CTkTextbox (Consolas, read-only)               │
    │   └── Buffer circular: max 5000 linhas O(1)     │
    │   └── Auto-scroll: sempre no final              │
    │                                                 │
    └─────────────────────────────────────────────────┘

Concorrência:
    - Pipeline roda em threading.Thread(daemon=True)
    - TkinterLogHandler usa widget.after(0, fn) para thread-safety
    - Nenhuma escrita no CTkTextbox ocorre fora do mainloop

Complexidade:
    - Inserção de log: O(1) amortizado
    - Trim de buffer: O(k) onde k = linhas excedentes (executado raramente)
    - Exportação: O(n) onde n = linhas no buffer
    - Memória: O(1) delimitada pelo MAX_BUFFER_LINES
"""
from __future__ import annotations

import logging
import threading
import time
from datetime import datetime, timezone
from pathlib import Path
from tkinter import filedialog
from typing import Any, Callable

import customtkinter as ctk

# ---------------------------------------------------------------------------
# Constantes do design system (importadas do main_gui se disponível)
# ---------------------------------------------------------------------------
try:
    from src.desktop.main_gui import C
except ImportError:
    C = {
        "bg": "#0a0e1a", "surface": "#111827", "border": "#1e293b",
        "accent": "#4fc3f7", "green": "#66bb6a", "red": "#ef5350",
        "orange": "#ffa726", "match": "#76ff03", "text": "#e0e0e0",
        "muted": "#9e9e9e",
    }

# Buffer máximo: 5000 linhas → memória constante O(1)
MAX_BUFFER_LINES = 5000

# Mapa de cores por nível de log
LOG_LEVEL_COLORS = {
    "DEBUG":    "#616161",   # cinza
    "INFO":     "#4fc3f7",   # azul claro (accent)
    "WARNING":  "#ffa726",   # laranja
    "ERROR":    "#ef5350",   # vermelho
    "CRITICAL": "#ff1744",   # vermelho intenso
    "SUCCESS":  "#66bb6a",   # verde
}

# Tag de cor para prefixo de timestamp
TIMESTAMP_COLOR = "#757575"


# ============================================================================
# LOGGING HANDLER — Thread-safe via widget.after()
# ============================================================================
class TkinterLogHandler(logging.Handler):
    """
    Handler que redireciona logs do módulo `logging` do Python para
    um CTkTextbox de forma thread-safe.

    O handler NUNCA escreve diretamente no widget. Toda operação de
    UI é agendada via `widget.after(0, callback)`, garantindo que a
    renderização ocorra exclusivamente na Main Thread do Tk.

    Uso:
        handler = TkinterLogHandler(terminal_frame)
        logging.getLogger().addHandler(handler)
    """

    def __init__(self, terminal: "LogTerminalFrame") -> None:
        super().__init__()
        self._terminal = terminal
        self.setFormatter(logging.Formatter("%(message)s"))

    def emit(self, record: logging.LogRecord) -> None:
        """Intercepta log e agenda inserção no mainloop."""
        try:
            msg = self.format(record)
            level = record.levelname
            timestamp = datetime.now().strftime("%H:%M:%S.%f")[:-3]
            formatted = f"[{timestamp}] [{level:<8}] {msg}"

            # Agenda no mainloop — NUNCA escreve direto
            self._terminal.textbox.after(
                0, lambda: self._terminal.append_log(formatted, level)
            )
        except Exception:
            self.handleError(record)


# ============================================================================
# LOG TERMINAL FRAME — Widget de observabilidade
# ============================================================================
class LogTerminalFrame(ctk.CTkFrame):
    """
    Terminal de log embutido com controles de execução e exportação.

    Características:
        - Font monoespaçada (Consolas) para estética de terminal
        - CTkTextbox read-only (state=disabled)
        - Buffer circular: max 5000 linhas, trim automático O(k)
        - Auto-scroll: rolagem automática para o final
        - Export to Markdown com metadados temporais

    Args:
        master: widget pai (CTkFrame, CTk, etc.)
        pipeline_callback: função a executar no botão "Start Pipeline".
            Se None, usa mock interno para demonstração.
        show_start_button: se True, exibe o botão "Start Pipeline".
    """

    def __init__(
        self,
        master: Any,
        pipeline_callback: Callable[[], None] | None = None,
        show_start_button: bool = True,
        **kwargs: Any,
    ) -> None:
        super().__init__(master, fg_color=C["surface"], corner_radius=8, **kwargs)

        self._pipeline_callback = pipeline_callback
        self._show_start = show_start_button
        self._line_count = 0
        self._is_running = False

        # Logger dedicado para este terminal
        self._logger = logging.getLogger("fluid.terminal")
        self._logger.setLevel(logging.DEBUG)
        self._logger.propagate = False

        self._build()
        self._setup_handler()
        self._write_banner()

    # ----------------------------------------------------------------
    # Construção da interface
    # ----------------------------------------------------------------
    def _build(self) -> None:
        """Monta top bar + área do terminal."""
        self.grid_rowconfigure(1, weight=1)
        self.grid_columnconfigure(0, weight=1)

        # --- Top Bar ---
        top_bar = ctk.CTkFrame(self, fg_color="transparent", height=40)
        top_bar.grid(row=0, column=0, sticky="ew", padx=8, pady=(8, 4))
        top_bar.grid_columnconfigure(1, weight=1)

        # Título
        ctk.CTkLabel(
            top_bar, text="📟  Terminal de Observabilidade",
            font=("Segoe UI", 13, "bold"), text_color=C["accent"],
        ).grid(row=0, column=0, sticky="w", padx=(4, 12))

        # Botão Start Pipeline (opcional)
        if self._show_start:
            self.start_btn = ctk.CTkButton(
                top_bar, text="▶  Start Pipeline", width=160, height=32,
                font=("Segoe UI", 12, "bold"),
                fg_color="#2e7d32", hover_color="#1b5e20",
                command=self._on_start_pipeline,
            )
            self.start_btn.grid(row=0, column=2, padx=(4, 4))

        # Botão Export
        self.export_btn = ctk.CTkButton(
            top_bar, text="📄 Export to Markdown", width=180, height=32,
            font=("Segoe UI", 11),
            fg_color="transparent", border_width=1, border_color=C["border"],
            hover_color="#1e293b",
            command=self._on_export,
        )
        self.export_btn.grid(row=0, column=3, padx=(4, 4))

        # Botão Clear
        self.clear_btn = ctk.CTkButton(
            top_bar, text="🗑️ Clear", width=80, height=32,
            font=("Segoe UI", 11),
            fg_color="transparent", border_width=1, border_color=C["border"],
            hover_color="#1e293b",
            command=self._on_clear,
        )
        self.clear_btn.grid(row=0, column=4, padx=(4, 4))

        # Contador de linhas
        self._line_label = ctk.CTkLabel(
            top_bar, text="0 linhas",
            font=("Consolas", 10), text_color=C["muted"],
        )
        self._line_label.grid(row=0, column=5, padx=(8, 4))

        # Separador
        sep = ctk.CTkFrame(self, height=1, fg_color=C["border"])
        sep.grid(row=0, column=0, sticky="sew", padx=8)

        # --- Área do Terminal ---
        self.textbox = ctk.CTkTextbox(
            self,
            font=("Consolas", 11),
            fg_color="#0d1117",
            text_color=C["text"],
            corner_radius=6,
            state="disabled",
            wrap="none",
            activate_scrollbars=True,
        )
        self.textbox.grid(row=1, column=0, sticky="nsew", padx=8, pady=(4, 8))

        # Configura tags de cor para cada nível de log
        for level, color in LOG_LEVEL_COLORS.items():
            self.textbox._textbox.tag_configure(f"level_{level}", foreground=color)
        self.textbox._textbox.tag_configure("timestamp", foreground=TIMESTAMP_COLOR)
        self.textbox._textbox.tag_configure("banner", foreground=C["accent"])

    def _setup_handler(self) -> None:
        """Registra o TkinterLogHandler no logger dedicado."""
        # Remove handlers anteriores para evitar duplicatas
        for h in self._logger.handlers[:]:
            self._logger.removeHandler(h)
        handler = TkinterLogHandler(self)
        handler.setLevel(logging.DEBUG)
        self._logger.addHandler(handler)

    def _write_banner(self) -> None:
        """Escreve o banner inicial no terminal."""
        banner = (
            "╔══════════════════════════════════════════════════════╗\n"
            "║  F.L.U.I.D — Pipeline Observation Terminal          ║\n"
            "║  Framework for Ligand-target Unified In-silico Discovery ║\n"
            "╚══════════════════════════════════════════════════════╝\n"
            ""
        )
        self.textbox.configure(state="normal")
        self.textbox.insert("end", banner)
        self.textbox.configure(state="disabled")
        self._line_count = banner.count("\n")
        self._update_line_label()

    # ----------------------------------------------------------------
    # API Pública: inserção de logs
    # ----------------------------------------------------------------
    def append_log(self, text: str, level: str = "INFO") -> None:
        """
        Insere texto no terminal. DEVE ser chamado na Main Thread.

        Para uso de background threads, utilize o logger:
            self._logger.info("mensagem")

        Args:
            text: linha de texto para inserir.
            level: nível do log (DEBUG, INFO, WARNING, ERROR, CRITICAL).
        """
        self.textbox.configure(state="normal")
        self.textbox.insert("end", text + "\n")
        self._line_count += 1

        # --- Trim de buffer: mantém memória O(1) ---
        if self._line_count > MAX_BUFFER_LINES:
            excess = self._line_count - MAX_BUFFER_LINES
            self.textbox.delete("1.0", f"{excess + 1}.0")
            self._line_count = MAX_BUFFER_LINES

        # Auto-scroll para o final
        self.textbox.see("end")
        self.textbox.configure(state="disabled")
        self._update_line_label()

    def get_logger(self) -> logging.Logger:
        """Retorna o logger dedicado deste terminal."""
        return self._logger

    # ----------------------------------------------------------------
    # Ações dos botões
    # ----------------------------------------------------------------
    def _on_start_pipeline(self) -> None:
        """Dispara pipeline em thread secundária."""
        if self._is_running:
            self._logger.warning("Pipeline já está em execução!")
            return

        self._is_running = True
        if self._show_start:
            self.start_btn.configure(
                state="disabled", text="⏳ Executando...",
                fg_color="#37474f",
            )

        callback = self._pipeline_callback or self._mock_pipeline
        thread = threading.Thread(target=self._run_pipeline, args=(callback,), daemon=True)
        thread.start()

    def _run_pipeline(self, callback: Callable[[], None]) -> None:
        """Wrapper que executa o callback e restaura o estado ao fim."""
        try:
            callback()
        except Exception as e:
            self._logger.error(f"Erro fatal no pipeline: {e}")
        finally:
            # Restaura botão via mainloop
            self.textbox.after(0, self._on_pipeline_finished)

    def _on_pipeline_finished(self) -> None:
        """Restaura o estado do botão após fim do pipeline."""
        self._is_running = False
        if self._show_start:
            self.start_btn.configure(
                state="normal", text="▶  Start Pipeline",
                fg_color="#2e7d32",
            )
        self._logger.info("Pipeline finalizado. Pronto para nova execução.")

    def _on_export(self) -> None:
        """Exporta conteúdo do terminal para Markdown."""
        # Leitura O(n) do buffer
        self.textbox.configure(state="normal")
        content = self.textbox.get("1.0", "end-1c")
        self.textbox.configure(state="disabled")

        if not content.strip():
            self._logger.warning("Terminal vazio — nada para exportar.")
            return

        path = filedialog.asksaveasfilename(
            title="Salvar relatório de log",
            defaultextension=".md",
            filetypes=[("Markdown", "*.md"), ("Texto", "*.txt")],
            initialfile=f"fluid_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.md",
        )

        if not path:
            return

        now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
        md_content = (
            f"# F.L.U.I.D — Relatório de Execução\n\n"
            f"**Data de exportação:** {now}  \n"
            f"**Linhas no buffer:** {self._line_count}  \n"
            f"**Buffer máximo:** {MAX_BUFFER_LINES}  \n\n"
            f"---\n\n"
            f"## Log de Execução\n\n"
            f"```log\n{content}\n```\n"
        )

        try:
            Path(path).write_text(md_content, encoding="utf-8")
            self._logger.info(f"Relatório exportado: {path}")
        except OSError as e:
            self._logger.error(f"Falha ao salvar relatório: {e}")

    def _on_clear(self) -> None:
        """Limpa o terminal e reseta o contador."""
        self.textbox.configure(state="normal")
        self.textbox.delete("1.0", "end")
        self.textbox.configure(state="disabled")
        self._line_count = 0
        self._update_line_label()
        self._write_banner()

    def _update_line_label(self) -> None:
        """Atualiza o contador de linhas na top bar."""
        self._line_label.configure(text=f"{self._line_count} linhas")

    # ----------------------------------------------------------------
    # Mock pipeline para demonstração standalone
    # ----------------------------------------------------------------
    def _mock_pipeline(self) -> None:
        """
        Pipeline simulado que emite logs assíncronos para demonstração.

        Reproduz as etapas reais do F.L.U.I.D com delays para simular
        processamento computacionalmente intensivo.
        """
        log = self._logger

        log.info("=" * 56)
        log.info("PIPELINE INICIADO — Modo demonstração")
        log.info("=" * 56)
        time.sleep(0.5)

        # --- Etapa 1: DataIngestor ---
        log.info("")
        log.info("[ETAPA 1/5] DataIngestor — Ingestão de sequências")
        log.info("-" * 48)
        time.sleep(0.3)
        log.info("  Carregando 1.247 sequências codificantes...")
        time.sleep(0.8)
        log.info("  Extraindo k-mers (k=7) para pré-filtro Jaccard...")
        time.sleep(0.6)
        log.info("  Pré-filtro: 1.247 → 842 candidatos (32.5% descartados)")
        time.sleep(0.4)
        log.info("  Greedy Centroid Clustering (threshold=0.85)...")
        time.sleep(1.0)
        log.info("  Clustering: 842 sequências → 23 centroides representativos")
        log.info("  Complexidade efetiva: O(N × L) = O(1247 × 892)")
        log.info("  ✓ DataIngestor concluído em 3.1s")

        # --- Etapa 2: RAGEngine ---
        log.info("")
        log.info("[ETAPA 2/5] RAGEngine — Busca semântica na literatura")
        log.info("-" * 48)
        time.sleep(0.3)
        log.info("  Query: 'catalytic residues SARS-CoV-2 3CLpro active site'")
        time.sleep(0.5)
        log.info("  Europe PMC: 47 artigos recuperados (Open Access)")
        time.sleep(0.4)
        log.info("  ChromaDB: indexação HNSW com 312 chunks")
        time.sleep(0.6)
        log.info("  Top-3 chunks recuperados (cosine similarity > 0.82)")
        time.sleep(0.8)
        log.info("  LLM: extraindo coordenadas estruturadas via Pydantic...")
        time.sleep(1.2)
        log.info("  Grid Box:  center=(-10.7, 12.4, 68.9) size=(22, 22, 22)")
        log.info("  Resíduos:  His41, Cys145, Met49, Gly143")
        log.info("  Confiança: 0.87 (acima do threshold 0.3)")
        log.info("  ✓ RAGEngine concluído em 3.8s")

        # --- Etapa 3: TargetScorer ---
        log.info("")
        log.info("[ETAPA 3/5] TargetScorer — Ranqueamento evolutivo")
        log.info("-" * 48)
        time.sleep(0.3)
        log.info("  Calculando dN/dS para 23 centroides...")
        time.sleep(0.7)
        for i, (prot, omega) in enumerate([
            ("3CLpro (P0DTD1)", 0.048), ("RdRp (P0DTD1)", 0.062),
            ("Spike (P0DTC2)", 0.891), ("NSP3 (P0DTD1)", 0.234),
        ]):
            time.sleep(0.3)
            log.info(f"  [{i+1}/23] {prot}: omega={omega:.3f}")

        time.sleep(0.4)
        log.info("  Normalização adaptativa (mu=1e-4, classe=rna_virus)")
        log.info("  Multiplicador metálico: Phi=2.0 (Zn2+ detectado)")
        log.info("  Penalização IDR: Psi=0.92 (4.2% regiões desordenadas)")
        time.sleep(0.5)
        log.info("  Score final (min-max normalizado no batch):")
        log.info("    #1  3CLpro       T=0.9412  ★ Alvo prioritário")
        log.info("    #2  RdRp         T=0.8731")
        log.info("    #3  NSP3         T=0.5124")
        log.info("    #23 Spike        T=0.0891  (seleção positiva — risco)")
        log.info("  ✓ TargetScorer concluído em 2.4s")

        # --- Etapa 4: DockingManager ---
        log.info("")
        log.info("[ETAPA 4/5] DockingManager — AutoDock Vina")
        log.info("-" * 48)
        time.sleep(0.3)
        log.info("  Receptor: 6LU7.pdb (preparado via Meeko)")
        time.sleep(0.4)

        ligands = [
            ("Nirmatrelvir",  "CC(C)(C)NC(=O)...", -8.7),
            ("Aspirina",      "CC(=O)Oc1ccccc1...", -5.2),
            ("Remdesivir",    "CCC(CC)COC(=O)...",  -9.1),
            ("Ibuprofeno",    "CC(C)Cc1ccc(cc1)...", -4.8),
            ("GC-376",        "CC(C)CC(NC(=O)...",  -7.8),
        ]

        for i, (name, smiles, dg) in enumerate(ligands):
            log.info(f"  Ligante {i+1}/5: {name}")
            log.info(f"    SMILES: {smiles}")
            time.sleep(0.3)
            log.info(f"    RDKit: conformer 3D gerado (ETKDG v3 + MMFF94)")
            time.sleep(0.2)
            log.info(f"    Meeko: PDBQT gerado ({name.lower()}.pdbqt)")
            time.sleep(0.8)

            if dg < -7.0:
                log.info(f"    Vina: dG = {dg:.1f} kcal/mol  ★ MATCH (< -7.0)")
            else:
                log.info(f"    Vina: dG = {dg:.1f} kcal/mol")
            time.sleep(0.2)

        log.info("  Garbage Collection: 5 dirs temporários removidos")
        log.info("  ✓ DockingManager concluído em 8.2s")

        # --- Etapa 5: Resultados ---
        log.info("")
        log.info("[ETAPA 5/5] Consolidação de resultados")
        log.info("-" * 48)
        time.sleep(0.5)
        log.info("  Matches (dG < -7.0 kcal/mol):")
        log.info("    ✅ Remdesivir    dG = -9.1 kcal/mol")
        log.info("    ✅ Nirmatrelvir  dG = -8.7 kcal/mol")
        log.info("    ✅ GC-376        dG = -7.8 kcal/mol")
        log.info("  Sem afinidade:")
        log.info("    — Aspirina      dG = -5.2 kcal/mol")
        log.info("    — Ibuprofeno    dG = -4.8 kcal/mol")
        time.sleep(0.3)

        log.info("")
        log.info("=" * 56)
        log.info("PIPELINE CONCLUÍDO — Tempo total: 17.5s")
        log.info(f"  Data: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        log.info(f"  Alvos analisados: 23")
        log.info(f"  Ligantes testados: 5")
        log.info(f"  Matches: 3/5 (60%)")
        log.info("=" * 56)


# ============================================================================
# ENTRY POINT STANDALONE — para demonstração independente
# ============================================================================
def main() -> None:
    """Abre o terminal como janela standalone para demonstração."""
    ctk.set_appearance_mode("dark")
    ctk.set_default_color_theme("dark-blue")

    root = ctk.CTk()
    root.title("F.L.U.I.D — Pipeline Terminal")
    root.geometry("960x620")
    root.minsize(800, 500)
    root.configure(fg_color=C["bg"])

    root.grid_rowconfigure(0, weight=1)
    root.grid_columnconfigure(0, weight=1)

    terminal = LogTerminalFrame(root, show_start_button=True)
    terminal.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)

    root.mainloop()


if __name__ == "__main__":
    main()
