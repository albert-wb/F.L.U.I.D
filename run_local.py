# =============================================================================
# F.L.U.I.D — Script de Inicialização Local (Sem Docker)
# =============================================================================
"""
Inicia o F.L.U.I.D em modo desenvolvimento local:
    1. Verifica dependências
    2. Inicia a API FastAPI (uvicorn)
    3. Inicia o Dashboard Dash

Modo: Standalone (sem Redis/Celery — tasks executam inline)

Uso:
    python run_local.py          # Inicia API + GUI Desktop
    python run_local.py --api    # Só a API
    python run_local.py --gui    # Só a GUI Desktop
    python run_local.py --dash   # Só o Dashboard web
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import os
import time
from pathlib import Path

# Força UTF-8 no stdout (Windows usa cp1252 por padrão)
if sys.stdout and hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


def check_dependencies() -> bool:
    """Verifica se as dependências mínimas estão instaladas."""
    missing = []
    required = ["fastapi", "uvicorn", "dash", "numpy", "pandas", "loguru", "pydantic"]

    for pkg in required:
        try:
            __import__(pkg)
        except ImportError:
            missing.append(pkg)

    if missing:
        print(f"\n❌ Dependências faltando: {', '.join(missing)}")
        print(f"\n   Instale com:")
        print(f"   pip install -r requirements.txt\n")
        return False

    print("✅ Dependências verificadas")
    return True


def _free_port(port: int) -> None:
    """Verifica se a porta está em uso e mata o processo ocupante."""
    import socket
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        if s.connect_ex(("127.0.0.1", port)) == 0:
            # Porta ocupada — mata via taskkill
            try:
                r = subprocess.run(
                    ["powershell", "-Command",
                     f"(Get-NetTCPConnection -LocalPort {port} "
                     f"-ErrorAction SilentlyContinue).OwningProcess"],
                    capture_output=True, text=True, timeout=5,
                )
                for pid in set(r.stdout.strip().split()):
                    if pid.isdigit() and int(pid) > 0:
                        subprocess.run(
                            ["taskkill", "/F", "/PID", pid],
                            capture_output=True, timeout=5,
                        )
            except Exception:
                pass
            time.sleep(1.5)  # Aguarda o SO liberar o socket


def start_api(port: int = 8000) -> subprocess.Popen:
    """Inicia a API FastAPI com uvicorn."""
    _free_port(port)

    print(f"\n🚀 Iniciando API em http://localhost:{port}")
    print(f"   📖 Docs: http://localhost:{port}/docs")
    print(f"   ❤️  Health: http://localhost:{port}/api/v1/health\n")

    return subprocess.Popen(
        [
            sys.executable, "-m", "uvicorn",
            "src.api.main:app",
            "--host", "127.0.0.1",
            "--port", str(port),
        ],
        cwd=str(Path(__file__).parent),
    )


def start_dashboard(port: int = 8050) -> subprocess.Popen:
    """Inicia o Dashboard Dash."""
    print(f"\n🌊 Iniciando Dashboard em http://localhost:{port}\n")

    return subprocess.Popen(
        [sys.executable, "-c", """
from src.dashboard.visualizer import Visualizer
viz = Visualizer()
viz.run(host='127.0.0.1', port={port}, debug=True)
""".replace("{port}", str(port))],
        cwd=str(Path(__file__).parent),
    )


def start_gui() -> subprocess.Popen:
    """Inicia o cliente desktop CustomTkinter."""
    print("\n🖥️  Iniciando GUI Desktop...\n")
    return subprocess.Popen(
        [sys.executable, "-c",
         "from src.desktop.main_gui import main; main()"],
        cwd=str(Path(__file__).parent),
    )


def main():
    parser = argparse.ArgumentParser(description="F.L.U.I.D — Launcher")
    parser.add_argument("--api", action="store_true", help="Inicia apenas a API")
    parser.add_argument("--gui", action="store_true", help="Inicia apenas a GUI Desktop")
    parser.add_argument("--dash", action="store_true", help="Inicia apenas o Dashboard web")
    parser.add_argument("--api-port", type=int, default=8000)
    parser.add_argument("--dash-port", type=int, default=8050)
    args = parser.parse_args()

    print("""
    +===========================================================+
    |                                                           |
    |   F.L.U.I.D -- Drug Discovery Platform                   |
    |   Framework for Ligand-target Unified In-silico Discovery |
    |                                                           |
    +===========================================================+
    """)

    if not check_dependencies():
        sys.exit(1)

    # Se nenhum flag, inicia API + GUI Desktop
    no_flags = not args.api and not args.gui and not args.dash

    processes = []

    try:
        if args.api or no_flags:
            processes.append(start_api(args.api_port))
            time.sleep(2)  # Espera API subir

        if args.gui or no_flags:
            processes.append(start_gui())

        if args.dash:
            processes.append(start_dashboard(args.dash_port))

        print("\n" + "=" * 50)
        print("  Pressione Ctrl+C para parar todos os serviços")
        print("=" * 50 + "\n")

        # Mantém rodando
        for p in processes:
            p.wait()

    except KeyboardInterrupt:
        print("\n\n🛑 Parando serviços...")
        for p in processes:
            p.terminate()
        print("✅ Todos os serviços parados.\n")


if __name__ == "__main__":
    main()
