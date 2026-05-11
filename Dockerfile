# =============================================================================
# F.L.U.I.D — Dockerfile Multi-stage
# Base Linux com RDKit, AutoDock Vina, e dependências Python
# =============================================================================

# --- Stage 1: Builder ---
FROM python:3.12-slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY requirements.txt .
RUN pip install --no-cache-dir --prefix=/install -r requirements.txt

# --- Stage 2: Runtime ---
FROM python:3.12-slim

LABEL maintainer="F.L.U.I.D Team"
LABEL description="Drug Discovery Pipeline — API + Celery Worker"

# Dependências de sistema para RDKit e Vina
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxrender1 \
    libxext6 \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Instala AutoDock Vina binário (Linux x86_64)
RUN wget -q https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64 \
    -O /usr/local/bin/vina \
    && chmod +x /usr/local/bin/vina

# Copia dependências Python do builder
COPY --from=builder /install /usr/local

# Cria diretório da aplicação
WORKDIR /app

# Copia código-fonte
COPY src/ ./src/
COPY requirements.txt .

# Cria workspace
RUN mkdir -p /app/workspace/structures /app/workspace/results /app/workspace/chroma_db

# Variáveis de ambiente
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    WORKSPACE_DIR=/app/workspace \
    VINA_BINARY_PATH=/usr/local/bin/vina \
    REDIS_URL=redis://redis:6379/0 \
    LOG_LEVEL=INFO

# Porta da API
EXPOSE 8000

# Default: inicia a API (override no docker-compose para worker)
CMD ["uvicorn", "src.api.main:app", "--host", "0.0.0.0", "--port", "8000"]
