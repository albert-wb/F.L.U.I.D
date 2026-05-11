# =============================================================================
# F.L.U.I.D — Celery App Configuration
# =============================================================================
"""
Configuração central da aplicação Celery para processamento assíncrono.

O broker padrão é Redis (configurável via REDIS_URL no .env).
Workers executam tarefas CPU-bound de docking e scoring em background,
evitando timeout na API/dashboard.

Para iniciar um worker:
    celery -A src.workers.celery_app worker --loglevel=info --concurrency=4
"""
from __future__ import annotations

import os

from celery import Celery
from dotenv import load_dotenv

load_dotenv()

BROKER_URL = os.getenv("REDIS_URL", "redis://localhost:6379/0")
RESULT_BACKEND = os.getenv("REDIS_URL", "redis://localhost:6379/0")

celery_app = Celery(
    "fluid",
    broker=BROKER_URL,
    backend=RESULT_BACKEND,
    include=[
        "src.workers.docking_tasks",
        "src.workers.scoring_tasks",
    ],
)

celery_app.conf.update(
    # Serialização
    task_serializer="json",
    result_serializer="json",
    accept_content=["json"],

    # Timezone
    timezone="America/Sao_Paulo",
    enable_utc=True,

    # Limites de segurança
    task_soft_time_limit=3600,     # 1h soft limit
    task_time_limit=7200,          # 2h hard kill
    task_acks_late=True,           # ACK após conclusão (retry-safe)
    worker_prefetch_multiplier=1,  # 1 task por vez (CPU-bound)

    # Resultados expiram em 24h
    result_expires=86400,

    # Retry
    task_default_retry_delay=60,
    task_max_retries=3,
)
