# =============================================================================
# F.L.U.I.D — Celery Stub (Modo Local sem Redis)
# =============================================================================
"""
Quando Celery/Redis não estão disponíveis, este módulo fornece um
stub que executa tasks de forma síncrona (inline), permitindo rodar
o F.L.U.I.D localmente sem infraestrutura de message broker.

Em produção com Docker, o celery_app.py real é usado.
"""
from __future__ import annotations
from typing import Any, Callable


class _InlineResult:
    """Simula AsyncResult do Celery."""
    def __init__(self, result: Any) -> None:
        self.result = result
        self.id = "inline-task"
        self.status = "SUCCESS"

    def get(self, timeout: int = None) -> Any:
        return self.result


class _InlineTask:
    """Stub de task que executa a função inline (síncrona)."""
    def __init__(self, func: Callable) -> None:
        self._func = func
        self.name = getattr(func, '__name__', 'task')

    def delay(self, *args: Any, **kwargs: Any) -> _InlineResult:
        """Executa inline em vez de despachar para broker."""
        try:
            result = self._func(None, *args, **kwargs)  # None = self (bind)
            return _InlineResult(result)
        except Exception:
            return _InlineResult(None)

    def __call__(self, *args: Any, **kwargs: Any) -> Any:
        return self._func(*args, **kwargs)


class _StubCeleryApp:
    """
    Stub que imita a API do Celery para desenvolvimento local.
    Tasks decoradas com @celery_app.task(...) são registradas
    e executadas inline quando .delay() é chamado.
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        self._tasks: dict[str, _InlineTask] = {}

    def task(self, *args: Any, **kwargs: Any) -> Callable:
        """Decorador que registra a task para execução inline."""
        def decorator(func: Callable) -> _InlineTask:
            task = _InlineTask(func)
            task.name = kwargs.get("name", func.__name__)
            self._tasks[task.name] = task
            return task
        # Suporta @celery_app.task e @celery_app.task(bind=True, ...)
        if args and callable(args[0]):
            return _InlineTask(args[0])
        return decorator

    @property
    def conf(self) -> Any:
        class _Conf:
            def update(self, *a: Any, **kw: Any) -> None:
                pass
        return _Conf()


# Tenta importar Celery real; fallback para stub
try:
    from celery import Celery as _RealCelery
    _HAS_CELERY = True
except ImportError:
    _HAS_CELERY = False


def create_celery_app() -> Any:
    """Factory que retorna Celery real ou stub conforme disponibilidade."""
    if _HAS_CELERY:
        import os
        from dotenv import load_dotenv
        load_dotenv()

        BROKER_URL = os.getenv("REDIS_URL", "redis://localhost:6379/0")
        app = _RealCelery(
            "fluid",
            broker=BROKER_URL,
            backend=BROKER_URL,
            include=[
                "src.workers.docking_tasks",
                "src.workers.scoring_tasks",
            ],
        )
        app.conf.update(
            task_serializer="json",
            result_serializer="json",
            accept_content=["json"],
            timezone="America/Sao_Paulo",
            enable_utc=True,
            task_soft_time_limit=3600,
            task_time_limit=7200,
            task_acks_late=True,
            worker_prefetch_multiplier=1,
            result_expires=86400,
        )
        return app
    else:
        import warnings
        warnings.warn(
            "⚠️  Celery/Redis não encontrados — usando modo INLINE (dev). "
            "Tasks serão executadas sincronamente.",
            RuntimeWarning,
            stacklevel=2,
        )
        return _StubCeleryApp()


celery_app = create_celery_app()
