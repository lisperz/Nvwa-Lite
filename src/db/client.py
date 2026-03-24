"""PostgreSQL connection pool — singleton, shared across the process."""

from __future__ import annotations

import logging
import os
from typing import TYPE_CHECKING

logger = logging.getLogger(__name__)

try:
    import psycopg2
    from psycopg2 import pool as pg_pool
    PSYCOPG2_AVAILABLE = True
except ImportError:
    PSYCOPG2_AVAILABLE = False

if TYPE_CHECKING:
    import psycopg2.extensions

_pool: pg_pool.ThreadedConnectionPool | None = None


def get_pool() -> pg_pool.ThreadedConnectionPool | None:
    """Return the shared connection pool, initialising it on first call.

    Returns None (with a logged warning) if DATABASE_URL is not set or
    psycopg2 is not installed — the rest of the app degrades gracefully.
    """
    global _pool  # noqa: PLW0603

    if _pool is not None:
        return _pool

    if not PSYCOPG2_AVAILABLE:
        logger.warning("psycopg2 not installed — DB logging disabled")
        return None

    dsn = os.getenv("DATABASE_URL")
    if not dsn:
        logger.warning("DATABASE_URL not set — DB logging disabled")
        return None

    try:
        _pool = pg_pool.ThreadedConnectionPool(minconn=1, maxconn=5, dsn=dsn)
        logger.info("PostgreSQL connection pool initialised")
    except Exception as e:
        logger.error("Failed to create PostgreSQL pool: %s", e)
        _pool = None

    return _pool


def get_conn():
    """Context manager: borrow a connection from the pool and return it."""
    import contextlib

    @contextlib.contextmanager
    def _ctx():
        pool = get_pool()
        if pool is None:
            yield None
            return
        conn = pool.getconn()
        try:
            yield conn
            conn.commit()
        except Exception:
            conn.rollback()
            raise
        finally:
            pool.putconn(conn)

    return _ctx()
