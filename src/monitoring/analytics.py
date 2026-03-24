"""Analytics service — queries RDS PostgreSQL for dashboard metrics.

All public methods return plain Python dicts/lists so the dashboard layer
has no SQL dependency.  Every method catches DB errors and returns empty
results so the dashboard degrades gracefully when the DB is unavailable.
"""

from __future__ import annotations

import logging
from typing import Any

from src.db.client import get_conn

logger = logging.getLogger(__name__)

# gpt-4o-mini pricing (USD per 1 000 tokens, as of 2025)
_INPUT_COST_PER_1K = 0.000150
_OUTPUT_COST_PER_1K = 0.000600


def _cost(input_tokens: int, output_tokens: int) -> float:
    return (
        input_tokens / 1000 * _INPUT_COST_PER_1K
        + output_tokens / 1000 * _OUTPUT_COST_PER_1K
    )


def _query(sql: str, params: tuple = ()) -> list[dict[str, Any]]:
    """Run a SELECT and return rows as dicts.  Returns [] on any error."""
    try:
        with get_conn() as conn:
            if conn is None:
                return []
            with conn.cursor() as cur:
                cur.execute(sql, params)
                cols = [d[0] for d in cur.description]
                return [dict(zip(cols, row)) for row in cur.fetchall()]
    except Exception as exc:
        logger.error("Analytics query failed: %s", exc)
        return []


# ---------------------------------------------------------------------------
# Overview KPIs
# ---------------------------------------------------------------------------

def get_overview(hours: int = 24) -> dict[str, Any]:
    """Top-level numbers for the KPI tiles."""
    rows = _query(
        """
        SELECT
            COUNT(DISTINCT cm.user_id)                          AS active_users,
            COUNT(DISTINCT cm.session_id)                       AS active_sessions,
            COUNT(*)                                            AS total_messages,
            COALESCE(SUM(te.cnt), 0)                           AS total_tool_calls,
            COALESCE(SUM(te.errors), 0)                        AS total_errors,
            COALESCE(SUM(tu.input_tokens), 0)                  AS input_tokens,
            COALESCE(SUM(tu.output_tokens), 0)                 AS output_tokens
        FROM chat_messages cm
        LEFT JOIN (
            SELECT session_id,
                   COUNT(*) AS cnt,
                   SUM(CASE WHEN status = 'error' THEN 1 ELSE 0 END) AS errors
            FROM tool_executions
            WHERE created_at >= NOW() - INTERVAL '%s hours'
            GROUP BY session_id
        ) te ON te.session_id = cm.session_id
        LEFT JOIN (
            SELECT session_id,
                   SUM(input_tokens) AS input_tokens,
                   SUM(output_tokens) AS output_tokens
            FROM token_usage
            WHERE created_at >= NOW() - INTERVAL '%s hours'
            GROUP BY session_id
        ) tu ON tu.session_id = cm.session_id
        WHERE cm.created_at >= NOW() - INTERVAL '%s hours'
        """,
        (hours, hours, hours),
    )
    if not rows:
        return {
            "active_users": 0, "active_sessions": 0, "total_messages": 0,
            "total_tool_calls": 0, "total_errors": 0,
            "input_tokens": 0, "output_tokens": 0, "estimated_cost_usd": 0.0,
        }
    r = rows[0]
    r["estimated_cost_usd"] = _cost(r["input_tokens"], r["output_tokens"])
    return r


# ---------------------------------------------------------------------------
# Hourly activity (for sparkline / bar chart)
# ---------------------------------------------------------------------------

def get_hourly_activity(hours: int = 24) -> list[dict[str, Any]]:
    """Message count grouped by hour, ordered ascending."""
    return _query(
        """
        SELECT
            date_trunc('hour', created_at) AS hour,
            COUNT(*) AS messages
        FROM chat_messages
        WHERE created_at >= NOW() - INTERVAL '%s hours'
          AND role = 'user'
        GROUP BY 1
        ORDER BY 1
        """,
        (hours,),
    )


# ---------------------------------------------------------------------------
# Per-user breakdown
# ---------------------------------------------------------------------------

def get_user_breakdown(hours: int = 24) -> list[dict[str, Any]]:
    """One row per user with aggregated stats, sorted by messages desc."""
    return _query(
        """
        SELECT
            cm.user_id,
            COUNT(DISTINCT cm.session_id)                       AS sessions,
            SUM(CASE WHEN cm.role = 'user' THEN 1 ELSE 0 END)  AS messages,
            COALESCE(SUM(te.tool_calls), 0)                    AS tool_calls,
            COALESCE(SUM(tu.total_tokens), 0)                  AS total_tokens,
            COALESCE(SUM(tu.input_tokens), 0)                  AS input_tokens,
            COALESCE(SUM(tu.output_tokens), 0)                 AS output_tokens,
            MAX(cm.created_at)                                  AS last_seen
        FROM chat_messages cm
        LEFT JOIN (
            SELECT user_id, session_id, COUNT(*) AS tool_calls
            FROM tool_executions
            WHERE created_at >= NOW() - INTERVAL '%s hours'
            GROUP BY user_id, session_id
        ) te ON te.user_id = cm.user_id AND te.session_id = cm.session_id
        LEFT JOIN (
            SELECT user_id, session_id,
                   SUM(input_tokens + output_tokens) AS total_tokens,
                   SUM(input_tokens) AS input_tokens,
                   SUM(output_tokens) AS output_tokens
            FROM token_usage
            WHERE created_at >= NOW() - INTERVAL '%s hours'
            GROUP BY user_id, session_id
        ) tu ON tu.user_id = cm.user_id AND tu.session_id = cm.session_id
        WHERE cm.created_at >= NOW() - INTERVAL '%s hours'
        GROUP BY cm.user_id
        ORDER BY messages DESC
        """,
        (hours, hours, hours),
    )


# ---------------------------------------------------------------------------
# Sessions
# ---------------------------------------------------------------------------

def get_recent_sessions(hours: int = 24, limit: int = 50) -> list[dict[str, Any]]:
    """Recent analysis sessions with metadata."""
    return _query(
        """
        SELECT
            s.session_id,
            s.user_id,
            s.filename,
            s.status,
            s.message_count,
            s.started_at,
            s.ended_at,
            EXTRACT(EPOCH FROM (COALESCE(s.ended_at, NOW()) - s.started_at))::int
                AS duration_seconds
        FROM analysis_sessions s
        WHERE s.started_at >= NOW() - INTERVAL '%s hours'
        ORDER BY s.started_at DESC
        LIMIT %s
        """,
        (hours, limit),
    )


def get_session_messages(session_id: str) -> list[dict[str, Any]]:
    """All chat messages for a specific session."""
    return _query(
        """
        SELECT role, content, tool_called, response_ms, created_at
        FROM chat_messages
        WHERE session_id = %s
        ORDER BY created_at
        """,
        (session_id,),
    )


# ---------------------------------------------------------------------------
# Tool analytics
# ---------------------------------------------------------------------------

def get_tool_stats(hours: int = 24) -> list[dict[str, Any]]:
    """Per-tool usage count, error count, avg duration — sorted by calls desc."""
    return _query(
        """
        SELECT
            tool_name,
            COUNT(*)                                                    AS total_calls,
            SUM(CASE WHEN status = 'error' THEN 1 ELSE 0 END)         AS errors,
            ROUND(AVG(duration_ms)::numeric, 1)                        AS avg_ms,
            ROUND(MAX(duration_ms)::numeric, 1)                        AS max_ms
        FROM tool_executions
        WHERE created_at >= NOW() - INTERVAL '%s hours'
        GROUP BY tool_name
        ORDER BY total_calls DESC
        """,
        (hours,),
    )


# ---------------------------------------------------------------------------
# Recent errors
# ---------------------------------------------------------------------------

def get_recent_errors(hours: int = 24, limit: int = 20) -> list[dict[str, Any]]:
    """Most recent failed tool executions."""
    return _query(
        """
        SELECT
            created_at,
            user_id,
            session_id,
            tool_name,
            result_preview AS error_preview
        FROM tool_executions
        WHERE status = 'error'
          AND created_at >= NOW() - INTERVAL '%s hours'
        ORDER BY created_at DESC
        LIMIT %s
        """,
        (hours, limit),
    )
