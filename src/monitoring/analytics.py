"""Analytics service — queries RDS PostgreSQL for dashboard metrics.

All public methods return plain Python dicts/lists so the dashboard layer
has no SQL dependency.  Every method catches DB errors and returns empty
results so the dashboard degrades gracefully when the DB is unavailable.
"""

from __future__ import annotations

import decimal
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
                rows = []
                for row in cur.fetchall():
                    rows.append({
                        k: float(v) if isinstance(v, decimal.Decimal) else v
                        for k, v in zip(cols, row)
                    })
                return rows
    except Exception as exc:
        logger.error("Analytics query failed: %s", exc)
        return []


# ---------------------------------------------------------------------------
# Overview KPIs
# ---------------------------------------------------------------------------

def get_overview(hours: int = 24) -> dict[str, Any]:
    """Top-level numbers for the KPI tiles.

    Uses separate queries to avoid fan-out JOINs that inflate counts.
    """
    msg_rows = _query(
        """
        SELECT
            COUNT(DISTINCT user_id)      AS active_users,
            COUNT(DISTINCT session_id)   AS active_sessions,
            COUNT(*)                     AS total_messages
        FROM chat_messages
        WHERE created_at >= NOW() - INTERVAL '%s hours'
        """,
        (hours,),
    )
    tool_rows = _query(
        """
        SELECT
            COUNT(*)                                             AS total_tool_calls,
            SUM(CASE WHEN status = 'error' THEN 1 ELSE 0 END)  AS total_errors
        FROM tool_executions
        WHERE created_at >= NOW() - INTERVAL '%s hours'
        """,
        (hours,),
    )
    token_rows = _query(
        """
        SELECT
            COALESCE(SUM(input_tokens), 0)  AS input_tokens,
            COALESCE(SUM(output_tokens), 0) AS output_tokens
        FROM token_usage
        WHERE created_at >= NOW() - INTERVAL '%s hours'
        """,
        (hours,),
    )

    m = msg_rows[0] if msg_rows else {"active_users": 0, "active_sessions": 0, "total_messages": 0}
    t = tool_rows[0] if tool_rows else {"total_tool_calls": 0, "total_errors": 0}
    tk = token_rows[0] if token_rows else {"input_tokens": 0, "output_tokens": 0}

    return {
        **m, **t, **tk,
        "estimated_cost_usd": _cost(tk["input_tokens"], tk["output_tokens"]),
    }


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
    """One row per user with aggregated stats, sorted by messages desc.

    Uses separate subqueries joined on user_id to avoid fan-out multiplication.
    """
    return _query(
        """
        SELECT
            cm.user_id,
            cm.sessions,
            cm.messages,
            COALESCE(te.tool_calls, 0)                          AS tool_calls,
            COALESCE(tu.input_tokens, 0)                        AS input_tokens,
            COALESCE(tu.output_tokens, 0)                       AS output_tokens,
            COALESCE(tu.input_tokens, 0) + COALESCE(tu.output_tokens, 0) AS total_tokens,
            cm.last_seen
        FROM (
            SELECT user_id,
                   COUNT(DISTINCT session_id) AS sessions,
                   SUM(CASE WHEN role = 'user' THEN 1 ELSE 0 END) AS messages,
                   MAX(created_at) AS last_seen
            FROM chat_messages
            WHERE created_at >= NOW() - INTERVAL '%s hours'
            GROUP BY user_id
        ) cm
        LEFT JOIN (
            SELECT user_id, COUNT(*) AS tool_calls
            FROM tool_executions
            WHERE created_at >= NOW() - INTERVAL '%s hours'
            GROUP BY user_id
        ) te ON te.user_id = cm.user_id
        LEFT JOIN (
            SELECT user_id,
                   SUM(input_tokens)  AS input_tokens,
                   SUM(output_tokens) AS output_tokens
            FROM token_usage
            WHERE created_at >= NOW() - INTERVAL '%s hours'
            GROUP BY user_id
        ) tu ON tu.user_id = cm.user_id
        ORDER BY cm.messages DESC
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
            s.created_at,
            s.ended_at,
            EXTRACT(EPOCH FROM (COALESCE(s.ended_at, NOW()) - s.created_at))::int
                AS duration_seconds
        FROM analysis_sessions s
        WHERE s.created_at >= NOW() - INTERVAL '%s hours'
        ORDER BY s.created_at DESC
        LIMIT %s
        """,
        (hours, limit),
    )


def get_session_messages(session_id: str) -> list[dict[str, Any]]:
    """All chat messages for a specific session, with artifacts."""
    messages = _query(
        """
        SELECT id, role, content, tool_called, response_ms, created_at
        FROM chat_messages
        WHERE session_id = %s
        ORDER BY created_at
        """,
        (session_id,),
    )
    if not messages:
        return []

    # Attach artifacts to each message
    msg_ids = [m["id"] for m in messages]
    placeholders = ",".join(["%s"] * len(msg_ids))
    artifacts = _query(
        f"""
        SELECT message_id, artifact_type, title, image_b64, csv_data, display_df, code
        FROM message_artifacts
        WHERE message_id IN ({placeholders})
        ORDER BY id
        """,
        tuple(msg_ids),
    )

    # Group artifacts by message_id
    from collections import defaultdict
    art_by_msg: dict[int, list] = defaultdict(list)
    for a in artifacts:
        art_by_msg[a["message_id"]].append(a)

    for m in messages:
        m["artifacts"] = art_by_msg.get(m["id"], [])

    return messages


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


# ---------------------------------------------------------------------------
# Feedback responses
# ---------------------------------------------------------------------------

def get_feedback_responses(hours: int = 24) -> list[dict[str, Any]]:
    """Get all feedback responses within time window."""
    return _query(
        """
        SELECT
            response_id,
            session_id,
            user_id,
            q1_score,
            q2_time_saved,
            q3_open_text,
            created_at
        FROM feedback_responses
        WHERE created_at >= NOW() - INTERVAL '%s hours'
        ORDER BY created_at DESC
        """,
        (hours,),
    )


def get_feedback_stats(hours: int = 24) -> dict[str, Any]:
    """Get feedback statistics."""
    rows = _query(
        """
        SELECT
            COUNT(*) as total_responses,
            AVG(q1_score) as avg_score,
            COUNT(CASE WHEN q1_score >= 4 THEN 1 END) as positive_count
        FROM feedback_responses
        WHERE created_at >= NOW() - INTERVAL '%s hours'
        """,
        (hours,),
    )
    return rows[0] if rows else {"total_responses": 0, "avg_score": 0, "positive_count": 0}
