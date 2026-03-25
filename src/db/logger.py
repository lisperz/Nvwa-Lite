"""DatabaseLogger — writes structured events to PostgreSQL.

Mirrors the EventLogger interface so call sites in core.py need no changes.
All methods are fire-and-forget: failures are logged but never raised.
"""

from __future__ import annotations

import json
import logging
from datetime import datetime
from typing import Any

from src.db.client import get_conn

logger = logging.getLogger(__name__)


class DatabaseLogger:
    """Persist agent/user events to RDS PostgreSQL."""

    # ------------------------------------------------------------------
    # Session bootstrap
    # ------------------------------------------------------------------

    def ensure_session(self, user_id: str, session_id: str, filename: str) -> None:
        """Insert analysis_session row if it doesn't exist yet.

        Called once per agent turn before any other log method.
        Safe to call multiple times — uses INSERT … ON CONFLICT DO NOTHING.
        """
        with get_conn() as conn:
            if conn is None:
                return
            with conn.cursor() as cur:
                cur.execute(
                    """
                    INSERT INTO analysis_sessions (session_id, user_id, filename, status)
                    VALUES (%s, %s, %s, 'active')
                    ON CONFLICT (session_id) DO NOTHING
                    """,
                    (session_id, user_id, filename),
                )

    def end_session(self, session_id: str) -> None:
        """Mark session as completed and record end time."""
        with get_conn() as conn:
            if conn is None:
                return
            with conn.cursor() as cur:
                cur.execute(
                    """
                    UPDATE analysis_sessions
                    SET status = 'completed', ended_at = NOW()
                    WHERE session_id = %s
                    """,
                    (session_id,),
                )

    # ------------------------------------------------------------------
    # Tool executions
    # ------------------------------------------------------------------

    def log_tool_execution(
        self,
        user_id: str,
        session_id: str,
        tool_name: str,
        args: dict[str, Any],
        result: str,
        duration_ms: float,
        status: str,
    ) -> None:
        """Persist a tool call record."""
        try:
            with get_conn() as conn:
                if conn is None:
                    return
                with conn.cursor() as cur:
                    cur.execute(
                        """
                        INSERT INTO tool_executions
                            (session_id, user_id, tool_name, args,
                             result_preview, duration_ms, status)
                        VALUES (%s, %s, %s, %s, %s, %s, %s)
                        """,
                        (
                            session_id,
                            user_id,
                            tool_name,
                            json.dumps(args),
                            result[:500] if len(result) > 500 else result,
                            round(duration_ms, 2),
                            status,
                        ),
                    )
        except Exception as e:
            logger.error("DB log_tool_execution failed: %s", e)

    # ------------------------------------------------------------------
    # Chat messages
    # ------------------------------------------------------------------

    def log_user_message(
        self,
        user_id: str,
        session_id: str,
        message: str,
        response_time_ms: float,
        tool_called: bool = False,
    ) -> None:
        """Persist user message + agent turn metadata."""
        try:
            with get_conn() as conn:
                if conn is None:
                    return
                with conn.cursor() as cur:
                    # User message row
                    cur.execute(
                        """
                        INSERT INTO chat_messages
                            (session_id, user_id, role, content, tool_called, response_ms)
                        VALUES (%s, %s, 'user', %s, %s, %s)
                        """,
                        (
                            session_id,
                            user_id,
                            message[:2000] if len(message) > 2000 else message,
                            tool_called,
                            round(response_time_ms, 2),
                        ),
                    )
                    # Bump message counter on the session
                    cur.execute(
                        """
                        UPDATE analysis_sessions
                        SET message_count = message_count + 1
                        WHERE session_id = %s
                        """,
                        (session_id,),
                    )
        except Exception as e:
            logger.error("DB log_user_message failed: %s", e)

    def log_assistant_message(
        self,
        user_id: str,
        session_id: str,
        message: str,
    ) -> int | None:
        """Persist assistant response to chat_messages. Returns the new row id."""
        try:
            with get_conn() as conn:
                if conn is None:
                    return None
                with conn.cursor() as cur:
                    cur.execute(
                        """
                        INSERT INTO chat_messages
                            (session_id, user_id, role, content, tool_called, response_ms)
                        VALUES (%s, %s, 'assistant', %s, false, null)
                        RETURNING id
                        """,
                        (
                            session_id,
                            user_id,
                            message[:2000] if len(message) > 2000 else message,
                        ),
                    )
                    return cur.fetchone()[0]
        except Exception as e:
            logger.error("DB log_assistant_message failed: %s", e)
            return None

    def log_artifacts(
        self,
        message_id: int,
        session_id: str,
        user_id: str,
        plot_results: list,
        table_results: list,
    ) -> None:
        """Persist plots and tables linked to an assistant message."""
        import base64
        try:
            with get_conn() as conn:
                if conn is None:
                    return
                with conn.cursor() as cur:
                    for pr in plot_results:
                        image_b64 = base64.b64encode(pr.image).decode("utf-8")
                        cur.execute(
                            """
                            INSERT INTO message_artifacts
                                (message_id, session_id, user_id, artifact_type,
                                 title, image_b64, code)
                            VALUES (%s, %s, %s, 'plot', %s, %s, %s)
                            """,
                            (message_id, session_id, user_id,
                             pr.message[:200], image_b64, pr.code),
                        )
                    for tr in table_results:
                        cur.execute(
                            """
                            INSERT INTO message_artifacts
                                (message_id, session_id, user_id, artifact_type,
                                 title, csv_data, display_df, code)
                            VALUES (%s, %s, %s, 'table', %s, %s, %s, %s)
                            """,
                            (message_id, session_id, user_id,
                             tr.message[:200], tr.csv_data, tr.display_df, tr.code),
                        )
        except Exception as e:
            logger.error("DB log_artifacts failed: %s", e)

    # ------------------------------------------------------------------
    # Token usage
    # ------------------------------------------------------------------

    def log_token_usage(
        self,
        user_id: str,
        session_id: str,
        input_tokens: int,
        output_tokens: int,
        model: str,
    ) -> None:
        """Persist LLM token consumption for cost tracking."""
        try:
            with get_conn() as conn:
                if conn is None:
                    return
                with conn.cursor() as cur:
                    cur.execute(
                        """
                        INSERT INTO token_usage
                            (session_id, user_id, model, input_tokens, output_tokens)
                        VALUES (%s, %s, %s, %s, %s)
                        """,
                        (session_id, user_id, model, input_tokens, output_tokens),
                    )
        except Exception as e:
            logger.error("DB log_token_usage failed: %s", e)
