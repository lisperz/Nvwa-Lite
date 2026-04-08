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

    def ensure_session(
        self,
        user_id: str,
        session_id: str,
        filename: str,
        dataset_metadata: dict[str, Any] | None = None,
    ) -> None:
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
                    INSERT INTO analysis_sessions (session_id, user_id, filename, status, dataset_metadata)
                    VALUES (%s, %s, %s, 'active', %s)
                    ON CONFLICT (session_id) DO NOTHING
                    """,
                    (session_id, user_id, filename, json.dumps(dataset_metadata) if dataset_metadata else None),
                )

    def end_session(self, session_id: str, end_reason: str = "normal") -> None:
        """Mark session as completed and record end time and reason."""
        with get_conn() as conn:
            if conn is None:
                return
            with conn.cursor() as cur:
                cur.execute(
                    """
                    UPDATE analysis_sessions
                    SET status = 'completed', ended_at = NOW(), end_reason = %s
                    WHERE session_id = %s
                    """,
                    (end_reason, session_id),
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
        error_stacktrace: str | None = None,
        turn_id: str | None = None,
        call_index: int | None = None,
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
                             result, duration_ms, status,
                             error_stacktrace, turn_id, call_index)
                        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                        """,
                        (
                            session_id,
                            user_id,
                            tool_name,
                            json.dumps(args),
                            result,
                            round(duration_ms, 2),
                            status,
                            error_stacktrace,
                            turn_id,
                            call_index,
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
        turn_id: str | None = None,
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
                            (session_id, user_id, role, content, tool_called, response_ms, turn_id)
                        VALUES (%s, %s, 'user', %s, %s, %s, %s)
                        """,
                        (
                            session_id,
                            user_id,
                            message,
                            tool_called,
                            round(response_time_ms, 2),
                            turn_id,
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
                            message,
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
        import os

        # Initialize S3 service if configured
        s3_service = None
        if os.getenv("S3_BUCKET_NAME"):
            try:
                from src.storage.service import S3StorageService
                s3_service = S3StorageService(
                    bucket_name=os.getenv("S3_BUCKET_NAME"),
                    region=os.getenv("AWS_REGION", "us-east-2")
                )
            except Exception as e:
                logger.warning(f"S3 service unavailable: {e}")

        try:
            with get_conn() as conn:
                if conn is None:
                    return
                with conn.cursor() as cur:
                    for idx, pr in enumerate(plot_results):
                        image_b64 = base64.b64encode(pr.image).decode("utf-8")

                        # Upload to S3 if available
                        s3_key = None
                        if s3_service:
                            try:
                                filename = f"plot_{message_id}_{idx}.png"
                                s3_key = s3_service.upload_result(user_id, session_id, "plot", pr.image, filename)
                            except Exception as e:
                                logger.warning(f"S3 plot upload failed: {e}")

                        cur.execute(
                            """
                            INSERT INTO message_artifacts
                                (message_id, session_id, user_id, artifact_type,
                                 title, image_b64, code, s3_key)
                            VALUES (%s, %s, %s, 'plot', %s, %s, %s, %s)
                            """,
                            (message_id, session_id, user_id,
                             pr.message[:200], image_b64, pr.code, s3_key),
                        )

                    for idx, tr in enumerate(table_results):
                        # Upload to S3 if available
                        s3_key = None
                        if s3_service:
                            try:
                                filename = f"table_{message_id}_{idx}.csv"
                                s3_key = s3_service.upload_result(user_id, session_id, "table", tr.csv_data.encode('utf-8'), filename)
                            except Exception as e:
                                logger.warning(f"S3 table upload failed: {e}")

                        cur.execute(
                            """
                            INSERT INTO message_artifacts
                                (message_id, session_id, user_id, artifact_type,
                                 title, csv_data, display_df, code, s3_key)
                            VALUES (%s, %s, %s, 'table', %s, %s, %s, %s, %s)
                            """,
                            (message_id, session_id, user_id,
                             tr.message[:200], tr.csv_data, tr.display_df, tr.code, s3_key),
                        )
        except Exception as e:
            logger.error("DB log_artifacts failed: %s", e)

    # ------------------------------------------------------------------
    # Feedback
    # ------------------------------------------------------------------

    def log_feedback(
        self,
        session_id: str,
        user_id: str,
        q1_score: int,
        q2_time_saved: str,
        q3_open_text: str | None = None,
    ) -> None:
        """Persist user feedback response."""
        try:
            with get_conn() as conn:
                if conn is None:
                    return
                with conn.cursor() as cur:
                    cur.execute(
                        """
                        INSERT INTO feedback_responses
                            (session_id, user_id, q1_score, q2_time_saved, q3_open_text)
                        VALUES (%s, %s, %s, %s, %s)
                        """,
                        (session_id, user_id, q1_score, q2_time_saved, q3_open_text),
                    )
        except Exception as e:
            logger.error("DB log_feedback failed: %s", e)


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
