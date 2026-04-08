"""Structured logging service for monitoring and analytics.

Implements the logging schema required by the Monitoring Protocol:
- Every log entry includes timestamp, user_id, session_id, task_type, event, payload
- Separate log files for different event types
- JSON format for easy parsing and analysis
"""

from __future__ import annotations

import json
import logging
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any


@dataclass
class LogEntry:
    """Structured log entry following the Monitoring Protocol schema."""

    timestamp: str
    user_id: str
    session_id: str
    task_type: str
    event: str
    payload: dict[str, Any]

    def to_json(self) -> str:
        """Convert log entry to JSON string."""
        return json.dumps(asdict(self), default=str)


class EventLogger:
    """Structured logging for monitoring and analytics.

    Creates separate log files for different event types:
    - tool_execution.log: Tool calls and their results
    - user_interaction.log: User messages and responses
    - system_metrics.log: Token usage, performance metrics
    """

    def __init__(self, log_dir: Path = Path("logs")):
        """Initialize the event logger.

        Args:
            log_dir: Directory to store log files.
        """
        self.log_dir = log_dir
        self.log_dir.mkdir(exist_ok=True)

        # Create separate loggers for different event types
        self.tool_logger = self._create_logger("tool_execution")
        self.interaction_logger = self._create_logger("user_interaction")
        self.metrics_logger = self._create_logger("system_metrics")

    def _create_logger(self, name: str) -> logging.Logger:
        """Create a logger with JSON formatting.

        Args:
            name: Logger name (used for log file name).

        Returns:
            Configured logger instance.
        """
        logger = logging.getLogger(f"nvwa.{name}")
        logger.setLevel(logging.INFO)

        # Avoid duplicate handlers
        if logger.handlers:
            return logger

        # File handler
        log_file = self.log_dir / f"{name}.log"
        handler = logging.FileHandler(log_file)
        handler.setLevel(logging.INFO)

        # No formatting - we'll write JSON directly
        formatter = logging.Formatter("%(message)s")
        handler.setFormatter(formatter)

        logger.addHandler(handler)
        return logger

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
        """Log a tool execution event.

        Args:
            user_id: User identifier.
            session_id: Session identifier.
            tool_name: Name of the tool that was executed.
            args: Tool arguments.
            result: Tool result (full content, no truncation).
            duration_ms: Execution duration in milliseconds.
            status: "success" or "error".
            error_stacktrace: Full traceback on error, None on success.
            turn_id: UUID for the current turn (links tool calls to the user message).
            call_index: 1-based index of this tool call within the turn.
        """
        payload: dict[str, Any] = {
            "args": args,
            "result": result,
            "duration_ms": round(duration_ms, 2),
            "status": status,
        }
        if turn_id is not None:
            payload["turn_id"] = turn_id
        if call_index is not None:
            payload["call_index"] = call_index
        if error_stacktrace is not None:
            payload["error_stacktrace"] = error_stacktrace
        entry = LogEntry(
            timestamp=datetime.utcnow().isoformat() + "Z",
            user_id=user_id,
            session_id=session_id,
            task_type=tool_name,
            event="tool_execution",
            payload=payload,
        )
        self.tool_logger.info(entry.to_json())

    def log_user_message(
        self,
        user_id: str,
        session_id: str,
        message: str,
        response_time_ms: float,
        tool_called: bool = False,
        turn_id: str | None = None,
    ) -> None:
        """Log a user interaction event.

        Args:
            user_id: User identifier.
            session_id: Session identifier.
            message: User's message (full content, no truncation).
            response_time_ms: Total response time in milliseconds.
            tool_called: Whether any tools were called.
            turn_id: UUID for the current turn.
        """
        payload: dict[str, Any] = {
            "message": message,
            "response_time_ms": round(response_time_ms, 2),
            "tool_called": tool_called,
        }
        if turn_id is not None:
            payload["turn_id"] = turn_id
        entry = LogEntry(
            timestamp=datetime.utcnow().isoformat() + "Z",
            user_id=user_id,
            session_id=session_id,
            task_type="chat_message",
            event="user_interaction",
            payload=payload,
        )
        self.interaction_logger.info(entry.to_json())

    def log_token_usage(
        self,
        user_id: str,
        session_id: str,
        input_tokens: int,
        output_tokens: int,
        model: str,
    ) -> None:
        """Log LLM token usage for cost tracking.

        Args:
            user_id: User identifier.
            session_id: Session identifier.
            input_tokens: Number of input tokens.
            output_tokens: Number of output tokens.
            model: Model identifier (e.g., "gpt-4o-mini").
        """
        entry = LogEntry(
            timestamp=datetime.utcnow().isoformat() + "Z",
            user_id=user_id,
            session_id=session_id,
            task_type="llm_call",
            event="token_usage",
            payload={
                "input_tokens": input_tokens,
                "output_tokens": output_tokens,
                "total_tokens": input_tokens + output_tokens,
                "model": model,
            },
        )
        self.metrics_logger.info(entry.to_json())

    def log_session_event(
        self,
        user_id: str,
        session_id: str,
        event_type: str,
        metadata: dict[str, Any] | None = None,
    ) -> None:
        """Log a session lifecycle event.

        Args:
            user_id: User identifier.
            session_id: Session identifier.
            event_type: Type of event ("session_start", "session_end", etc.).
            metadata: Additional event metadata.
        """
        entry = LogEntry(
            timestamp=datetime.utcnow().isoformat() + "Z",
            user_id=user_id,
            session_id=session_id,
            task_type="session",
            event=event_type,
            payload=metadata or {},
        )
        self.interaction_logger.info(entry.to_json())

    def log_error(
        self,
        user_id: str,
        session_id: str,
        error_type: str,
        error_message: str,
        context: dict[str, Any] | None = None,
    ) -> None:
        """Log an error event.

        Args:
            user_id: User identifier.
            session_id: Session identifier.
            error_type: Type of error.
            error_message: Error message.
            context: Additional error context.
        """
        entry = LogEntry(
            timestamp=datetime.utcnow().isoformat() + "Z",
            user_id=user_id,
            session_id=session_id,
            task_type="error",
            event=error_type,
            payload={
                "error_message": error_message,
                "context": context or {},
            },
        )
        self.metrics_logger.info(entry.to_json())
