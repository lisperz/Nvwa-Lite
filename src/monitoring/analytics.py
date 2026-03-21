"""Analytics service for parsing and aggregating structured logs.

Provides real-time metrics for the monitoring dashboard by reading
and analyzing JSON log files.
"""

from __future__ import annotations

import json
from collections import defaultdict
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any


class AnalyticsService:
    """Parse and aggregate structured logs for dashboard display."""

    def __init__(self, log_dir: Path | str = "logs") -> None:
        """Initialize analytics service.

        Args:
            log_dir: Directory containing JSON log files.
        """
        self.log_dir = Path(log_dir)

    def _read_log_file(self, filename: str) -> list[dict[str, Any]]:
        """Read and parse a JSON log file."""
        log_path = self.log_dir / filename
        if not log_path.exists():
            return []

        entries = []
        with open(log_path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    entries.append(json.loads(line))
                except json.JSONDecodeError:
                    continue
        return entries

    def get_active_users(self, hours: int = 24) -> set[str]:
        """Get unique users active in the last N hours."""
        cutoff = datetime.now(timezone.utc) - timedelta(hours=hours)
        users = set()

        for log_file in ["user_interaction.log", "tool_execution.log"]:
            entries = self._read_log_file(log_file)
            for entry in entries:
                try:
                    timestamp = datetime.fromisoformat(entry["timestamp"].replace("Z", "+00:00"))
                    if timestamp >= cutoff:
                        users.add(entry["user_id"])
                except (KeyError, ValueError):
                    continue

        return users

    def get_total_sessions(self, hours: int = 24) -> int:
        """Get total unique sessions in the last N hours."""
        cutoff = datetime.now(timezone.utc) - timedelta(hours=hours)
        sessions = set()

        entries = self._read_log_file("user_interaction.log")
        for entry in entries:
            try:
                timestamp = datetime.fromisoformat(entry["timestamp"].replace("Z", "+00:00"))
                if timestamp >= cutoff:
                    sessions.add(entry["session_id"])
            except (KeyError, ValueError):
                continue

        return len(sessions)

    def get_tool_usage_stats(self, hours: int = 24) -> dict[str, int]:
        """Get tool usage counts in the last N hours."""
        cutoff = datetime.now(timezone.utc) - timedelta(hours=hours)
        tool_counts = defaultdict(int)

        entries = self._read_log_file("tool_execution.log")
        for entry in entries:
            try:
                timestamp = datetime.fromisoformat(entry["timestamp"].replace("Z", "+00:00"))
                if timestamp >= cutoff and entry["event"] == "tool_execution":
                    tool_name = entry["payload"].get("tool_name", "unknown")
                    tool_counts[tool_name] += 1
            except (KeyError, ValueError):
                continue

        return dict(tool_counts)

    def get_error_count(self, hours: int = 24) -> int:
        """Get total error count in the last N hours."""
        cutoff = datetime.now(timezone.utc) - timedelta(hours=hours)
        error_count = 0

        entries = self._read_log_file("tool_execution.log")
        for entry in entries:
            try:
                timestamp = datetime.fromisoformat(entry["timestamp"].replace("Z", "+00:00"))
                if timestamp >= cutoff:
                    status = entry["payload"].get("status")
                    if status == "error":
                        error_count += 1
            except (KeyError, ValueError):
                continue

        return error_count

    def get_total_token_usage(self, hours: int = 24) -> dict[str, int]:
        """Get aggregated token usage in the last N hours."""
        cutoff = datetime.now(timezone.utc) - timedelta(hours=hours)
        total_tokens = 0
        prompt_tokens = 0
        completion_tokens = 0

        entries = self._read_log_file("system_metrics.log")
        for entry in entries:
            try:
                timestamp = datetime.fromisoformat(entry["timestamp"].replace("Z", "+00:00"))
                if timestamp >= cutoff and entry["event"] == "token_usage":
                    payload = entry["payload"]
                    total_tokens += payload.get("total_tokens", 0)
                    prompt_tokens += payload.get("input_tokens", 0)
                    completion_tokens += payload.get("output_tokens", 0)
            except (KeyError, ValueError):
                continue

        return {
            "total": total_tokens,
            "prompt": prompt_tokens,
            "completion": completion_tokens,
        }

    def get_average_response_time(self, hours: int = 24) -> float:
        """Get average tool execution time in the last N hours (in seconds)."""
        cutoff = datetime.now(timezone.utc) - timedelta(hours=hours)
        durations = []

        entries = self._read_log_file("tool_execution.log")
        for entry in entries:
            try:
                timestamp = datetime.fromisoformat(entry["timestamp"].replace("Z", "+00:00"))
                if timestamp >= cutoff and entry["event"] == "tool_execution":
                    duration_ms = entry["payload"].get("duration_ms", 0)
                    if duration_ms > 0:
                        durations.append(duration_ms / 1000.0)
            except (KeyError, ValueError):
                continue

        return sum(durations) / len(durations) if durations else 0.0

    def get_recent_errors(self, limit: int = 10) -> list[dict[str, Any]]:
        """Get most recent errors."""
        entries = self._read_log_file("tool_execution.log")
        errors = []

        for entry in entries:
            try:
                if entry["payload"].get("status") == "error":
                    errors.append({
                        "timestamp": entry["timestamp"],
                        "user_id": entry["user_id"],
                        "session_id": entry["session_id"],
                        "tool_name": entry["payload"].get("tool_name", "unknown"),
                        "error": entry["payload"].get("result", "Unknown error"),
                    })
            except KeyError:
                continue

        # Sort by timestamp descending and return most recent
        errors.sort(key=lambda x: x["timestamp"], reverse=True)
        return errors[:limit]

    def get_user_activity(self, hours: int = 24) -> list[dict[str, Any]]:
        """Get per-user activity summary."""
        cutoff = datetime.now(timezone.utc) - timedelta(hours=hours)
        user_stats = defaultdict(lambda: {"messages": 0, "tools": 0, "sessions": set()})

        # Count messages
        entries = self._read_log_file("user_interaction.log")
        for entry in entries:
            try:
                timestamp = datetime.fromisoformat(entry["timestamp"].replace("Z", "+00:00"))
                if timestamp >= cutoff:
                    user_id = entry["user_id"]
                    user_stats[user_id]["messages"] += 1
                    user_stats[user_id]["sessions"].add(entry["session_id"])
            except (KeyError, ValueError):
                continue

        # Count tool executions
        entries = self._read_log_file("tool_execution.log")
        for entry in entries:
            try:
                timestamp = datetime.fromisoformat(entry["timestamp"].replace("Z", "+00:00"))
                if timestamp >= cutoff and entry["event"] == "tool_execution":
                    user_id = entry["user_id"]
                    user_stats[user_id]["tools"] += 1
            except (KeyError, ValueError):
                continue

        # Convert to list format
        result = []
        for user_id, stats in user_stats.items():
            result.append({
                "user_id": user_id,
                "messages": stats["messages"],
                "tools": stats["tools"],
                "sessions": len(stats["sessions"]),
            })

        result.sort(key=lambda x: x["messages"], reverse=True)
        return result

