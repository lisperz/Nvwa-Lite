"""Unit tests for logging changes (Changes 1–6, logging audit 2026-04-07).

Tests call EventLogger and DatabaseLogger methods directly — no live Postgres,
no agent invocation needed. DatabaseLogger DB calls are mocked at the
connection level. EventLogger writes to a temp directory.

Run with: pytest tests/unit/test_logging.py -v
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from src.logging.service import EventLogger


@pytest.fixture(autouse=True)
def reset_event_loggers():
    """Clear EventLogger file handlers between tests.

    EventLogger uses module-level named loggers (nvwa.*) which are singletons.
    Without this, the handler guard in _create_logger returns the old handler
    from the previous test's tmp_path, causing FileNotFoundError.
    """
    for name in ["nvwa.tool_execution", "nvwa.user_interaction", "nvwa.system_metrics"]:
        lg = logging.getLogger(name)
        for h in lg.handlers[:]:
            h.close()
            lg.removeHandler(h)
    yield


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_event_logger(tmp_path: Path) -> EventLogger:
    return EventLogger(log_dir=tmp_path)


def last_log_line(log_file: Path) -> dict:
    lines = [l for l in log_file.read_text().strip().splitlines() if l]
    return json.loads(lines[-1])


# ---------------------------------------------------------------------------
# Change 1 — No truncation
# ---------------------------------------------------------------------------

class TestNoTruncation:
    def test_tool_result_full_content(self, tmp_path):
        el = make_event_logger(tmp_path)
        long_result = "x" * 3000
        el.log_tool_execution(
            user_id="u1", session_id="s1", tool_name="tool_a",
            args={}, result=long_result, duration_ms=100.0, status="success",
        )
        entry = last_log_line(tmp_path / "tool_execution.log")
        assert entry["payload"]["result"] == long_result
        assert len(entry["payload"]["result"]) == 3000

    def test_tool_result_no_result_preview_key(self, tmp_path):
        el = make_event_logger(tmp_path)
        el.log_tool_execution(
            user_id="u1", session_id="s1", tool_name="tool_a",
            args={}, result="short", duration_ms=10.0, status="success",
        )
        entry = last_log_line(tmp_path / "tool_execution.log")
        assert "result_preview" not in entry["payload"]
        assert "result" in entry["payload"]

    def test_user_message_full_content(self, tmp_path):
        el = make_event_logger(tmp_path)
        long_msg = "q" * 3000
        el.log_user_message(
            user_id="u1", session_id="s1", message=long_msg,
            response_time_ms=500.0,
        )
        entry = last_log_line(tmp_path / "user_interaction.log")
        assert entry["payload"]["message"] == long_msg
        assert len(entry["payload"]["message"]) == 3000

    def test_user_message_no_preview_key(self, tmp_path):
        el = make_event_logger(tmp_path)
        el.log_user_message(
            user_id="u1", session_id="s1", message="hello",
            response_time_ms=100.0,
        )
        entry = last_log_line(tmp_path / "user_interaction.log")
        assert "message_preview" not in entry["payload"]
        assert "message" in entry["payload"]


# ---------------------------------------------------------------------------
# Change 2 — Error stacktrace
# ---------------------------------------------------------------------------

class TestErrorStacktrace:
    def test_stacktrace_present_on_error(self, tmp_path):
        el = make_event_logger(tmp_path)
        tb = "Traceback (most recent call last):\n  File 'tools.py', line 42\nKeyError: 'x'"
        el.log_tool_execution(
            user_id="u1", session_id="s1", tool_name="tool_a",
            args={}, result="err", duration_ms=50.0, status="error",
            error_stacktrace=tb,
        )
        entry = last_log_line(tmp_path / "tool_execution.log")
        assert entry["payload"]["error_stacktrace"] == tb

    def test_stacktrace_absent_on_success(self, tmp_path):
        el = make_event_logger(tmp_path)
        el.log_tool_execution(
            user_id="u1", session_id="s1", tool_name="tool_a",
            args={}, result="ok", duration_ms=50.0, status="success",
        )
        entry = last_log_line(tmp_path / "tool_execution.log")
        assert "error_stacktrace" not in entry["payload"]

    def test_stacktrace_captured_from_real_exception(self, tmp_path):
        import traceback
        try:
            raise ValueError("cluster 999 not found")
        except ValueError:
            tb = traceback.format_exc()

        el = make_event_logger(tmp_path)
        el.log_tool_execution(
            user_id="u1", session_id="s1", tool_name="de_tool",
            args={}, result="err", duration_ms=30.0, status="error",
            error_stacktrace=tb,
        )
        entry = last_log_line(tmp_path / "tool_execution.log")
        assert "ValueError" in entry["payload"]["error_stacktrace"]
        assert "cluster 999 not found" in entry["payload"]["error_stacktrace"]


# ---------------------------------------------------------------------------
# Change 3 — Dataset metadata at session start (DatabaseLogger)
# ---------------------------------------------------------------------------

class TestDatasetMetadata:
    def _make_db_logger(self):
        from src.db.logger import DatabaseLogger
        return DatabaseLogger()

    def _run_ensure_session(self, metadata):
        db = self._make_db_logger()
        captured = {}

        def fake_execute(sql, params):
            captured["sql"] = sql
            captured["params"] = params

        mock_cur = MagicMock()
        mock_cur.execute.side_effect = fake_execute
        mock_conn = MagicMock()
        mock_conn.__enter__ = MagicMock(return_value=mock_conn)
        mock_conn.__exit__ = MagicMock(return_value=False)
        mock_conn.cursor.return_value.__enter__ = MagicMock(return_value=mock_cur)
        mock_conn.cursor.return_value.__exit__ = MagicMock(return_value=False)

        with patch("src.db.logger.get_conn", return_value=mock_conn):
            db.ensure_session("u1", "s1", "test.h5ad", dataset_metadata=metadata)

        return captured

    def test_full_metadata_stored(self):
        meta = {
            "n_cells": 2638, "n_genes": 1838,
            "has_umap": True, "has_clustering": True, "clustering_key": "leiden",
        }
        captured = self._run_ensure_session(meta)
        params = captured["params"]
        stored = json.loads(params[3])
        assert stored["n_cells"] == 2638
        assert stored["clustering_key"] == "leiden"
        assert stored["has_umap"] is True

    def test_has_umap_false(self):
        # Real production state: GSE223414 lost UMAP in h5ad conversion
        meta = {
            "n_cells": 4832, "n_genes": 22150,
            "has_umap": False, "has_clustering": True, "clustering_key": "leiden",
        }
        captured = self._run_ensure_session(meta)
        stored = json.loads(captured["params"][3])
        assert stored["has_umap"] is False

    def test_has_clustering_false(self):
        meta = {
            "n_cells": 1000, "n_genes": 500,
            "has_umap": False, "has_clustering": False, "clustering_key": None,
        }
        captured = self._run_ensure_session(meta)
        stored = json.loads(captured["params"][3])
        assert stored["has_clustering"] is False
        assert stored["clustering_key"] is None

    def test_no_metadata_stores_null(self):
        captured = self._run_ensure_session(None)
        assert captured["params"][3] is None


# ---------------------------------------------------------------------------
# Change 4 — end_reason
# ---------------------------------------------------------------------------

class TestEndReason:
    def _run_end_session(self, reason):
        from src.db.logger import DatabaseLogger
        db = DatabaseLogger()
        captured = {}

        def fake_execute(sql, params):
            captured["params"] = params

        mock_cur = MagicMock()
        mock_cur.execute.side_effect = fake_execute
        mock_conn = MagicMock()
        mock_conn.__enter__ = MagicMock(return_value=mock_conn)
        mock_conn.__exit__ = MagicMock(return_value=False)
        mock_conn.cursor.return_value.__enter__ = MagicMock(return_value=mock_cur)
        mock_conn.cursor.return_value.__exit__ = MagicMock(return_value=False)

        with patch("src.db.logger.get_conn", return_value=mock_conn):
            db.end_session("s1", end_reason=reason)

        return captured

    @pytest.mark.parametrize("reason", ["normal", "error", "max_iterations"])
    def test_end_reason_written(self, reason):
        captured = self._run_end_session(reason)
        assert captured["params"][0] == reason

    def test_default_end_reason_is_normal(self):
        from src.db.logger import DatabaseLogger
        db = DatabaseLogger()
        captured = {}

        def fake_execute(sql, params):
            captured["params"] = params

        mock_cur = MagicMock()
        mock_cur.execute.side_effect = fake_execute
        mock_conn = MagicMock()
        mock_conn.__enter__ = MagicMock(return_value=mock_conn)
        mock_conn.__exit__ = MagicMock(return_value=False)
        mock_conn.cursor.return_value.__enter__ = MagicMock(return_value=mock_cur)
        mock_conn.cursor.return_value.__exit__ = MagicMock(return_value=False)

        with patch("src.db.logger.get_conn", return_value=mock_conn):
            db.end_session("s1")  # no reason arg

        assert captured["params"][0] == "normal"


# ---------------------------------------------------------------------------
# Change 5 — turn_id and call_index
# ---------------------------------------------------------------------------

class TestTurnIdAndCallIndex:
    def test_turn_id_in_tool_log(self, tmp_path):
        el = make_event_logger(tmp_path)
        el.log_tool_execution(
            user_id="u1", session_id="s1", tool_name="tool_a",
            args={}, result="ok", duration_ms=10.0, status="success",
            turn_id="turn-abc", call_index=1,
        )
        entry = last_log_line(tmp_path / "tool_execution.log")
        assert entry["payload"]["turn_id"] == "turn-abc"

    def test_call_index_increments(self, tmp_path):
        el = make_event_logger(tmp_path)
        el.log_tool_execution(
            user_id="u1", session_id="s1", tool_name="violin_plot",
            args={}, result="ok", duration_ms=10.0, status="success",
            turn_id="turn-abc", call_index=1,
        )
        el.log_tool_execution(
            user_id="u1", session_id="s1", tool_name="feature_plot",
            args={}, result="ok", duration_ms=10.0, status="success",
            turn_id="turn-abc", call_index=2,
        )
        lines = [json.loads(l) for l in (tmp_path / "tool_execution.log").read_text().strip().splitlines()]
        assert lines[0]["payload"]["call_index"] == 1
        assert lines[1]["payload"]["call_index"] == 2
        assert lines[0]["payload"]["turn_id"] == lines[1]["payload"]["turn_id"]

    def test_turn_id_in_user_message_log(self, tmp_path):
        el = make_event_logger(tmp_path)
        el.log_user_message(
            user_id="u1", session_id="s1", message="hello",
            response_time_ms=100.0, turn_id="turn-abc",
        )
        entry = last_log_line(tmp_path / "user_interaction.log")
        assert entry["payload"]["turn_id"] == "turn-abc"

    def test_different_turns_have_different_turn_ids(self, tmp_path):
        el = make_event_logger(tmp_path)
        el.log_tool_execution(
            user_id="u1", session_id="s1", tool_name="tool_a",
            args={}, result="ok", duration_ms=10.0, status="success",
            turn_id="turn-001", call_index=1,
        )
        el.log_tool_execution(
            user_id="u1", session_id="s1", tool_name="tool_a",
            args={}, result="ok", duration_ms=10.0, status="success",
            turn_id="turn-002", call_index=1,
        )
        lines = [json.loads(l) for l in (tmp_path / "tool_execution.log").read_text().strip().splitlines()]
        assert lines[0]["payload"]["turn_id"] != lines[1]["payload"]["turn_id"]

    def test_turn_id_absent_when_not_passed(self, tmp_path):
        el = make_event_logger(tmp_path)
        el.log_tool_execution(
            user_id="u1", session_id="s1", tool_name="tool_a",
            args={}, result="ok", duration_ms=10.0, status="success",
        )
        entry = last_log_line(tmp_path / "tool_execution.log")
        assert "turn_id" not in entry["payload"]
        assert "call_index" not in entry["payload"]


# ---------------------------------------------------------------------------
# Change 6 — Rename result_preview → result
# ---------------------------------------------------------------------------

class TestResultColumnRename:
    def test_result_key_present(self, tmp_path):
        el = make_event_logger(tmp_path)
        el.log_tool_execution(
            user_id="u1", session_id="s1", tool_name="tool_a",
            args={}, result="output", duration_ms=10.0, status="success",
        )
        entry = last_log_line(tmp_path / "tool_execution.log")
        assert "result" in entry["payload"]

    def test_result_preview_key_absent(self, tmp_path):
        el = make_event_logger(tmp_path)
        el.log_tool_execution(
            user_id="u1", session_id="s1", tool_name="tool_a",
            args={}, result="output", duration_ms=10.0, status="success",
        )
        entry = last_log_line(tmp_path / "tool_execution.log")
        assert "result_preview" not in entry["payload"]
