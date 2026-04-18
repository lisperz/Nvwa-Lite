"""Unit tests for `src/logging/uns_snapshot.py` (T-035.b).

Covers the two pure functions:
- `snapshot_uns(adata, tool_name, prior_keys) -> dict`
- `emit_uns_snapshot(logger, user_id, session_id, tool_name, adata, prior_keys) -> None`

Run with: pytest tests/unit/test_uns_snapshot.py -v
"""

from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import MagicMock

import numpy as np
import pandas as pd

from src.logging.uns_snapshot import emit_uns_snapshot, snapshot_uns


def _fake_adata(uns: dict) -> SimpleNamespace:
    """Minimal stand-in for AnnData — we only touch `.uns`."""
    return SimpleNamespace(uns=uns)


# ---------------------------------------------------------------------------
# T1 — shape of snapshot_uns return
# ---------------------------------------------------------------------------

class TestSnapshotShape:
    def test_returns_expected_keys(self):
        adata = _fake_adata({"rank_genes_groups": {"names": None}, "pairwise_de_result": {}})
        snap = snapshot_uns(adata, "compare_groups_de", prior_keys=["rank_genes_groups"])
        assert set(snap) == {"tool_name", "keys_after", "new_keys", "payload_preview"}
        assert snap["tool_name"] == "compare_groups_de"
        assert sorted(snap["keys_after"]) == ["pairwise_de_result", "rank_genes_groups"]
        assert snap["new_keys"] == ["pairwise_de_result"]

    def test_no_new_keys_when_prior_matches(self):
        adata = _fake_adata({"rank_genes_groups": {"names": None}})
        snap = snapshot_uns(adata, "marker_genes", prior_keys=["rank_genes_groups"])
        assert snap["new_keys"] == []

    def test_empty_prior_marks_all_as_new(self):
        adata = _fake_adata({"a": 1, "b": 2})
        snap = snapshot_uns(adata, "t", prior_keys=[])
        assert sorted(snap["new_keys"]) == ["a", "b"]


# ---------------------------------------------------------------------------
# T2 — dict values show one level of sub-keys, not recursive
# ---------------------------------------------------------------------------

class TestDictPreview:
    def test_dict_preview_shows_one_level(self):
        adata = _fake_adata(
            {"pairwise_de_result": {"group1": "A", "group2": "B", "results_df": pd.DataFrame({"x": [1]})}}
        )
        snap = snapshot_uns(adata, "compare_groups_de", prior_keys=[])
        preview = snap["payload_preview"]["pairwise_de_result"]
        assert preview["type"] == "dict"
        assert sorted(preview["subkeys"]) == ["group1", "group2", "results_df"]
        # No recursive descent — results_df is not expanded
        assert "results_df" not in preview or isinstance(preview.get("results_df"), str) is False
        assert "type" in preview and "subkeys" in preview
        assert len(preview) == 2  # type + subkeys only

    def test_nested_dict_not_recursed(self):
        adata = _fake_adata({"outer": {"inner": {"deep": 1}}})
        snap = snapshot_uns(adata, "t", prior_keys=[])
        preview = snap["payload_preview"]["outer"]
        assert preview["subkeys"] == ["inner"]
        # inner dict contents must not leak
        assert "deep" not in str(preview)


# ---------------------------------------------------------------------------
# T3 — DataFrame/ndarray show type + shape, not contents
# ---------------------------------------------------------------------------

class TestArrayLikePreview:
    def test_dataframe_shows_type_and_shape(self):
        df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
        adata = _fake_adata({"results": df})
        snap = snapshot_uns(adata, "t", prior_keys=[])
        preview = snap["payload_preview"]["results"]
        assert preview["type"] == "DataFrame"
        assert preview["shape"] == (3, 2)
        # Contents must not leak — preview has only type + shape
        assert set(preview.keys()) == {"type", "shape"}

    def test_ndarray_shows_type_and_shape(self):
        arr = np.zeros((4, 5))
        adata = _fake_adata({"matrix": arr})
        snap = snapshot_uns(adata, "t", prior_keys=[])
        preview = snap["payload_preview"]["matrix"]
        assert preview["type"] == "ndarray"
        assert preview["shape"] == (4, 5)

    def test_scalar_shows_type_only(self):
        adata = _fake_adata({"n": 42, "label": "hello"})
        snap = snapshot_uns(adata, "t", prior_keys=[])
        assert snap["payload_preview"]["n"] == {"type": "int"}
        assert snap["payload_preview"]["label"] == {"type": "str"}


# ---------------------------------------------------------------------------
# T4 — emit_uns_snapshot calls log_session_event correctly
# ---------------------------------------------------------------------------

class TestEmit:
    def test_calls_log_session_event_with_uns_snapshot_type(self):
        logger = MagicMock()
        adata = _fake_adata({"rank_genes_groups": {"names": None}, "pairwise_de_result": {"group1": "A"}})
        emit_uns_snapshot(
            logger=logger,
            user_id="u1",
            session_id="s1",
            tool_name="compare_groups_de",
            adata=adata,
            prior_keys=["rank_genes_groups"],
        )
        logger.log_session_event.assert_called_once()
        kwargs = logger.log_session_event.call_args.kwargs
        assert kwargs["user_id"] == "u1"
        assert kwargs["session_id"] == "s1"
        assert kwargs["event_type"] == "uns_snapshot"
        md = kwargs["metadata"]
        assert md["tool_name"] == "compare_groups_de"
        assert md["new_keys"] == ["pairwise_de_result"]
        assert "rank_genes_groups" in md["keys_after"]
        assert md["payload_preview"]["pairwise_de_result"]["type"] == "dict"

    def test_no_raise_when_logger_is_none(self):
        """Guard: caller may pass None when binding is absent."""
        adata = _fake_adata({"a": 1})
        # Should simply no-op; no exception.
        emit_uns_snapshot(
            logger=None,
            user_id="u1",
            session_id="s1",
            tool_name="t",
            adata=adata,
            prior_keys=[],
        )
