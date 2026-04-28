"""Unit tests for src/core/results.py — type construction and invariants.

Covers TextResult, ArtifactResult, and ToolExecutionError. No LLM, no adata —
these are pure dataclass + exception tests.
"""

from __future__ import annotations

import pytest

from src.core.results import (
    ArtifactResult,
    TextResult,
    ToolExecutionError,
)


# ---------------------------------------------------------------------------
# TextResult
# ---------------------------------------------------------------------------


class TestTextResult:
    def test_default_status_success(self):
        r = TextResult(text="hello")
        assert r.status == "success"
        assert r.error_message is None

    def test_explicit_error_status(self):
        r = TextResult(text="oops", status="error", error_message="boom")
        assert r.status == "error"
        assert r.error_message == "boom"

    def test_warning_status(self):
        r = TextResult(text="careful", status="warning")
        assert r.status == "warning"

    def test_frozen_cannot_mutate(self):
        r = TextResult(text="immutable")
        with pytest.raises(Exception):
            r.text = "changed"  # type: ignore[misc]


# ---------------------------------------------------------------------------
# ArtifactResult
# ---------------------------------------------------------------------------


class TestArtifactResultImage:
    def _build(self, **kwargs) -> ArtifactResult:
        defaults = dict(
            text="plot ready",
            artifact_kind="image",
            tool_name="some_plot",
            params_used={"x": 1},
            image_bytes=b"\x89PNG-fake",
        )
        defaults.update(kwargs)
        return ArtifactResult(**defaults)  # type: ignore[arg-type]

    def test_image_bytes_set(self):
        r = self._build()
        assert r.image_bytes == b"\x89PNG-fake"
        assert r.artifact_kind == "image"

    def test_csv_fields_default_none(self):
        r = self._build()
        assert r.csv_data is None
        assert r.display_df is None

    def test_default_entities_acted_on_is_empty_list(self):
        r = self._build()
        assert r.entities_acted_on == []

    def test_default_status_success(self):
        r = self._build()
        assert r.status == "success"
        assert r.error_message is None

    def test_code_optional(self):
        r = self._build(code="sc.pl.umap(adata)")
        assert r.code == "sc.pl.umap(adata)"


class TestArtifactResultCsv:
    def _build(self, **kwargs) -> ArtifactResult:
        defaults = dict(
            text="table ready",
            artifact_kind="csv",
            tool_name="some_table",
            params_used={"threshold": 0.05},
            csv_data="a,b,c\n1,2,3",
        )
        defaults.update(kwargs)
        return ArtifactResult(**defaults)  # type: ignore[arg-type]

    def test_csv_data_set(self):
        r = self._build()
        assert r.csv_data == "a,b,c\n1,2,3"
        assert r.artifact_kind == "csv"

    def test_image_bytes_default_none(self):
        r = self._build()
        assert r.image_bytes is None

    def test_display_df_optional(self):
        r = self._build(display_df="| a | b |\n|---|---|\n| 1 | 2 |")
        assert r.display_df == "| a | b |\n|---|---|\n| 1 | 2 |"

    def test_entities_acted_on_carries_through(self):
        r = self._build(entities_acted_on=["pct_counts_mt", "n_genes_by_counts"])
        assert r.entities_acted_on == ["pct_counts_mt", "n_genes_by_counts"]


class TestArtifactResultRequiredFields:
    def test_text_required(self):
        with pytest.raises(TypeError):
            ArtifactResult(  # type: ignore[call-arg]
                artifact_kind="image",
                tool_name="x",
                params_used={},
            )

    def test_artifact_kind_required(self):
        with pytest.raises(TypeError):
            ArtifactResult(  # type: ignore[call-arg]
                text="x",
                tool_name="x",
                params_used={},
            )

    def test_tool_name_required(self):
        with pytest.raises(TypeError):
            ArtifactResult(  # type: ignore[call-arg]
                text="x",
                artifact_kind="image",
                params_used={},
            )

    def test_params_used_required(self):
        with pytest.raises(TypeError):
            ArtifactResult(  # type: ignore[call-arg]
                text="x",
                artifact_kind="image",
                tool_name="x",
            )

    def test_frozen_cannot_mutate(self):
        r = ArtifactResult(
            text="x", artifact_kind="image", tool_name="x", params_used={},
        )
        with pytest.raises(Exception):
            r.text = "changed"  # type: ignore[misc]


# ---------------------------------------------------------------------------
# ToolExecutionError
# ---------------------------------------------------------------------------


class TestToolExecutionError:
    def test_message_carried(self):
        e = ToolExecutionError("disk full")
        assert str(e) == "disk full"

    def test_tool_name_default_none(self):
        e = ToolExecutionError("disk full")
        assert e.tool_name is None

    def test_tool_name_set(self):
        e = ToolExecutionError("disk full", tool_name="qc_summary_table")
        assert e.tool_name == "qc_summary_table"

    def test_is_exception(self):
        with pytest.raises(ToolExecutionError):
            raise ToolExecutionError("boom")
