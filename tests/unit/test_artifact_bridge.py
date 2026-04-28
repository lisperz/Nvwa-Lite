"""Unit tests for src/agent/artifacts.consume_artifact_result — the bridge.

Verifies ArtifactResult → PlotResult / TableResult translation populates the
new artifacts singleton correctly so UI consumption (which aggregates this
plus the legacy src.agent.tools singletons) sees the right content.

No LLM, no adata — these are pure bridge-translation tests.
"""

from __future__ import annotations

import pytest

from src.agent import artifacts
from src.core.results import ArtifactResult


@pytest.fixture(autouse=True)
def clear_buffer():
    """Isolation: clear new-path buffer before and after each test."""
    artifacts.clear()
    yield
    artifacts.clear()


# ---------------------------------------------------------------------------
# Image bridge
# ---------------------------------------------------------------------------


class TestImageBridge:
    def _build(self, **kwargs) -> ArtifactResult:
        defaults = dict(
            text="plot generated",
            artifact_kind="image",
            tool_name="some_plot",
            params_used={"color": "leiden"},
            image_bytes=b"\x89PNG-fake-bytes",
            code="sc.pl.umap(adata, color='leiden')",
        )
        defaults.update(kwargs)
        return ArtifactResult(**defaults)  # type: ignore[arg-type]

    def test_image_pushes_to_plot_buffer(self):
        artifacts.consume_artifact_result(self._build())
        assert len(artifacts.get_plot_results()) == 1
        assert len(artifacts.get_table_results()) == 0

    def test_image_preserves_bytes(self):
        artifacts.consume_artifact_result(self._build(image_bytes=b"specific-bytes"))
        plot = artifacts.get_plot_results()[0]
        assert plot.image == b"specific-bytes"

    def test_image_preserves_code(self):
        artifacts.consume_artifact_result(
            self._build(code="my_custom_code()"),
        )
        plot = artifacts.get_plot_results()[0]
        assert plot.code == "my_custom_code()"

    def test_image_text_becomes_message(self):
        artifacts.consume_artifact_result(
            self._build(text="UMAP rendered for 3 conditions."),
        )
        plot = artifacts.get_plot_results()[0]
        assert plot.message == "UMAP rendered for 3 conditions."

    def test_image_missing_bytes_falls_back_to_empty(self):
        """ArtifactResult with image_bytes=None → bridge writes b'' (gatekeeper artifact_empty would catch)."""
        r = ArtifactResult(
            text="x", artifact_kind="image", tool_name="x", params_used={},
            image_bytes=None,
        )
        artifacts.consume_artifact_result(r)
        plot = artifacts.get_plot_results()[0]
        assert plot.image == b""

    def test_image_missing_code_falls_back_to_empty_string(self):
        r = ArtifactResult(
            text="x", artifact_kind="image", tool_name="x", params_used={},
            image_bytes=b"data", code=None,
        )
        artifacts.consume_artifact_result(r)
        plot = artifacts.get_plot_results()[0]
        assert plot.code == ""

    def test_multiple_image_pushes_accumulate(self):
        artifacts.consume_artifact_result(self._build(text="first"))
        artifacts.consume_artifact_result(self._build(text="second"))
        plots = artifacts.get_plot_results()
        assert len(plots) == 2
        assert plots[0].message == "first"
        assert plots[1].message == "second"


# ---------------------------------------------------------------------------
# CSV bridge
# ---------------------------------------------------------------------------


class TestCsvBridge:
    def _build(self, **kwargs) -> ArtifactResult:
        defaults = dict(
            text="table ready",
            artifact_kind="csv",
            tool_name="qc_summary_table",
            params_used={},
            csv_data="metric,median\npct_counts_mt,5.2",
            display_df="| metric | median |\n|---|---|\n| pct_counts_mt | 5.2 |",
            code="qc_summary_table(adata)",
        )
        defaults.update(kwargs)
        return ArtifactResult(**defaults)  # type: ignore[arg-type]

    def test_csv_pushes_to_table_buffer(self):
        artifacts.consume_artifact_result(self._build())
        assert len(artifacts.get_table_results()) == 1
        assert len(artifacts.get_plot_results()) == 0

    def test_csv_preserves_csv_data(self):
        artifacts.consume_artifact_result(self._build(csv_data="a,b\n1,2"))
        table = artifacts.get_table_results()[0]
        assert table.csv_data == "a,b\n1,2"

    def test_csv_preserves_display_df(self):
        artifacts.consume_artifact_result(
            self._build(display_df="| custom | markdown |\n|---|---|\n| x | y |"),
        )
        table = artifacts.get_table_results()[0]
        assert table.display_df == "| custom | markdown |\n|---|---|\n| x | y |"

    def test_csv_text_becomes_message(self):
        artifacts.consume_artifact_result(
            self._build(text="3 metrics flagged."),
        )
        table = artifacts.get_table_results()[0]
        assert table.message == "3 metrics flagged."

    def test_csv_preserves_code(self):
        artifacts.consume_artifact_result(self._build(code="custom_call()"))
        table = artifacts.get_table_results()[0]
        assert table.code == "custom_call()"

    def test_csv_missing_csv_data_falls_back_to_empty(self):
        r = ArtifactResult(
            text="x", artifact_kind="csv", tool_name="x", params_used={},
            csv_data=None,
        )
        artifacts.consume_artifact_result(r)
        table = artifacts.get_table_results()[0]
        assert table.csv_data == ""

    def test_csv_missing_display_df_falls_back_to_empty_string(self):
        r = ArtifactResult(
            text="x", artifact_kind="csv", tool_name="x", params_used={},
            csv_data="a,b\n1,2", display_df=None,
        )
        artifacts.consume_artifact_result(r)
        table = artifacts.get_table_results()[0]
        assert table.display_df == ""


# ---------------------------------------------------------------------------
# Buffer isolation
# ---------------------------------------------------------------------------


class TestBufferLifecycle:
    def test_clear_resets_both_buffers(self):
        artifacts.consume_artifact_result(ArtifactResult(
            text="img", artifact_kind="image", tool_name="x", params_used={},
            image_bytes=b"png",
        ))
        artifacts.consume_artifact_result(ArtifactResult(
            text="csv", artifact_kind="csv", tool_name="x", params_used={},
            csv_data="a,b",
        ))
        assert len(artifacts.get_plot_results()) == 1
        assert len(artifacts.get_table_results()) == 1

        artifacts.clear()
        assert len(artifacts.get_plot_results()) == 0
        assert len(artifacts.get_table_results()) == 0

    def test_get_plot_results_returns_copy(self):
        """Mutating the returned list must not mutate the buffer."""
        artifacts.consume_artifact_result(ArtifactResult(
            text="x", artifact_kind="image", tool_name="x", params_used={},
            image_bytes=b"data",
        ))
        plots = artifacts.get_plot_results()
        plots.clear()
        # Buffer should still have the entry
        assert len(artifacts.get_plot_results()) == 1

    def test_get_table_results_returns_copy(self):
        artifacts.consume_artifact_result(ArtifactResult(
            text="x", artifact_kind="csv", tool_name="x", params_used={},
            csv_data="a,b",
        ))
        tables = artifacts.get_table_results()
        tables.clear()
        assert len(artifacts.get_table_results()) == 1
