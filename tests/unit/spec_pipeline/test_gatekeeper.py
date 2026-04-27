"""Exhaustive-branch unit tests for src/agent/gatekeeper.py.

Covers every branch of the public ``check()`` function:
- all 5 block reasons: tool_mismatch, artifact_claim_false, artifact_empty,
  error_masked, entity_missing
- first-block-wins priority order
- retry_hint populated vs None (only artifact_empty returns None)
- pass path when all checks clear
- error-indicator keyword list behaviour
- artifact-claim phrase matching
- entity_missing only for field_type params

Registry isolation: ``_clean_registry`` autouse fixture snapshots and restores
REGISTRY so no test leaks state.

No adata needed — gatekeeper is a pure function over Spec + strings + result lists.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

import pytest

import src.tools.registry as reg_module
from src.agent.gatekeeper import (
    GatekeeperResult,
    _ARTIFACT_CLAIM_PHRASES,
    _ERROR_INDICATORS,
    check,
)
from src.core.spec import Spec
from src.tools.registry import register


# ---------------------------------------------------------------------------
# Registry isolation
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _clean_registry():
    """Save and restore REGISTRY around every test."""
    snapshot = dict(reg_module.REGISTRY)
    yield
    reg_module.REGISTRY.clear()
    reg_module.REGISTRY.update(snapshot)


# ---------------------------------------------------------------------------
# Minimal result-shaped objects (no MagicMock)
# ---------------------------------------------------------------------------


@dataclass
class PlotResult:
    """Minimal plot artifact — mirrors the attributes gatekeeper walks."""
    img_bytes: Optional[bytes] = None


@dataclass
class PlotResultAlt:
    """Same shape but uses the 'image_bytes' fallback attribute name."""
    image_bytes: Optional[bytes] = None


@dataclass
class PlotResultPng:
    """Uses the 'png_bytes' fallback attribute name."""
    png_bytes: Optional[bytes] = None


@dataclass
class PlotResultBytes:
    """Uses the 'bytes' fallback attribute name."""
    bytes: Optional[bytes] = None  # noqa: A003


class _FakeDF:
    """Minimal dataframe-like with a .shape attribute."""

    def __init__(self, n_rows: int) -> None:
        self.shape = (n_rows, 3)


@dataclass
class TableResult:
    """Minimal table artifact."""
    dataframe: Optional[Any] = None


@dataclass
class TableResultDF:
    """Uses 'df' fallback attribute."""
    df: Optional[Any] = None


@dataclass
class TableResultData:
    """Uses 'data' fallback attribute."""
    data: Optional[Any] = None


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_spec(tool_name: str = "my_tool", params: Optional[dict] = None) -> Spec:
    return Spec(
        scenario_id="test-scenario",
        tool_name=tool_name,
        params=params or {},
    )


def _register_tool(name: str, field_params: Optional[dict[str, str]] = None) -> None:
    """Register a throwaway tool with the given name via @register.

    ``field_params`` maps param-name -> field_type string.
    Non-field_type params can also be included as a plain function signature.
    """
    field_params = field_params or {}

    # Build a function whose signature matches the requested params.
    # We create it dynamically via exec so the __name__ matches ``name``.
    param_names = list(field_params.keys())

    if param_names:
        sig_str = ", ".join(param_names)
        func_src = f"def {name}({sig_str}): pass"
    else:
        func_src = f"def {name}(): pass"

    ns: dict = {}
    exec(func_src, ns)  # noqa: S102
    fn = ns[name]

    extras = {p: {"field_type": ft} for p, ft in field_params.items()}
    register(description=f"test tool {name}", params=extras)(fn)


def _register_tool_with_mixed_params(
    name: str,
    field_params: dict[str, str],
    mechanical_params: list[str],
) -> None:
    """Register a tool with both field_type params and mechanical (no field_type) params."""
    param_names = list(field_params.keys()) + mechanical_params
    sig_str = ", ".join(param_names)
    func_src = f"def {name}({sig_str}): pass"
    ns: dict = {}
    exec(func_src, ns)  # noqa: S102
    fn = ns[name]

    extras = {p: {"field_type": ft} for p, ft in field_params.items()}
    # mechanical_params get no extras — field_type stays None
    register(description=f"test tool {name}", params=extras)(fn)


# ---------------------------------------------------------------------------
# Pass path
# ---------------------------------------------------------------------------


class TestPassPath:
    def test_all_clear_returns_pass_status(self):
        spec = _make_spec()
        result = check(spec, "tool ran OK", "The analysis completed.", [], [])
        assert result.status == "pass"

    def test_all_clear_returns_pass_reason(self):
        spec = _make_spec()
        result = check(spec, "tool ran OK", "The analysis completed.", [], [])
        assert result.reason == "pass"

    def test_all_clear_retry_hint_is_none(self):
        spec = _make_spec()
        result = check(spec, "tool ran OK", "The analysis completed.", [], [])
        assert result.retry_hint is None

    def test_all_clear_returns_gatekeeper_result(self):
        spec = _make_spec()
        result = check(spec, "tool ran OK", "The analysis completed.", [], [])
        assert isinstance(result, GatekeeperResult)

    def test_pass_explanation_nonempty(self):
        spec = _make_spec()
        result = check(spec, "tool ran OK", "The analysis completed.", [], [])
        assert result.explanation


# ---------------------------------------------------------------------------
# tool_mismatch
# ---------------------------------------------------------------------------


class TestToolMismatch:
    def test_fires_when_other_tool_named_spec_not_named(self):
        _register_tool("spec_tool")
        _register_tool("other_tool")
        spec = _make_spec("spec_tool")
        # Responder mentions other_tool but NOT spec_tool
        result = check(spec, "ok", "The other_tool was used here.", [], [])
        assert result.status == "block"
        assert result.reason == "tool_mismatch"

    def test_retry_hint_populated_for_tool_mismatch(self):
        _register_tool("spec_tool")
        _register_tool("other_tool")
        spec = _make_spec("spec_tool")
        result = check(spec, "ok", "The other_tool ran.", [], [])
        assert result.retry_hint is not None
        assert len(result.retry_hint) > 0

    def test_suppressed_when_spec_tool_also_named(self):
        """If responder names both spec tool AND another tool, no mismatch."""
        _register_tool("spec_tool")
        _register_tool("other_tool")
        spec = _make_spec("spec_tool")
        result = check(spec, "ok", "The other_tool and spec_tool were used.", [], [])
        # tool_mismatch should NOT fire (spec_mentioned is True)
        assert result.reason != "tool_mismatch"

    def test_suppressed_when_only_spec_tool_named(self):
        _register_tool("spec_tool")
        spec = _make_spec("spec_tool")
        result = check(spec, "ok", "The spec_tool ran successfully.", [], [])
        assert result.reason != "tool_mismatch"

    def test_suppressed_when_no_tools_named_at_all(self):
        _register_tool("spec_tool")
        _register_tool("other_tool")
        spec = _make_spec("spec_tool")
        result = check(spec, "ok", "The analysis completed without issues.", [], [])
        assert result.reason != "tool_mismatch"

    def test_case_insensitive_matching(self):
        _register_tool("spec_tool")
        _register_tool("other_tool")
        spec = _make_spec("spec_tool")
        # other_tool in uppercase — should still match because we lowercase both
        result = check(spec, "ok", "The OTHER_TOOL was called.", [], [])
        assert result.status == "block"
        assert result.reason == "tool_mismatch"

    def test_tool_not_in_registry_no_mismatch(self):
        """If spec.tool_name has no registration, get_tool_names may not include it;
        no others_mentioned -> no mismatch even if text is arbitrary."""
        reg_module.REGISTRY.clear()
        spec = _make_spec("nonexistent_tool")
        result = check(spec, "ok", "Nothing relevant here.", [], [])
        assert result.reason != "tool_mismatch"

    def test_explanation_contains_spec_tool_name(self):
        _register_tool("spec_tool")
        _register_tool("other_tool")
        spec = _make_spec("spec_tool")
        result = check(spec, "ok", "other_tool was used.", [], [])
        assert "spec_tool" in result.explanation

    def test_multiple_other_tools_mentioned(self):
        _register_tool("spec_tool")
        _register_tool("other_tool_a")
        _register_tool("other_tool_b")
        spec = _make_spec("spec_tool")
        result = check(spec, "ok", "other_tool_a and other_tool_b were invoked.", [], [])
        assert result.status == "block"
        assert result.reason == "tool_mismatch"


# ---------------------------------------------------------------------------
# artifact_claim_false
# ---------------------------------------------------------------------------


class TestArtifactClaimFalse:
    @pytest.mark.parametrize("phrase", _ARTIFACT_CLAIM_PHRASES)
    def test_fires_for_each_claim_phrase(self, phrase: str):
        spec = _make_spec()
        result = check(spec, "ok", f"The chart {phrase} for you.", [], [])
        assert result.status == "block"
        assert result.reason == "artifact_claim_false"

    def test_retry_hint_populated(self):
        spec = _make_spec()
        result = check(spec, "ok", "has been generated.", [], [])
        assert result.retry_hint is not None

    def test_suppressed_when_plot_results_present(self):
        spec = _make_spec()
        plot = PlotResult(img_bytes=b"data")
        result = check(spec, "ok", "has been generated.", [plot], [])
        assert result.reason != "artifact_claim_false"

    def test_suppressed_when_table_results_present(self):
        spec = _make_spec()
        table = TableResult(dataframe=_FakeDF(5))
        result = check(spec, "ok", "has been generated.", [], [table])
        assert result.reason != "artifact_claim_false"

    def test_suppressed_when_both_artifacts_present(self):
        spec = _make_spec()
        plot = PlotResult(img_bytes=b"data")
        table = TableResult(dataframe=_FakeDF(5))
        result = check(spec, "ok", "has been generated.", [plot], [table])
        assert result.reason != "artifact_claim_false"

    def test_not_fired_when_no_claim_phrase(self):
        spec = _make_spec()
        result = check(spec, "ok", "The analysis ran. No plot was produced.", [], [])
        assert result.reason != "artifact_claim_false"

    def test_case_insensitive_phrase_matching(self):
        """Phrases must match case-insensitively since responder text is lowered."""
        spec = _make_spec()
        # Uppercase original should still match after lowering
        result = check(spec, "ok", "THE CHART HAS BEEN GENERATED.", [], [])
        assert result.status == "block"
        assert result.reason == "artifact_claim_false"

    def test_explanation_contains_matched_phrase(self):
        spec = _make_spec()
        result = check(spec, "ok", "has been plotted for you.", [], [])
        assert "has been plotted" in result.explanation


# ---------------------------------------------------------------------------
# artifact_empty
# ---------------------------------------------------------------------------


class TestArtifactEmpty:
    def test_fires_for_zero_byte_plot_img_bytes(self):
        spec = _make_spec()
        plot = PlotResult(img_bytes=b"")
        result = check(spec, "ok", "All good.", [plot], [])
        assert result.status == "block"
        assert result.reason == "artifact_empty"

    def test_fires_for_zero_byte_plot_image_bytes_attr(self):
        spec = _make_spec()
        plot = PlotResultAlt(image_bytes=b"")
        result = check(spec, "ok", "All good.", [plot], [])
        assert result.status == "block"
        assert result.reason == "artifact_empty"

    def test_fires_for_zero_byte_plot_png_bytes_attr(self):
        spec = _make_spec()
        plot = PlotResultPng(png_bytes=b"")
        result = check(spec, "ok", "All good.", [plot], [])
        assert result.status == "block"
        assert result.reason == "artifact_empty"

    def test_fires_for_zero_byte_plot_bytes_attr(self):
        spec = _make_spec()
        plot = PlotResultBytes(bytes=b"")
        result = check(spec, "ok", "All good.", [plot], [])
        assert result.status == "block"
        assert result.reason == "artifact_empty"

    def test_fires_for_zero_row_table_dataframe_attr(self):
        spec = _make_spec()
        table = TableResult(dataframe=_FakeDF(0))
        result = check(spec, "ok", "All good.", [], [table])
        assert result.status == "block"
        assert result.reason == "artifact_empty"

    def test_fires_for_zero_row_table_df_attr(self):
        spec = _make_spec()
        table = TableResultDF(df=_FakeDF(0))
        result = check(spec, "ok", "All good.", [], [table])
        assert result.status == "block"
        assert result.reason == "artifact_empty"

    def test_fires_for_zero_row_table_data_attr(self):
        spec = _make_spec()
        table = TableResultData(data=_FakeDF(0))
        result = check(spec, "ok", "All good.", [], [table])
        assert result.status == "block"
        assert result.reason == "artifact_empty"

    def test_retry_hint_is_none_for_empty_plot(self):
        spec = _make_spec()
        plot = PlotResult(img_bytes=b"")
        result = check(spec, "ok", "All good.", [plot], [])
        assert result.retry_hint is None

    def test_retry_hint_is_none_for_empty_table(self):
        spec = _make_spec()
        table = TableResult(dataframe=_FakeDF(0))
        result = check(spec, "ok", "All good.", [], [table])
        assert result.retry_hint is None

    def test_nonempty_plot_does_not_fire(self):
        spec = _make_spec()
        plot = PlotResult(img_bytes=b"PNG_DATA_HERE")
        result = check(spec, "ok", "All good.", [plot], [])
        assert result.reason != "artifact_empty"

    def test_nonempty_table_does_not_fire(self):
        spec = _make_spec()
        table = TableResult(dataframe=_FakeDF(10))
        result = check(spec, "ok", "All good.", [], [table])
        assert result.reason != "artifact_empty"

    def test_plot_with_no_known_attr_skips_check(self):
        """A plot object with none of the known image-byte attrs must not block."""
        class UnknownPlot:
            other_attr = b""

        spec = _make_spec()
        result = check(spec, "ok", "All good.", [UnknownPlot()], [])
        assert result.reason != "artifact_empty"

    def test_table_with_no_known_attr_skips_check(self):
        """A table object with none of the known df attrs must not block."""
        class UnknownTable:
            other_attr = None

        spec = _make_spec()
        result = check(spec, "ok", "All good.", [], [UnknownTable()])
        assert result.reason != "artifact_empty"

    def test_table_shape_raises_attribute_error_skips(self):
        """If shape[0] raises, the check must be skipped (no false block)."""
        class BadDF:
            @property
            def shape(self):
                raise AttributeError("no shape")

        spec = _make_spec()
        table = TableResult(dataframe=BadDF())
        result = check(spec, "ok", "All good.", [], [table])
        assert result.reason != "artifact_empty"

    def test_empty_plot_list_does_not_fire(self):
        spec = _make_spec()
        result = check(spec, "ok", "All good.", [], [])
        assert result.reason != "artifact_empty"


# ---------------------------------------------------------------------------
# error_masked
# ---------------------------------------------------------------------------


class TestErrorMasked:
    def test_fires_when_error_prefix_and_no_indicator(self):
        spec = _make_spec()
        result = check(spec, "Error: something went wrong", "All looks great!", [], [])
        assert result.status == "block"
        assert result.reason == "error_masked"

    def test_retry_hint_populated(self):
        spec = _make_spec()
        result = check(spec, "Error: something went wrong", "All looks great!", [], [])
        assert result.retry_hint is not None

    def test_suppressed_when_no_error_prefix(self):
        spec = _make_spec()
        result = check(spec, "Success: analysis done", "All good.", [], [])
        assert result.reason != "error_masked"

    def test_suppressed_when_tool_output_empty_string(self):
        spec = _make_spec()
        result = check(spec, "", "All good.", [], [])
        assert result.reason != "error_masked"

    def test_leading_whitespace_before_error_prefix(self):
        """lstrip() means leading whitespace before 'Error:' is ignored."""
        spec = _make_spec()
        result = check(spec, "   Error: oops", "All looks great!", [], [])
        assert result.status == "block"
        assert result.reason == "error_masked"

    @pytest.mark.parametrize("keyword", _ERROR_INDICATORS)
    def test_suppressed_when_each_error_keyword_present(self, keyword: str):
        spec = _make_spec()
        responder = f"Unfortunately there was an {keyword} with the analysis."
        result = check(spec, "Error: bad", responder, [], [])
        assert result.reason != "error_masked", (
            f"error_masked fired even though '{keyword}' appeared in responder"
        )

    def test_case_insensitive_indicator_matching(self):
        """Keyword matching is on lowered_responder, so uppercase in original is fine."""
        spec = _make_spec()
        result = check(spec, "Error: bad", "An ERROR occurred.", [], [])
        assert result.reason != "error_masked"

    def test_explanation_mentions_error_colon(self):
        spec = _make_spec()
        result = check(spec, "Error: bad", "Looks fine!", [], [])
        assert "Error:" in result.explanation

    def test_fires_only_for_error_colon_not_other_error_words(self):
        """tool_output containing the word 'error' but not starting with 'Error:' must not block."""
        spec = _make_spec()
        result = check(
            spec,
            "The tool encountered an internal error during processing",
            "Looks fine!",
            [],
            [],
        )
        assert result.reason != "error_masked"


# ---------------------------------------------------------------------------
# entity_missing
# ---------------------------------------------------------------------------


class TestEntityMissing:
    """Covers the no-artifact branch of `_check_entity_missing`.

    Every case here passes empty `plot_results` and `table_results` lists, so
    `tool_output` is the only signal the check has. That is the correct
    branch to test today — the only @register tool (`dataset_overview`) has
    zero field_type params, so the live system never trips this loop.

    When the first @register plot/table tool lands (T-040), entities-acted-on
    will move to a structured field on PlotResult/TableResult and the check
    will read from there when artifacts are present. The cases below stay
    valid as the no-artifact-fallback branch; new cases for the
    artifacts-present branch will be added alongside the contract change.
    """

    def test_fires_when_field_type_param_absent_from_output(self):
        _register_tool("plot_gene", {"gene": "gene"})
        spec = _make_spec("plot_gene", {"gene": "BRCA1"})
        # BRCA1 is not in tool_output
        result = check(spec, "Plotting complete.", "The plot is ready.", [], [])
        assert result.status == "block"
        assert result.reason == "entity_missing"

    def test_retry_hint_populated(self):
        _register_tool("plot_gene", {"gene": "gene"})
        spec = _make_spec("plot_gene", {"gene": "BRCA1"})
        result = check(spec, "Plotting complete.", "The plot is ready.", [], [])
        assert result.retry_hint is not None

    def test_suppressed_when_entity_present_in_output(self):
        _register_tool("plot_gene", {"gene": "gene"})
        spec = _make_spec("plot_gene", {"gene": "BRCA1"})
        result = check(spec, "Plotting BRCA1 complete.", "The plot is ready.", [], [])
        assert result.reason != "entity_missing"

    def test_case_insensitive_entity_matching(self):
        _register_tool("plot_gene", {"gene": "gene"})
        spec = _make_spec("plot_gene", {"gene": "BRCA1"})
        # tool_output has lowercase version
        result = check(spec, "plotting brca1 complete.", "The plot is ready.", [], [])
        assert result.reason != "entity_missing"

    def test_not_fired_for_non_field_type_param(self):
        """Mechanical params (title, ncols, show_legend) must never trigger entity_missing."""
        _register_tool_with_mixed_params(
            "plot_umap",
            field_params={},  # no field_type params
            mechanical_params=["title", "ncols"],
        )
        spec = _make_spec("plot_umap", {"title": "My Title", "ncols": 2})
        result = check(spec, "UMAP plotted.", "The plot is ready.", [], [])
        assert result.reason != "entity_missing"

    def test_not_fired_when_tool_not_in_registry(self):
        """If get_tool returns None, entity_missing must be skipped."""
        reg_module.REGISTRY.clear()
        spec = _make_spec("nonexistent_tool", {"gene": "TP53"})
        result = check(spec, "output here", "All good.", [], [])
        assert result.reason != "entity_missing"

    def test_list_param_all_present(self):
        _register_tool("multi_gene_plot", {"genes": "gene"})
        spec = _make_spec("multi_gene_plot", {"genes": ["TP53", "MYC"]})
        result = check(spec, "Plotting TP53 and MYC.", "Ready.", [], [])
        assert result.reason != "entity_missing"

    def test_list_param_some_missing(self):
        _register_tool("multi_gene_plot", {"genes": "gene"})
        spec = _make_spec("multi_gene_plot", {"genes": ["TP53", "MYC"]})
        # Only TP53 present, MYC absent
        result = check(spec, "Plotting TP53.", "Ready.", [], [])
        assert result.status == "block"
        assert result.reason == "entity_missing"

    def test_list_param_all_missing(self):
        _register_tool("multi_gene_plot", {"genes": "gene"})
        spec = _make_spec("multi_gene_plot", {"genes": ["TP53", "MYC"]})
        result = check(spec, "No genes found.", "Ready.", [], [])
        assert result.status == "block"
        assert result.reason == "entity_missing"

    def test_param_value_none_skipped(self):
        """If spec.params[field_type_param] is None, skip that param (no block)."""
        _register_tool("plot_gene", {"gene": "gene"})
        spec = _make_spec("plot_gene", {"gene": None})
        result = check(spec, "Plotting complete.", "Ready.", [], [])
        assert result.reason != "entity_missing"

    def test_param_not_in_spec_params_skipped(self):
        """If the tool has a field_type param but spec.params doesn't include it, skip."""
        _register_tool("plot_gene", {"gene": "gene"})
        spec = _make_spec("plot_gene", {})  # no 'gene' key at all
        result = check(spec, "Plotting complete.", "Ready.", [], [])
        assert result.reason != "entity_missing"

    def test_explanation_contains_missing_entity(self):
        _register_tool("plot_gene", {"gene": "gene"})
        spec = _make_spec("plot_gene", {"gene": "EGFR"})
        result = check(spec, "Plotting done.", "Ready.", [], [])
        assert "EGFR" in result.explanation

    def test_mixed_field_and_mechanical_only_field_checked(self):
        """Only gene (field_type) triggers; title (no field_type) is ignored."""
        _register_tool_with_mixed_params(
            "fancy_plot",
            field_params={"gene": "gene"},
            mechanical_params=["title"],
        )
        spec = _make_spec("fancy_plot", {"gene": "PTEN", "title": "My Title"})
        # 'PTEN' absent from output; 'My Title' also absent — only PTEN should trigger
        result = check(spec, "Plot done without entities.", "Ready.", [], [])
        assert result.status == "block"
        assert result.reason == "entity_missing"
        assert "PTEN" in result.explanation
        # 'My Title' should NOT appear in the missing list
        assert "My Title" not in result.explanation

    def test_cell_type_field_type_fires(self):
        _register_tool("filter_cells", {"cell_type": "cell_type"})
        spec = _make_spec("filter_cells", {"cell_type": "T cell"})
        result = check(spec, "Filtered cells.", "Done.", [], [])
        assert result.status == "block"
        assert result.reason == "entity_missing"

    def test_cell_type_field_type_suppressed_when_present(self):
        _register_tool("filter_cells", {"cell_type": "cell_type"})
        spec = _make_spec("filter_cells", {"cell_type": "T cell"})
        result = check(spec, "Filtered T cell population.", "Done.", [], [])
        assert result.reason != "entity_missing"


# ---------------------------------------------------------------------------
# First-block-wins priority order
# ---------------------------------------------------------------------------


class TestFirstBlockWins:
    def test_tool_mismatch_beats_artifact_claim_false(self):
        """tool_mismatch (priority 1) must fire before artifact_claim_false (priority 2)."""
        _register_tool("spec_tool")
        _register_tool("other_tool")
        spec = _make_spec("spec_tool")
        # Trigger tool_mismatch: mention other_tool without spec_tool
        # Trigger artifact_claim_false: include a claim phrase with no artifacts
        responder = "The other_tool has been generated the plot."
        result = check(spec, "ok", responder, [], [])
        assert result.reason == "tool_mismatch"

    def test_tool_mismatch_beats_artifact_empty(self):
        _register_tool("spec_tool")
        _register_tool("other_tool")
        spec = _make_spec("spec_tool")
        plot = PlotResult(img_bytes=b"")  # triggers artifact_empty
        result = check(spec, "ok", "The other_tool ran.", [plot], [])
        assert result.reason == "tool_mismatch"

    def test_tool_mismatch_beats_error_masked(self):
        _register_tool("spec_tool")
        _register_tool("other_tool")
        spec = _make_spec("spec_tool")
        result = check(spec, "Error: failed", "The other_tool ran.", [], [])
        assert result.reason == "tool_mismatch"

    def test_tool_mismatch_beats_entity_missing(self):
        _register_tool_with_mixed_params(
            "spec_tool", field_params={"gene": "gene"}, mechanical_params=[]
        )
        _register_tool("other_tool")
        spec = _make_spec("spec_tool", {"gene": "BRCA1"})
        # other_tool in responder triggers tool_mismatch; BRCA1 absent would trigger entity_missing
        result = check(spec, "no entity here", "The other_tool did stuff.", [], [])
        assert result.reason == "tool_mismatch"

    def test_artifact_claim_false_beats_artifact_empty(self):
        """artifact_claim_false (priority 2) must fire before artifact_empty (priority 3).

        artifact_claim_false fires when no artifacts and a claim phrase exists.
        artifact_empty fires when an artifact IS present but empty.
        These are mutually exclusive by design: artifact_claim_false checks
        (plot_results or table_results) == empty, while artifact_empty checks
        the contents of non-empty lists. They cannot both fire simultaneously.

        We verify the conceptual priority by confirming artifact_claim_false
        fires with no artifacts + a claim phrase (before artifact_empty would
        even be reached), and returns the correct reason.
        """
        spec = _make_spec()
        result = check(spec, "ok", "has been generated.", [], [])
        assert result.reason == "artifact_claim_false"

    def test_artifact_empty_beats_error_masked(self):
        """artifact_empty (priority 3) must fire before error_masked (priority 4)."""
        spec = _make_spec()
        plot = PlotResult(img_bytes=b"")
        # Trigger error_masked: output starts with Error: and responder doesn't ack
        result = check(spec, "Error: failed", "All good.", [plot], [])
        assert result.reason == "artifact_empty"

    def test_error_masked_beats_entity_missing(self):
        """error_masked (priority 4) must fire before entity_missing (priority 5)."""
        _register_tool("plot_gene", {"gene": "gene"})
        spec = _make_spec("plot_gene", {"gene": "EGFR"})
        # EGFR absent from output triggers entity_missing;
        # Error: prefix + no indicator triggers error_masked
        result = check(spec, "Error: failed", "All good.", [], [])
        assert result.reason == "error_masked"

    def test_artifact_empty_beats_entity_missing(self):
        """artifact_empty (priority 3) must fire before entity_missing (priority 5)."""
        _register_tool("plot_gene", {"gene": "gene"})
        spec = _make_spec("plot_gene", {"gene": "EGFR"})
        plot = PlotResult(img_bytes=b"")
        result = check(spec, "plot done", "Done.", [plot], [])
        assert result.reason == "artifact_empty"

    def test_full_priority_chain_tool_mismatch_wins_over_all(self):
        """When all five triggers are possible, tool_mismatch wins."""
        _register_tool_with_mixed_params(
            "spec_tool", field_params={"gene": "gene"}, mechanical_params=[]
        )
        _register_tool("other_tool")
        spec = _make_spec("spec_tool", {"gene": "TP53"})
        plot = PlotResult(img_bytes=b"")  # artifact_empty candidate
        # Responder: names other_tool (tool_mismatch), claims artifact (artifact_claim_false
        # suppressed because plot_results non-empty), no error ack (error_masked candidate),
        # TP53 absent from output (entity_missing candidate)
        result = check(
            spec,
            "Error: failed",
            "The other_tool has been generated.",
            [plot],
            [],
        )
        assert result.reason == "tool_mismatch"


# ---------------------------------------------------------------------------
# GatekeeperResult dataclass contract
# ---------------------------------------------------------------------------


class TestGatekeeperResultContract:
    def test_status_literal_pass(self):
        r = GatekeeperResult(status="pass", reason="pass", explanation="ok")
        assert r.status == "pass"

    def test_status_literal_block(self):
        r = GatekeeperResult(status="block", reason="tool_mismatch", explanation="x")
        assert r.status == "block"

    def test_retry_hint_defaults_none(self):
        r = GatekeeperResult(status="pass", reason="pass", explanation="ok")
        assert r.retry_hint is None

    def test_retry_hint_settable(self):
        r = GatekeeperResult(
            status="block", reason="error_masked", explanation="x", retry_hint="fix it"
        )
        assert r.retry_hint == "fix it"

    def test_all_block_reasons_valid(self):
        reasons = [
            "pass",
            "tool_mismatch",
            "artifact_claim_false",
            "artifact_empty",
            "error_masked",
            "entity_missing",
        ]
        for reason in reasons:
            r = GatekeeperResult(
                status="pass" if reason == "pass" else "block",
                reason=reason,
                explanation="test",
            )
            assert r.reason == reason


# ---------------------------------------------------------------------------
# Artifact-claim phrase completeness
# ---------------------------------------------------------------------------


class TestArtifactClaimPhraseList:
    def test_phrase_list_is_nonempty(self):
        assert len(_ARTIFACT_CLAIM_PHRASES) > 0

    def test_all_phrases_are_lowercase(self):
        """All phrases must be lowercase (responder text is lowered before matching)."""
        for phrase in _ARTIFACT_CLAIM_PHRASES:
            assert phrase == phrase.lower(), f"Phrase not lowercase: {phrase!r}"

    def test_phrase_list_contains_has_been_generated(self):
        assert "has been generated" in _ARTIFACT_CLAIM_PHRASES

    def test_phrase_list_contains_ive_variants(self):
        assert any("i've" in p for p in _ARTIFACT_CLAIM_PHRASES)

    def test_phrase_list_contains_i_have_variants(self):
        assert any("i have" in p for p in _ARTIFACT_CLAIM_PHRASES)


# ---------------------------------------------------------------------------
# Error-indicator keyword list
# ---------------------------------------------------------------------------


class TestErrorIndicatorList:
    def test_error_indicators_nonempty(self):
        assert len(_ERROR_INDICATORS) > 0

    def test_all_indicators_lowercase(self):
        for kw in _ERROR_INDICATORS:
            assert kw == kw.lower(), f"Indicator not lowercase: {kw!r}"

    def test_error_keyword_present(self):
        assert "error" in _ERROR_INDICATORS

    def test_failed_keyword_present(self):
        assert "failed" in _ERROR_INDICATORS

    def test_unable_keyword_present(self):
        assert "unable" in _ERROR_INDICATORS
