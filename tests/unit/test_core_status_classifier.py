"""Unit tests for _classify_tool_status in src/agent/core.py.

T-033 verification: classifier infers (status, error_msg) from a tool return
value. Error contract is capital-E ``Error:`` prefix on the first line.

Run with: pytest tests/unit/test_core_status_classifier.py -v
"""

import pytest

from src.agent.core import _classify_tool_status


# ---------------------------------------------------------------------------
# Error path — canonicalized prefix forms
# Sample of the 9 prefix forms, post-canonicalization (all start with `Error:`).
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "result",
    [
        "Error: Unknown tool 'foo'.",
        "Error: Unexpected — division by zero",
        "Error: Preprocessing — adata.X is None",
        "Error: Marker gene analysis — group 'X' not found",
        "Error: Cluster DEG analysis — rank_genes_groups failed",
        "Error: Pairwise DE — group 'late' has 0 cells",
        "Error: Inspecting metadata — obs column missing",
        "Error: Calculating mitochondrial percentage — no MT- genes",
        # Real customer-session incidents (2026-04-15 demo logs)
        "Error: Gene 'FOOBAR' not found in dataset.",
        "Error: Could not find category value 'late_stage'.",
    ],
)
def test_error_prefix_classifies_as_error(result):
    status, error_msg = _classify_tool_status(result)
    assert status == "error"
    assert error_msg is not None
    assert error_msg.startswith("Error:")


def test_error_with_leading_whitespace_still_classifies():
    # Classifier lstrips before checking prefix.
    status, error_msg = _classify_tool_status("   Error: leading spaces")
    assert status == "error"
    assert error_msg == "Error: leading spaces"


def test_error_with_multiline_body_truncates_to_first_line():
    result = "Error: first line failure\ntraceback line 1\ntraceback line 2"
    status, error_msg = _classify_tool_status(result)
    assert status == "error"
    assert error_msg == "Error: first line failure"
    assert "traceback" not in error_msg


def test_error_msg_capped_at_200_chars():
    long_tail = "x" * 500
    result = f"Error: {long_tail}"
    status, error_msg = _classify_tool_status(result)
    assert status == "error"
    assert len(error_msg) <= 200


# ---------------------------------------------------------------------------
# Success path — plain tool return strings
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "result",
    [
        "Plot saved to plots/umap_001.png",
        "Differential expression complete: 142 genes significant at p<0.05.",
        "Found 3 clusters matching 'T cell'.",
        "PlotResult(figure_path='plots/foo.png', caption='...', ...)",
        "Summary: n_cells=12345, n_genes=20000, obs_columns=['celltype','condition']",
    ],
)
def test_success_strings_classify_as_success(result):
    status, error_msg = _classify_tool_status(result)
    assert status == "success"
    assert error_msg is None


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

def test_empty_string_is_success():
    assert _classify_tool_status("") == ("success", None)


def test_none_is_success():
    # Non-str tool returns fall through to success — tool didn't raise.
    assert _classify_tool_status(None) == ("success", None)


def test_non_string_return_is_success():
    # Tools currently return str, but classifier must be safe for future
    # typed returns (PlotResult dict, etc.).
    assert _classify_tool_status({"figure_path": "..."}) == ("success", None)
    assert _classify_tool_status(42) == ("success", None)


def test_lowercase_error_prefix_is_success():
    # Contract is capital-E `Error:` — lowercase is not an error marker.
    status, error_msg = _classify_tool_status("error: something weird")
    assert status == "success"
    assert error_msg is None


def test_error_substring_mid_body_is_success():
    # "No error" in body must not flip status.
    result = "Gene count: 12345. No error: all good."
    status, error_msg = _classify_tool_status(result)
    assert status == "success"
    assert error_msg is None
