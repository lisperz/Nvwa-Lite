"""Unit + integration tests for the output honesty guard (src/agent/output_guard.py).

T-029 hotfix. The guard detects fabricated artifact-creation claims from the LLM
(e.g. "the feature plot has been generated") when no tool was called, and drives
a one-shot retry. If the LLM re-hallucinates, the agent returns a user-facing
rephrase message instead of continuing to loop.

Run with: pytest tests/unit/test_output_guard.py -v
"""

from __future__ import annotations

from unittest.mock import MagicMock

import pytest

from src.agent.output_guard import (
    GuardResult,
    REPHRASE_FALLBACK_MESSAGE,
    build_corrective_system_message,
    check_output_honesty,
)


# ---------------------------------------------------------------------------
# Guard helper — triggered cases
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("text", [
    # Passive creation-claims
    "The feature plot of MKI67 expression has been generated, split by condition.",
    "The heatmap has been plotted across cell types.",
    "The UMAP has been created for your review.",
    # Active creation-claims
    "I've generated the UMAP plot for you.",
    "I've plotted the violin plot for MKI67.",
    "I have created the heatmap.",
    "I've produced the feature plot.",
])
def test_guard_triggers_on_plot_claim_without_tool(text):
    result = check_output_honesty(text, [])
    assert result.triggered is True, f"expected trigger for: {text!r}"
    assert result.matched_phrase is not None


@pytest.mark.parametrize("text", [
    "I've saved the CSV for download.",
    "I've created the table of DE results.",
    "I have saved the results file.",
    "The CSV has been generated.",
    "The table has been created.",
])
def test_guard_triggers_on_table_or_csv_claim_without_tool(text):
    result = check_output_honesty(text, [])
    assert result.triggered is True, f"expected trigger for: {text!r}"
    assert result.matched_phrase is not None


# ---------------------------------------------------------------------------
# Guard helper — not triggered
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("text", [
    "Here's your feature plot of MKI67.",
    "I've saved the CSV.",
    "The heatmap has been generated.",
])
def test_guard_passes_when_tool_was_called(text):
    fake_tool_call = {"name": "feature_plot", "args": {"gene": "MKI67"}, "id": "c1"}
    result = check_output_honesty(text, [fake_tool_call])
    assert result.triggered is False
    assert result.matched_phrase is None


@pytest.mark.parametrize("text", [
    "",
    "Please specify which cell type you want to analyze.",
    "The dataset has 45460 cells and 25639 genes.",
    "Which gene would you like to visualize?",
    "I cannot complete that request because the gene is not in the dataset.",
    # Commentary on existing artifacts — must NOT fire (these words appeared
    # in an earlier phrase list and caused false positives).
    "The plot shows the distribution of MKI67 across clusters.",
    "The visualization shows a clear separation of cell types.",
    "The table shows the top markers per cluster.",
    "Here's your answer: the dataset has 45k cells.",
    "Here is the list of available tools.",
])
def test_guard_passes_on_neutral_text(text):
    result = check_output_honesty(text, [])
    assert result.triggered is False
    assert result.matched_phrase is None


def test_guard_is_case_insensitive():
    result = check_output_honesty("THE PLOT HAS BEEN GENERATED", [])
    assert result.triggered is True


def test_guard_result_has_expected_shape():
    result = check_output_honesty("Here's your plot", [])
    assert isinstance(result, GuardResult)
    assert hasattr(result, "triggered")
    assert hasattr(result, "matched_phrase")
    assert hasattr(result, "reason")
    assert isinstance(result.reason, str)


# ---------------------------------------------------------------------------
# Corrective SystemMessage
# ---------------------------------------------------------------------------


def test_corrective_message_interpolates_matched_phrase():
    msg = build_corrective_system_message("the plot has been generated")
    assert "the plot has been generated" in msg.lower()


def test_corrective_message_has_anti_cache_instruction():
    # The MKI67 bug transcript suggests GPT-4o-mini may repeat a cached reply
    # when the corrective message is too generic. Must explicitly instruct
    # the model not to repeat the prior reply.
    msg = build_corrective_system_message("the plot has been generated")
    lowered = msg.lower()
    assert "do not repeat" in lowered or "not repeat" in lowered


def test_corrective_message_offers_honest_fallback():
    # The model must be allowed to admit it cannot, not forced into a tool call.
    msg = build_corrective_system_message("the plot has been generated")
    lowered = msg.lower()
    assert "cannot" in lowered or "plainly" in lowered or "unable" in lowered


# ---------------------------------------------------------------------------
# User-facing rephrase fallback
# ---------------------------------------------------------------------------


def test_rephrase_fallback_message_is_defined():
    assert isinstance(REPHRASE_FALLBACK_MESSAGE, str)
    assert len(REPHRASE_FALLBACK_MESSAGE) > 20


def test_rephrase_fallback_suggests_rephrasing():
    assert "rephrase" in REPHRASE_FALLBACK_MESSAGE.lower()


# ---------------------------------------------------------------------------
# AgentRunner integration — retry state machine
# ---------------------------------------------------------------------------


def _make_llm_response(content: str = "", tool_calls: list | None = None):
    """Build a MagicMock that mimics the LangChain AIMessage interface."""
    r = MagicMock()
    r.content = content
    r.tool_calls = tool_calls or []
    r.response_metadata = {}
    return r


class _FakeTool:
    """Minimal tool stand-in. MagicMock.name is reserved, so use a real class."""

    def __init__(self, name: str, result: str = "ok"):
        self.name = name
        self._result = result

    def invoke(self, args):
        return self._result


def test_agent_retries_after_first_guard_trigger_and_succeeds():
    """LLM hallucinates → guard fires → corrective SystemMessage appended →
    LLM calls tool on retry → guard passes → success."""
    from src.agent.core import AgentRunner

    hallucinated = _make_llm_response("The feature plot has been generated.")
    with_tool = _make_llm_response(
        "",
        [{"name": "feature_plot", "args": {"gene": "MKI67"}, "id": "c1"}],
    )
    final = _make_llm_response("Here is the plot you requested.")

    mock_llm = MagicMock()
    mock_llm.invoke = MagicMock(side_effect=[hallucinated, with_tool, final])

    runner = AgentRunner(
        llm_with_tools=mock_llm,
        tools=[_FakeTool("feature_plot", "plot_ok")],
        system_prompt="test-prompt",
    )
    response = runner.invoke("draw a feature plot of MKI67 by condition")

    assert mock_llm.invoke.call_count == 3
    assert response.tool_called is True
    assert response.text == "Here is the plot you requested."


def test_agent_exhausts_after_second_guard_trigger_and_returns_rephrase_message():
    """LLM hallucinates → corrective msg → LLM STILL hallucinates →
    guard exhausted → user-facing rephrase message returned, loop exits."""
    from src.agent.core import AgentRunner

    hallucinated_1 = _make_llm_response("The feature plot has been generated.")
    hallucinated_2 = _make_llm_response("I've generated the feature plot, split by condition.")

    mock_llm = MagicMock()
    mock_llm.invoke = MagicMock(side_effect=[hallucinated_1, hallucinated_2])

    runner = AgentRunner(
        llm_with_tools=mock_llm,
        tools=[],
        system_prompt="test-prompt",
    )
    response = runner.invoke("draw a feature plot of MKI67 by condition")

    assert mock_llm.invoke.call_count == 2
    assert response.text == REPHRASE_FALLBACK_MESSAGE
    assert response.tool_called is False


def test_agent_no_retry_when_first_response_is_clean():
    """LLM produces non-artifact text → guard passes → no retry, single turn."""
    from src.agent.core import AgentRunner

    clean = _make_llm_response("The dataset has 45460 cells and 25639 genes.")

    mock_llm = MagicMock()
    mock_llm.invoke = MagicMock(side_effect=[clean])

    runner = AgentRunner(
        llm_with_tools=mock_llm,
        tools=[],
        system_prompt="test-prompt",
    )
    response = runner.invoke("how many cells are in the dataset?")

    assert mock_llm.invoke.call_count == 1
    assert response.text == "The dataset has 45460 cells and 25639 genes."
    assert response.tool_called is False


def test_agent_no_retry_when_tool_is_called_with_artifact_claim():
    """LLM calls a tool AND claims artifact creation → guard passes (tool_calls non-empty)."""
    from src.agent.core import AgentRunner

    with_tool = _make_llm_response(
        "",
        [{"name": "feature_plot", "args": {"gene": "MKI67"}, "id": "c1"}],
    )
    final = _make_llm_response("The feature plot has been generated.")

    mock_llm = MagicMock()
    mock_llm.invoke = MagicMock(side_effect=[with_tool, final])

    runner = AgentRunner(
        llm_with_tools=mock_llm,
        tools=[_FakeTool("feature_plot", "plot_ok")],
        system_prompt="test-prompt",
    )
    response = runner.invoke("draw a feature plot of MKI67")

    assert mock_llm.invoke.call_count == 2
    assert response.tool_called is True
    assert response.text == "The feature plot has been generated."
