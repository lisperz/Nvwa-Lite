"""Unit tests for src/agent/extractor.py.

Coverage strategy:
- `extract` (public API) is tested by patching `_llm_call_json`, the provider
  seam below our wrapper. The subject under test is our Spec construction +
  schema enforcement, not the LLM's output.
- `_format_tool_catalog` is a pure function over the REGISTRY — tested
  directly with controlled tool entries.
- `_llm_call_json` edge cases (missing API key, empty content) require an
  empty SDK response. The OpenAI client constructor is patched at the
  module boundary; subject under test is our Python branch on empty / bad
  content, not the LLM model.
"""

from __future__ import annotations

from typing import Any
from unittest.mock import MagicMock

import pytest

import src.agent.extractor as extractor_module
import src.core.registry as reg_module
from src.agent.extractor import (
    _format_tool_catalog,
    extract,
)
from src.core.spec import Spec
from src.core.registry import register


# ---------------------------------------------------------------------------
# Registry isolation
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _clean_registry():
    snapshot = dict(reg_module.REGISTRY)
    yield
    reg_module.REGISTRY.clear()
    reg_module.REGISTRY.update(snapshot)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _patch_llm(monkeypatch, response_dict: dict) -> dict[str, Any]:
    """Replace extractor._llm_call_json with a fake; capture inputs."""
    captured: dict[str, Any] = {}

    def fake(system, user, history, model):
        captured["system"] = system
        captured["user"] = user
        captured["history"] = history
        captured["model"] = model
        return dict(response_dict)

    monkeypatch.setattr(extractor_module, "_llm_call_json", fake)
    return captured


def _fake_openai_response(content: str | None):
    msg = MagicMock()
    msg.content = content
    choice = MagicMock()
    choice.message = msg
    resp = MagicMock()
    resp.choices = [choice]
    return resp


def _patch_openai_client(monkeypatch, response) -> MagicMock:
    client = MagicMock()
    client.chat.completions.create.return_value = response
    monkeypatch.setattr(extractor_module, "OpenAI", lambda api_key: client)
    monkeypatch.setenv("OPENAI_API_KEY", "test-key")
    return client


# ---------------------------------------------------------------------------
# extract — happy path
# ---------------------------------------------------------------------------


class TestExtractHappyPath:
    def test_returns_spec_instance(self, monkeypatch):
        _patch_llm(monkeypatch, {
            "scenario_id": "plot",
            "tool_name": "plot_umap",
            "pre_canonical_params": {"gene": "BRCA1"},
        })
        result = extract("plot BRCA1")
        assert isinstance(result, Spec)

    def test_scenario_id_passed_through(self, monkeypatch):
        _patch_llm(monkeypatch, {
            "scenario_id": "qc",
            "tool_name": "dataset_overview",
            "pre_canonical_params": {},
        })
        result = extract("what's in this data")
        assert result.scenario_id == "qc"

    def test_tool_name_passed_through(self, monkeypatch):
        _patch_llm(monkeypatch, {
            "scenario_id": "plot",
            "tool_name": "plot_umap",
            "pre_canonical_params": {},
        })
        result = extract("plot it")
        assert result.tool_name == "plot_umap"

    def test_pre_canonical_params_passed_through(self, monkeypatch):
        _patch_llm(monkeypatch, {
            "scenario_id": "plot",
            "tool_name": "plot_gene",
            "pre_canonical_params": {"gene": "gapdh", "ncols": 2},
        })
        result = extract("plot gapdh")
        assert result.pre_canonical_params == {"gene": "gapdh", "ncols": 2}

    def test_decline_sentinel_tool_name_none(self, monkeypatch):
        """tool_name='none' is the explicit decline sentinel — must round-trip."""
        _patch_llm(monkeypatch, {
            "scenario_id": "unknown",
            "tool_name": "none",
            "pre_canonical_params": {},
        })
        result = extract("run trajectory analysis")
        assert result.tool_name == "none"

    def test_user_prompt_passed_to_seam(self, monkeypatch):
        captured = _patch_llm(monkeypatch, {
            "scenario_id": "plot",
            "tool_name": "plot_umap",
            "pre_canonical_params": {},
        })
        extract("USER_PROMPT_MARKER")
        assert captured["user"] == "USER_PROMPT_MARKER"

    def test_chat_history_none_becomes_empty_list(self, monkeypatch):
        captured = _patch_llm(monkeypatch, {
            "scenario_id": "plot",
            "tool_name": "plot_umap",
            "pre_canonical_params": {},
        })
        extract("plot it", chat_history=None)
        assert captured["history"] == []

    def test_chat_history_passed_through(self, monkeypatch):
        captured = _patch_llm(monkeypatch, {
            "scenario_id": "plot",
            "tool_name": "plot_umap",
            "pre_canonical_params": {},
        })
        history = [("user", "earlier"), ("assistant", "ok")]
        extract("now plot", chat_history=history)
        assert captured["history"] == history

    def test_default_model_is_gpt_4o_mini(self, monkeypatch):
        captured = _patch_llm(monkeypatch, {
            "scenario_id": "plot",
            "tool_name": "plot_umap",
            "pre_canonical_params": {},
        })
        extract("plot it")
        assert captured["model"] == "gpt-4o-mini"


# ---------------------------------------------------------------------------
# extract — schema enforcement (NB1 fix)
# ---------------------------------------------------------------------------


class TestExtractSchemaEnforcement:
    def test_missing_tool_name_raises_explicit_runtime_error(self, monkeypatch):
        """NB1 fix: previously raised bare KeyError, now an explicit RuntimeError."""
        _patch_llm(monkeypatch, {
            "scenario_id": "plot",
            # tool_name missing
            "pre_canonical_params": {},
        })
        with pytest.raises(RuntimeError, match="tool_name"):
            extract("plot it")

    def test_missing_tool_name_error_message_actionable(self, monkeypatch):
        """The error message must surface what was missing, not a bare 'tool_name'."""
        _patch_llm(monkeypatch, {"scenario_id": "plot"})
        with pytest.raises(RuntimeError) as exc_info:
            extract("plot it")
        msg = str(exc_info.value)
        assert "Extractor" in msg
        assert "tool_name" in msg
        assert "missing" in msg.lower()

    def test_missing_scenario_id_defaults_to_unknown(self, monkeypatch):
        _patch_llm(monkeypatch, {
            "tool_name": "plot_umap",
            "pre_canonical_params": {},
        })
        result = extract("plot it")
        assert result.scenario_id == "unknown"

    def test_missing_pre_canonical_params_defaults_to_empty(self, monkeypatch):
        _patch_llm(monkeypatch, {
            "scenario_id": "plot",
            "tool_name": "plot_umap",
        })
        result = extract("plot it")
        assert result.pre_canonical_params == {}

    def test_tool_name_coerced_to_string(self, monkeypatch):
        """If the LLM returns a non-string tool_name (rare), it's str()'d."""
        _patch_llm(monkeypatch, {
            "scenario_id": "plot",
            "tool_name": 42,  # malformed but recoverable
            "pre_canonical_params": {},
        })
        result = extract("plot it")
        assert result.tool_name == "42"


# ---------------------------------------------------------------------------
# _format_tool_catalog
# ---------------------------------------------------------------------------


class TestFormatToolCatalog:
    def test_empty_registry_returns_marker(self):
        reg_module.REGISTRY.clear()
        assert _format_tool_catalog() == "(no tools registered)"

    def test_single_tool_with_no_params(self):
        reg_module.REGISTRY.clear()

        @register(description="Simple tool description.")
        def simple_tool():
            pass

        result = _format_tool_catalog()
        assert "simple_tool" in result
        assert "Simple tool description." in result

    def test_required_param_marked_required(self):
        reg_module.REGISTRY.clear()

        @register(description="Plot a gene.")
        def plot_gene(gene: str):
            pass

        result = _format_tool_catalog()
        assert "gene" in result
        assert "required" in result

    def test_optional_param_shows_default(self):
        reg_module.REGISTRY.clear()

        @register(description="Plot UMAP.")
        def plot_umap(ncols: int = 2):
            pass

        result = _format_tool_catalog()
        assert "ncols" in result
        assert "default=2" in result

    def test_enum_param_shows_options(self):
        reg_module.REGISTRY.clear()

        @register(
            description="Plot with style.",
            params={"style": {"enum": ["dot", "violin"]}},
        )
        def plot_style(style: str = "dot"):
            pass

        result = _format_tool_catalog()
        assert "[one of" in result
        assert "dot" in result
        assert "violin" in result

    def test_description_param_appended(self):
        reg_module.REGISTRY.clear()

        @register(
            description="Plot a gene.",
            params={"gene": {"description": "Gene symbol to plot."}},
        )
        def plot_gene(gene: str):
            pass

        result = _format_tool_catalog()
        assert "Gene symbol to plot." in result

    def test_multiple_tools_separated(self):
        reg_module.REGISTRY.clear()

        @register(description="Tool A description.")
        def tool_a():
            pass

        @register(description="Tool B description.")
        def tool_b():
            pass

        result = _format_tool_catalog()
        assert "tool_a" in result
        assert "tool_b" in result

    def test_adata_param_filtered_from_catalog(self):
        """adata is injected by the dispatcher; must not appear in the LLM-facing catalog."""
        reg_module.REGISTRY.clear()

        @register(description="Reads adata.")
        def reads_adata(adata, gene: str):
            pass

        result = _format_tool_catalog()
        assert "gene" in result
        # "adata" must not appear as a param name. (It may appear inside the
        # tool description if author included it, but the synthetic tool above
        # doesn't.)
        # Check the line-by-line: no line should be "    adata: ..."
        for line in result.splitlines():
            assert not line.strip().startswith("adata:")


# ---------------------------------------------------------------------------
# _llm_call_json — auth + content edge cases
# ---------------------------------------------------------------------------


class TestLLMCallJsonAuth:
    def test_missing_api_key_raises(self, monkeypatch):
        monkeypatch.delenv("OPENAI_API_KEY", raising=False)
        with pytest.raises(RuntimeError, match="OPENAI_API_KEY"):
            extractor_module._llm_call_json(
                system="sys", user="usr", history=[], model="gpt-4o-mini",
            )

    def test_empty_content_raises(self, monkeypatch):
        """Extractor (unlike responder) raises on empty — caller falls back to LangChain."""
        _patch_openai_client(monkeypatch, _fake_openai_response(""))
        with pytest.raises(RuntimeError, match="empty"):
            extractor_module._llm_call_json(
                system="sys", user="usr", history=[], model="gpt-4o-mini",
            )

    def test_none_content_raises(self, monkeypatch):
        _patch_openai_client(monkeypatch, _fake_openai_response(None))
        with pytest.raises(RuntimeError, match="empty"):
            extractor_module._llm_call_json(
                system="sys", user="usr", history=[], model="gpt-4o-mini",
            )

    def test_invalid_json_raises(self, monkeypatch):
        """Bad JSON → json.loads raises; caller catches and falls back."""
        _patch_openai_client(monkeypatch, _fake_openai_response("not json {{{"))
        with pytest.raises(Exception):  # json.JSONDecodeError is a ValueError subclass
            extractor_module._llm_call_json(
                system="sys", user="usr", history=[], model="gpt-4o-mini",
            )

    def test_valid_json_parsed_to_dict(self, monkeypatch):
        _patch_openai_client(
            monkeypatch,
            _fake_openai_response('{"tool_name": "x", "scenario_id": "y"}'),
        )
        result = extractor_module._llm_call_json(
            system="sys", user="usr", history=[], model="gpt-4o-mini",
        )
        assert result == {"tool_name": "x", "scenario_id": "y"}

    def test_response_format_json_object_passed(self, monkeypatch):
        """Extractor must request structured JSON output from the API."""
        client = _patch_openai_client(
            monkeypatch, _fake_openai_response('{"tool_name": "x"}'),
        )
        extractor_module._llm_call_json(
            system="sys", user="usr", history=[], model="gpt-4o-mini",
        )
        kwargs = client.chat.completions.create.call_args[1]
        assert kwargs["response_format"] == {"type": "json_object"}

    def test_history_translated_to_chat_messages(self, monkeypatch):
        client = _patch_openai_client(
            monkeypatch, _fake_openai_response('{"tool_name": "x"}'),
        )
        extractor_module._llm_call_json(
            system="sys",
            user="usr",
            history=[("user", "hi"), ("assistant", "hello"), ("garbage", "x")],
            model="gpt-4o-mini",
        )
        sent = client.chat.completions.create.call_args[1]["messages"]
        roles = [m["role"] for m in sent]
        # system + (user, assistant, garbage→user) + final user
        assert roles == ["system", "user", "assistant", "user", "user"]

    def test_temperature_zero_passed(self, monkeypatch):
        client = _patch_openai_client(
            monkeypatch, _fake_openai_response('{"tool_name": "x"}'),
        )
        extractor_module._llm_call_json(
            system="sys", user="usr", history=[], model="gpt-4o-mini",
        )
        assert client.chat.completions.create.call_args[1]["temperature"] == 0
