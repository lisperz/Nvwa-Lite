"""Unit tests for src/agent/responder.py.

Coverage strategy:
- Pure helpers (`_compute_defaults_applied`, `_format_canonicalizations`,
  `_safe_repr`) are tested directly — no LLM, no mocks.
- `draft_response` (the public API) is tested by patching `_llm_call_text`,
  the provider seam below our wrapper. The subject under test is our
  prompt-construction + retry-hint logic, not the LLM's output.
- `_llm_call_text`'s empty-content fallback (the B1 fix) requires an empty
  raw response from the SDK. We patch `responder._llm_call_text`'s OpenAI
  client construction at the module boundary; the subject under test is our
  Python branch on empty content, not the LLM model.

All registry-touching tests use a `_clean_registry` autouse fixture to snapshot
and restore REGISTRY around each test.
"""

from __future__ import annotations

from typing import Any
from unittest.mock import MagicMock

import pytest

import src.agent.responder as responder_module
import src.tools.registry as reg_module
from src.agent.responder import (
    ResponderResult,
    _compute_defaults_applied,
    _format_canonicalizations,
    _safe_repr,
    draft_response,
)
from src.core.spec import Canonicalization, Spec
from src.tools.registry import ParamSpec, ToolEntry, register


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


def _make_spec(
    tool_name: str = "my_tool",
    params: dict | None = None,
    pre_canonical: dict | None = None,
    canonicalizations: list | None = None,
) -> Spec:
    return Spec(
        scenario_id="test",
        tool_name=tool_name,
        params=params or {},
        pre_canonical_params=pre_canonical or {},
        canonicalizations_applied=canonicalizations or [],
    )


def _make_entry(name: str, params: list[ParamSpec]) -> ToolEntry:
    return ToolEntry(
        name=name,
        kind="atomic",
        description=f"test tool {name}",
        params=params,
        callable=lambda *a, **k: None,
    )


def _patch_llm(
    monkeypatch,
    text: str = "narrative",
    in_tok: int = 100,
    out_tok: int = 50,
) -> dict[str, Any]:
    """Replace responder._llm_call_text with a fake; capture inputs for assertions."""
    captured: dict[str, Any] = {}

    def fake(system, user, history, model):
        captured["system"] = system
        captured["user"] = user
        captured["history"] = history
        captured["model"] = model
        return ResponderResult(text=text, input_tokens=in_tok, output_tokens=out_tok)

    monkeypatch.setattr(responder_module, "_llm_call_text", fake)
    return captured


def _fake_openai_response(
    content: str | None,
    prompt_tok: int = 10,
    comp_tok: int = 0,
):
    """Build a stub openai client response with the given content + usage."""
    msg = MagicMock()
    msg.content = content
    choice = MagicMock()
    choice.message = msg
    usage = MagicMock()
    usage.prompt_tokens = prompt_tok
    usage.completion_tokens = comp_tok
    resp = MagicMock()
    resp.choices = [choice]
    resp.usage = usage
    return resp


def _patch_openai_client(monkeypatch, response) -> MagicMock:
    """Patch responder.OpenAI so client.chat.completions.create returns `response`."""
    client = MagicMock()
    client.chat.completions.create.return_value = response
    monkeypatch.setattr(responder_module, "OpenAI", lambda api_key: client)
    monkeypatch.setenv("OPENAI_API_KEY", "test-key")
    return client


# ---------------------------------------------------------------------------
# draft_response — happy path + prompt construction
# ---------------------------------------------------------------------------


class TestDraftResponseHappyPath:
    def test_returns_responder_result(self, monkeypatch):
        _patch_llm(monkeypatch)
        result = draft_response(_make_spec(), "tool ran", chat_history=None)
        assert isinstance(result, ResponderResult)

    def test_text_propagated_from_llm(self, monkeypatch):
        _patch_llm(monkeypatch, text="hi there")
        result = draft_response(_make_spec(), "tool ran", chat_history=None)
        assert result.text == "hi there"

    def test_token_usage_propagated(self, monkeypatch):
        _patch_llm(monkeypatch, in_tok=500, out_tok=42)
        result = draft_response(_make_spec(), "tool ran", chat_history=None)
        assert result.input_tokens == 500
        assert result.output_tokens == 42

    def test_prompt_contains_tool_name(self, monkeypatch):
        captured = _patch_llm(monkeypatch)
        draft_response(_make_spec(tool_name="plot_umap"), "ok", chat_history=None)
        assert "plot_umap" in captured["system"]

    def test_prompt_contains_tool_output(self, monkeypatch):
        captured = _patch_llm(monkeypatch)
        draft_response(_make_spec(), "TOOL_OUTPUT_MARKER_xyz", chat_history=None)
        assert "TOOL_OUTPUT_MARKER_xyz" in captured["system"]

    def test_chat_history_none_becomes_empty_list(self, monkeypatch):
        captured = _patch_llm(monkeypatch)
        draft_response(_make_spec(), "ok", chat_history=None)
        assert captured["history"] == []

    def test_chat_history_passed_through(self, monkeypatch):
        captured = _patch_llm(monkeypatch)
        history = [("user", "hello"), ("assistant", "hi")]
        draft_response(_make_spec(), "ok", chat_history=history)
        assert captured["history"] == history

    def test_default_model_is_gpt_4o_mini(self, monkeypatch):
        captured = _patch_llm(monkeypatch)
        draft_response(_make_spec(), "ok")
        assert captured["model"] == "gpt-4o-mini"

    def test_model_override_propagates(self, monkeypatch):
        captured = _patch_llm(monkeypatch)
        draft_response(_make_spec(), "ok", model="gpt-5-pro")
        assert captured["model"] == "gpt-5-pro"


# ---------------------------------------------------------------------------
# draft_response — retry_hint
# ---------------------------------------------------------------------------


class TestDraftResponseRetryHint:
    def test_retry_hint_prepended_to_prompt(self, monkeypatch):
        captured = _patch_llm(monkeypatch)
        draft_response(
            _make_spec(),
            "ok",
            retry_hint="The narrative claimed an artifact that doesn't exist.",
        )
        assert "RETRY" in captured["system"]
        assert "doesn't exist" in captured["system"]

    def test_no_retry_hint_means_no_retry_block(self, monkeypatch):
        captured = _patch_llm(monkeypatch)
        draft_response(_make_spec(), "ok")
        assert "RETRY" not in captured["system"]

    def test_retry_hint_includes_anti_hallucination_framing(self, monkeypatch):
        captured = _patch_llm(monkeypatch)
        draft_response(_make_spec(), "ok", retry_hint="something")
        # Must frame as a constraint ("revise only") to discourage new fabrications
        assert "Revise only" in captured["system"]


# ---------------------------------------------------------------------------
# draft_response — tool description lookup
# ---------------------------------------------------------------------------


class TestDraftResponseToolLookup:
    def test_unknown_tool_uses_unknown_marker(self, monkeypatch):
        reg_module.REGISTRY.clear()
        captured = _patch_llm(monkeypatch)
        draft_response(_make_spec(tool_name="not_in_registry"), "ok")
        assert "(unknown tool)" in captured["system"]

    def test_known_tool_description_in_prompt(self, monkeypatch):
        reg_module.REGISTRY.clear()

        @register(description="DESC_MARKER_for_test")
        def my_test_tool():
            pass

        captured = _patch_llm(monkeypatch)
        draft_response(_make_spec(tool_name="my_test_tool"), "ok")
        assert "DESC_MARKER_for_test" in captured["system"]


# ---------------------------------------------------------------------------
# _compute_defaults_applied
# ---------------------------------------------------------------------------


class TestComputeDefaultsApplied:
    def test_none_entry_returns_empty(self):
        assert _compute_defaults_applied(_make_spec(), None) == []

    def test_param_with_default_not_in_spec_listed(self):
        entry = _make_entry("t", [
            ParamSpec(name="ncols", type="int", required=False, default=2),
        ])
        spec = _make_spec(params={})
        assert _compute_defaults_applied(spec, entry) == ["ncols"]

    def test_param_with_default_in_spec_not_listed(self):
        entry = _make_entry("t", [
            ParamSpec(name="ncols", type="int", required=False, default=2),
        ])
        spec = _make_spec(params={"ncols": 4})
        assert _compute_defaults_applied(spec, entry) == []

    def test_param_with_no_default_not_listed(self):
        entry = _make_entry("t", [
            ParamSpec(name="gene", type="str", required=True, default=None),
        ])
        spec = _make_spec(params={})
        assert _compute_defaults_applied(spec, entry) == []

    def test_mixed_params_only_unsupplied_with_default_listed(self):
        entry = _make_entry("t", [
            ParamSpec(name="gene", type="str", required=True, default=None),
            ParamSpec(name="ncols", type="int", required=False, default=2),
            ParamSpec(name="title", type="str", required=False, default="Plot"),
        ])
        spec = _make_spec(params={"gene": "BRCA1", "title": "My Plot"})
        # Only ncols: not supplied + has default → defaulted
        assert _compute_defaults_applied(spec, entry) == ["ncols"]

    def test_empty_params_list_returns_empty(self):
        entry = _make_entry("t", [])
        assert _compute_defaults_applied(_make_spec(), entry) == []


# ---------------------------------------------------------------------------
# _format_canonicalizations
# ---------------------------------------------------------------------------


class TestFormatCanonicalizations:
    def test_empty_returns_none_marker(self):
        assert _format_canonicalizations(_make_spec()) == "(none)"

    def test_single_canonicalization_no_context(self):
        c = Canonicalization(field="gene", raw="gapdh", canonical="GAPDH")
        spec = _make_spec(canonicalizations=[c])
        assert _format_canonicalizations(spec) == "gene: 'gapdh' → 'GAPDH'"

    def test_canonicalization_with_context_includes_marker(self):
        c = Canonicalization(
            field="gene",
            raw="brca",
            canonical="BRCA1",
            context={"reasoning": "exact match"},
        )
        spec = _make_spec(canonicalizations=[c])
        result = _format_canonicalizations(spec)
        assert "gene: 'brca' → 'BRCA1'" in result
        assert "context:" in result

    def test_multiple_canonicalizations_separated_by_semicolon(self):
        spec = _make_spec(canonicalizations=[
            Canonicalization(field="gene", raw="gapdh", canonical="GAPDH"),
            Canonicalization(field="cell_type", raw="t", canonical="T cell"),
        ])
        result = _format_canonicalizations(spec)
        assert "gene: 'gapdh' → 'GAPDH'" in result
        assert "cell_type: 't' → 'T cell'" in result
        assert "; " in result


# ---------------------------------------------------------------------------
# _safe_repr
# ---------------------------------------------------------------------------


class TestSafeRepr:
    def test_simple_dict(self):
        assert _safe_repr({"a": 1}) == "{'a': 1}"

    def test_simple_string(self):
        assert _safe_repr("hi") == "'hi'"

    def test_repr_raising_object_falls_back_to_str(self):
        class WeirdObj:
            def __repr__(self):
                raise RuntimeError("no repr")

            def __str__(self):
                return "STR_FALLBACK"

        assert _safe_repr(WeirdObj()) == "STR_FALLBACK"

    def test_none_value(self):
        assert _safe_repr(None) == "None"


# ---------------------------------------------------------------------------
# _llm_call_text — empty-content fallback (B1 fix)
# ---------------------------------------------------------------------------


class TestLLMCallTextEmptyContent:
    """Covers the B1 bug fix: empty LLM content must not crash the turn.

    The OpenAI client constructor is patched at the module boundary so the
    client returns a deterministic response. This is patching the SDK seam,
    not faking LLM intelligence — the subject under test is our Python
    branch when content is empty.
    """

    def test_empty_string_does_not_raise(self, monkeypatch):
        _patch_openai_client(monkeypatch, _fake_openai_response(""))
        result = responder_module._llm_call_text(
            system="sys", user="usr", history=[], model="gpt-4o-mini",
        )
        assert isinstance(result, ResponderResult)

    def test_none_content_does_not_raise(self, monkeypatch):
        _patch_openai_client(monkeypatch, _fake_openai_response(None))
        result = responder_module._llm_call_text(
            system="sys", user="usr", history=[], model="gpt-4o-mini",
        )
        assert isinstance(result, ResponderResult)

    def test_empty_returns_nonempty_fallback_text(self, monkeypatch):
        _patch_openai_client(monkeypatch, _fake_openai_response(""))
        result = responder_module._llm_call_text(
            system="sys", user="usr", history=[], model="gpt-4o-mini",
        )
        assert result.text

    def test_fallback_does_not_contain_artifact_claim_phrases(self, monkeypatch):
        """Fallback narrative must not trip the gatekeeper's artifact_claim_false check."""
        from src.agent.gatekeeper import _ARTIFACT_CLAIM_PHRASES

        _patch_openai_client(monkeypatch, _fake_openai_response(""))
        result = responder_module._llm_call_text(
            system="sys", user="usr", history=[], model="gpt-4o-mini",
        )
        lowered = result.text.lower()
        for phrase in _ARTIFACT_CLAIM_PHRASES:
            assert phrase not in lowered, (
                f"Fallback contains artifact-claim phrase: {phrase!r}"
            )

    def test_empty_propagates_token_usage_from_call(self, monkeypatch):
        _patch_openai_client(
            monkeypatch, _fake_openai_response("", prompt_tok=123, comp_tok=0),
        )
        result = responder_module._llm_call_text(
            system="sys", user="usr", history=[], model="gpt-4o-mini",
        )
        assert result.input_tokens == 123
        assert result.output_tokens == 0

    def test_whitespace_only_content_uses_fallback(self, monkeypatch):
        _patch_openai_client(monkeypatch, _fake_openai_response("   \n  \t  "))
        result = responder_module._llm_call_text(
            system="sys", user="usr", history=[], model="gpt-4o-mini",
        )
        # Stripping the original yields "" → fallback fires
        assert "  " not in result.text  # not the raw whitespace
        assert result.text


# ---------------------------------------------------------------------------
# _llm_call_text — happy path + auth + history
# ---------------------------------------------------------------------------


class TestLLMCallTextHappyPath:
    def test_missing_api_key_raises(self, monkeypatch):
        monkeypatch.delenv("OPENAI_API_KEY", raising=False)
        with pytest.raises(RuntimeError, match="OPENAI_API_KEY"):
            responder_module._llm_call_text(
                system="sys", user="usr", history=[], model="gpt-4o-mini",
            )

    def test_normal_content_returned_stripped(self, monkeypatch):
        _patch_openai_client(monkeypatch, _fake_openai_response("  hello  ", 10, 5))
        result = responder_module._llm_call_text(
            system="sys", user="usr", history=[], model="gpt-4o-mini",
        )
        assert result.text == "hello"
        assert result.input_tokens == 10
        assert result.output_tokens == 5

    def test_history_translated_to_chat_messages(self, monkeypatch):
        client = _patch_openai_client(
            monkeypatch, _fake_openai_response("ok", 1, 1),
        )
        responder_module._llm_call_text(
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
            monkeypatch, _fake_openai_response("ok", 1, 1),
        )
        responder_module._llm_call_text(
            system="sys", user="usr", history=[], model="gpt-4o-mini",
        )
        kwargs = client.chat.completions.create.call_args[1]
        assert kwargs["temperature"] == 0


# ---------------------------------------------------------------------------
# ResponderResult dataclass
# ---------------------------------------------------------------------------


class TestResponderResultContract:
    def test_fields_settable(self):
        r = ResponderResult(text="x", input_tokens=1, output_tokens=2)
        assert r.text == "x"
        assert r.input_tokens == 1
        assert r.output_tokens == 2
