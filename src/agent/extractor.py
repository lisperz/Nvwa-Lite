"""LLM-based spec extractor — converts a user prompt into a Spec via OpenAI JSON output.

Inline prompt (coupled to the Spec + tool catalog). prompts.py owns the
conversational/system prompt; this file owns the spec-extraction prompt only.

Provider portability: the only OpenAI-specific code is `_llm_call_json` below.
Schema build, prompt rendering, and Spec construction are provider-agnostic.
Switching providers means swapping the body of that one function.
"""

from __future__ import annotations

import json
import logging
import os

from openai import OpenAI

# Importing src.tools triggers @register side-effects from family modules so REGISTRY
# is populated before extractor reads it. Safe to keep even if app already imported it.
import src.tools  # noqa: F401

from src.core.spec import Spec
from src.tools.registry import REGISTRY, get_tool_names

logger = logging.getLogger(__name__)


_DEFAULT_MODEL = "gpt-4o-mini"


_EXTRACTOR_PROMPT_TMPL = """You convert a user query about a single-cell RNA-seq dataset into a structured Spec.

Output a JSON object with EXACTLY these three keys:
- "scenario_id": a short string label for the task type. Use "plot", "compare", "explore", "qc", or "unknown".
- "tool_name": one of the tool names listed below, OR the literal string "none" if no listed tool fits the user's intent.
- "pre_canonical_params": an object mapping the chosen tool's parameter names to the user's values. Use the user's words verbatim (e.g. "gapdh" not "GAPDH"); a downstream resolver canonicalizes them. Use {{}} when tool_name is "none".

Available tools:
{tool_catalog}

When to use tool_name="none":
- The user's intent cannot be fulfilled by any tool listed above.
- Examples: trajectory/velocity analysis, dataset integration, plot types not listed, unrelated visualization requests.
- Do NOT force-pick the closest-named tool. Returning "none" routes the request to the legacy fallback path; force-picking would mis-dispatch and confuse the user.

Do not invent parameter values not implied by the user query. If a required parameter is missing from the user query, leave it out of pre_canonical_params — the validator will surface it as needing input.

Output ONLY the JSON object, no commentary."""


def extract(
    user_prompt: str,
    chat_history: list[tuple[str, str]] | None = None,
    model: str = _DEFAULT_MODEL,
) -> Spec:
    """Extract a Spec from the user prompt via structured-output LLM call.

    The extractor populates spec.pre_canonical_params with the user's verbatim
    values; the resolver downstream builds spec.params (additive pattern).
    """
    system = _EXTRACTOR_PROMPT_TMPL.format(tool_catalog=_format_tool_catalog())
    response_dict = _llm_call_json(
        system=system,
        user=user_prompt,
        history=chat_history or [],
        model=model,
    )
    if "tool_name" not in response_dict:
        raise RuntimeError("Extractor: response missing 'tool_name' key.")
    return Spec(
        scenario_id=str(response_dict.get("scenario_id", "unknown")),
        tool_name=str(response_dict["tool_name"]),
        pre_canonical_params=dict(response_dict.get("pre_canonical_params", {})),
    )


def _format_tool_catalog() -> str:
    """Render registered tools + their LLM-facing params for the extractor prompt."""
    if not REGISTRY:
        return "(no tools registered)"
    lines: list[str] = []
    for entry in REGISTRY.values():
        lines.append(f"- {entry.name}: {entry.description}")
        for p in entry.params:
            req = "required" if p.required else f"default={p.default!r}"
            enum = f" [one of {p.enum}]" if p.enum else ""
            desc = f" — {p.description}" if p.description else ""
            lines.append(f"    {p.name}: {p.type}, {req}{enum}{desc}")
    return "\n".join(lines)


def _llm_call_json(
    system: str,
    user: str,
    history: list[tuple[str, str]],
    model: str,
) -> dict:
    """Single seam for the LLM provider. Returns parsed JSON dict.

    Switching to a different provider (Anthropic, local model, etc.) means
    swapping the body of this function only.
    """
    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        raise RuntimeError("Extractor: OPENAI_API_KEY not set in environment.")

    client = OpenAI(api_key=api_key)

    messages: list[dict] = [{"role": "system", "content": system}]
    for role, content in history:
        api_role = "assistant" if role == "assistant" else "user"
        messages.append({"role": api_role, "content": content})
    messages.append({"role": "user", "content": user})

    response = client.chat.completions.create(
        model=model,
        messages=messages,
        response_format={"type": "json_object"},
        temperature=0,
    )
    content = response.choices[0].message.content
    if not content:
        raise RuntimeError("Extractor: empty response from LLM.")
    return json.loads(content)
