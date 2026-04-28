"""Post-dispatch response drafter — pairs tool output with a short user-facing narrative.

The third LLM stage in the spec pipeline (router → extractor → resolver → validator →
dispatcher → **responder**). Inline prompt; provider seam mirrors `extractor.py`.

Responder is a leaf: input text → output text. It never makes tool-calling decisions.
For workflow tools that internally call multiple atomic tools, the workflow consolidates
results into a single `tool_output` string before responder sees it.

Logger-free by design. Returns `(narrative, token_usage)` so core.py owns logging.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from typing import TYPE_CHECKING

from openai import OpenAI

from src.core.spec import Spec
from src.core.registry import get_tool

if TYPE_CHECKING:
    pass

logger = logging.getLogger(__name__)


_DEFAULT_MODEL = "gpt-4o-mini"


_RESPONDER_PROMPT_TMPL = """You write a short, accurate user-facing narrative for the result of a tool call on a single-cell RNA-seq dataset.

The user asked something. A tool ran. You receive the spec describing what the tool did and the tool's raw output. Your job is to narrate concisely what was done and what the user is seeing — pairing with any artifact (plot, table) the tool already produced.

Rules:
1. **Stick to what actually happened.** Do not describe data or steps the tool did not perform. Do not invent observations from the artifact.
2. **Disclose canonicalizations.** If the spec shows a value was canonicalized (raw → canonical), surface it briefly: "I used GAPDH (you wrote 'gapdh', case-corrected)".
3. **Disclose defaults.** If parameters were filled by tool defaults (not provided by the user), mention them briefly: "Using ncols=2 (default)".
4. **Disclose resolver context.** If a canonicalization carries a `context` note (e.g., disambiguation reasoning), surface it.
5. **Always state the tool name used** by its canonical name (e.g., `plot_umap`, `dataset_overview`). Without this, the user can't verify which tool actually ran. The downstream gatekeeper relies on this disclosure to validate your narrative.
6. **If tool_output starts with "Error:", do not pretend it succeeded.** Rephrase the error in user-friendly terms and suggest a next step.
7. Keep the narrative short — typically 1-3 sentences. The artifact (plot/table) carries the visual; you carry the explanation.
8. Do not output markdown headers, bullet lists for the disclosure parts, or metadata. Plain prose.

Spec context:
- tool: {tool_name}
- tool description: {tool_description}
- raw user inputs (pre-canonicalization): {pre_canonical_params}
- canonical params actually used: {params}
- canonicalizations applied: {canonicalizations}
- defaults applied (not provided by user): {defaults_applied}

Tool output (raw):
---
{tool_output}
---

Write the user-facing narrative now."""


@dataclass
class ResponderResult:
    """Output of `draft_response`. Token counts let core.py log usage uniformly."""

    text: str
    input_tokens: int
    output_tokens: int


def draft_response(
    spec: Spec,
    tool_output: str,
    chat_history: list[tuple[str, str]] | None = None,
    model: str = _DEFAULT_MODEL,
    retry_hint: str | None = None,
) -> ResponderResult:
    """Draft a short user-facing narrative pairing the tool output with the spec.

    `retry_hint` is set by core.py on a gatekeeper-block retry. It prepends an
    anti-hallucination corrective block to the prompt, framed as a constraint
    ("revise only to fix this") to minimize the risk of the hint priming new
    inventions.

    Returns a ResponderResult with the narrative text and token usage. core.py
    is responsible for logging via EventLogger.
    """
    entry = get_tool(spec.tool_name)
    tool_description = entry.description if entry else "(unknown tool)"
    defaults_applied = _compute_defaults_applied(spec, entry)

    prompt = _RESPONDER_PROMPT_TMPL.format(
        tool_name=spec.tool_name,
        tool_description=tool_description,
        pre_canonical_params=_safe_repr(spec.pre_canonical_params),
        params=_safe_repr(spec.params),
        canonicalizations=_format_canonicalizations(spec),
        defaults_applied=defaults_applied if defaults_applied else "(none)",
        tool_output=tool_output,
    )

    if retry_hint:
        prompt = (
            f"IMPORTANT — RETRY. Your previous draft was blocked by the deterministic "
            f"gatekeeper for this reason: {retry_hint}\n"
            f"Revise only to fix this. Do not introduce information beyond the tool "
            f"output below.\n\n"
        ) + prompt

    return _llm_call_text(
        system=prompt,
        user="Generate the narrative.",
        history=chat_history or [],
        model=model,
    )


def _compute_defaults_applied(spec: Spec, entry) -> list[str]:
    """List param names whose values come from tool-signature defaults (not user-provided).

    Computed fresh on demand: a param is "defaulted" if it has a non-None default in the
    registry AND was not present in spec.params (the resolver's output, what the tool sees).
    Per step 6 amendment, this is no longer carried on Spec.
    """
    if entry is None:
        return []
    return [
        p.name
        for p in entry.params
        if p.name not in spec.params and p.default is not None
    ]


def _format_canonicalizations(spec: Spec) -> str:
    """Render canonicalizations_applied compactly for the prompt."""
    if not spec.canonicalizations_applied:
        return "(none)"
    parts: list[str] = []
    for c in spec.canonicalizations_applied:
        ctx = f" [context: {c.context}]" if c.context else ""
        parts.append(f"{c.field}: '{c.raw}' → '{c.canonical}'{ctx}")
    return "; ".join(parts)


def _safe_repr(value) -> str:
    """Compact repr for prompt embedding; falls back to str() on weird types."""
    try:
        return repr(value)
    except Exception:
        return str(value)


def _llm_call_text(
    system: str,
    user: str,
    history: list[tuple[str, str]],
    model: str,
) -> ResponderResult:
    """Single seam for the LLM provider. Returns narrative text + token usage.

    Switching to a different provider (Anthropic, local model, etc.) means
    swapping the body of this function only.
    """
    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        raise RuntimeError("Responder: OPENAI_API_KEY not set in environment.")

    client = OpenAI(api_key=api_key)

    messages: list[dict] = [{"role": "system", "content": system}]
    for role, content in history:
        api_role = "assistant" if role == "assistant" else "user"
        messages.append({"role": api_role, "content": content})
    messages.append({"role": "user", "content": user})

    response = client.chat.completions.create(
        model=model,
        messages=messages,
        temperature=0,
    )
    text = (response.choices[0].message.content or "").strip()
    usage = response.usage
    if not text:
        # Empty LLM content is rare but real (content-filter, transient API
        # glitch). The tool already ran; don't crash the turn — return a
        # neutral fallback so the gatekeeper still gets a chance to validate.
        logger.warning("Responder: empty response from LLM; returning fallback narrative.")
        text = (
            "The tool ran, but I couldn't generate a narrative for the result. "
            "Please check the output above or rephrase your request."
        )
    return ResponderResult(
        text=text,
        input_tokens=usage.prompt_tokens if usage else 0,
        output_tokens=usage.completion_tokens if usage else 0,
    )
