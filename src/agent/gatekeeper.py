"""Deterministic gatekeeper — validates responder's narrative against ground truth.

Pair-with-responder analog to validator-with-extractor: every LLM stage is followed
by a deterministic check. Gatekeeper runs after the responder; catches severe +
silent failures without any LLM call.

Scope (MVP):
  - tool_mismatch:        responder names a different tool than was dispatched
  - artifact_claim_false: responder claims artifact creation but no PlotResult/TableResult
  - artifact_empty:       artifact produced but zero bytes / zero rows
  - error_masked:         tool returned "Error:" but responder doesn't acknowledge
  - entity_missing:       a field_type param's value is absent from tool_output

Out of scope for deterministic gatekeeper (deferred to T-004 LLM-as-judge layer):
  - Responder narrating paraphrased inventions (responder prompt's job)
  - Output quality / scientific correctness (tool's / Layer-4's job)
  - Silent substitution between legacy and @register tools (needs LLM intent match)

Block priority (first-block-wins, most-severe first):
  tool_mismatch → artifact_claim_false → artifact_empty → error_masked → entity_missing
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal, Optional

from src.core.spec import Spec
from src.tools.registry import get_tool, get_tool_names


# Artifact-claim phrases (ported from src/agent/output_guard.py)
_ARTIFACT_CLAIM_PHRASES: tuple[str, ...] = (
    "has been generated",
    "has been plotted",
    "has been created",
    "has been saved",
    "has been produced",
    "i've generated",
    "i have generated",
    "i've plotted",
    "i have plotted",
    "i've created",
    "i have created",
    "i've saved",
    "i have saved",
    "i've produced",
    "i have produced",
)

# Error-acknowledgment keywords (Q9 option ii — broader list)
_ERROR_INDICATORS: tuple[str, ...] = (
    "error",
    "couldn't",
    "didn't work",
    "doesn't work",
    "unable",
    "failed",
    "problem",
    "issue",
    "wrong",
    "no result",
    "empty",
    "not able",
    "cannot",
)


@dataclass
class GatekeeperResult:
    """Deterministic verdict on the responder's narrative.

    `retry_hint` is populated only for reasons where a responder re-draft could
    plausibly fix the issue. For reasons where the root cause is upstream
    (e.g. artifact_empty — tool produced nothing), retry_hint is None and
    core.py should surface to the user directly.
    """

    status: Literal["pass", "block"]
    reason: Literal[
        "pass",
        "tool_mismatch",
        "artifact_claim_false",
        "artifact_empty",
        "error_masked",
        "entity_missing",
    ]
    explanation: str
    retry_hint: Optional[str] = None


def check(
    spec: Spec,
    tool_output: str,
    responder_text: str,
    plot_results: list,
    table_results: list,
) -> GatekeeperResult:
    """Run checks in priority order; return first block or pass."""
    lowered_responder = responder_text.lower()

    for check_fn in (
        _check_tool_mismatch,
        _check_artifact_claim_false,
        _check_artifact_empty,
        _check_error_masked,
        _check_entity_missing,
    ):
        result = check_fn(
            spec=spec,
            tool_output=tool_output,
            lowered_responder=lowered_responder,
            plot_results=plot_results,
            table_results=table_results,
        )
        if result is not None:
            return result

    return GatekeeperResult(
        status="pass",
        reason="pass",
        explanation="All deterministic checks passed.",
        retry_hint=None,
    )


def _check_tool_mismatch(
    spec: Spec,
    lowered_responder: str,
    **_: Any,
) -> Optional[GatekeeperResult]:
    """Block when responder names a registered tool different from spec.tool_name.

    Uses substring match — imperfect when tool names are substrings of each
    other (e.g. "plot" vs "plot_umap"); acceptable for MVP since registered
    tool names are distinct multi-word identifiers.
    """
    spec_tool_lower = spec.tool_name.lower()
    registered = get_tool_names()

    spec_mentioned = spec_tool_lower in lowered_responder
    others_mentioned = [
        name for name in registered
        if name != spec.tool_name and name.lower() in lowered_responder
    ]

    if others_mentioned and not spec_mentioned:
        return GatekeeperResult(
            status="block",
            reason="tool_mismatch",
            explanation=(
                f"The response mentions tool(s) {others_mentioned!r} but the tool that actually "
                f"ran was {spec.tool_name!r}. The narrative is describing a different tool than "
                f"was dispatched."
            ),
            retry_hint=(
                f"Your previous draft mentioned {others_mentioned!r}, but the tool that actually "
                f"ran was '{spec.tool_name}'. Revise your narrative to name the correct tool."
            ),
        )
    return None


def _check_artifact_claim_false(
    lowered_responder: str,
    plot_results: list,
    table_results: list,
    **_: Any,
) -> Optional[GatekeeperResult]:
    """Block when responder claims artifact creation but no artifact exists."""
    if plot_results or table_results:
        return None

    matched = [p for p in _ARTIFACT_CLAIM_PHRASES if p in lowered_responder]
    if not matched:
        return None

    return GatekeeperResult(
        status="block",
        reason="artifact_claim_false",
        explanation=(
            f"The response claims artifact creation ({matched[0]!r}) but no plot or table was "
            f"actually produced. The narrative is describing output that doesn't exist."
        ),
        retry_hint=(
            f"Your previous draft included {matched[0]!r}, implying a plot or table was created. "
            f"No artifact was produced. Revise to avoid claiming artifact creation, OR clearly "
            f"state that the tool didn't produce one."
        ),
    )


def _check_artifact_empty(
    plot_results: list,
    table_results: list,
    **_: Any,
) -> Optional[GatekeeperResult]:
    """Block when a produced artifact is zero-size / zero-rows.

    Uses getattr with a small fallback list for attribute-name variation.
    If no known attribute is found, skips the check (conservative — no false block).
    """
    for plot in plot_results:
        img_bytes: Any = None
        for attr in ("img_bytes", "image_bytes", "png_bytes", "bytes"):
            img_bytes = getattr(plot, attr, None)
            if img_bytes is not None:
                break
        if img_bytes is not None and hasattr(img_bytes, "__len__") and len(img_bytes) == 0:
            return GatekeeperResult(
                status="block",
                reason="artifact_empty",
                explanation=(
                    "A plot result was produced but its image bytes are empty. The tool "
                    "reported success but the artifact is zero-size."
                ),
                retry_hint=None,  # not retry-able — responder redraft won't fix empty bytes
            )

    for table in table_results:
        df: Any = None
        for attr in ("dataframe", "df", "data"):
            df = getattr(table, attr, None)
            if df is not None:
                break
        if df is not None:
            try:
                n_rows = int(df.shape[0])
            except (AttributeError, TypeError, ValueError):
                continue
            if n_rows == 0:
                return GatekeeperResult(
                    status="block",
                    reason="artifact_empty",
                    explanation=(
                        "A table result was produced but has zero rows. The tool reported "
                        "success but the table is empty."
                    ),
                    retry_hint=None,
                )

    return None


def _check_error_masked(
    tool_output: str,
    lowered_responder: str,
    **_: Any,
) -> Optional[GatekeeperResult]:
    """Block when tool_output starts with 'Error:' but responder doesn't acknowledge."""
    if not tool_output.lstrip().startswith("Error:"):
        return None

    if any(w in lowered_responder for w in _ERROR_INDICATORS):
        return None

    return GatekeeperResult(
        status="block",
        reason="error_masked",
        explanation=(
            "The tool returned an error (output starts with 'Error:') but the response doesn't "
            "acknowledge it. The user would be misled about what happened."
        ),
        retry_hint=(
            "The tool's output is an error message (starts with 'Error:'). Your previous draft "
            "didn't acknowledge this. Rewrite to clearly state the tool errored and suggest a "
            "next step to the user."
        ),
    )


def _check_entity_missing(
    spec: Spec,
    tool_output: str,
    **_: Any,
) -> Optional[GatekeeperResult]:
    """Block when a field_type param's value is absent (case-insensitive) from tool_output.

    Only checks params with a field_type set (gene, cell_type, condition, obs_key, …).
    Non-field_type params (title, ncols, show_legend) are mechanical config — not
    expected to echo in tool_output.

    Deferred contract (T-040 plot-tool migration):
        Plot/table tools whose primary output is artifact bytes (PNG, dataframe)
        return only a terse status string in tool_output; the entity values
        live in the artifact, not the text. Substring-matching tool_output for
        those tools would false-fire even on correct calls. The fix is a
        structured `entities_acted_on` field on PlotResult/TableResult that
        this check reads when artifacts are present, falling back to the
        tool_output substring scan only when no artifact exists.

        Today the only @register tool is `dataset_overview`, which has zero
        field_type params, so this loop is a trivial no-op for every live
        call. The contract above will be locked when the first @register plot
        tool migrates in T-040; the existing test_gatekeeper.py
        `TestEntityMissing` cases all use empty plot_results / table_results
        and remain valid as the no-artifact branch.
    """
    entry = get_tool(spec.tool_name)
    if entry is None:
        return None

    lowered_output = tool_output.lower()
    missing: list[str] = []

    for p in entry.params:
        if p.field_type is None:
            continue
        value = spec.params.get(p.name)
        if value is None:
            continue
        values = value if isinstance(value, list) else [value]
        for v in values:
            v_str = str(v)
            if v_str and v_str.lower() not in lowered_output:
                missing.append(v_str)

    if not missing:
        return None

    return GatekeeperResult(
        status="block",
        reason="entity_missing",
        explanation=(
            f"The user asked about {missing!r}, but the tool output doesn't reference any of "
            f"these entities. The tool may have ignored the parameters or produced output too "
            f"terse to confirm they were acted on."
        ),
        retry_hint=(
            f"The tool output doesn't reference the user's requested entities: {missing!r}. "
            f"Either the tool ignored them or its output is too terse to confirm. Revise the "
            f"narrative to either explicitly name the entities from tool_output OR acknowledge "
            f"that the tool output doesn't confirm they were used."
        ),
    )
