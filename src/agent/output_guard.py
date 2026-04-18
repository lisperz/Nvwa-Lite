"""Output-honesty guard for fabricated artifact-creation claims.

T-029 hotfix. Detects assistant replies that claim a plot, table, or CSV was
produced when no tool was called in the turn, and enables a one-shot retry
with a corrective SystemMessage. If the LLM re-hallucinates on retry, the
agent surfaces a user-facing rephrase message.

# THROWAWAY: replaced by gatekeeper.check_output_honesty() in T-004.
# The signature `check_output_honesty(assistant_text, tool_calls) -> GuardResult`
# is pre-shaped for the T-004 gatekeeper so the swap at core.py is mechanical
# (phrase-match swaps out for LLM-as-judge or structured-response check).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass
class GuardResult:
    """Result of an output-honesty check.

    Attributes:
        triggered: True iff the assistant text claims artifact creation
            without a matching tool call.
        matched_phrase: The specific phrase that matched, or None.
        reason: Human-readable explanation (for logs + debugging).
    """

    triggered: bool
    matched_phrase: str | None
    reason: str


# Phrases that indicate a claim of artifact (plot / table / CSV) CREATION.
# Kept lowercase; matching is case-insensitive. Scoped tightly to creation
# verbs ("generated", "plotted", "created", "saved", "produced") to avoid
# false positives on legitimate commentary (e.g. "the plot shows X",
# "here's your answer"). Covers both passive ("has been ...") and active
# ("I've ...") forms. These substrings match for any artifact type — the
# noun (plot/table/CSV) need not appear in the phrase because the phrase
# alone is a claim of creation.
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


REPHRASE_FALLBACK_MESSAGE = (
    "I wasn't able to complete that action. Could you please rephrase your "
    "request? Try specifying the gene name, plot type, and grouping variable "
    "more explicitly (for example: \"feature plot of MKI67 grouped by condition\")."
)


def check_output_honesty(
    assistant_text: str,
    tool_calls: list[Any] | None,
) -> GuardResult:
    """Check whether the assistant claimed artifact creation without a tool call.

    Args:
        assistant_text: The LLM's assistant-text content for this iteration.
        tool_calls: The LLM's tool-call list for this iteration (empty/None if
            the model produced text only).

    Returns:
        GuardResult with triggered=True when a claim phrase appears in the
        text AND tool_calls is empty. Any non-empty tool_calls list short-
        circuits to triggered=False — the honest-with-tool case.
    """
    if tool_calls:
        return GuardResult(
            triggered=False,
            matched_phrase=None,
            reason="tool_calls non-empty; honest tool-call path",
        )

    if not assistant_text:
        return GuardResult(
            triggered=False,
            matched_phrase=None,
            reason="assistant_text empty",
        )

    lowered = assistant_text.lower()
    for phrase in _ARTIFACT_CLAIM_PHRASES:
        if phrase in lowered:
            return GuardResult(
                triggered=True,
                matched_phrase=phrase,
                reason=(
                    f"assistant_text contains artifact-claim phrase {phrase!r} "
                    "but tool_calls is empty"
                ),
            )

    return GuardResult(
        triggered=False,
        matched_phrase=None,
        reason="no artifact-claim phrase matched",
    )


def build_corrective_system_message(matched_phrase: str) -> str:
    """Build the corrective SystemMessage injected on guard trigger.

    Written to disrupt the pattern-match/cache behavior observed in the MKI67
    bug (10:58:46 hallucination vs 10:59:39 success at near-identical wording).
    Two load-bearing properties:

    - Interpolates the matched phrase so the model is forced to reason about
      its own prior text (novel token sequence, breaks deterministic repeat).
    - Explicitly instructs "do not repeat" to act as anti-cache signal.
    """
    return (
        f"Your previous reply claimed artifact creation ({matched_phrase!r}) "
        "but you did not call any tool in this turn. This is not allowed. "
        "You must now either (a) call the appropriate tool to fulfill the "
        "user's request, or (b) state plainly that you cannot and explain "
        "why. Do not repeat the previous reply."
    )
