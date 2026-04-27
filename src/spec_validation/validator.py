"""Deterministic post-resolver validator — final gate before dispatch.

Role: verify the incoming spec is actionable by tools. Tools dispatched after a
passing validation should not need to check or flag warnings about their inputs.
Separate pipeline stage from the resolver: validator only runs on resolver-clean
specs (core.py short-circuits on resolver-emitted issues and never invokes the
validator in that case).

Validator is a pure judge — never mutates the input spec. Three checks:
  1. unknown_tool — spec.tool_name is in the registry
  2. missing      — every required param is present in spec.params
  3. wrong_type   — present params match the expected Python type (shallow)

Missing-optional-with-default is NOT an issue; Python applies the tool's own
signature default at dispatch time. Response-layer disclosure of which defaults
were applied is computed on demand (not carried on Spec).
"""

from __future__ import annotations

from typing import Any

from src.core.spec import Spec
from src.spec_validation.result import Issue, ValidationResult
from src.tools.registry import get_tool, get_tool_names


def validate(spec: Spec) -> tuple[Spec, ValidationResult]:
    """Check spec against the tool registry. Return the unchanged spec + result.

    Status:
      - "ok" if no issues
      - "needs_input" if any issue

    ("defaults_disclosed" is never emitted here — defaults are not validator's
    concern; response layer reports them separately if needed.)
    """
    entry = get_tool(spec.tool_name)
    if entry is None:
        return spec, ValidationResult(
            status="needs_input",
            stage="validator",
            issues=[Issue(
                field="tool_name",
                reason="unknown_tool",
                suggestions=get_tool_names(),
            )],
        )

    issues: list[Issue] = []

    # missing — required params absent
    for p in entry.params:
        if p.required and p.name not in spec.params:
            issues.append(Issue(field=f"params.{p.name}", reason="missing"))

    # wrong_type — present params match expected shallow type
    for p in entry.params:
        if p.name not in spec.params:
            continue
        if not _type_matches(spec.params[p.name], p.type):
            issues.append(Issue(
                field=f"params.{p.name}",
                reason="wrong_type",
                suggestions=[p.type],
            ))

    status = "needs_input" if issues else "ok"
    return spec, ValidationResult(status=status, stage="validator", issues=issues)


def _type_matches(value: Any, type_tag: str) -> bool:
    """Shallow type check against the registry's string type tag.

    Only rejects obvious mismatches on common Python builtins; defers complex
    or unknown tags to the extractor's JSON schema (the primary type enforcer).
    Intentionally forgiving to avoid false-positive wrong_type issues.
    """
    if type_tag == "Any":
        return True

    # Unwrap Optional[X] → X (or accept None)
    if type_tag.startswith("Optional[") and type_tag.endswith("]"):
        if value is None:
            return True
        return _type_matches(value, type_tag[len("Optional["):-1])

    # Strip generic subscript for shallow check: "list[str]" → "list"
    base = type_tag.split("[", 1)[0].strip()

    # Asymmetric on purpose: int params accept bool (subclass), but bool params
    # reject plain int. The belt-and-suspenders role (see module docstring) catches
    # callers who pass 0/1 where True/False was meant; the JSON schema is the
    # primary enforcer and emits bool-only at the extractor layer already.
    simple_types: dict[str, tuple[type, ...]] = {
        "str": (str,),
        "int": (int,),              # bool is an int subclass; accept intentionally
        "float": (int, float),      # ints are accepted where floats are expected
        "bool": (bool,),            # strict: int 0/1 is NOT accepted (asymmetric, intentional)
        "list": (list,),
        "dict": (dict,),
    }
    py_types = simple_types.get(base)
    if py_types is None:
        return True  # unknown tag — don't false-positive
    return isinstance(value, py_types)
