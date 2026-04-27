"""Entity canonicalization orchestrator — walks spec.pre_canonical_params, dispatches by field_type, builds spec.params additively."""

from __future__ import annotations

from typing import Any, Optional, TYPE_CHECKING

from src.core.spec import Canonicalization, Spec
from src.domain.resolver.cell_type_lookup import lookup_cell_type_name
from src.domain.resolver.condition_lookup import lookup_condition_name
from src.domain.resolver.gene_lookup import lookup_gene_name
from src.spec_validation.result import Issue
from src.tools.registry import get_tool

if TYPE_CHECKING:
    from anndata import AnnData


def resolve(spec: Spec, adata: "AnnData") -> tuple[Spec, list[Issue]]:
    """Build canonical spec.params from spec.pre_canonical_params; return new Spec + unresolved issues.

    Additive pattern: pre_canonical_params is the extractor's output (untouched here);
    params is built fresh — pass-through for fields with no field_type, canonicalized
    for fields with a registered field_type. Failed canonicalizations leave the field
    out of params and emit an Issue.
    """
    entry = get_tool(spec.tool_name)
    if entry is None:
        return spec, []  # unknown tool — validator catches this

    field_types = {p.name: p.field_type for p in entry.params if p.field_type}

    new_params: dict[str, Any] = {}
    canonicalizations: list[Canonicalization] = []
    issues: list[Issue] = []

    for param_name, value in spec.pre_canonical_params.items():
        field_type = field_types.get(param_name)
        if field_type is None:
            new_params[param_name] = value
            continue

        if isinstance(value, list):
            resolved_list, cans, list_issues = _resolve_list(
                param_name, field_type, value, adata
            )
            new_params[param_name] = resolved_list
            canonicalizations.extend(cans)
            issues.extend(list_issues)
        else:
            canonical, suggestions, reason, context = _dispatch(
                field_type, str(value), adata
            )
            if reason is None and canonical is not None:
                new_params[param_name] = canonical
                if canonical != value or context is not None:
                    canonicalizations.append(Canonicalization(
                        field=f"params.{param_name}",
                        raw=str(value),
                        canonical=str(canonical),
                        context=context,
                    ))
            else:
                if reason == "not_found":
                    other = _cross_field_check(str(value), adata, excluding=field_type)
                    if other:
                        issues.append(Issue(
                            field=f"params.{param_name}",
                            reason="wrong_field",
                            suggestions=[f"{other}={value}"],
                        ))
                        continue
                issues.append(Issue(
                    field=f"params.{param_name}",
                    reason=reason or "not_found",
                    suggestions=suggestions,
                ))

    updated = spec.model_copy(update={
        "params": new_params,
        "canonicalizations_applied": spec.canonicalizations_applied + canonicalizations,
    })
    return updated, issues


def _dispatch(
    field_type: str, raw: str, adata: "AnnData"
) -> tuple[Optional[str], list[str], Optional[str], Optional[dict]]:
    """Scalar dispatch to the lookup primitive for this field_type."""
    if field_type == "gene":
        r = lookup_gene_name(adata, raw)
        if r.matched:
            return r.resolved_name, [], None, None
        return None, r.candidates, "not_found", None

    if field_type == "cell_type":
        r = lookup_cell_type_name(adata, raw)
        if r.matched:
            return r.resolved_name, [], None, None
        reason = "ambiguous" if r.candidates else "not_found"
        return None, r.candidates, reason, None

    if field_type == "condition":
        r = lookup_condition_name(adata, raw)
        if r.matched:
            return r.resolved_name, [], None, {"obs_column": r.obs_column}
        reason = "ambiguous" if r.candidates else "not_found"
        return None, r.candidates, reason, None

    return None, [], "wrong_type", None


def _resolve_list(
    param_name: str,
    field_type: str,
    values: list[Any],
    adata: "AnnData",
) -> tuple[list[Any], list[Canonicalization], list[Issue]]:
    """Resolve each list element; collect per-element canonicalizations + issues."""
    resolved_list: list[Any] = list(values)
    cans: list[Canonicalization] = []
    issues: list[Issue] = []
    for idx, element in enumerate(values):
        canonical, suggestions, reason, context = _dispatch(
            field_type, str(element), adata
        )
        if reason is None and canonical is not None:
            if canonical != element or context is not None:
                resolved_list[idx] = canonical
                cans.append(Canonicalization(
                    field=f"params.{param_name}[{idx}]",
                    raw=str(element),
                    canonical=str(canonical),
                    context=context,
                ))
        else:
            issues.append(Issue(
                field=f"params.{param_name}[{idx}]",
                reason=reason or "not_found",
                suggestions=suggestions,
            ))
    return resolved_list, cans, issues


def _cross_field_check(
    raw: str, adata: "AnnData", excluding: str
) -> Optional[str]:
    """Return the field_type where raw resolves, if not in the excluded category."""
    if excluding != "gene" and lookup_gene_name(adata, raw).matched:
        return "gene"
    if excluding != "cell_type" and lookup_cell_type_name(adata, raw).matched:
        return "cell_type"
    if excluding != "condition" and lookup_condition_name(adata, raw).matched:
        return "condition"
    return None
