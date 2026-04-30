"""Entity canonicalization orchestrator — walks spec.pre_canonical_params, dispatches by field_type, builds spec.params additively."""

from __future__ import annotations

from typing import Any, Optional, TYPE_CHECKING

from src.core.spec import Canonicalization, Spec
from src.domain.resolver.cell_type_lookup import lookup_cell_type_name
from src.domain.resolver.column_classifier import Role, get_column_role
from src.domain.resolver.condition_lookup import lookup_condition_name
from src.domain.resolver.gene_lookup import lookup_gene_name
from src.domain.resolver.obs_column_lookup import lookup_obs_column
from src.spec_validation.result import Issue
from src.core.registry import get_tool

if TYPE_CHECKING:
    from anndata import AnnData


# <dependent_param>: <sibling_param_whose_resolved_value_names_the_column>.
# When resolving the dependent param, its value is looked up against
# adata.obs[<resolved sibling>].unique() via _dispatch_subset_value.
# Generalize to per-tool @register declarations only when 3+ pair shapes appear.
_CROSS_FIELD_SIBLINGS: dict[str, str] = {
    "subset_value": "subset_key",
    "group1": "groupby",
    "group2": "groupby",
}

# Pair-aware sibling correction registry: when both dependent params in the
# tuple-key fail to resolve in the column named by the sibling param (value),
# but BOTH succeed in some other obs column, swap the sibling param to that
# column. Closes the LLM-picks-wrong-groupby class of failures (Yalu §5
# DAI-18 prompt: "DE of EC between WT and CKO" → LLM extracts groupby=
# cell_type but WT/CKO live in the condition column).
_PAIR_AWARE_CORRECTIONS: dict[tuple[str, str], str] = {
    ("group1", "group2"): "groupby",
}


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

    # Cross-field dependency table: <param_name> resolves against the column
    # named in <sibling_param_name>. Sibling-target params resolve first so
    # the dependent param's lookup can route by the sibling column's role.
    # Yalu scenarios served:
    #   subset_value ↔ subset_key  — §3B.6 condition subset, §3D.3 cell-type subset, ...
    #   group1/group2 ↔ groupby    — §5.1/§5.2 run_de + volcano_plot
    sibling_targets = set(_CROSS_FIELD_SIBLINGS.values())
    sorted_items = sorted(
        spec.pre_canonical_params.items(),
        key=lambda kv: 0 if kv[0] in sibling_targets else 1,
    )

    correction_applied = False

    for param_name, value in sorted_items:
        # At the transition from sibling-target params to dependent params,
        # apply pair-aware corrections. This swaps a wrong-groupby into a
        # right-groupby BEFORE the dependent params (group1/group2) dispatch,
        # so they resolve against the corrected column.
        if not correction_applied and param_name not in sibling_targets:
            _apply_pair_aware_corrections(
                spec.pre_canonical_params, new_params, canonicalizations, adata,
            )
            correction_applied = True

        field_type = field_types.get(param_name)
        if field_type is None:
            new_params[param_name] = value
            continue

        # Empty string is the LLM's way of saying "not provided" — skip
        # canonicalization so the tool body's default fires, rather than
        # forcing a needs_input clarification turn over an empty value.
        if isinstance(value, str) and not value.strip():
            continue

        # Cross-field dispatch: dependent param's lookup uses the sibling
        # param's resolved column (caught only after the sibling resolved
        # into new_params via the sibling_targets ordering above).
        subset_key_col: Optional[str] = None
        sibling_param = _CROSS_FIELD_SIBLINGS.get(param_name)
        if sibling_param is not None:
            sib_value = new_params.get(sibling_param)
            if isinstance(sib_value, str) and sib_value:
                subset_key_col = sib_value

        if isinstance(value, list):
            resolved_list, cans, list_issues = _resolve_list(
                param_name, field_type, value, adata,
                subset_key_col=subset_key_col,
            )
            new_params[param_name] = resolved_list
            canonicalizations.extend(cans)
            issues.extend(list_issues)
        else:
            canonical, suggestions, reason, context = _dispatch(
                field_type, str(value), adata,
                subset_key_col=subset_key_col,
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
    field_type: str, raw: str, adata: "AnnData",
    *, subset_key_col: Optional[str] = None,
) -> tuple[Optional[str], list[str], Optional[str], Optional[dict]]:
    """Scalar dispatch to the lookup primitive for this field_type.

    ``subset_key_col`` is the resolved obs column that an earlier-resolved
    ``subset_key`` param canonicalized to. Set only when dispatching the
    sibling ``subset_value`` param; routes by the column's role so a
    ``condition``-column subset resolves against condition values, etc.
    Yalu §3B.6 (subset by condition column) depends on this.
    """
    if subset_key_col is not None:
        return _dispatch_subset_value(raw, adata, subset_key_col)

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

    if field_type == "obs_column":
        r = lookup_obs_column(adata, raw)
        if r.matched:
            return r.resolved_name, [], None, None
        reason = "ambiguous" if r.strategy == "substring_ambiguous" else "not_found"
        return None, r.candidates, reason, None

    return None, [], "wrong_type", None


def _dispatch_subset_value(
    raw: str, adata: "AnnData", subset_key_col: str,
) -> tuple[Optional[str], list[str], Optional[str], Optional[dict]]:
    """Resolve a subset_value against the column subset_key was resolved to.

    Routes by the column's classified role: CELL_TYPE → cell_type_lookup,
    CONDITION → condition_lookup (both with obs_col passed). Anything else
    (CLUSTERING / SAMPLE_ID / BATCH / OTHER / AMBIGUOUS / None) falls back
    to a direct case-insensitive match against ``adata.obs[col].unique()``,
    since no fuzzy lookup primitive applies.
    """
    role = get_column_role(adata, subset_key_col)

    if role == Role.CELL_TYPE:
        r = lookup_cell_type_name(adata, raw, obs_col=subset_key_col)
        if r.matched:
            return r.resolved_name, [], None, None
        reason = "ambiguous" if r.candidates else "not_found"
        return None, r.candidates, reason, None

    if role == Role.CONDITION:
        r = lookup_condition_name(adata, raw, obs_col=subset_key_col)
        if r.matched:
            return r.resolved_name, [], None, {"obs_column": r.obs_column}
        reason = "ambiguous" if r.candidates else "not_found"
        return None, r.candidates, reason, None

    # Direct unique-value match (CLUSTERING / SAMPLE_ID / BATCH / OTHER /
    # AMBIGUOUS / None / role unknown). Yalu §3B.6 with orig.ident
    # (classified SAMPLE_ID) lands here.
    if subset_key_col not in adata.obs.columns:
        return None, [], "not_found", None
    uniques = [str(v) for v in adata.obs[subset_key_col].unique()]
    if raw in uniques:
        return raw, [], None, None
    raw_lower = raw.lower()
    for v in uniques:
        if v.lower() == raw_lower:
            return v, [], None, None
    return None, uniques[:20], "not_found", None


def _resolve_list(
    param_name: str,
    field_type: str,
    values: list[Any],
    adata: "AnnData",
    *, subset_key_col: Optional[str] = None,
) -> tuple[list[Any], list[Canonicalization], list[Issue]]:
    """Resolve each list element; collect per-element canonicalizations + issues."""
    resolved_list: list[Any] = list(values)
    cans: list[Canonicalization] = []
    issues: list[Issue] = []
    for idx, element in enumerate(values):
        canonical, suggestions, reason, context = _dispatch(
            field_type, str(element), adata,
            subset_key_col=subset_key_col,
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


def _find_column_containing_pair(
    g1_raw: str, g2_raw: str, adata: "AnnData", *, exclude: Optional[str] = None,
) -> Optional[str]:
    """Find the unique obs column where BOTH raw values appear. Returns None
    if 0 columns match OR 2+ columns match (ambiguous → no silent correction).
    Used by pair-aware groupby correction.
    """
    candidates: list[str] = []
    for col in adata.obs.columns:
        if col == exclude:
            continue
        col_vals = {str(v) for v in adata.obs[col].unique()}
        if g1_raw in col_vals and g2_raw in col_vals:
            candidates.append(col)
    return candidates[0] if len(candidates) == 1 else None


def _apply_pair_aware_corrections(
    pre_canonical: dict, new_params: dict,
    canonicalizations: list[Canonicalization], adata: "AnnData",
) -> None:
    """Swap a sibling-column param when its dependent-pair's values don't
    fit the LLM-extracted column but both fit another column unambiguously.
    Mutates ``new_params`` and appends to ``canonicalizations`` in-place.

    Yalu §5 DAI-18 case: "DE of EC between WT and CKO" → LLM extracts
    groupby=cell_type instead of the condition column. WT and CKO both
    fail in cell_type but both exist in orig.ident → swap groupby. Surfaces
    via Canonicalization so the responder discloses the correction (Q1=b
    visible).

    Asymmetric failures (one value found, the other not — likely a typo)
    and ambiguous multi-column matches are left to the standard error path.
    Pair-aware correction only fires on confident pair-symmetric mismatches.
    """
    for (dep1, dep2), sibling in _PAIR_AWARE_CORRECTIONS.items():
        g1_raw = pre_canonical.get(dep1)
        g2_raw = pre_canonical.get(dep2)
        sibling_col = new_params.get(sibling)
        if not (
            isinstance(g1_raw, str) and g1_raw
            and isinstance(g2_raw, str) and g2_raw
            and isinstance(sibling_col, str) and sibling_col
            and sibling_col in adata.obs.columns
        ):
            continue
        sibling_vals = {str(v) for v in adata.obs[sibling_col].unique()}
        g1_in_sibling = g1_raw in sibling_vals
        g2_in_sibling = g2_raw in sibling_vals
        if g1_in_sibling and g2_in_sibling:
            continue  # both fit current sibling — no correction needed
        if g1_in_sibling != g2_in_sibling:
            continue  # asymmetric failure — likely typo, not wrong-sibling
        better_col = _find_column_containing_pair(
            g1_raw, g2_raw, adata, exclude=sibling_col,
        )
        if better_col is None:
            continue  # ambiguous or no match — let standard error path run
        new_params[sibling] = better_col
        canonicalizations.append(Canonicalization(
            field=f"params.{sibling}",
            raw=sibling_col,
            canonical=better_col,
            context={
                "reason": "pair_aware_correction",
                "dependents": [dep1, dep2],
                "explanation": (
                    f"{dep1}={g1_raw!r} and {dep2}={g2_raw!r} both exist in "
                    f"{better_col!r} but not in {sibling_col!r} — swapped "
                    f"{sibling} to the column containing both values."
                ),
            },
        ))
