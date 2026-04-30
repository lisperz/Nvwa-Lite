"""Condition resolution: scan declared condition columns in adata.uns."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData


@dataclass
class ConditionLookupResult:
    matched: bool
    input_name: str
    resolved_name: str | None
    obs_column: str | None
    strategy: str | None
    candidates: list[str] = field(default_factory=list)  # "column:value" when ambiguous
    message: str = ""


def lookup_condition_name(
    adata: "AnnData", raw: str, obs_col: str | None = None,
) -> ConditionLookupResult:
    """Resolve a condition value against adata.obs.

    If ``obs_col`` is provided AND present in adata.obs, search only that
    column (used by the resolver's subset_value cross-field dispatch when
    subset_key has been resolved to a CONDITION-role column). Otherwise
    scan all declared condition columns from adata.uns["nvwa_meta"].
    """
    if obs_col is not None and obs_col in adata.obs.columns:
        condition_cols = [obs_col]
    else:
        condition_cols = _get_condition_cols(adata)

    exact_matches: list[tuple[str, str]] = []
    for col in condition_cols:
        if col not in adata.obs.columns:
            continue
        if raw in [str(x) for x in adata.obs[col].unique()]:
            exact_matches.append((col, raw))

    if len(exact_matches) == 1:
        col, val = exact_matches[0]
        return ConditionLookupResult(
            matched=True, input_name=raw, resolved_name=val,
            obs_column=col, strategy="exact",
            message=f"Condition '{raw}' found in column '{col}' (exact match).",
        )
    if len(exact_matches) > 1:
        return ConditionLookupResult(
            matched=False, input_name=raw, resolved_name=None,
            obs_column=None, strategy="exact",
            candidates=[f"{c}:{v}" for c, v in exact_matches],
            message=f"'{raw}' is ambiguous — matches multiple condition columns.",
        )

    raw_n = _normalize(raw)
    norm_matches: list[tuple[str, str]] = []
    for col in condition_cols:
        if col not in adata.obs.columns:
            continue
        for v in [str(x) for x in adata.obs[col].unique()]:
            if _normalize(v) == raw_n:
                norm_matches.append((col, v))
                break

    if len(norm_matches) == 1:
        col, val = norm_matches[0]
        return ConditionLookupResult(
            matched=True, input_name=raw, resolved_name=val,
            obs_column=col, strategy="normalized",
            message=f"Condition '{raw}' resolved to '{col}:{val}' (normalized match).",
        )
    if len(norm_matches) > 1:
        return ConditionLookupResult(
            matched=False, input_name=raw, resolved_name=None,
            obs_column=None, strategy="normalized",
            candidates=[f"{c}:{v}" for c, v in norm_matches],
            message=f"'{raw}' is ambiguous — matches multiple condition columns.",
        )

    return ConditionLookupResult(
        matched=False, input_name=raw, resolved_name=None,
        obs_column=None, strategy=None,
        message=f"Condition '{raw}' not found in any of: {condition_cols}.",
    )


def _get_condition_cols(adata: "AnnData") -> list[str]:
    uns = getattr(adata, "uns", {}) or {}
    meta = uns.get("nvwa_meta", {}) if hasattr(uns, "get") else {}
    return meta.get("condition_cols", ["condition"])


def _normalize(s: str) -> str:
    return re.sub(r"[\s_\-]+", "_", s.strip().lower())
