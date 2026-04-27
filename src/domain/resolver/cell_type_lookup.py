"""Cell-type resolution: normalize + exact, with initials matching."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData


@dataclass
class CellTypeLookupResult:
    matched: bool
    input_name: str
    resolved_name: str | None
    strategy: str | None
    candidates: list[str] = field(default_factory=list)
    message: str = ""


def lookup_cell_type_name(
    adata: "AnnData",
    raw: str,
    obs_col: str = "cell_type",
) -> CellTypeLookupResult:
    """Resolve a cell-type string against adata.obs[obs_col] unique values."""
    if obs_col not in adata.obs.columns:
        return CellTypeLookupResult(
            matched=False, input_name=raw, resolved_name=None, strategy=None,
            message=f"Column '{obs_col}' not found in adata.obs.",
        )

    # Guard against degenerate queries: empty / whitespace-only raw can otherwise
    # exact-match a whitespace-only label stored in obs, which is semantically wrong.
    if not (raw or "").strip():
        return CellTypeLookupResult(
            matched=False, input_name=raw, resolved_name=None, strategy=None,
            message="Cell-type query is empty or whitespace-only.",
        )

    all_names = [str(x) for x in adata.obs[obs_col].unique()]

    if raw in all_names:
        return CellTypeLookupResult(
            matched=True, input_name=raw, resolved_name=raw, strategy="exact",
            message=f"Cell type '{raw}' found (exact match).",
        )

    raw_n = _normalize(raw)
    norm_map = {_normalize(n): n for n in all_names}
    if raw_n in norm_map:
        resolved = norm_map[raw_n]
        return CellTypeLookupResult(
            matched=True, input_name=raw, resolved_name=resolved, strategy="normalized",
            message=f"Cell type '{raw}' resolved to '{resolved}' (normalized match).",
        )

    if len(raw) <= 3 and raw.isupper():
        upper_matches = [n for n in all_names if _initials_uppercase(n) == raw]
        if len(upper_matches) == 1:
            return CellTypeLookupResult(
                matched=True, input_name=raw, resolved_name=upper_matches[0],
                strategy="initials_uppercase",
                message=f"Cell type '{raw}' resolved to '{upper_matches[0]}' (uppercase initials).",
            )
        if len(upper_matches) > 1:
            return CellTypeLookupResult(
                matched=False, input_name=raw, resolved_name=None,
                strategy="initials_uppercase", candidates=upper_matches,
                message=f"'{raw}' matches multiple cell types by uppercase initials.",
            )

        fl_matches = [n for n in all_names if _initials_first_letters(n) == raw]
        if len(fl_matches) == 1:
            return CellTypeLookupResult(
                matched=True, input_name=raw, resolved_name=fl_matches[0],
                strategy="initials_firstletter",
                message=f"Cell type '{raw}' resolved to '{fl_matches[0]}' (first-letter initials).",
            )
        if len(fl_matches) > 1:
            return CellTypeLookupResult(
                matched=False, input_name=raw, resolved_name=None,
                strategy="initials_firstletter", candidates=fl_matches,
                message=f"'{raw}' matches multiple cell types by first-letter initials.",
            )

    return CellTypeLookupResult(
        matched=False, input_name=raw, resolved_name=None, strategy=None,
        message=f"Cell type '{raw}' not found in '{obs_col}' ({len(all_names)} available).",
    )


def _normalize(s: str) -> str:
    return re.sub(r"[\s_\-]+", "_", s.strip().lower())


def _initials_uppercase(name: str) -> str:
    return "".join(c for c in name if c.isupper())


def _initials_first_letters(name: str) -> str:
    tokens = re.split(r"[\s_\-]+", name)
    return "".join(t[0].upper() for t in tokens if t)
