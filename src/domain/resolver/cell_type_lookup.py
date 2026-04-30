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


_CELL_TYPE_COL_CANDIDATES: tuple[str, ...] = (
    "cell_type", "celltype", "annotation", "label",
)


def _resolve_cell_type_col(adata: "AnnData", obs_col: str | None) -> str | None:
    """Resolve which obs column holds cell-type labels.

    Caller-supplied ``obs_col`` wins when present in adata.obs. Otherwise
    walks the candidate list (matches the obs_column resolver's semantic
    map for cell-type columns).
    """
    if obs_col is not None and obs_col in adata.obs.columns:
        return obs_col
    for cand in _CELL_TYPE_COL_CANDIDATES:
        if cand in adata.obs.columns:
            return cand
    return None


def lookup_cell_type_name(
    adata: "AnnData",
    raw: str,
    obs_col: str | None = None,
) -> CellTypeLookupResult:
    """Resolve a cell-type string against the cell-type obs column's values.

    If ``obs_col`` is None or absent, auto-discovers from the candidate chain
    (``cell_type``, ``celltype``, ``annotation``, ``label``). This matches the
    obs_column resolver's semantic map and unblocks subset_value resolution
    on datasets where the cell-type column isn't literally named ``cell_type``.
    """
    resolved_col = _resolve_cell_type_col(adata, obs_col)
    if resolved_col is None:
        return CellTypeLookupResult(
            matched=False, input_name=raw, resolved_name=None, strategy=None,
            message=(
                f"No cell-type column found in adata.obs. Looked for: "
                f"{', '.join(_CELL_TYPE_COL_CANDIDATES)}."
            ),
        )
    obs_col = resolved_col

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

    # Substring + simple plural-stem match for broad terms (e.g.
    # "cardiomyocytes" → ["Early cardiomyocyte", "Ventricular cardiomyocyte"]).
    # Min length 3 avoids junk matches from very short queries.
    raw_lower = raw.strip().lower()
    if len(raw_lower) >= 3:
        raw_stem = raw_lower[:-1] if raw_lower.endswith("s") else raw_lower
        substring_hits = sorted({
            n for n in all_names
            if raw_lower in n.lower() or (raw_stem != raw_lower and raw_stem in n.lower())
        })
        if len(substring_hits) == 1:
            return CellTypeLookupResult(
                matched=True, input_name=raw, resolved_name=substring_hits[0],
                strategy="substring",
                message=f"Cell type '{raw}' resolved to '{substring_hits[0]}' (substring match).",
            )
        if len(substring_hits) > 1:
            return CellTypeLookupResult(
                matched=False, input_name=raw, resolved_name=None,
                strategy="broad_match_ambiguous", candidates=substring_hits,
                message=(
                    f"'{raw}' matches {len(substring_hits)} cell types: "
                    f"{', '.join(substring_hits)}."
                ),
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
