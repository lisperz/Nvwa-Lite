"""Obs-column-name resolution: map natural-language column references to actual
``adata.obs`` columns.

Closes the gap where the spec-pipeline extractor passes the user's literal
phrase (e.g. "leiden cluster", "cell type column") into params like
``color_by`` / ``split_by`` / ``subset_key`` and no resolver branch
canonicalizes it. The legacy LangChain path solved this at the prompt layer
(prompts.py:26-31, 109); the spec pipeline solves it here as a structured
field-type lookup.

Strategy chain: exact → normalized → semantic → substring. Multi-match
returns ``matched=False`` with ``candidates``, which the orchestrating
resolver translates into an ``Issue(reason="ambiguous")``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData


# Semantic map of dimension keywords to ordered candidate column lists.
# Mirrors src/agent/prompts.py:109 ("Dimension Resolution Rules"). More-specific
# keywords are listed first so substring iteration finds them before generic ones.
_SEMANTIC_MAP: list[tuple[str, list[str]]] = [
    ("seurat_clusters", ["seurat_clusters"]),
    ("seurat clusters", ["seurat_clusters"]),
    ("leiden", ["leiden"]),
    ("louvain", ["louvain"]),
    ("cluster", ["leiden", "louvain", "seurat_clusters", "cluster"]),
    ("cell type", ["cell_type", "celltype", "annotation", "label"]),
    ("celltype", ["cell_type", "celltype", "annotation", "label"]),
    ("cell identity", ["cell_type", "celltype", "annotation", "label", "ident"]),
    ("identity", ["cell_type", "celltype", "annotation", "label", "ident"]),
    ("annotation", ["cell_type", "celltype", "annotation", "label"]),
    ("label", ["cell_type", "celltype", "annotation", "label"]),
    ("condition", ["orig.ident", "condition", "sample", "batch"]),
    ("sample", ["orig.ident", "condition", "sample", "batch"]),
    ("treatment", ["orig.ident", "condition", "sample", "batch"]),
    ("batch", ["orig.ident", "condition", "sample", "batch"]),
]


@dataclass
class ObsColumnLookupResult:
    matched: bool
    input_name: str
    resolved_name: str | None
    strategy: str | None
    candidates: list[str] = field(default_factory=list)
    message: str = ""


def lookup_obs_column(adata: "AnnData", raw: str) -> ObsColumnLookupResult:
    """Resolve ``raw`` against ``adata.obs.columns`` via exact / normalized /
    semantic / substring strategies.
    """
    all_cols = [str(c) for c in adata.obs.columns]

    if not (raw or "").strip():
        return ObsColumnLookupResult(
            matched=False, input_name=raw, resolved_name=None, strategy=None,
            candidates=all_cols,
            message="obs column query is empty or whitespace-only.",
        )

    if raw in all_cols:
        return ObsColumnLookupResult(
            matched=True, input_name=raw, resolved_name=raw, strategy="exact",
            message=f"obs column '{raw}' found (exact match).",
        )

    raw_norm = _normalize(raw)
    norm_map = {_normalize(c): c for c in all_cols}
    if raw_norm in norm_map:
        resolved = norm_map[raw_norm]
        return ObsColumnLookupResult(
            matched=True, input_name=raw, resolved_name=resolved,
            strategy="normalized",
            message=f"obs column '{raw}' resolved to '{resolved}' (normalized match).",
        )

    raw_lower = raw.strip().lower()
    for keyword, candidates in _SEMANTIC_MAP:
        if keyword in raw_lower or _normalize(keyword) in raw_norm:
            for cand in candidates:
                if cand in all_cols:
                    return ObsColumnLookupResult(
                        matched=True, input_name=raw, resolved_name=cand,
                        strategy="semantic",
                        message=(
                            f"obs column '{raw}' resolved to '{cand}' "
                            f"(semantic match for '{keyword}')."
                        ),
                    )

    if len(raw_lower) >= 3:
        substring_hits = sorted({
            c for c in all_cols
            if raw_lower in c.lower() or raw_norm in _normalize(c)
        })
        if len(substring_hits) == 1:
            return ObsColumnLookupResult(
                matched=True, input_name=raw, resolved_name=substring_hits[0],
                strategy="substring",
                message=(
                    f"obs column '{raw}' resolved to '{substring_hits[0]}' "
                    f"(substring match)."
                ),
            )
        if len(substring_hits) > 1:
            return ObsColumnLookupResult(
                matched=False, input_name=raw, resolved_name=None,
                strategy="substring_ambiguous", candidates=substring_hits,
                message=(
                    f"'{raw}' matches {len(substring_hits)} obs columns: "
                    f"{', '.join(substring_hits)}."
                ),
            )

    return ObsColumnLookupResult(
        matched=False, input_name=raw, resolved_name=None, strategy=None,
        candidates=all_cols,
        message=(
            f"'{raw}' not found in obs columns. Available: {', '.join(all_cols)}."
        ),
    )


def _normalize(s: str) -> str:
    """Lowercase, strip ' column'/'col' suffix, drop all non-alphanumerics.

    Stripping non-alphanumerics (rather than collapsing to ``_``) bridges
    underscore-vs-no-separator (``cell_type`` ≡ ``celltype``), camelCase
    vs snake_case, and dot.notation variations.
    """
    s = s.strip().lower()
    for suffix in (" column", " col", "_column", "_col"):
        if s.endswith(suffix):
            s = s[: -len(suffix)]
            break
    return re.sub(r"[^a-z0-9]+", "", s)
