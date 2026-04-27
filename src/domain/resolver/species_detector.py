"""Upload-time species detector for single-cell datasets.

Writes adata.uns["nvwa_meta"]["species"] with a value (human | mouse | rat | None),
confidence, source, and optional candidates. Heuristic-only: no LLM. Ambiguous
outcomes surface via low confidence + candidates — dataset_overview() presents
them to the user.

Preserves any keys already present in adata.uns["nvwa_meta"] (e.g. schema written
by column_classifier.py running earlier or later in the load_dataset path).
"""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass
from typing import TYPE_CHECKING, Any, Optional

if TYPE_CHECKING:
    from anndata import AnnData

logger = logging.getLogger(__name__)


@dataclass
class SpeciesResult:
    value: Optional[str]                         # "human" | "mouse" | "rat" | None
    confidence: str                               # "high" | "low"
    source: str                                   # "declared" | "ensembl_prefix" | "gene_case" | "unknown"
    candidates: Optional[list[str]] = None


_SAMPLE_SIZE: int = 100

_MAJORITY_THRESHOLD: float = 0.7

# Ensembl ID prefix → canonical species
_ENSEMBL_PREFIXES: dict[str, str] = {
    "ENSG": "human",
    "ENSMUSG": "mouse",
    "ENSRNOG": "rat",
}

# Top-level adata.uns keys checked for an existing species declaration
_DECLARATION_KEYS: tuple[str, ...] = ("species", "organism")

# Canonical species set for fallthrough / candidates
_ALL_SPECIES: list[str] = ["human", "mouse", "rat"]

# Normalize known declared values (case-insensitive) → canonical species.
# Covers scientific names, common names, and common genome builds.
_DECLARATION_NORMALIZATION: dict[str, str] = {
    # Human
    "human": "human",
    "homo sapiens": "human",
    "hg38": "human",
    "hg19": "human",
    "grch38": "human",
    "grch37": "human",
    # Mouse
    "mouse": "mouse",
    "mus musculus": "mouse",
    "mm10": "mouse",
    "mm39": "mouse",
    "grcm38": "mouse",
    "grcm39": "mouse",
    # Rat
    "rat": "rat",
    "rattus norvegicus": "rat",
    "rn6": "rat",
    "rn7": "rat",
    "mratbn7.2": "rat",
}


def detect_species(adata: "AnnData") -> None:
    """Detect species from adata; merge result into adata.uns['nvwa_meta']['species'].

    Idempotent. Preserves existing nvwa_meta keys outside this module's scope
    (schema, condition_cols, ambiguous_cols, user_resolutions).
    """
    result = (
        _from_declaration(adata)
        or _from_ensembl_prefix(adata)
        or _from_gene_case(adata)
        or SpeciesResult(
            value=None,
            confidence="low",
            source="unknown",
            candidates=list(_ALL_SPECIES),
        )
    )

    merged = _existing_meta(adata)
    merged["species"] = asdict(result)
    adata.uns["nvwa_meta"] = merged
    logger.info(
        "species_detector: value=%s, confidence=%s, source=%s",
        result.value, result.confidence, result.source,
    )


def _from_declaration(adata: "AnnData") -> Optional[SpeciesResult]:
    """Use adata.uns['species'] or adata.uns['organism'] if set to a recognized value."""
    uns: Any = getattr(adata, "uns", {}) or {}
    if not hasattr(uns, "get"):
        return None
    for key in _DECLARATION_KEYS:
        raw = uns.get(key)
        if raw is None:
            continue
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8", errors="ignore")
        if not isinstance(raw, str):
            continue
        normalized = _DECLARATION_NORMALIZATION.get(raw.strip().lower())
        if normalized:
            return SpeciesResult(value=normalized, confidence="high", source="declared")
    return None


def _from_ensembl_prefix(adata: "AnnData") -> Optional[SpeciesResult]:
    """Count Ensembl prefix matches across a sample of var_names; majority wins."""
    var_names = _sample_var_names(adata)
    if not var_names:
        return None
    counts: dict[str, int] = {species: 0 for species in _ENSEMBL_PREFIXES.values()}
    # Iterate prefixes longest-first so ENSMUSG wins before ENSG on mouse IDs.
    ordered = sorted(_ENSEMBL_PREFIXES.items(), key=lambda kv: -len(kv[0]))
    for name in var_names:
        for prefix, species in ordered:
            if name.startswith(prefix):
                counts[species] += 1
                break
    total = len(var_names)
    for species, hits in counts.items():
        if hits / total >= _MAJORITY_THRESHOLD:
            return SpeciesResult(value=species, confidence="high", source="ensembl_prefix")
    return None


def _from_gene_case(adata: "AnnData") -> Optional[SpeciesResult]:
    """Majority-vote by gene-symbol case convention. Mouse/rat share Title case."""
    var_names = _sample_var_names(adata)
    if not var_names:
        return None
    all_upper = sum(1 for n in var_names if _is_symbol_all_upper(n))
    title_case = sum(1 for n in var_names if _is_symbol_title_case(n))
    total = len(var_names)
    if all_upper / total >= _MAJORITY_THRESHOLD:
        return SpeciesResult(value="human", confidence="high", source="gene_case")
    if title_case / total >= _MAJORITY_THRESHOLD:
        return SpeciesResult(
            value=None,
            confidence="low",
            source="gene_case",
            candidates=["mouse", "rat"],
        )
    return None


def _sample_var_names(adata: "AnnData") -> list[str]:
    """Return up to _SAMPLE_SIZE var_names as non-empty strings."""
    var_names = getattr(adata, "var_names", None)
    if var_names is None:
        return []
    try:
        raw = list(var_names[:_SAMPLE_SIZE])
    except Exception:
        return []
    return [str(n) for n in raw if n]


def _is_symbol_all_upper(name: str) -> bool:
    """All alphabetic chars are uppercase (e.g. GAPDH, IL2RA, MT-CO1). Digits/hyphens ignored."""
    has_letter = False
    for ch in name:
        if ch.isalpha():
            has_letter = True
            if not ch.isupper():
                return False
    return has_letter


def _is_symbol_title_case(name: str) -> bool:
    """Leading uppercase letter followed by at least one lowercase letter (e.g. Gapdh, Trp53)."""
    has_upper_first = False
    has_lower = False
    for ch in name:
        if not ch.isalpha():
            continue
        if not has_upper_first:
            if not ch.isupper():
                return False
            has_upper_first = True
        elif ch.islower():
            has_lower = True
    return has_upper_first and has_lower


def _existing_meta(adata: "AnnData") -> dict[str, Any]:
    """Return a mutable copy of adata.uns['nvwa_meta'] if it's a dict, else empty."""
    uns: Any = getattr(adata, "uns", {}) or {}
    meta = uns.get("nvwa_meta", {}) if hasattr(uns, "get") else {}
    return dict(meta) if isinstance(meta, dict) else {}
