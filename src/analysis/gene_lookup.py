"""Gene name lookup and resolution against a loaded AnnData dataset.

Provides four ordered strategies to resolve ambiguous gene names:
1. Exact match
2. Case-insensitive match
3. Species prefix strip/swap (Mm_/Hs_ prefix, UPPER/Title case toggle)
4. Fuzzy match via difflib
"""

from __future__ import annotations

import difflib
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData

_SPECIES_PREFIXES = ("Mm_", "Hs_", "mm_", "hs_")


@dataclass
class GeneLookupResult:
    matched: bool
    input_name: str
    resolved_name: str | None  # exact name in dataset, or None
    strategy: str | None       # "exact" | "case_insensitive" | "species_prefix" | "fuzzy"
    candidates: list[str] = field(default_factory=list)  # top fuzzy candidates when not matched
    message: str = ""


def _all_var_names(adata: AnnData) -> list[str]:
    # Restrict to adata.var_names only — downstream tools never read from raw,
    # so including raw.var_names would let lookup_gene return a name that the
    # next tool call cannot find, producing a confusing error after a "success".
    return list(adata.var_names)


def _try_exact(gene: str, names_set: set[str]) -> str | None:
    return gene if gene in names_set else None


def _try_case_insensitive(gene: str, names: list[str]) -> str | None:
    gene_lower = gene.lower()
    for name in names:
        if name.lower() == gene_lower:
            return name
    return None


def _try_species_prefix(gene: str, names_set: set[str]) -> str | None:
    """Try stripping species prefixes and toggling case on the stripped part.

    Only acts on inputs that carry a known prefix (Mm_/Hs_/mm_/hs_).
    Bare ALLCAPS inputs (e.g. IL2) are intentionally left to fall through to
    fuzzy matching — otherwise IL2 → Il2 could silently resolve to the wrong
    gene when the user actually meant IL2RA.
    """
    stripped: str | None = None
    for prefix in _SPECIES_PREFIXES:
        if gene.startswith(prefix):
            stripped = gene[len(prefix):]
            break

    if stripped is None:
        return None

    candidates: list[str] = [stripped]
    if stripped == stripped.upper() and len(stripped) > 1:
        # e.g. Mm_TNNT2 → stripped=TNNT2 → also try Tnnt2
        candidates.append(stripped[0].upper() + stripped[1:].lower())
    elif stripped[0].isupper() and stripped[1:] == stripped[1:].lower() and len(stripped) > 1:
        # e.g. Mm_Tnnt2 → stripped=Tnnt2 → also try TNNT2
        candidates.append(stripped.upper())

    for candidate in candidates:
        if candidate in names_set:
            return candidate
    return None


def _try_fuzzy(gene: str, names: list[str], cutoff: float = 0.6, n: int = 5) -> list[str]:
    return difflib.get_close_matches(gene, names, n=n, cutoff=cutoff)


def lookup_gene_name(adata: AnnData, gene: str) -> GeneLookupResult:
    """Resolve a gene name against the loaded dataset using four ordered strategies.

    Args:
        adata: The annotated data matrix.
        gene: Gene name to look up (may be a typo, wrong case, or alias).

    Returns:
        GeneLookupResult with matched=True and resolved_name if found,
        or matched=False with candidates list if unresolved.
    """
    all_names = _all_var_names(adata)
    names_set = set(all_names)

    # Strategy 1: exact
    if _try_exact(gene, names_set):
        return GeneLookupResult(
            matched=True,
            input_name=gene,
            resolved_name=gene,
            strategy="exact",
            message=f"Gene '{gene}' found in dataset (exact match).",
        )

    # Strategy 2: case-insensitive
    resolved = _try_case_insensitive(gene, all_names)
    if resolved:
        return GeneLookupResult(
            matched=True,
            input_name=gene,
            resolved_name=resolved,
            strategy="case_insensitive",
            message=f"Gene '{gene}' resolved to '{resolved}' (case-insensitive match). Use '{resolved}' for analysis.",
        )

    # Strategy 3: species prefix / case convention toggle
    resolved = _try_species_prefix(gene, names_set)
    if resolved:
        return GeneLookupResult(
            matched=True,
            input_name=gene,
            resolved_name=resolved,
            strategy="species_prefix",
            message=f"Gene '{gene}' resolved to '{resolved}' (species/case convention). Use '{resolved}' for analysis.",
        )

    # Strategy 4: fuzzy
    candidates = _try_fuzzy(gene, all_names)
    if candidates:
        candidate_str = ", ".join(candidates)
        return GeneLookupResult(
            matched=False,
            input_name=gene,
            resolved_name=None,
            strategy="fuzzy",
            candidates=candidates,
            message=(
                f"Gene '{gene}' not found. Closest matches in dataset: {candidate_str}. "
                f"Please confirm which gene you meant."
            ),
        )

    return GeneLookupResult(
        matched=False,
        input_name=gene,
        resolved_name=None,
        strategy=None,
        candidates=[],
        message=(
            f"Gene '{gene}' not found in this dataset ({len(all_names)} genes available). "
            f"No close matches found. Please check the gene name."
        ),
    )
