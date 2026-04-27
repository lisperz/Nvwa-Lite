"""Inspection tools — summarize / describe / check the loaded dataset.

First file populated under the new tool registry (`src/tools/registry.py`).
T-040 will migrate legacy inspection tools (dataset_info, check_data_status,
inspect_metadata, DE tables, cluster_mapping) into this file later.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from src.tools.registry import register

if TYPE_CHECKING:
    from anndata import AnnData


# User-facing labels for Role values surfaced in ambiguity questions
_ROLE_HUMAN: dict[str, str] = {
    "condition": "experimental condition",
    "sample_id": "sample ID",
    "batch": "batch label",
    "cell_type": "cell type",
    "clustering": "cluster assignment",
    "qc_metric": "QC metric",
}


@register(
    description=(
        "Top-level dataset summary only: overall cells, genes, species, condition "
        "columns, and cell-type column. Surface any columns where the heuristic is "
        "unsure so the user can clarify. Use this for prompts like \"what's in this "
        "data\" or \"give me a data overview\". Do NOT use for per-group / per-cluster "
        "/ per-condition counts or aggregations — those require a different tool."
    ),
)
def dataset_overview(adata: "AnnData") -> str:
    """Render a human-readable overview of the loaded dataset.

    Reads `adata.uns["nvwa_meta"]` (populated by column_classifier + species_detector
    at upload time) and `adata.n_obs / n_vars`. Produces a single string combining
    the summary with any ambiguity questions, per the combined-single-turn UX
    decision for the column-classification pipeline.
    """
    meta = _get_meta(adata)
    schema = meta.get("schema", {}) if isinstance(meta, dict) else {}

    lines: list[str] = [f"Dataset: {adata.n_obs:,} cells × {adata.n_vars:,} genes."]

    species_line = _species_line(meta.get("species"))
    if species_line:
        lines.append(species_line)

    condition_line = _conditions_line(meta.get("condition_cols", []), schema)
    if condition_line:
        lines.append(condition_line)

    cell_types_line = _cell_types_line(schema)
    if cell_types_line:
        lines.append(cell_types_line)

    questions = _ambiguity_questions(meta, schema)
    if questions:
        lines.append("")
        lines.append("I need your help on a couple of things before we dive in:")
        lines.extend(questions)

    return "\n".join(lines)


def _get_meta(adata: "AnnData") -> dict[str, Any]:
    uns: Any = getattr(adata, "uns", {}) or {}
    meta = uns.get("nvwa_meta", {}) if hasattr(uns, "get") else {}
    return meta if isinstance(meta, dict) else {}


def _species_line(species: Any) -> str:
    if not isinstance(species, dict):
        return ""
    value = species.get("value")
    if not value or species.get("confidence") != "high":
        return ""  # low-confidence species surfaces as an ambiguity question
    source_note = {
        "declared": "declared in dataset metadata",
        "ensembl_prefix": "detected from Ensembl prefixes",
        "gene_case": "detected from gene-symbol case",
    }.get(species.get("source", ""), "detected")
    return f"Species: {value} ({source_note})."


def _conditions_line(condition_cols: Any, schema: dict[str, dict]) -> str:
    if not isinstance(condition_cols, list) or not condition_cols:
        return ""
    parts: list[str] = []
    for col in condition_cols:
        entry = schema.get(col, {})
        values = entry.get("values")
        if values:
            parts.append(f"{col} ({', '.join(values)})")
        else:
            parts.append(f"{col} ({entry.get('n_unique', '?')} values)")
    return f"Conditions detected: {'; '.join(parts)}."


def _cell_types_line(schema: dict[str, dict]) -> str:
    ct_cols = [col for col, entry in schema.items() if entry.get("role") == "cell_type"]
    if not ct_cols:
        return ""
    col = ct_cols[0]
    entry = schema[col]
    n_unique = entry.get("n_unique", 0)
    values = entry.get("values")
    if values:
        if n_unique > len(values):
            preview = ", ".join(values[:5])
            return f"Cell types annotated: {n_unique} types ({preview}, …)."
        return f"Cell types annotated: {n_unique} types ({', '.join(values)})."
    return f"Cell types annotated: {n_unique} types in column '{col}'."


def _ambiguity_questions(meta: dict[str, Any], schema: dict[str, dict]) -> list[str]:
    questions: list[str] = []

    # Column-level ambiguities from column_classifier
    for col in meta.get("ambiguous_cols", []) or []:
        entry = schema.get(col, {})
        values = entry.get("values") or []
        candidates = entry.get("candidates") or []
        n_unique = entry.get("n_unique", len(values))

        values_str = ", ".join(values[:6])
        if len(values) > 6:
            values_str += ", …"

        humanized = [_ROLE_HUMAN.get(c, c) for c in candidates]
        cand_str = _join_candidates(humanized)

        detail = f" ({values_str})" if values_str else ""
        questions.append(
            f"- Column '{col}' has {n_unique} unique values{detail}. Is this a {cand_str}?"
        )

    # Species ambiguity from species_detector
    species = meta.get("species")
    if isinstance(species, dict) and species.get("confidence") == "low":
        candidates = species.get("candidates") or []
        cand_str = _join_candidates(candidates)
        if species.get("source") == "gene_case":
            questions.append(
                f"- Species: gene symbols are Title-case, common to {cand_str}. Which is this dataset?"
            )
        else:
            questions.append(
                f"- Species: couldn't detect from the data. Is this {cand_str}?"
            )

    return questions


def _join_candidates(items: list[str]) -> str:
    if not items:
        return "something"
    if len(items) == 1:
        return items[0]
    if len(items) == 2:
        return f"{items[0]} or {items[1]}"
    return f"{', '.join(items[:-1])}, or {items[-1]}"
