"""Tools for cell-type-restricted (subset) analysis.

Provides tools that filter the dataset to a specific cell type or group
before generating plots, enabling targeted comparisons like
"Show TNNT2 in cardiomyocytes across conditions".
"""

from __future__ import annotations

import difflib
import logging
from typing import TYPE_CHECKING

from langchain_core.tools import tool

from src.agent.viz_state import update_viz_state
from src.plotting.executor import PlotResult, plot_feature, plot_violin

if TYPE_CHECKING:
    from anndata import AnnData

logger = logging.getLogger(__name__)

_adata: AnnData | None = None
_plot_results: list[PlotResult] = []
_plot_generated_callback = None


def _get_adata() -> AnnData:
    if _adata is None:
        raise RuntimeError("No dataset loaded.")
    return _adata


def _get_cluster_key() -> str:
    adata = _get_adata()
    for key in ("leiden", "louvain", "seurat_clusters", "cluster", "clusters"):
        if key in adata.obs.columns:
            return key
    return "leiden"


def _store_and_return(result: PlotResult) -> str:
    _plot_results.append(result)
    if _plot_generated_callback:
        _plot_generated_callback()
    return f"Plot generated successfully.\nCode: {result.code}\n{result.message}"


def _resolve_subset_values(adata: AnnData, key: str, value: str) -> list[str]:
    """Resolve a subset value to one or more matching values.

    For broad terms like "cardiomyocytes", returns ALL matching values
    (e.g., ["Early cardiomyocyte", "Ventricular cardiomyocyte"]).

    Tries in order: exact match, case-insensitive exact, substring (all matches),
    then difflib fuzzy.
    """
    available = [str(v) for v in adata.obs[key].unique()]

    # Exact match — single value
    if value in available:
        return [value]

    # Case-insensitive exact match
    lower_map: dict[str, str] = {}
    for v in available:
        lower_map.setdefault(v.lower(), v)
    if value.lower() in lower_map:
        return [lower_map[value.lower()]]

    # Substring match — collect ALL values that contain the query (or vice versa)
    query_lower = value.lower()
    # Normalize common plurals: "cardiomyocytes" → "cardiomyocyte"
    query_stem = query_lower.rstrip("s")
    matches = []
    for av in available:
        av_lower = av.lower()
        if query_lower in av_lower or query_stem in av_lower:
            matches.append(av)
    if matches:
        return sorted(matches)

    # Fuzzy match via difflib — single best match
    fuzzy = difflib.get_close_matches(value, available, n=1, cutoff=0.6)
    if fuzzy:
        return [fuzzy[0]]

    return []


def _build_subset(adata: AnnData, key: str, values: list[str]):
    """Build a subset of adata matching any of the given values.

    Returns (adata_subset, label_str) where label_str describes what was matched.
    """
    mask = adata.obs[key].astype(str).isin(values)
    adata_sub = adata[mask].copy()
    if len(values) == 1:
        label = values[0]
    else:
        label = f"{', '.join(values)}"
    return adata_sub, label


def _report_condition_coverage(
    adata_full: AnnData, adata_sub: AnnData, groupby: str, subset_label: str,
) -> str:
    """Report which groups are present/missing after subsetting."""
    all_groups = set(str(v) for v in adata_full.obs[groupby].unique())
    present_groups = set(str(v) for v in adata_sub.obs[groupby].unique())
    missing_groups = sorted(all_groups - present_groups)

    if not missing_groups:
        return ""

    return (
        f"\n\nNote: {len(missing_groups)} of {len(all_groups)} {groupby} groups "
        f"have no {subset_label} cells and are not shown: {', '.join(missing_groups)}."
    )


def _compute_group_summary(
    adata_sub: AnnData, genes: str | list[str], groupby: str,
) -> str:
    """Compute per-group mean/median expression for quantitative comparison."""
    import numpy as np

    gene_name = genes if isinstance(genes, str) else genes[0]

    # Get expression values
    if gene_name in adata_sub.var_names:
        expr = adata_sub[:, gene_name].X
    elif adata_sub.raw is not None and gene_name in adata_sub.raw.var_names:
        expr = adata_sub.raw[:, gene_name].X
    else:
        return f"(Could not compute summary for {gene_name})"

    if hasattr(expr, "toarray"):
        expr = expr.toarray().flatten()
    else:
        expr = np.asarray(expr).flatten()

    groups = adata_sub.obs[groupby].values
    lines = []
    for group in sorted(adata_sub.obs[groupby].unique()):
        mask = groups == group
        vals = expr[mask]
        n = len(vals)
        mean_val = float(np.mean(vals))
        median_val = float(np.median(vals))
        lines.append(f"  {group}: n={n:,}, mean={mean_val:.3f}, median={median_val:.3f}")

    return "\n".join(lines)


@tool
def subset_violin_plot(
    genes: str,
    subset_key: str,
    subset_value: str,
    groupby: str = "",
) -> str:
    """Generate a violin plot for a SUBSET of cells filtered by a metadata column.

    Use this when the user wants to compare gene expression WITHIN a specific
    cell type across conditions, or within a specific condition across cell types.

    IMPORTANT — Broad cell type matching:
    The subset_value supports broad matching. For example:
    - "cardiomyocytes" will match ALL cardiomyocyte subtypes in the dataset
      (e.g., "Early cardiomyocyte", "Ventricular cardiomyocyte", etc.)
    - "T cells" will match all T cell subtypes
    This ensures broad biological categories include all relevant subtypes.

    Examples:
    - "Show TNNT2 in cardiomyocytes across conditions"
      → subset_key="cell_type", subset_value="cardiomyocyte", groupby="orig.ident"
    - "Compare CD3E between PA-IVS and Control in T cells"
      → subset_key="cell_type", subset_value="T cell", groupby="orig.ident"

    Args:
        genes: Gene name(s) to plot. Comma-separated for multiple (e.g. 'TNNT2,MYH7').
        subset_key: The metadata column to filter on (e.g., 'cell_type', 'leiden').
        subset_value: The value to filter for (e.g., 'cardiomyocyte', '0').
                      Broad terms match ALL related subtypes automatically.
        groupby: The observation key to group the subset by (e.g., 'orig.ident').
                 If empty, auto-detects clustering key.
    """
    adata = _get_adata()

    # Validate subset_key
    if subset_key not in adata.obs.columns:
        available = ", ".join(sorted(adata.obs.columns))
        return f"Error: '{subset_key}' not found in observations. Available: {available}"

    # Resolve subset_value — may match multiple subtypes
    resolved = _resolve_subset_values(adata, subset_key, subset_value)
    if not resolved:
        available = ", ".join(sorted(str(v) for v in adata.obs[subset_key].unique()))
        return (
            f"Error: '{subset_value}' not found in column '{subset_key}'.\n"
            f"Available values: {available}"
        )

    # Auto-detect groupby
    if not groupby:
        groupby = _get_cluster_key()

    # Create subset from all matching values
    adata_sub, label = _build_subset(adata, subset_key, resolved)
    n_cells = adata_sub.n_obs

    if n_cells == 0:
        return f"Error: No cells found for {subset_key} matching '{subset_value}'."

    # Validate groupby on subset
    if groupby not in adata_sub.obs.columns:
        return f"Error: '{groupby}' not found in subset observations."

    # Parse genes
    if "," in genes:
        gene_list = [g.strip() for g in genes.split(",")]
    else:
        gene_list = genes.strip()

    try:
        result = plot_violin(adata_sub, genes=gene_list, groupby=groupby)

        # Report condition coverage
        coverage_note = _report_condition_coverage(adata, adata_sub, groupby, label)

        # Compute per-group summary statistics for the agent to report
        summary_lines = _compute_group_summary(adata_sub, gene_list, groupby)

        # Build subset description
        if len(resolved) == 1:
            subset_desc = f"{resolved[0]}"
        else:
            subset_desc = f"{len(resolved)} matching types: {', '.join(resolved)}"

        # Build realistic code that reflects the actual subset operation
        subset_code = (
            f"# Subset to {subset_key} matching: {resolved}\n"
            f"mask = adata.obs[\"{subset_key}\"].isin({resolved})\n"
            f"adata_sub = adata[mask].copy()  # {n_cells:,} cells\n"
            f"{result.code.replace('adata', 'adata_sub')}"
        )

        enhanced = PlotResult(
            image=result.image,
            code=subset_code,
            message=(
                f"Violin plot of {genes} in {subset_desc} ({n_cells:,} cells), "
                f"grouped by {groupby}.{coverage_note}"
            ),
        )

        result_str = _store_and_return(enhanced)

        # Append quantitative summary so the agent can answer comparison questions
        result_str += f"\n\n--- Expression Summary ({genes} in {subset_desc}) ---\n"
        result_str += summary_lines
        result_str += (
            "\n\nIMPORTANT: Use the summary above to DIRECTLY ANSWER any comparison "
            "question. State which group is higher/lower and by how much."
        )
        update_viz_state(
            "violin", groupby=groupby,
            genes=gene_list if isinstance(gene_list, list) else [gene_list],
        )
        return result_str

    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Subset violin plot failed")
        return f"Unexpected error: {e}"


@tool
def subset_feature_plot(
    gene: str,
    subset_key: str,
    subset_value: str,
    split_by: str = "",
) -> str:
    """Generate a feature plot (UMAP colored by gene expression) for a SUBSET of cells.

    Use this when the user wants to see spatial gene expression on UMAP for just
    one cell type or group, optionally split by another variable (e.g., condition).

    IMPORTANT — Broad cell type matching:
    The subset_value supports broad matching. For example:
    - "cardiomyocytes" will match ALL cardiomyocyte subtypes in the dataset
      (e.g., "Early cardiomyocyte", "Ventricular cardiomyocyte", etc.)
    This ensures broad biological categories include all relevant subtypes.

    If some conditions have no cells of that type after filtering, those
    conditions are reported as missing rather than silently dropped.

    Examples:
    - "Show TNNT2 expression on UMAP for cardiomyocytes only"
      → subset_key="cell_type", subset_value="cardiomyocyte"
    - "Show TNNT2 in cardiomyocytes split by condition"
      → subset_key="cell_type", subset_value="cardiomyocyte", split_by="orig.ident"

    Args:
        gene: Gene name to visualize.
        subset_key: The metadata column to filter on (e.g., 'cell_type').
        subset_value: The value to filter for (e.g., 'cardiomyocyte').
                      Broad terms match ALL related subtypes automatically.
        split_by: Optional observation key to split the subset into panels.
    """
    adata = _get_adata()

    # Validate subset_key
    if subset_key not in adata.obs.columns:
        available = ", ".join(sorted(adata.obs.columns))
        return f"Error: '{subset_key}' not found in observations. Available: {available}"

    # Resolve subset_value — may match multiple subtypes
    resolved = _resolve_subset_values(adata, subset_key, subset_value)
    if not resolved:
        available = ", ".join(sorted(str(v) for v in adata.obs[subset_key].unique()))
        return (
            f"Error: '{subset_value}' not found in column '{subset_key}'.\n"
            f"Available values: {available}"
        )

    # Create subset from all matching values
    adata_sub, label = _build_subset(adata, subset_key, resolved)
    n_cells = adata_sub.n_obs

    if n_cells == 0:
        return f"Error: No cells found for {subset_key} matching '{subset_value}'."

    try:
        result = plot_feature(
            adata_sub, gene=gene,
            split_by=split_by if split_by else None,
        )

        # Report condition coverage if split_by is used
        coverage_note = ""
        if split_by and split_by in adata.obs.columns:
            coverage_note = _report_condition_coverage(adata, adata_sub, split_by, label)

        # Build subset description
        if len(resolved) == 1:
            subset_desc = f"{resolved[0]}"
        else:
            subset_desc = f"{len(resolved)} matching types: {', '.join(resolved)}"

        # Build realistic code that reflects the actual subset operation
        subset_code = (
            f"# Subset to {subset_key} matching: {resolved}\n"
            f"mask = adata.obs[\"{subset_key}\"].isin({resolved})\n"
            f"adata_sub = adata[mask].copy()  # {n_cells:,} cells\n"
            f"{result.code.replace('adata', 'adata_sub')}"
        )

        enhanced = PlotResult(
            image=result.image,
            code=subset_code,
            message=(
                f"Feature plot of {gene} in {subset_desc} ({n_cells:,} cells)"
                + (f", split by {split_by}." if split_by else ".")
                + coverage_note
            ),
        )

        result_str = _store_and_return(enhanced)
        update_viz_state("feature", genes=[gene], split_by=split_by or None)
        return result_str

    except ValueError as e:
        return f"Error: {e}"
    except Exception as e:
        logger.exception("Subset feature plot failed")
        return f"Unexpected error: {e}"


def get_subset_tools() -> list:
    """Return all subset analysis tools."""
    return [subset_violin_plot, subset_feature_plot]
