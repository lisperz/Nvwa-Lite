"""Marker gene survey — Yalu Layer 1 §4.1.

This module hosts ``find_all_markers``, the multi-group marker survey
tool. Structurally distinct from ``src/domain/de.py``'s single-target-
group tools:

- Different scanpy call shape: no ``groups`` arg → scanpy iterates all
  unique values in the groupby column.
- Own uns slot ``nvwa_all_markers`` (separate from de.py's pairwise and
  one-vs-rest slots — one writer per slot, name-by-content).
- Output is a long-form DataFrame with a ``group`` column identifying
  which cell type each row belongs to, vs de.py's single-group result.

Future §4.3 (dot_plot_top_markers) will consume this slot.

UX contract:
- ``find_all_markers.text`` does NOT include a volcano offer (Yalu §4
  doesn't describe one) and does NOT include the statistical disclaimer
  (Yalu §4 frames marker discovery as characterization, not comparison).
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import scanpy as sc

from src.agent.viz_state import update_viz_state
from src.core.registry import register
from src.core.results import ArtifactResult, ToolExecutionError

if TYPE_CHECKING:
    from anndata import AnnData

logger = logging.getLogger(__name__)


_ALL_MARKERS_UNS_KEY = "nvwa_all_markers"


@register(
    description=(
        "Marker gene survey across ALL cell types (Yalu §4.1). For each "
        "value in the groupby column, finds genes specifically expressed "
        "in that group vs all other cells. Writes scanpy-native survey "
        f"results to adata.uns[{_ALL_MARKERS_UNS_KEY!r}]. Use find_markers "
        "for a single specific cell type; use run_de for pairwise group A "
        "vs group B."
    ),
    params={
        "groupby": {
            "description": (
                "obs column containing cell type labels (typically the "
                "cell-type annotation column). The survey iterates each "
                "unique value in this column."
            ),
            "field_type": "obs_column",
        },
        "method": {
            "description": (
                "Statistical method: 'wilcoxon' (default), 't-test', or "
                "'logreg'."
            ),
        },
        "pvals_adj_threshold": {
            "description": (
                "Adj p-value threshold for filtering. Default 0.05."
            ),
        },
        "logfc_threshold": {
            "description": (
                "log2FC threshold for filtering (only-positive markers "
                "when threshold>=0). Default 0 — Yalu §4.1 default returns "
                "positive markers."
            ),
        },
        "n_top_genes": {
            "description": (
                "Top N markers shown PER cell type in the summary text. "
                "Default 10."
            ),
        },
    },
)
def find_all_markers(
    adata: "AnnData",
    groupby: str,
    method: str = "wilcoxon",
    pvals_adj_threshold: float = 0.05,
    logfc_threshold: float = 0.0,
    n_top_genes: int = 10,
) -> ArtifactResult:
    """Multi-group marker survey across all groups in groupby; covers Yalu §4.1."""
    if groupby not in adata.obs.columns:
        raise ToolExecutionError(
            f"groupby column '{groupby}' not in adata.obs.",
            tool_name="find_all_markers",
        )

    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method,
        pts=True,
        key_added=_ALL_MARKERS_UNS_KEY,
    )
    markers_df = sc.get.rank_genes_groups_df(
        adata, group=None, key=_ALL_MARKERS_UNS_KEY,
    )

    filtered = markers_df[
        (markers_df["pvals_adj"] < pvals_adj_threshold)
        & (markers_df["logfoldchanges"] > logfc_threshold)
    ]

    n_groups = int(filtered["group"].nunique())
    n_total = int(len(filtered))

    text_lines = [
        f"Markers across {n_groups} {groupby} values "
        f"({adata.n_obs:,} cells, method={method}).",
        f"{n_total} significant markers total "
        f"(padj<{pvals_adj_threshold}, logFC>{logfc_threshold}).",
        "",
    ]
    for celltype, group_df in filtered.groupby("group", sort=False):
        top = group_df.nlargest(n_top_genes, "logfoldchanges")
        gene_list = ", ".join(top["names"].tolist())
        text_lines.append(f"{celltype}: {gene_list}")

    text = "\n".join(text_lines)

    update_viz_state(
        "find_all_markers",
        groupby=groupby,
    )

    code = (
        f"sc.tl.rank_genes_groups(adata, groupby={groupby!r}, "
        f"method={method!r}, key_added={_ALL_MARKERS_UNS_KEY!r})\n"
        f"markers_df = sc.get.rank_genes_groups_df(adata, group=None, "
        f"key={_ALL_MARKERS_UNS_KEY!r})"
    )

    return ArtifactResult(
        text=text,
        artifact_kind="csv",
        tool_name="find_all_markers",
        params_used={
            "groupby": groupby,
            "method": method,
            "pvals_adj_threshold": pvals_adj_threshold,
            "logfc_threshold": logfc_threshold,
            "n_top_genes": n_top_genes,
        },
        entities_acted_on=[groupby],
        csv_data=filtered.to_csv(index=False),
        display_df=filtered.head(20).to_markdown(index=False),
        code=code,
    )
