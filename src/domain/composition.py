"""Cell type composition — Yalu Layer 1 §6.

Single ``@register`` tool covering Yalu §6A.1–§6A.3 + §6B.1–§6B.3 via
internal branching on ``mode`` + optional ``celltypes`` subset (§6B).

Modes:
- ``'count'``: absolute cell counts; grouped bars per cell type, condition-colored.
- ``'proportion_by_celltype'``: proportion within each cell type; grouped bars,
  condition-colored.
- ``'proportion_by_condition'``: proportion within each condition; stacked bars
  colored by cell type.

The ``reference`` param anchors blue in the condition color palette for the
two grouped-bars modes; it is ignored in ``'proportion_by_condition'`` (stacks
are colored by cell type, not condition). See param description.

Yalu §6 explicit precondition: dataset must have ≥2 conditions. Tool body
validates and raises ``ToolExecutionError`` if the condition column has fewer
than 2 unique values.

Subset (§6B) uses a ``celltypes: list[str]`` param (NOT the generic
subset_key/subset_value pair from heatmap/dot/violin) — §6B always subsets
the celltype axis; the generic machinery would be unused flexibility.

Naming divergence from heatmap/dot convention: this tool uses domain names
``celltype_col`` / ``condition_col`` instead of ``groupby`` / ``split_by``.
§6's two columns are equally required and roles flip per mode; asymmetric
naming would lie about the structure. Convention vs. correctness — the
asymmetric idiom doesn't fit here.
"""

from __future__ import annotations

import io
import logging
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

from src.agent.viz_state import update_viz_state
from src.core.registry import register
from src.core.results import ArtifactResult, ToolExecutionError
from src.domain.plot_style import condition_color_map

if TYPE_CHECKING:
    from anndata import AnnData

logger = logging.getLogger(__name__)


_VALID_MODES = ("count", "proportion_by_celltype", "proportion_by_condition")


@register(
    description=(
        "Cell type composition bar plot across conditions (Yalu §6). Three "
        "modes: 'count' (absolute cell numbers, grouped bars per cell type), "
        "'proportion_by_celltype' (proportion within each cell type, grouped "
        "bars), 'proportion_by_condition' (proportion within each condition, "
        "stacked bars colored by cell type). Optional celltypes list "
        "subsets to specific cell types (Yalu §6B). Requires at least 2 "
        "conditions in the dataset."
    ),
    params={
        "celltype_col": {
            "description": (
                "obs column containing cell type labels (typically the "
                "cell-type annotation column)."
            ),
            "field_type": "obs_column",
        },
        "condition_col": {
            "description": (
                "obs column containing condition labels (treatment/control, "
                "sample identifier, etc.). Yalu §6 requires ≥2 unique values."
            ),
            "field_type": "obs_column",
        },
        "mode": {
            "description": (
                "Composition rendering mode: 'count' (absolute cell counts), "
                "'proportion_by_celltype' (proportion within each cell type), "
                "or 'proportion_by_condition' (proportion within each "
                "condition, stacked bars)."
            ),
            "enum": ["count", "proportion_by_celltype", "proportion_by_condition"],
        },
        "celltypes": {
            "description": (
                "Optional list of cell-type values to subset before plotting "
                "(Yalu §6B). When omitted, all cell types are included (§6A)."
            ),
            "field_type": "cell_type",
        },
        "reference": {
            "description": (
                "Optional reference / control condition value. Used in 'count' "
                "and 'proportion_by_celltype' modes to anchor the reference "
                "condition to blue in the condition color palette. Ignored "
                "in 'proportion_by_condition' mode (stacks are colored by "
                "cell type, not condition)."
            ),
            "field_type": "condition",
        },
    },
)
def composition_barplot(
    adata: "AnnData",
    celltype_col: str,
    condition_col: str,
    mode: str = "count",
    celltypes: list[str] | None = None,
    reference: str = "",
) -> ArtifactResult:
    """Cell type composition bar plot; covers Yalu §6A.1–§6A.3 + §6B.1–§6B.3."""
    if mode not in _VALID_MODES:
        raise ToolExecutionError(
            f"mode must be one of {_VALID_MODES}; got {mode!r}.",
            tool_name="composition_barplot",
        )

    if celltype_col not in adata.obs.columns:
        raise ToolExecutionError(
            f"celltype_col '{celltype_col}' not in adata.obs.",
            tool_name="composition_barplot",
        )

    if condition_col not in adata.obs.columns:
        raise ToolExecutionError(
            f"condition_col '{condition_col}' not in adata.obs.",
            tool_name="composition_barplot",
        )

    n_conditions = int(adata.obs[condition_col].nunique())
    if n_conditions < 2:
        raise ToolExecutionError(
            f"Composition analysis requires at least 2 conditions in "
            f"adata.obs[{condition_col!r}]; found {n_conditions} (Yalu §6).",
            tool_name="composition_barplot",
        )

    cts_subset = list(celltypes) if celltypes else []
    if cts_subset:
        mask = adata.obs[celltype_col].astype(str).isin(cts_subset)
        if not mask.any():
            raise ToolExecutionError(
                f"None of the requested cell types {cts_subset} found in "
                f"adata.obs[{celltype_col!r}].",
                tool_name="composition_barplot",
            )
        adata = adata[mask].copy()

    count_df = (
        adata.obs.groupby([celltype_col, condition_col], observed=True)
        .size()
        .reset_index(name="count")
    )

    if hasattr(adata.obs[celltype_col], "cat"):
        celltypes_list = [
            str(c) for c in adata.obs[celltype_col].cat.categories
            if (adata.obs[celltype_col] == c).any()
        ]
    else:
        celltypes_list = sorted(adata.obs[celltype_col].astype(str).unique())
    conditions = sorted({str(c) for c in adata.obs[condition_col].dropna().unique()})

    if mode in ("count", "proportion_by_celltype"):
        color_map = condition_color_map(conditions, reference=reference or None)
        image_bytes, code = _composition_grouped_bars(
            count_df, celltypes_list, conditions, color_map,
            celltype_col, condition_col, mode,
        )
    else:
        image_bytes, code = _composition_stacked_bars(
            count_df, celltypes_list, conditions,
            celltype_col, condition_col,
        )

    update_viz_state(
        "composition_barplot",
        groupby=celltype_col,
        split_by=condition_col,
        subset_value=cts_subset or None,
        reference=reference or None,
    )

    parts = [
        f"Cell type composition ({mode}): {len(celltypes_list)} cell types "
        f"× {n_conditions} conditions"
    ]
    if cts_subset:
        parts.append(f"; subset to {cts_subset}")
    if reference and mode != "proportion_by_condition":
        parts.append(f"; reference condition: {reference}")
    parts.append(".")

    return ArtifactResult(
        text="".join(parts),
        artifact_kind="image",
        tool_name="composition_barplot",
        params_used={
            "celltype_col": celltype_col,
            "condition_col": condition_col,
            "mode": mode,
            "celltypes": cts_subset or None,
            "reference": reference or None,
        },
        entities_acted_on=list(conditions) + cts_subset,
        image_bytes=image_bytes,
        code=code,
    )


def _composition_grouped_bars(
    count_df,
    celltypes_list: list[str],
    conditions: list[str],
    color_map: dict[str, str],
    celltype_col: str,
    condition_col: str,
    mode: str,
) -> tuple[bytes, str]:
    """Yalu §6A.1 / §6A.2: grouped bars per cell type, condition-colored."""
    df = count_df.copy()
    if mode == "proportion_by_celltype":
        df["value"] = df.groupby(celltype_col, observed=True)["count"].transform(
            lambda x: x / x.sum() if x.sum() > 0 else 0
        )
        ylabel = "Proportion"
        title = "Proportion of each cell type across conditions"
    else:
        df["value"] = df["count"]
        ylabel = "Cell Count"
        title = "Cell number by cell type and condition"

    fig, ax = plt.subplots(figsize=(12, 5))
    x = np.arange(len(celltypes_list))
    n_cond = len(conditions)
    width = 0.8 / max(n_cond, 1)

    for i, cond in enumerate(conditions):
        values = []
        for ct in celltypes_list:
            row = df[
                (df[celltype_col].astype(str) == ct)
                & (df[condition_col].astype(str) == cond)
            ]
            values.append(float(row["value"].values[0]) if len(row) else 0.0)
        ax.bar(
            x + i * width - (n_cond - 1) * width / 2,
            values, width,
            label=cond, color=color_map.get(cond, "#888888"), alpha=0.85,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(celltypes_list, rotation=45, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(title=condition_col)

    if mode == "proportion_by_celltype":
        ax.set_ylim(0, 1.1)
    else:
        cur_max = ax.get_ylim()[1]
        ax.set_ylim(0, cur_max * 1.1)

    plt.tight_layout()
    code = (
        f"# {title}\n"
        f"count_df = adata.obs.groupby([{celltype_col!r}, {condition_col!r}]).size()\n"
        f"# mode={mode!r}; grouped bars, condition-colored "
        f"(reference anchors blue if set)"
    )
    return _capture_png(), code


def _composition_stacked_bars(
    count_df,
    celltypes_list: list[str],
    conditions: list[str],
    celltype_col: str,
    condition_col: str,
) -> tuple[bytes, str]:
    """Yalu §6A.3: stacked bars per condition, colored by cell type."""
    df = count_df.copy()
    df["proportion"] = df.groupby(condition_col, observed=True)["count"].transform(
        lambda x: x / x.sum() if x.sum() > 0 else 0
    )

    n_celltypes = len(celltypes_list)
    palettes = (
        sc.pl.palettes.default_20,
        sc.pl.palettes.default_28,
        sc.pl.palettes.default_102,
    )
    palette = next(
        (p for p in palettes if n_celltypes <= len(p)),
        sc.pl.palettes.default_102,
    )
    celltype_colors = {
        ct: palette[i % len(palette)] for i, ct in enumerate(celltypes_list)
    }

    fig_height = max(4, n_celltypes * 0.3)
    fig, ax = plt.subplots(figsize=(4, fig_height))
    x = np.arange(len(conditions))
    bottom = np.zeros(len(conditions))

    for ct in celltypes_list:
        proportions = []
        for cond in conditions:
            row = df[
                (df[celltype_col].astype(str) == ct)
                & (df[condition_col].astype(str) == cond)
            ]
            proportions.append(float(row["proportion"].values[0]) if len(row) else 0.0)
        proportions_arr = np.array(proportions)
        ax.bar(
            x, proportions_arr, bottom=bottom,
            label=ct, color=celltype_colors[ct], alpha=0.85, width=0.5,
        )
        bottom += proportions_arr

    ax.set_xticks(x)
    ax.set_xticklabels(conditions)
    ax.set_xlim(-0.5, len(conditions) - 0.5)
    ax.set_ylabel("Proportion")
    ax.set_title("Cell type composition per condition")
    ax.legend(
        bbox_to_anchor=(1.05, 1), loc="upper left",
        fontsize=8, title=celltype_col,
    )
    ax.set_ylim(0, 1.05)

    plt.tight_layout()
    code = (
        f"# Cell type composition per condition (stacked, celltype-colored)\n"
        f"count_df = adata.obs.groupby([{celltype_col!r}, {condition_col!r}]).size()\n"
        f"# mode='proportion_by_condition'; reference param ignored in this mode"
    )
    return _capture_png(), code


def _capture_png() -> bytes:
    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=150, bbox_inches="tight", facecolor="white")
    plt.close("all")
    buf.seek(0)
    return buf.read()
