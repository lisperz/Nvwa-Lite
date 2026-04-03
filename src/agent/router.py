"""Task Router — classifies incoming prompts before agent invocation.

Implements a two-tier rule-based classifier:
  Tier 1: Detect Layer 2b (interpretive/conversational) linguistic markers.
  Tier 2: Match Layer 2a tool keywords.
  Fallback: Return 'ambiguous' if neither tier matches.

No LLM calls. Fully deterministic and testable without mocking.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class RouterResult:
    """Classification result from the task router."""

    layer: str              # "2a", "2b", or "ambiguous"
    task_type: str | None   # tool name if layer is "2a", None otherwise
    confidence: str         # "high" or "low"
    matched_on: str | None  # keyword or pattern that triggered the match


# ---------------------------------------------------------------------------
# Tier 1 — Layer 2b patterns
# Phrases that indicate interpretive or conversational intent.
# Checked before tool matching. Any match short-circuits to Layer 2b.
# ---------------------------------------------------------------------------
LAYER_2B_PATTERNS: list[str] = [
    "what does this",
    "what does that",
    "what do these",
    "what do those",
    "why is",
    "why are",
    "why does",
    "why do",
    "explain",
    "interpret",
    "biologically",
    "significance of",
    "what can you tell me",
    "what should i",
    "is this normal",
    "suggest",
    "what's interesting",
    "what is interesting",
    "make sense of",
    "give me insight",
    "what insight",
]

# ---------------------------------------------------------------------------
# Tier 2 — Tool intent map
# Maps each tool name to a set of keyword phrases.
# Checked in order — first match wins.
# More specific entries must come before more generic ones.
# All keywords matched as case-insensitive substrings.
# ---------------------------------------------------------------------------
TOOL_INTENT_MAP: list[tuple[str, set[str]]] = [
    # --- Visualization ---
    # highlight_cluster and feature_plot must come before umap_plot —
    # their prompts often contain "umap" as context, not as the task.
    ("highlight_cluster", {
        "highlight cluster", "focus on cluster", "show only cluster",
    }),
    ("feature_plot", {
        "feature plot", "gene on umap", "expression on umap",
        "where is expressed", "color umap by gene",
    }),
    ("umap_plot", {
        "umap", "embedding", "dimensionality reduction", "cell map", "split by",
    }),
    ("violin_plot", {
        "violin", "expression distribution", "gene distribution",
    }),
    ("dotplot", {
        "dot plot", "dotplot", "bubble plot",
    }),
    ("heatmap_plot", {
        "heatmap", "heat map", "expression matrix",
    }),
    ("scatter_plot", {
        "scatter plot", "gene vs gene", "correlation between genes",
    }),
    ("volcano_plot_tool", {
        "volcano", "fold change plot", "visualize de results",
    }),
    # --- Differential Expression (specific before generic) ---
    ("compare_groups_de", {
        "disease vs", "vs normal", "condition vs", "pairwise", "compare between",
        "compare groups",
    }),
    ("get_cluster_degs", {
        "markers for cluster", "what defines cluster", "cluster markers",
        "top genes for cluster", "degs for cluster",
    }),
    ("differential_expression", {
        "differential expression", "de analysis", "deg", "differentially expressed",
        "all markers", "find markers",
    }),
    # --- Marker Genes ---
    ("get_top_markers", {
        "top markers", "top genes", "best markers", "marker genes per cluster",
    }),
    # --- Preprocessing ---
    ("preprocess_data", {
        "preprocess", "normalize", "run umap", "process data",
        "prepare data", "run clustering", "cluster the data",
    }),
    # --- Dataset / Metadata ---
    # summarize_obs_column must come before dataset_info —
    # "how many cells in X" is more specific than "how many cells".
    ("summarize_obs_column", {
        "how many cells in", "how many cells are in", "how many cells per",
        "count by", "cells per", "distribution of cells", "summarize column",
    }),
    ("dataset_info", {
        "what's in my data", "dataset info", "how many cells",
        "how many genes", "data overview", "data summary",
    }),
    ("check_data_status", {
        "data status", "processing status", "what stage", "is data ready",
    }),
    ("inspect_metadata", {
        "obs columns", "cell attributes", "available columns",
        "what variables", "list columns",
    }),
    ("get_cluster_mapping", {
        "cluster mapping", "cluster names", "cluster ids",
        "what are the clusters", "cluster annotation", "list clusters",
    }),
    # --- QC ---
    ("calculate_mito_pct", {
        "mitochondrial", "mito percentage", "mito pct", "mt genes", "mito genes",
    }),
    ("summarize_qc_metrics_tool", {
        "qc metrics", "quality control", "qc summary",
        "read depth", "genes per cell",
    }),
    # --- Export ---
    ("get_de_results_table", {
        "de table", "export de", "download markers",
        "de results table", "get all de",
    }),
    # --- Advanced Analysis ---
    ("calculate_average_expression", {
        "average expression", "mean expression", "average gene expression",
    }),
    ("find_highest_expression", {
        "highest expression", "most expressed", "highest expressing",
        "expresses highest", "which cell type expresses",
    }),
    ("rename_cluster", {
        "rename cluster", "annotate cluster", "label cluster", "assign cell type",
    }),
]


def classify_intent(prompt: str) -> RouterResult:
    """Classify an incoming user prompt as Layer 2a, 2b, or ambiguous.

    Tier 1: Check for Layer 2b interpretive patterns first.
    Tier 2: Match against known tool keyword sets for Layer 2a.
    Fallback: Return ambiguous if neither tier matches.

    Args:
        prompt: The raw user input string.

    Returns:
        RouterResult with layer, task_type, confidence, and matched_on.
    """
    normalized = prompt.lower().strip()

    # Tier 1 — interpretive language signals Layer 2b
    for pattern in LAYER_2B_PATTERNS:
        if pattern in normalized:
            return RouterResult(
                layer="2b",
                task_type=None,
                confidence="high",
                matched_on=pattern,
            )

    # Tier 2 — keyword match to specific tool
    for tool_name, keywords in TOOL_INTENT_MAP:
        for keyword in keywords:
            if keyword in normalized:
                return RouterResult(
                    layer="2a",
                    task_type=tool_name,
                    confidence="high",
                    matched_on=keyword,
                )

    # No match
    return RouterResult(
        layer="ambiguous",
        task_type=None,
        confidence="low",
        matched_on=None,
    )