"""Unit tests for the task router (src/agent/router.py).

All tests are deterministic — no API calls, no mocking required.
Run with: pytest tests/test_router.py -v
"""

import pytest
from src.agent.router import RouterResult, classify_intent


# ---------------------------------------------------------------------------
# Layer 2a — known tool prompts
# One representative prompt per tool. Asserts correct tool name is returned.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("prompt,expected_tool", [
    # Visualization
    ("Show me a UMAP colored by cell type",             "umap_plot"),
    ("Plot a UMAP split by condition",                  "umap_plot"),
    ("Show me the embedding of my cells",               "umap_plot"),
    ("Plot a violin plot for MKI67",                    "violin_plot"),
    ("Show the expression distribution of ACTB",        "violin_plot"),
    ("Show a dot plot of top marker genes",             "dotplot"),
    ("Generate a dotplot across clusters",              "dotplot"),
    ("Show a feature plot for CD3E",                    "feature_plot"),
    ("Show CD3E expression on umap",                    "feature_plot"),
    ("Plot a heatmap of the top genes",                 "heatmap_plot"),
    ("Show a heat map across cell types",               "heatmap_plot"),
    ("Scatter plot of MALAT1 vs ACTB",                  "scatter_plot"),
    ("Show gene vs gene correlation between CD3E and CD4", "scatter_plot"),
    ("Show the volcano plot for cluster 2",             "volcano_plot_tool"),
    # Differential expression
    ("Run pairwise comparison disease vs normal",       "compare_groups_de"),
    ("Compare groups: early vs late cardiomyocytes",    "compare_groups_de"),
    ("Run pairwise DE analysis",                        "compare_groups_de"),
    ("What are the cluster markers for cluster 3",      "get_cluster_degs"),
    ("Show markers for cluster 1",                      "get_cluster_degs"),
    ("Run differential expression analysis",            "differential_expression"),
    ("Find differentially expressed genes",             "differential_expression"),
    ("Run DE analysis across all clusters",             "differential_expression"),
    ("Show the top markers per cluster",                "get_top_markers"),
    ("Get the best markers for each group",             "get_top_markers"),
    # Preprocessing
    ("Preprocess my data",                              "preprocess_data"),
    ("Normalize and cluster the data",                  "preprocess_data"),
    ("Run clustering on my dataset",                    "preprocess_data"),
    # Dataset / Metadata
    ("How many cells are in this dataset",              "summarize_obs_column"),
    ("Give me a data overview",                         "dataset_info"),
    ("How many genes are in my data",                   "dataset_info"),
    ("What's the data status",                          "check_data_status"),
    ("Is data ready for analysis",                      "check_data_status"),
    ("What obs columns are available",                  "inspect_metadata"),
    ("Show available columns",                          "inspect_metadata"),
    ("Show the cluster names",                          "get_cluster_mapping"),
    ("What are the clusters called",                    "get_cluster_mapping"),
    # QC
    ("Calculate mitochondrial percentage",              "calculate_mito_pct"),
    ("Show mito pct per cell",                          "calculate_mito_pct"),
    ("How many cells in each cell type",                "summarize_obs_column"),
    ("Count by condition",                              "summarize_obs_column"),
    ("Show qc metrics summary",                         "summarize_qc_metrics_tool"),
    ("Run quality control on the dataset",              "summarize_qc_metrics_tool"),
    # Export
    ("Get the DE results table",                        "get_de_results_table"),
    ("Export DE results",                               "get_de_results_table"),
    # Advanced analysis
    ("Calculate average expression of CD3E",            "calculate_average_expression"),
    ("Show mean expression per cluster",                "calculate_average_expression"),
    ("Which cluster has highest expression of MKI67",   "find_highest_expression"),
    ("Find the highest expressing cluster for ACTB",    "find_highest_expression"),
    ("Highlight cluster 3 on the UMAP",                 "highlight_cluster"),
    ("Rename cluster 2 to T cells",                     "rename_cluster"),
    ("Annotate cluster 5 as macrophages",               "rename_cluster"),
])
def test_layer_2a_tool_match(prompt, expected_tool):
    result = classify_intent(prompt)
    assert result.layer == "2a", (
        f"Expected layer '2a' for prompt: '{prompt}'\n"
        f"Got layer='{result.layer}', matched_on='{result.matched_on}'"
    )
    assert result.task_type == expected_tool, (
        f"Expected task_type='{expected_tool}' for prompt: '{prompt}'\n"
        f"Got task_type='{result.task_type}', matched_on='{result.matched_on}'"
    )
    assert result.confidence == "high"
    assert result.matched_on is not None


# ---------------------------------------------------------------------------
# Layer 2b — interpretive and conversational prompts
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("prompt", [
    # Explanation requests
    "Explain what these clusters represent",
    "Explain why cluster 3 is separate",
    "Interpret the DE results",
    "What does this mean biologically",
    "What does this result mean",
    "What does that pattern suggest",
    # Causal questions
    "Why is MKI67 high in cluster 2",
    "Why are these cells separated",
    "Why does cluster 1 have low expression",
    "Why do cardiomyocytes cluster together",
    # Judgment / guidance
    "Is this normal for a heart dataset",
    "Does this suggest disease progression",
    "What should I focus on next",
    "What can you tell me about these results",
    # Biological interpretation
    "What is the biological significance of this finding",
    "Does this make sense of the clustering pattern",
    "Give me insight into the cell type composition",
    "What's interesting about this dataset",
])
def test_layer_2b_interpretive(prompt):
    result = classify_intent(prompt)
    assert result.layer == "2b", (
        f"Expected layer '2b' for prompt: '{prompt}'\n"
        f"Got layer='{result.layer}', matched_on='{result.matched_on}'"
    )
    assert result.task_type is None
    assert result.confidence == "high"
    assert result.matched_on is not None


# ---------------------------------------------------------------------------
# Ambiguous — prompts that should not match either tier
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("prompt", [
    "Hello",
    "Thanks",
    "Can you help me",
    "What next",
    "OK",
    "I want to analyze my data",  # too vague to match a tool
])
def test_ambiguous(prompt):
    result = classify_intent(prompt)
    assert result.layer == "ambiguous", (
        f"Expected 'ambiguous' for prompt: '{prompt}'\n"
        f"Got layer='{result.layer}', matched_on='{result.matched_on}'"
    )
    assert result.task_type is None
    assert result.confidence == "low"
    assert result.matched_on is None


# ---------------------------------------------------------------------------
# RouterResult is immutable (frozen dataclass)
# ---------------------------------------------------------------------------

def test_router_result_is_immutable():
    result = classify_intent("show me a UMAP")
    with pytest.raises((AttributeError, TypeError)):
        result.layer = "2b"  # type: ignore[misc]


# ---------------------------------------------------------------------------
# Case insensitivity
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("prompt", [
    "SHOW ME A UMAP",
    "Show Me A UMAP",
    "show me a umap",
])
def test_case_insensitive(prompt):
    result = classify_intent(prompt)
    assert result.layer == "2a"
    assert result.task_type == "umap_plot"


# ---------------------------------------------------------------------------
# CEO test cases — prompts from the April 1 GSE223414 session
# These are the real failures we need to correctly classify.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("prompt,expected_layer,expected_tool", [
    # MKI67 query — asking which cell type expresses it highest
    # Root: wrong dimension returned (cluster IDs vs cell types)
    # This is a find_highest_expression request → 2a
    ("Which cell type expresses MKI67 highest",     "2a", "find_highest_expression"),
    # UMAP split + color by cell type — required 3 correction rounds
    # Clear visualization request → 2a
    ("Show me a UMAP split by condition colored by cell type", "2a", "umap_plot"),
    # DE pairwise: Early cardiomyocyte disease vs normal
    # Pairwise comparison → 2a
    ("Run pairwise DE for early cardiomyocyte disease vs normal", "2a", "compare_groups_de"),
    # Dot plot top 3 per cell type
    # Dot plot request → 2a
    ("Show a dot plot of top 3 genes per cell type", "2a", "dotplot"),
    # Cell count query — same query returned different counts
    # Dataset info / summarize → 2a
    ("How many cells are in the early cardiomyocyte condition", "2a", "summarize_obs_column"),
])
def test_ceo_test_cases(prompt, expected_layer, expected_tool):
    result = classify_intent(prompt)
    assert result.layer == expected_layer, (
        f"CEO test case failed for: '{prompt}'\n"
        f"Expected layer='{expected_layer}', got='{result.layer}', matched_on='{result.matched_on}'"
    )
    if expected_tool:
        assert result.task_type == expected_tool, (
            f"CEO test case failed for: '{prompt}'\n"
            f"Expected task_type='{expected_tool}', got='{result.task_type}'"
        )
