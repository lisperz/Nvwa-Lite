# Bug Fix: Heatmap Shows 9 Genes Instead of 10

**Date**: 2026-03-01
**Issue**: When user requests "top 10 marker genes of each cluster", heatmap only displays 9 genes
**Status**: Fixed

## Root Cause

When collecting top N marker genes from multiple clusters, duplicate genes appear across clusters. Scanpy's heatmap function automatically removes duplicate genes, resulting in fewer genes than expected.

**Example with 8 clusters, requesting 10 genes per cluster:**
- Total genes collected: 8 × 10 = 80 genes
- Duplicate genes: 14
- Unique genes displayed: 80 - 14 = 66 genes

The agent was directly passing gene names to the heatmap without handling duplicates, causing confusion when the displayed count didn't match the request.

## Solution

Created a new helper module and tool to intelligently extract marker genes:

### 1. New Module: `src/analysis/marker_genes.py`

Provides three helper functions:

- `get_top_marker_genes_per_cluster()`: Gets top N genes per cluster, removes duplicates while preserving order
- `get_top_marker_genes_overall()`: Gets top N genes ranked by significance across all clusters
- `get_top_marker_genes_per_cluster_exact()`: Gets exactly N genes per cluster (may have duplicates)

### 2. New Tool: `get_top_markers`

Added to `src/agent/tools.py`:

```python
@tool
def get_top_markers(n_genes_per_cluster: int = 10, groupby: str = "leiden") -> str:
    """Get top marker genes for each cluster from DE results."""
    genes = get_top_marker_genes_per_cluster(adata, n_genes=n_genes_per_cluster, groupby=groupby)
    return ",".join(genes)
```

Returns comma-separated list of unique genes with clear messaging about total count.

### 3. Updated System Prompt

Added workflow guidance in `src/agent/prompts.py`:

- Listed `get_top_markers` in core analysis tools section
- Added usage instructions: "Use this BEFORE creating heatmaps to get the gene list"
- Added example interaction showing the heatmap workflow
- Explained that ~60-80 unique genes are expected when requesting 10 per cluster

## Benefits

1. **Transparency**: Agent now explicitly tells users how many unique genes are displayed
2. **Correctness**: Properly handles duplicate genes instead of silently removing them
3. **Flexibility**: Three different extraction strategies for different use cases
4. **Guidance**: System prompt guides agent to use the correct workflow

## Testing

The fix should be tested with:
1. Request "show me top 10 marker genes of each cluster"
2. Verify agent uses `get_top_markers` before creating heatmap
3. Verify agent explains the unique gene count in response
4. Verify heatmap displays correctly with all unique genes

## Files Changed

- `src/analysis/marker_genes.py` (NEW): Helper functions for marker gene extraction
- `src/agent/tools.py`: Added `get_top_markers` tool and import
- `src/agent/prompts.py`: Added workflow guidance and example interaction
