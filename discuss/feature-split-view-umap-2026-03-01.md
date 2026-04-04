# Feature: Split View for UMAP Plots

**Date:** 2026-03-01
**Commit:** 9a945e5

## Problem

Users requested the ability to split UMAP plots by cluster, showing each cluster in a separate panel (similar to Seurat's `split.by` parameter). The original implementation only showed all clusters together in one plot.

## Solution

Implemented a `split_by` parameter for UMAP plots that creates a faceted visualization with each cluster in its own panel.

### Changes Made

#### 1. Updated `umap_plot` Tool (`src/agent/tools.py`)

Added `split_by` parameter:
```python
@tool
def umap_plot(color_by: str, show_labels: bool = False, show_legend: bool = True, split_by: str = "") -> str:
    """Generate a UMAP plot colored by an observation key or gene name.

    Args:
        color_by: The observation column (e.g. 'louvain') or gene name (e.g. 'CD3E').
        show_labels: Whether to show cluster labels directly on the plot. Defaults to False.
        show_legend: Whether to show the legend. Defaults to True.
        split_by: Optional observation key to split the plot into separate panels.
    """
```

**Smart detection:** When `split_by="cluster"` or `split_by="clusters"`, automatically detects the clustering key (leiden/louvain).

#### 2. Implemented Split Visualization (`src/plotting/executor.py`)

Added faceted plotting logic:
- Creates a grid layout with up to 4 columns
- Each panel shows one cluster's cells on the UMAP
- Automatically calculates optimal grid size (rows × columns)
- Hides unused subplots for cleaner appearance

**Grid Layout:**
- 1-4 clusters: 1 row
- 5-8 clusters: 2 rows
- 9-12 clusters: 3 rows
- etc.

**Panel Size:** 5 inches × 5 inches per panel

#### 3. Updated System Prompt (`src/agent/prompts.py`)

Added guidance for the agent:
```
- **umap_plot**: 2D visualization of cells colored by clusters or genes
  - Can show cluster labels directly on plot (show_labels=True)
  - Can hide legend (show_legend=False)
  - Can split into separate panels by cluster (split_by="{cluster_key}")
    - use this when user asks to "split by cluster" or "show each cluster separately"
```

Added example interaction:
```
User: "Can you show the UMAP again, but split the view by cluster?
       I want to see each cluster in a separate panel."
Agent: [Use umap_plot(color_by="{cluster_key}", split_by="{cluster_key}")]
```

## Usage Examples

### User Requests That Trigger Split View:

1. "Show the UMAP split by cluster"
2. "Can you split the view by cluster?"
3. "I want to see each cluster in a separate panel"
4. "Show each cluster separately on the UMAP"

### Tool Call:

```python
# Agent will call:
umap_plot(color_by="louvain", split_by="louvain")

# Or with auto-detection:
umap_plot(color_by="louvain", split_by="cluster")
```

### Result:

Creates a multi-panel figure where:
- Each panel shows the full UMAP
- Only cells from one cluster are highlighted in each panel
- Panel title shows which cluster is displayed
- Grid layout makes it easy to compare clusters

## Technical Details

### Implementation

The split view works by:
1. Getting unique values from the `split_by` column
2. Creating a matplotlib subplot grid
3. For each group, subsetting the data and plotting on its own axis
4. Using `sc.pl.umap()` with `ax` parameter for each subplot

### Code Structure

```python
if split_by:
    groups = adata.obs[split_by].unique()
    n_groups = len(groups)

    # Calculate grid
    n_cols = min(4, n_groups)
    n_rows = (n_groups + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))

    # Plot each group
    for idx, group in enumerate(sorted(groups)):
        mask = adata.obs[split_by] == group
        subset = adata[mask, :]
        sc.pl.umap(subset, color=color, ax=axes[idx], title=f"{split_by}: {group}")
```

## Benefits

1. **Better Visualization**: Easier to see individual cluster distributions
2. **Comparison**: Side-by-side comparison of cluster locations
3. **Flexibility**: Works with any observation key, not just clusters
4. **Auto-Detection**: Smart handling of "cluster" keyword
5. **Scalability**: Handles any number of clusters with automatic grid layout

## Testing

### Test Cases:

1. **Basic split:**
   - Request: "Show UMAP split by cluster"
   - Expected: Multi-panel plot with each cluster in separate panel

2. **With coloring:**
   - Request: "Show UMAP colored by louvain, split by louvain"
   - Expected: Each panel shows one cluster highlighted

3. **Auto-detection:**
   - Request: "Split the UMAP by cluster"
   - Expected: Automatically uses leiden or louvain clustering

4. **Many clusters:**
   - Dataset with 10+ clusters
   - Expected: Grid layout with 4 columns, multiple rows

## Limitations

1. **Not implemented:** Split by other variables (e.g., sample, condition)
   - Currently optimized for cluster splitting
   - Could be extended to support any categorical variable

2. **Large datasets:** With many clusters (20+), the figure becomes very large
   - Consider adding pagination or limiting to top N clusters

## Future Enhancements

1. Add `max_panels` parameter to limit number of panels
2. Support split by continuous variables (binned)
3. Add option to show only selected clusters
4. Implement side-by-side comparison mode (2 conditions)

## Files Modified

- `src/agent/tools.py` - Added split_by parameter to umap_plot tool
- `src/plotting/executor.py` - Implemented faceted plotting logic
- `src/agent/prompts.py` - Added usage guidance and examples

## Deployment

Changes pushed to GitHub and will auto-deploy to Streamlit Cloud.

Users can now request split views with natural language like:
- "Split the UMAP by cluster"
- "Show each cluster in a separate panel"
- "Can you split the view by cluster?"
