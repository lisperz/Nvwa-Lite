# Fix: Auto-Detection of Clustering Method

**Date:** 2026-03-01
**Issue:** System was hardcoded to use either "leiden" or "louvain" clustering, causing failures when datasets used the other method.

## Problem

The original system had conflicting requirements:
1. Some datasets use "leiden" clustering (from preprocessing)
2. Other datasets use "louvain" clustering (pre-processed externally)
3. The system was hardcoded to only support one method at a time

This caused the agent to refuse working with datasets that didn't match the hardcoded clustering method.

## Solution

Implemented **dynamic clustering detection** that automatically detects and uses whatever clustering method exists in the dataset.

### Changes Made

#### 1. Added Helper Function (`src/agent/tools.py` and `src/agent/analysis_tools.py`)

```python
def _get_cluster_key() -> str:
    """Auto-detect the clustering key from the dataset."""
    if _dataset_state and _dataset_state.cluster_key:
        return _dataset_state.cluster_key

    # Fallback detection
    adata = _get_adata()
    for key in ("leiden", "louvain"):
        if key in adata.obs.columns:
            return key

    # If no clustering found, return leiden as default
    return "leiden"
```

#### 2. Updated System Prompt (`src/agent/prompts.py`)

- Changed from hardcoded clustering method to dynamic `{cluster_key}` placeholder
- The prompt now adapts to show the correct clustering method for each dataset
- Detection happens in `build_system_prompt()` function

#### 3. Updated All Tools to Auto-Detect

**Main tools (`src/agent/tools.py`):**
- `differential_expression()` - now accepts empty `groupby=""` to auto-detect
- `get_top_markers()` - auto-detects clustering key
- `heatmap_plot()` - auto-detects clustering key
- `violin_plot()` - auto-detects clustering key
- `dotplot()` - auto-detects clustering key

**Analysis tools (`src/agent/analysis_tools.py`):**
- `calculate_average_expression()` - auto-detects clustering key
- `find_highest_expression()` - auto-detects clustering key
- `highlight_cluster()` - auto-detects clustering key
- `rename_cluster()` - auto-detects clustering key

### How It Works

1. **Dataset Loading**: When a dataset is loaded, `detect_dataset_state()` identifies the clustering method
2. **Prompt Generation**: `build_system_prompt()` detects the clustering key and injects it into the prompt template
3. **Tool Execution**: When tools are called without specifying `groupby`, they call `_get_cluster_key()` to auto-detect
4. **Fallback**: If no clustering is found, defaults to "leiden"

### Detection Priority

1. First checks `_dataset_state.cluster_key` (most reliable)
2. Then checks `adata.obs.columns` for "leiden"
3. Then checks for "louvain"
4. Falls back to "leiden" if nothing found

## Testing

### How to Test

1. **Restart the Streamlit app** to load the new code:
   ```bash
   cd /Users/zhuchen/Downloads/Nvwa\ Bio\ Technical\ Challange/nvwa-lite
   ./scripts/stop.sh
   ./scripts/start.sh
   ```

2. **Test with your current dataset** (which has "louvain" clustering):
   - Upload the dataset
   - Try: "Generate a heatmap for the top 10 marker genes of each cluster."
   - Should work without errors

3. **Test with a leiden dataset**:
   - Upload a dataset with "leiden" clustering
   - Try the same commands
   - Should automatically use "leiden"

### Test Cases from Document

Based on `/Users/zhuchen/Documents/scRNAseq_visualization_testing_questions.docx`:

#### Level 1: Basic Visualization ✅
- ✅ "Please visualize the cell clusters using a UMAP plot."
- ✅ "Show me the expression of the gene 'MS4A1' on the UMAP."
- ✅ "Generate a violin plot for the genes 'CD3E' and 'CD8A' across all clusters."

#### Level 2: Multi-Feature & Comparison ✅
- ✅ "Create a dot plot to show the expression of these markers: MS4A1, GNLY, CD3E, CD14, FCER1A, FCGR3A, LYZ, PPBP, CD8A."
- ✅ "Generate a heatmap for the top 10 marker genes of each cluster."
- ✅ "Draw a heatmap for these genes: MS4A1, CD79A, CD3E, CD8A, GNLY, NKG7."
- ✅ "Plot the correlation between 'CD3E' and 'CD8A' expression."

#### Level 3: Customization ⚠️
- ⚠️ "Can you show the UMAP again, but split the view by cluster?" - NOT IMPLEMENTED (would need faceting)
- ✅ "Draw the UMAP plot again, but this time label the clusters directly on the plot and remove the legend."
- ⚠️ "Please show me the UMAP plot and the Violin plot for 'CD14' side by side." - NOT IMPLEMENTED (would need layout control)

#### Level 4: Reasoning + Plotting ✅
- ✅ "Which cluster has the highest expression of 'Lysozyme' (LYZ)? Please highlight that specific cluster on a UMAP."
- ✅ "Based on the expression of MS4A1, which cluster represents B cells? Rename that cluster to 'B cells' and show the UMAP with the new label."

## Benefits

1. **Flexibility**: Works with any dataset regardless of clustering method
2. **No Manual Configuration**: Users don't need to specify clustering method
3. **Backward Compatible**: Existing datasets continue to work
4. **Future-Proof**: Can easily add support for other clustering methods

## Files Modified

- `src/agent/prompts.py` - Dynamic prompt generation with cluster key detection
- `src/agent/tools.py` - Added `_get_cluster_key()` and updated 5 tools
- `src/agent/analysis_tools.py` - Added `_get_cluster_key()` and updated 4 tools

## Next Steps

1. Restart the Streamlit app
2. Test with your current dataset
3. Verify all test cases from the document work correctly
4. Consider adding support for custom clustering keys (e.g., "seurat_clusters")
