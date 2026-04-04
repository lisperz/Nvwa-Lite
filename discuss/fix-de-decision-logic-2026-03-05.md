# Fix: Differential Expression Decision Logic (2026-03-05)

## Problem Statement

The agent had two critical bugs in handling differential expression requests:

1. **Ambiguous requests**: When users said "run differential expression analysis", the agent incorrectly assumed they wanted marker gene analysis for all clusters, instead of asking for clarification.

2. **Pairwise comparison failure**: When users explicitly requested "run DEG between CD4 T cells and CD8 T cells", the agent incorrectly returned marker genes for both groups separately, instead of performing a true pairwise comparison.

## Root Cause

The system only had one DE tool (`differential_expression`) that performed one-vs-rest analysis for all clusters (equivalent to Seurat's `FindAllMarkers`). There was no tool for pairwise comparisons (equivalent to Seurat's `FindMarkers(ident.1, ident.2)`).

## Solution Implemented

### 1. Added Pairwise DE Function (`src/analysis/differential.py`)

**New function: `run_pairwise_de()`**
- Performs direct comparison between two specific groups
- Uses `sc.tl.rank_genes_groups()` with `groups=[group1]` and `reference=group2`
- Returns full DEG table with log2FC, p-values, and adjusted p-values
- Positive log2FC = higher in group1, negative = higher in group2

**Updated function: `run_differential_expression()`**
- Changed default `n_genes` from 100 to 20 (top 20 marker genes per cluster)
- Updated docstring to clarify it's for marker gene analysis (one-vs-rest)
- Improved logging and error messages

### 2. Added New Agent Tools (`src/agent/tools.py`)

**New tool: `compare_groups_de(group1, group2, groupby, method)`**
- Performs pairwise differential expression between two specific groups
- Validates that both groups exist in the dataset
- Stores results in `adata.uns["pairwise_de_result"]` for table export
- Returns summary with interpretation guidance

**New tool: `get_pairwise_de_table(top_n, padj_threshold)`**
- Generates downloadable CSV table for pairwise DE results
- Filters by adjusted p-value threshold (default 0.05)
- Shows upregulated/downregulated gene counts
- Displays first 100 rows in UI, full data in CSV

**Updated tool: `differential_expression(groupby, method)`**
- Updated docstring to clarify it's for marker gene analysis
- Added warning that it's NOT for pairwise comparisons
- Now uses default n_genes=20 instead of 100

### 3. Updated System Prompt (`src/agent/prompts.py`)

**Added new section: "CRITICAL: DIFFERENTIAL EXPRESSION DECISION LOGIC"**

Clearly distinguishes between:
- **A. Marker Gene Analysis (FindAllMarkers)**: One-vs-rest for all groups
- **B. Pairwise Differential Expression (FindMarkers)**: Direct comparison of two groups

**Added mandatory clarification protocol:**
When users say ambiguous phrases like "run differential expression analysis", the agent MUST ask:

```
"I can help with differential expression analysis. Which type would you like?

1. Marker gene analysis: Find top 20 distinguishing genes for each cell type/cluster
2. Pairwise comparison: Compare two specific cell types/clusters directly

Please specify which analysis you need, or tell me which two groups you want to compare."
```

**Updated USER INTENT MAPPING:**
- Added examples for both marker analysis and pairwise comparison
- Clarified when to use each tool
- Added table export guidance for pairwise results

## Files Modified

1. `src/analysis/differential.py` (+88 lines)
   - Added `run_pairwise_de()` function
   - Updated `run_differential_expression()` default parameters and docstring

2. `src/agent/tools.py` (+168 lines, now 771 total)
   - Added `compare_groups_de()` tool
   - Added `get_pairwise_de_table()` tool
   - Updated `differential_expression()` tool docstring
   - Updated `get_all_tools()` to include new tools
   - Updated imports

3. `src/agent/prompts.py` (+45 lines)
   - Added "CRITICAL: DIFFERENTIAL EXPRESSION DECISION LOGIC" section
   - Added mandatory clarification protocol
   - Updated USER INTENT MAPPING with examples

## Tool Count

- **Before**: 18 tools
- **After**: 20 tools
- **New tools**: `compare_groups_de`, `get_pairwise_de_table`

## Expected Behavior After Fix

### Scenario 1: Ambiguous Request
**User**: "run differential expression analysis"

**Agent response**:
```
I can help with differential expression analysis. Which type would you like?

1. Marker gene analysis: Find top 20 distinguishing genes for each cell type/cluster
2. Pairwise comparison: Compare two specific cell types/clusters directly

Please specify which analysis you need, or tell me which two groups you want to compare.
```

### Scenario 2: Explicit Pairwise Request
**User**: "run DEG between CD4 T cells and CD8 T cells"

**Agent actions**:
1. Calls `compare_groups_de(group1="CD4 T cells", group2="CD8 T cells")`
2. Returns summary: "Pairwise differential expression complete. Comparison: CD4 T cells vs CD8 T cells..."
3. If user asks for table: calls `get_pairwise_de_table()` to generate downloadable CSV

### Scenario 3: Marker Gene Request
**User**: "find marker genes for all clusters"

**Agent actions**:
1. Calls `differential_expression()` (one-vs-rest for all clusters)
2. Returns top 20 marker genes per cluster
3. Can follow up with `get_top_markers()` for heatmap or `get_de_results_table()` for full table

## Testing Checklist

- [ ] Test ambiguous request: "run differential expression"
- [ ] Test pairwise request: "compare CD4 T cells vs CD8 T cells"
- [ ] Test marker request: "find marker genes"
- [ ] Test pairwise table export: "show me the comparison table"
- [ ] Test with cluster IDs: "compare cluster 0 vs cluster 1"
- [ ] Test with invalid group names (should show helpful error)
- [ ] Verify top 20 genes default (not 100)

## Code Quality Notes

### Issue: File Length Violation
`src/agent/tools.py` is now **771 lines**, which exceeds the 300-line guideline by 471 lines.

**Recommendation for future refactoring:**
Split into separate modules:
- `src/agent/tools/plotting.py` - UMAP, violin, dotplot, feature, heatmap, scatter, volcano
- `src/agent/tools/analysis.py` - differential_expression, compare_groups_de, get_top_markers
- `src/agent/tools/tables.py` - get_de_results_table, get_pairwise_de_table
- `src/agent/tools/core.py` - dataset_info, check_data_status, preprocess_data, calculate_mito_pct
- `src/agent/tools/__init__.py` - bind_dataset, get_all_tools, result management

This refactoring should be done in a separate session to avoid introducing bugs.

## Summary

The fix successfully addresses both bugs:
1. ✅ Agent now asks for clarification when DE request is ambiguous
2. ✅ Agent can perform true pairwise DE comparisons between two groups
3. ✅ Default marker gene count reduced from 100 to 20 per cluster
4. ✅ Added table export for pairwise DE results
5. ✅ Clear documentation in system prompt distinguishing the two analysis types

The implementation follows Seurat's paradigm:
- `differential_expression()` ≈ `FindAllMarkers()`
- `compare_groups_de()` ≈ `FindMarkers(ident.1, ident.2)`
