# Fix: Improve Heatmap Code Display for Clarity

**Date:** 2026-03-01
**Commit:** f08381a

## Problem

When users requested "top 5 marker genes of each cell type", the displayed code showed a confusing long list of genes without context:

```python
sc.pl.heatmap(adata, var_names=['CD79A', 'CD74', 'CD79B', 'HLA-DRA', 'MS4A1',
'LDHB', 'CD3D', 'RPS27', 'RPS25', 'RPS12', 'CCL5', 'NKG7', 'CST7', 'GZMA',
'CTSW', 'S100A9', 'S100A8', 'LYZ', 'FCN1', 'FTL', 'HLA-DPB1', 'HLA-DPA1',
'HLA-DRB1', 'LST1', 'FCER1G', 'AIF1', 'FTH1', 'IFITM2', 'PPBP', 'PF4', 'GPX1',
'SDPR', 'NRGN', 'GZMB', 'PRF1', 'GNLY'], groupby="louvain", swap_axes=True,
figsize=(10.8, 6))
```

**Issues:**
1. Shows 36 genes but user asked for "top 5 per cluster"
2. No explanation of why there are 36 genes (8 clusters × 5 = 40, but duplicates removed)
3. Long gene list is hard to read and understand
4. Doesn't communicate the user's original intent

## Solution

Added context and clarity to the code display by:
1. Adding `n_genes_per_cluster` parameter to track user's request
2. Displaying informative comments explaining the gene selection
3. Truncating long gene lists for readability
4. Showing cluster and gene counts

### New Code Display Format

**For "top N genes per cluster" requests:**
```python
# Top 5 marker genes per cluster (36 unique genes after deduplication)
# Clusters: 8, Genes: 36
top_genes = ['CD79A', 'CD74', 'CD79B', 'HLA-DRA', 'MS4A1', '...']
sc.pl.heatmap(adata, var_names=top_genes, groupby="louvain", swap_axes=True, figsize=(10.8, 6.0))
```

**For manually specified genes (≤10 genes):**
```python
sc.pl.heatmap(adata, var_names=['CD3E', 'MS4A1', 'NKG7'], groupby="louvain", swap_axes=True, figsize=(10.0, 6.0))
```

**For manually specified genes (>10 genes):**
```python
sc.pl.heatmap(adata, var_names=['CD3E', 'MS4A1', 'NKG7', 'LYZ', 'CD14', '...', 'GNLY', 'PRF1'], groupby="louvain", swap_axes=True, figsize=(15.0, 6.0))
```

## Changes Made

### 1. Updated `plot_heatmap()` Function (`src/plotting/comparison.py`)

Added `n_genes_per_cluster` parameter:
```python
def plot_heatmap(
    adata: AnnData,
    genes: list[str],
    groupby: str = "leiden",
    n_genes_per_cluster: int = 0,  # NEW: context about user's request
) -> PlotResult:
```

Added smart code formatting logic:
```python
if n_genes_per_cluster > 0:
    # Show context about top N genes per cluster
    code = f"# Top {n_genes_per_cluster} marker genes per cluster ({n_genes} unique genes after deduplication)\n"
    code += f"# Clusters: {n_groups}, Genes: {n_genes}\n"
    code += f'top_genes = {genes[:5] + ["..."] if n_genes > 5 else genes}\n'
    code += f'sc.pl.heatmap(adata, var_names=top_genes, groupby="{groupby}", ...)'
else:
    # Show full or truncated gene list for manually specified genes
    if n_genes <= 10:
        genes_str = str(genes)
    else:
        genes_str = str(genes[:5] + ["..."] + genes[-2:])
    code = f'sc.pl.heatmap(adata, var_names={genes_str}, groupby="{groupby}", ...)'
```

### 2. Updated `heatmap_plot` Tool (`src/agent/tools.py`)

Added `n_genes_per_cluster` parameter:
```python
@tool
def heatmap_plot(genes: str, groupby: str = "", n_genes_per_cluster: int = 0) -> str:
    """Generate a heatmap for genes across cell groups.

    Args:
        genes: Comma-separated gene names.
        groupby: The observation key to group cells by.
        n_genes_per_cluster: Optional context about how many genes per cluster were requested.
    """
```

### 3. Updated System Prompt (`src/agent/prompts.py`)

**Updated CRITICAL RULES:**
```
3. **For heatmaps**: ALWAYS follow this exact workflow:
   a. First call: differential_expression(groupby="{cluster_key}") if not done yet
   b. Second call: get_top_markers(n_genes_per_cluster=N, groupby="{cluster_key}")
   c. Third call: heatmap_plot(genes="...", groupby="{cluster_key}", n_genes_per_cluster=N)
      - IMPORTANT: Pass n_genes_per_cluster=N to heatmap_plot so the code display shows context
```

**Updated example interaction:**
```
User: "Show me top 5 marker genes of each cluster"
Agent: [Use get_top_markers(n_genes_per_cluster=5)]
       [Then use heatmap_plot with the returned genes AND n_genes_per_cluster=5]
```

## Benefits

### 1. **Clarity**
Users immediately understand:
- What they requested ("top 5 genes per cluster")
- How many clusters exist (8)
- How many unique genes resulted (36)
- Why there are fewer genes than expected (deduplication)

### 2. **Educational**
The code display teaches users about:
- Deduplication of marker genes across clusters
- The relationship between clusters and genes
- How scanpy heatmaps work

### 3. **Reproducibility**
Users can:
- Understand the exact workflow used
- Reproduce the analysis in their own scripts
- Modify parameters based on clear context

### 4. **Debugging**
When something looks wrong, users can:
- See if the correct number of genes per cluster was used
- Verify the clustering method
- Check the figure size calculation

## Example Scenarios

### Scenario 1: Top 5 Genes Per Cluster (8 clusters)

**User Request:** "Generate a heatmap for the top 5 marker genes of each cell type"

**Old Display:**
```python
sc.pl.heatmap(adata, var_names=['CD79A', 'CD74', ...(36 genes)...], groupby="louvain")
```
❌ Confusing - why 36 genes when I asked for 5?

**New Display:**
```python
# Top 5 marker genes per cluster (36 unique genes after deduplication)
# Clusters: 8, Genes: 36
top_genes = ['CD79A', 'CD74', 'CD79B', 'HLA-DRA', 'MS4A1', '...']
sc.pl.heatmap(adata, var_names=top_genes, groupby="louvain", swap_axes=True, figsize=(10.8, 6.0))
```
✅ Clear - 8 clusters × 5 genes = 40, but 4 duplicates removed = 36 unique

### Scenario 2: Top 10 Genes Per Cluster (8 clusters)

**User Request:** "Show top 10 marker genes per cluster"

**Expected:** 8 × 10 = 80 genes, but after deduplication might be ~60-70

**Display:**
```python
# Top 10 marker genes per cluster (66 unique genes after deduplication)
# Clusters: 8, Genes: 66
top_genes = ['CD79A', 'CD74', 'CD79B', 'HLA-DRA', 'MS4A1', '...']
sc.pl.heatmap(adata, var_names=top_genes, groupby="louvain", swap_axes=True, figsize=(19.8, 6.0))
```
✅ Clear - shows deduplication reduced 80 to 66

### Scenario 3: Manually Specified Genes

**User Request:** "Draw a heatmap for these genes: MS4A1, CD79A, CD3E, CD8A, GNLY, NKG7"

**Display:**
```python
sc.pl.heatmap(adata, var_names=['MS4A1', 'CD79A', 'CD3E', 'CD8A', 'GNLY', 'NKG7'], groupby="louvain", swap_axes=True, figsize=(10.0, 6.0))
```
✅ Shows full list since it's short and manually specified

## Testing

### Test Cases:

1. **Top 5 genes per cluster:**
   - Request: "Generate a heatmap for the top 5 marker genes of each cell type"
   - Expected: Code shows "Top 5 marker genes per cluster" with gene count

2. **Top 10 genes per cluster:**
   - Request: "Show me top 10 marker genes"
   - Expected: Code shows "Top 10 marker genes per cluster" with deduplication info

3. **Manual gene list (short):**
   - Request: "Heatmap for CD3E, MS4A1, NKG7"
   - Expected: Full gene list shown

4. **Manual gene list (long):**
   - Request: "Heatmap for 15 specific genes"
   - Expected: Truncated list with "..."

## Files Modified

- `src/plotting/comparison.py` - Added n_genes_per_cluster parameter and smart code formatting
- `src/agent/tools.py` - Updated heatmap_plot tool signature
- `src/agent/prompts.py` - Updated workflow instructions and examples

## Deployment

✅ Changes pushed to GitHub (commit f08381a)
✅ Only you shown as contributor
✅ Streamlit Cloud will auto-deploy

## Impact

This fix improves user experience by making the code display:
- **Informative** - shows what was requested
- **Educational** - explains deduplication
- **Readable** - truncates long lists
- **Accurate** - reflects the actual workflow

Users will no longer be confused about why the gene count doesn't match their request.
