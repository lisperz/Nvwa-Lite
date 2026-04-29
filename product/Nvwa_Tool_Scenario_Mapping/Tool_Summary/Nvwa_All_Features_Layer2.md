# Nvwa Bio — All Features
## Layer 2: Workflow Spec (Tool Layer)
*April 2026 | Updated after validation | For Chen*

---

## Global Setup

```python
# Run once before all visualizations
sc.settings.set_figure_params(dpi=150, facecolor='white')
# Do NOT set figsize globally — each tool manages its own figsize
```

## Data Ingestion (run once after user uploads h5ad)

```python
# Ensure adata.X contains logcounts
adata.X = adata.layers['logcounts']

# Ensure UMAP coordinates are in correct format
adata.obsm['X_umap'] = np.array(adata.obsm['X_umap'])

# Calculate QC metrics if not present
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
```

---

## Global Parameters

- `CELLTYPE_COL` — obs column for cell type annotation
- `CONDITION_COL` — obs column for condition/sample grouping

## Condition Color Standard

- **Reference/Control** → Blue `#2E86C1`
- **Treatment/Experimental** → Red `#E74C3C`
- **3rd condition** → Green `#27AE60`
- **4th condition** → Purple `#8E44AD`

---

## Shared Reusable Tools

| Tool | Parameters | Description |
|--|--|--|
| `subset_data` | `adata, col, value` | Subset adata to a single value |
| `subset_data` | `adata, col, values=[list]` | Subset adata to multiple values |
| `create_combined_column` | `adata, col1, col2` | Creates `col1_col2` obs column |

---

## 1. UMAP Visualization

**Default:** tab20 colormap, `frameon=False`

### Scenario 1.1 — All cell types · No condition split
- **Type:** Single tool call
- **Tool:** `umap_plot`
- **Parameters:** `adata, color=CELLTYPE_COL`

### Scenario 1.2 — All cell types · Split by condition
- **Type:** Single tool call
- **Tool:** `umap_plot_split`
- **Parameters:** `adata, color=CELLTYPE_COL, CONDITION_COL`
- **Internal logic:** `get_umap_range(adata)` → unified axis range, subplot merge with figsize `(5 * n_conditions, 5)`

### Scenario 1.3 — Subset specific cell types · No condition split
- **Type:** Workflow
- **Steps:**
  1. `subset_data(adata, CELLTYPE_COL, values=[list])`
  2. `umap_plot(adata_sub, color=CELLTYPE_COL)`

### Scenario 1.4 — Subset specific cell types · Split by condition
- **Type:** Workflow
- **Steps:**
  1. `subset_data(adata, CELLTYPE_COL, values=[list])`
  2. `umap_plot_split(adata_sub, color=CELLTYPE_COL, CONDITION_COL)`

---

## 2. QC Analysis & Visualization

**QC Metrics:** `['nFeature_RNA', 'nCount_RNA', 'pct_counts_mt']`

### Scenario 2A.1 — All cells · No condition split
- **Type:** Single tool call
- **Tool:** `qc_violin_plot`
- **Parameters:** `adata`
- **Internal:** `multi_panel=True`, x-axis labels `ha='right'`

### Scenario 2A.2 — All cells · Split by condition
- **Type:** Single tool call
- **Tool:** `qc_violin_plot_split`
- **Parameters:** `adata, CONDITION_COL`
- **Internal:** One multi-panel figure per condition, `fig.suptitle(condition)`

### Scenario 2B.1 — QC Summary Table
- **Type:** Single tool call
- **Tool:** `qc_summary_table`
- **Parameters:** `adata`
- **Output:** Table with min, max, mean, median per metric + automated flagging
- **Flagging thresholds (user-configurable):**
  - `pct_counts_mt` median > 10% → flag
  - `nFeature_RNA` median < 200 → flag

---

## 3. Gene Expression Visualization

### 3A. Feature Plot

**Default colormap:** `Reds`

#### Scenario 3A.1 — Single gene · All cells · No condition split
- **Type:** Single tool call
- **Tool:** `feature_plot`
- **Parameters:** `gene, adata`

#### Scenario 3A.2 — Single gene · All cells · Split by condition
- **Type:** Single tool call
- **Tool:** `feature_plot_split`
- **Parameters:** `gene, adata, CONDITION_COL`
- **Internal:** `get_expression_range` → unified colorbar, subplot merge `(5 * n, 5)`

#### Scenario 3A.3 — Single gene · Subset one cell type · No condition split
- **Type:** Workflow
- **Steps:**
  1. `subset_data(adata, CELLTYPE_COL, celltype)`
  2. `feature_plot(gene, adata_sub)`

#### Scenario 3A.4 — Single gene · Subset one cell type · Split by condition
- **Type:** Workflow
- **Steps:**
  1. `subset_data(adata, CELLTYPE_COL, celltype)`
  2. `feature_plot_split(gene, adata_sub, CONDITION_COL)`
- **Internal:** `get_expression_range` + `get_umap_range` → unified scale and axis range

---

### 3B. Violin Plot

**Default figsize:** `(8, 4)` single, `(15, 6)` split
**Default:** x-axis labels `rotation=45, ha='right'`, `set_xlabel('')`, title auto-generated, `set_ylim(bottom=0)`

#### Scenario 3B.1 — Single gene · X-axis: cell type · No condition split
- **Type:** Single tool call
- **Tool:** `violin_plot`
- **Parameters:** `gene, adata, groupby=CELLTYPE_COL`

#### Scenario 3B.2 — Single gene · X-axis: cell type · Split by condition
- **Type:** Single tool call
- **Tool:** `violin_plot_split`
- **Parameters:** `gene, adata, CELLTYPE_COL, CONDITION_COL`
- **Internal:** seaborn violinplot with `hue=condition`

#### Scenario 3B.3 — Single gene · X-axis: condition
- **Type:** Single tool call
- **Tool:** `violin_plot`
- **Parameters:** `gene, adata, groupby=CONDITION_COL`

#### Scenario 3B.4 — Single gene · Subset one cell type · X-axis: condition
- **Type:** Workflow
- **Steps:**
  1. `subset_data(adata, CELLTYPE_COL, celltype)`
  2. `violin_plot(gene, adata_sub, groupby=CONDITION_COL)`

#### Scenario 3B.5 — Single gene · Subset multiple cell types · Split by condition
- **Type:** Workflow
- **Steps:**
  1. `subset_data(adata, CELLTYPE_COL, values=[list])`
  2. `violin_plot_split(gene, adata_sub, CELLTYPE_COL, CONDITION_COL)`

#### Scenario 3B.6 — Single gene · Subset one condition · X-axis: cell type
- **Type:** Workflow
- **Steps:**
  1. `subset_data(adata, CONDITION_COL, condition)`
  2. `violin_plot(gene, adata_sub, groupby=CELLTYPE_COL)`

---

### 3C. Dot Plot

#### Scenario 3C.1 — Multiple genes · All cell types · No condition split
- **Type:** Single tool call
- **Tool:** `dot_plot`
- **Parameters:** `genes, adata, groupby=CELLTYPE_COL`

#### Scenario 3C.2 — Multiple genes · All cell types · Split by condition
- **Type:** Workflow
- **Steps:**
  1. `create_combined_column(adata, CELLTYPE_COL, CONDITION_COL)`
  2. `dot_plot(genes, adata, groupby=celltype_condition)`
  3. `color_yaxis_by_condition(ax, CONDITION_COL)`

#### Scenario 3C.3 — Multiple genes · Subset specific cell types · No condition split
- **Type:** Workflow
- **Steps:**
  1. `subset_data(adata, CELLTYPE_COL, values=[list])`
  2. `dot_plot(genes, adata_sub, groupby=CELLTYPE_COL)`

#### Scenario 3C.4 — Multiple genes · Subset specific cell types · Split by condition
- **Type:** Workflow
- **Steps:**
  1. `subset_data(adata, CELLTYPE_COL, values=[list])`
  2. `create_combined_column(adata_sub, CELLTYPE_COL, CONDITION_COL)`
  3. `dot_plot(genes, adata_sub, groupby=celltype_condition)`
  4. `color_yaxis_by_condition(ax, CONDITION_COL)`

---

### 3D. Heatmap

**Default colormap:** `viridis`
**Default:** `swap_axes=True`, fontsize=6

#### Scenario 3D.1 — Multiple genes · Grouped by cell type · No condition split
- **Type:** Single tool call
- **Tool:** `heatmap_plot`
- **Parameters:** `genes, adata, groupby=CELLTYPE_COL`

#### Scenario 3D.2 — Multiple genes · Grouped by cell type · Split by condition
- **Type:** Single tool call
- **Tool:** `heatmap_plot_split`
- **Parameters:** `genes, adata, CELLTYPE_COL, CONDITION_COL`
- **Internal:** `create_combined_column` → find groupby ax → `color_xaxis_by_condition`

#### Scenario 3D.3 — Multiple genes · Subset specific cell types · No condition split
- **Type:** Workflow
- **Steps:**
  1. `subset_data(adata, CELLTYPE_COL, values=[list])`
  2. `heatmap_plot(genes, adata_sub, groupby=CELLTYPE_COL)`

#### Scenario 3D.4 — Multiple genes · Subset specific cell types · Split by condition
- **Type:** Workflow
- **Steps:**
  1. `subset_data(adata, CELLTYPE_COL, values=[list])`
  2. `heatmap_plot_split(genes, adata_sub, CELLTYPE_COL, CONDITION_COL)`

#### Scenario 3D.5 — Multiple genes · Grouped by condition · Individual cell level
- **Type:** Single tool call
- **Tool:** `heatmap_plot`
- **Parameters:** `genes, adata, groupby=CONDITION_COL`
- **Note:** Individual cell level, not aggregated

---

## 4. Find Markers

**Default:** `only_pos=True`, `pvals_adj < 0.05`, ranked by `logfoldchanges`

### Scenario 4.1 — Find marker genes for all cell types
- **Type:** Single tool call
- **Tool:** `find_all_markers`
- **Parameters:** `adata, CELLTYPE_COL, pvals_adj_threshold=0.05, logfc_threshold=0, n_top_genes=10`
- **Output:** Print top n genes per cell type + CSV of all significant markers

### Scenario 4.2 — Find marker genes for a specific cell type
- **Type:** Single tool call
- **Tool:** `find_markers`
- **Parameters:** `adata, CELLTYPE_COL, celltype, pvals_adj_threshold=0.05, logfc_threshold=0, n_top_genes=10`
- **Output:** Print top n genes + CSV

### Scenario 4.3 — Show top marker genes as dot plot
- **Type:** Workflow (user-triggered)
- **Steps:**
  1. `find_all_markers(adata, CELLTYPE_COL)` — if not already run
  2. `dot_plot_top_markers(adata, markers_df, CELLTYPE_COL, n_genes=5)`
- **Note:** `n_genes` user-configurable, default 5. Genes ordered by cell type order in adata.

---

## 5. Differential Expression Analysis

**Default:** Wilcoxon rank-sum, `pvals_adj < 0.05`, `lfc_threshold=0.25`
**Statistical disclaimer:** Always displayed after DE results.
**Post-DE:** Agent automatically asks if user wants volcano plot.

### Scenario 5.1 — Cell type vs cell type
- **Type:** Single tool call
- **Tool:** `run_de`
- **Parameters:** `adata, CELLTYPE_COL, group1, group2, pvals_adj_threshold=0.05`
- **Output:** Summary (total, up in group1, up in group2) + top 10 genes per group + CSV

### Scenario 5.2 — Condition vs condition · All cell types
- **Type:** Single tool call
- **Tool:** `run_de`
- **Parameters:** `adata, CONDITION_COL, group1, group2, pvals_adj_threshold=0.05`

### Scenario 5.3 — Condition vs condition · Single cell type
- **Type:** Workflow
- **Steps:**
  1. `subset_data(adata, CELLTYPE_COL, celltype)`
  2. `create_combined_column(adata_sub, CELLTYPE_COL, CONDITION_COL)`
  3. `run_de(adata_sub, groupby=celltype_condition, group1, group2)`

### Scenario 5.4 — Condition vs condition · Complex grouping
- **Type:** Workflow
- **Note:** Not validated on ifnb — requires dataset with complex condition structure
- **Steps:**
  1. `create_group_column(adata, CONDITION_COL, group_map={group1: [list], group2: [list]})`
  2. `run_de(adata, groupby=de_group, group1, group2)`
- **Note:** Agent probes user for condition structure before running

### Volcano Plot
- **Type:** Single tool call
- **Tool:** `volcano_plot`
- **Parameters:** `de_df, group1, group2, pvals_adj_threshold=0.05, lfc_threshold=0.25, lfc_clip=6, n_label=10`
- **Triggered by:** Agent prompt after any DE analysis

---

## 6. Cell Type Composition Analysis

**Note:** Agent checks for condition column before running any composition analysis.

### Scenario 6A.1 — All cell types · Absolute cell number
- **Type:** Single tool call
- **Tool:** `composition_barplot`
- **Parameters:** `adata, CELLTYPE_COL, CONDITION_COL, mode='count'`
- **Figsize:** `(12, 5)`

### Scenario 6A.2 — All cell types · Proportion per cell type
- **Type:** Single tool call
- **Tool:** `composition_barplot`
- **Parameters:** `adata, CELLTYPE_COL, CONDITION_COL, mode='proportion_by_celltype'`
- **Figsize:** `(12, 5)`, `set_ylim(0, 1.1)`

### Scenario 6A.3 — All cell types · Proportion per condition
- **Type:** Single tool call
- **Tool:** `composition_barplot`
- **Parameters:** `adata, CELLTYPE_COL, CONDITION_COL, mode='proportion_by_condition'`
- **Figsize:** dynamic `(4, max(4, n_celltypes * 0.3))`, `set_xlim(-0.5, n_conditions - 0.5)`

### Scenario 6B.1 — Subset cell types · Absolute cell number
- **Type:** Workflow
- **Steps:**
  1. `subset_data(adata, CELLTYPE_COL, values=[list])`
  2. `composition_barplot(adata_sub, CELLTYPE_COL, CONDITION_COL, mode='count')`

### Scenario 6B.2 — Subset cell types · Proportion per cell type
- **Type:** Workflow
- **Steps:**
  1. `subset_data(adata, CELLTYPE_COL, values=[list])`
  2. `composition_barplot(adata_sub, CELLTYPE_COL, CONDITION_COL, mode='proportion_by_celltype')`

### Scenario 6B.3 — Subset cell types · Proportion per condition
- **Type:** Workflow
- **Steps:**
  1. `subset_data(adata, CELLTYPE_COL, values=[list])`
  2. `composition_barplot(adata_sub, CELLTYPE_COL, CONDITION_COL, mode='proportion_by_condition')`

---
