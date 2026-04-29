# Nvwa Bio — All Features
## Layer 1: User Guide — Scenarios & Standard Prompts
*April 2026*

---

This guide shows you how to ask Nvwa for the analysis you want. Each scenario includes a standard prompt you can use or adapt. Bracketed terms like [Gene A], [Cell Type], or [Condition] should be replaced with your actual labels.

---

## 1. UMAP Visualization

UMAP plots show the overall distribution of cells in your dataset, colored by cell type. Use UMAP to visualize cell population structure, assess clustering quality, and compare cell distributions across conditions.

**Scenario 1.1** — All cell types · No condition split
> Show the UMAP colored by cell type.

**Scenario 1.2** — All cell types · Split by condition
> Show the UMAP colored by cell type, split by condition.
> *(Note: "Condition" refers to the grouping variable in your dataset, such as treatment/control, disease/healthy, etc.)*

**Scenario 1.3** — Subset specific cell types · No condition split
> Show the UMAP for [Cell Type A] and [Cell Type B] only.

**Scenario 1.4** — Subset specific cell types · Split by condition
> Show the UMAP for [Cell Type A] and [Cell Type B] only, split by condition.
> *(Note: "Condition" refers to the grouping variable in your dataset, such as treatment/control, disease/healthy, etc.)*

---

## 2. QC Analysis & Visualization

### 2A. QC Violin Plot

QC violin plots display the distribution of three key quality metrics across your cells: nFeature_RNA (genes detected per cell), nCount_RNA (total UMI counts per cell), and percent.mt (mitochondrial gene fraction). Use these plots to assess data quality and verify that filtering thresholds were appropriate.

**Scenario 2A.1** — All cells · No condition split
> Show the QC metrics for all cells.

**Scenario 2A.2** — All cells · Split by condition
> Show the QC metrics for all cells, split by condition.
> *(Note: "Condition" refers to the grouping variable in your dataset, such as treatment/control, disease/healthy, etc.)*

### 2B. QC Summary Table

The QC summary table provides a statistical overview of your dataset's quality metrics, including range, mean, and median for nFeature_RNA, nCount_RNA, and percent.mt. Nvwa will flag any metrics that fall outside typical quality ranges to help you quickly identify potential issues.

**Scenario 2B.1** — Dataset-level QC summary with automated flagging
> Give me a QC summary of my dataset.
> *(Note: Nvwa will flag elevated mitochondrial fraction (percent.mt median > 10–15%) and low gene detection (nFeature_RNA median < 200) as potential quality concerns. These are observations, not pass/fail judgments — you decide how to interpret them in the context of your experiment.)*

---

## 3. Gene Expression Visualization

### 3A. Feature Plot

Feature plots display gene expression on your UMAP, where color intensity indicates expression level. Use Feature Plot to visualize where a gene is expressed across your cell populations.

**Scenario 3A.1** — Single gene · All cells · No condition split
> Show the expression of [Gene A] on the UMAP.

**Scenario 3A.2** — Single gene · All cells · Split by condition
> Show the expression of [Gene A] on the UMAP, split by condition.
> *(Note: "Condition" refers to the grouping variable in your dataset, such as treatment/control, disease/healthy, etc.)*

**Scenario 3A.3** — Single gene · Subset one cell type · No condition split
> Show the expression of [Gene A] in [Cell Type] only.

**Scenario 3A.4** — Single gene · Subset one cell type · Split by condition
> Show the expression of [Gene A] in [Cell Type] only, split by condition.
> *(Note: "Condition" refers to the grouping variable in your dataset, such as treatment/control, disease/healthy, etc.)*

### 3B. Violin Plot

Violin plots show the distribution of gene expression across groups. Use Violin Plot to compare expression levels between cell types or conditions, and to assess expression variability within a population.

**Scenario 3B.1** — Single gene · X-axis: cell type · No condition split · No subset
> Show the expression of [Gene A] across all cell types.

**Scenario 3B.2** — Single gene · X-axis: cell type · Split by condition · No subset
> Show the expression of [Gene A] across all cell types, split by condition.
> *(Note: "Condition" refers to the grouping variable in your dataset, such as treatment/control, disease/healthy, etc.)*

**Scenario 3B.3** — Single gene · X-axis: condition · No subset
> Show the expression of [Gene A] across all conditions.

**Scenario 3B.4** — Single gene · Subset one cell type · X-axis: condition
> Show the expression of [Gene A] in [Cell Type] only, grouped by condition.

**Scenario 3B.5** — Single gene · Subset multiple cell types · X-axis: cell type · Split by condition
> Show the expression of [Gene A] in [Cell Type A], [Cell Type B], and [Cell Type C], split by condition.
> *(Note: "Condition" refers to the grouping variable in your dataset, such as treatment/control, disease/healthy, etc.)*

**Scenario 3B.6** — Single gene · Subset one condition · X-axis: cell type
> Show the expression of [Gene A] across all cell types in [Condition A] only.

### 3C. Dot Plot

Dot plots display multiple genes across cell types simultaneously. Dot size represents the percentage of cells expressing the gene; dot color represents average expression level. Use Dot Plot for an efficient multi-gene overview.

**Scenario 3C.1** — Multiple genes · All cell types · No condition split
> Show a dot plot of [top marker genes / Gene A, Gene B, Gene C] across all cell types.

**Scenario 3C.2** — Multiple genes · All cell types · Split by condition
> Show a dot plot of [top marker genes / Gene A, Gene B, Gene C] across all cell types, split by condition.
> *(Note: "Condition" refers to the grouping variable in your dataset, such as treatment/control, disease/healthy, etc.)*

**Scenario 3C.3** — Multiple genes · Subset specific cell types · No condition split
> Show a dot plot of [Gene A, Gene B, Gene C] in [Cell Type A] and [Cell Type B] only.

**Scenario 3C.4** — Multiple genes · Subset specific cell types · Split by condition
> Show a dot plot of [Gene A, Gene B, Gene C] in [Cell Type A] and [Cell Type B] only, split by condition.
> *(Note: "Condition" refers to the grouping variable in your dataset, such as treatment/control, disease/healthy, etc.)*

### 3D. Heatmap

Heatmaps display expression of multiple genes across cells or groups using color intensity. Unlike dot plots, heatmaps can show individual cell-level variation. Use Heatmap to visualize expression patterns and identify co-expressed gene modules.

**Scenario 3D.1** — Multiple genes · Grouped by cell type · No condition split · No subset
> Show a heatmap of [top marker genes / Gene A, Gene B, Gene C] across all cell types.

**Scenario 3D.2** — Multiple genes · Grouped by cell type · Split by condition · No subset
> Show a heatmap of [top marker genes / Gene A, Gene B, Gene C] across all cell types, split by condition.
> *(Note: "Condition" refers to the grouping variable in your dataset, such as treatment/control, disease/healthy, etc.)*

**Scenario 3D.3** — Multiple genes · Subset specific cell types · No condition split
> Show a heatmap of [Gene A, Gene B, Gene C] in [Cell Type A] and [Cell Type B] only.

**Scenario 3D.4** — Multiple genes · Subset specific cell types · Split by condition
> Show a heatmap of [Gene A, Gene B, Gene C] in [Cell Type A] and [Cell Type B] only, split by condition.
> *(Note: "Condition" refers to the grouping variable in your dataset, such as treatment/control, disease/healthy, etc.)*

**Scenario 3D.5** — Multiple genes · Grouped by condition · Individual cell level
> Show a heatmap of [Gene A, Gene B, Gene C] grouped by condition.
> *(Note: This scenario displays expression at the individual cell level rather than aggregated averages, allowing you to visualize expression heterogeneity within each condition.)*

---

## 4. Find Markers

Marker gene identification finds genes that are specifically expressed in a cell type compared to all other cells. Use Find Markers to characterize cell populations, validate cell type annotations, or identify candidate genes for further investigation.

**Scenario 4.1** — Find marker genes for all cell types
> Find marker genes for all cell types.
> *(Note: Nvwa will return the top positive marker genes for each cell type, ranked by log2 fold change.)*

**Scenario 4.2** — Find marker genes for a specific cell type
> Find marker genes for [Cell Type].
> *(Note: Nvwa will return positive marker genes for the specified cell type compared to all other cells.)*

**Scenario 4.3** — Show top marker genes as dot plot
> Show a dot plot of the top [5] marker genes for each cell type.
> *(Note: Default is top 5 per cell type. You can specify a different number.)*

---

## 5. Differential Expression Analysis

Differential expression (DE) analysis identifies genes that are significantly upregulated or downregulated between two groups. Nvwa supports comparisons between cell types, between conditions, and between conditions within a specific cell type.

> ⚠️ **Statistical note:** DE results in Nvwa are based on single-cell level testing (Wilcoxon rank-sum). This treats each cell as an independent observation, which can produce overly significant p-values when multiple samples are present. For publication, consider pseudo-bulk analysis with DESeq2, especially for condition vs condition comparisons with multiple samples per condition.

**Scenario 5.1** — Cell type vs cell type · All cells
> Find differentially expressed genes between [Cell Type A] and [Cell Type B].

**Scenario 5.2** — Condition vs condition · All cell types
> Find differentially expressed genes between [Condition A] and [Condition B].
> *(Note: "Condition" refers to the grouping variable in your dataset, such as treatment/control, disease/healthy, etc.)*

**Scenario 5.3** — Condition vs condition · Single cell type
> Find differentially expressed genes between [Condition A] and [Condition B] in [Cell Type] only.

**Scenario 5.4** — Condition vs condition · Complex grouping
> Find differentially expressed genes between [Group 1] and [Group 2], where [Group 1] includes [Condition A, Condition B] and [Group 2] includes [Condition C, Condition D].
> *(Note: Nvwa will ask you to confirm the grouping structure before running the analysis.)*

### Volcano Plot

A volcano plot visualizes DE results by plotting log2 fold change (x-axis) against statistical significance (y-axis). After completing any DE analysis, Nvwa will automatically ask:
> "Would you like me to generate a volcano plot for these results?"

Simply reply **yes** to generate the plot, or **no** to skip it.

---

## 6. Cell Type Composition Analysis

Cell type composition analysis compares how cell populations are distributed across conditions. All scenarios in this module require at least two conditions in your dataset.

### 6A. All Cell Types

**Scenario 6A.1** — All cell types · Absolute cell number · Grouped by condition
> Show the cell number for each cell type across conditions.

**Scenario 6A.2** — All cell types · Proportion per cell type · Grouped by condition
> Show the proportion of each cell type across conditions.

**Scenario 6A.3** — All cell types · Proportion per condition
> Show the cell type composition for each condition.

### 6B. Subset Cell Types

**Scenario 6B.1** — Subset cell types · Absolute cell number · Grouped by condition
> Show the cell number for [Cell Type A], [Cell Type B], and [Cell Type C] across conditions.

**Scenario 6B.2** — Subset cell types · Proportion per cell type · Grouped by condition
> Show the proportion of [Cell Type A], [Cell Type B], and [Cell Type C] across conditions.

**Scenario 6B.3** — Subset cell types · Proportion per condition
> Show the cell type composition of [Cell Type A], [Cell Type B], and [Cell Type C] for each condition.

---
