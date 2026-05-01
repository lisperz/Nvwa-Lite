# nvwa-mvp Regression Report

**Generated:** 2026-05-01 17:40:22 UTC  
**Dataset:** `/Users/yxchen/GitHub/nvwa-mvp/local/data/heart_snRNAseq_slim_new.h5ad`  
**Model:** `gpt-4o-mini`  
**Duration:** 207.2s  

---

## Summary

| Result | Count |
|--------|-------|
| ✅ Pass  | 36 |
| ❌ Fail  | 0 |
| ⚠️ Warn  | 1 |
| 💥 Error | 0 |
| ⏭ Skip  | 0 |
| **Total** | **37** |

---

## Per-case Results

| case_id | category | status | routed_via | duration | artifacts | failure reason |
|---------|----------|--------|------------|----------|-----------|----------------|
| L1-1.1 | data_analysis | ✅ PASS | spec | 3.3s | 1 plot | — |
| L1-1.2 | data_analysis | ✅ PASS | spec | 3.5s | 1 plot | — |
| L1-1.3 | data_analysis | ✅ PASS | spec | 2.9s | 1 plot | — |
| L1-1.4 | data_analysis | ✅ PASS | spec | 4.1s | 1 plot | — |
| L1-2A.1 | data_analysis | ✅ PASS | spec | 3.1s | 1 plot | — |
| L1-2A.2 | data_analysis | ✅ PASS | spec | 4.9s | 1 plot | — |
| L1-2B.1 | data_analysis | ✅ PASS | spec | 2.2s | 1 table | — |
| L1-3A.1 | data_analysis | ✅ PASS | spec | 2.4s | 1 plot | — |
| L1-3A.2 | data_analysis | ✅ PASS | spec | 3.3s | 1 plot | — |
| L1-3A.3 | data_analysis | ✅ PASS | spec | 4.1s | 1 plot | — |
| L1-3A.4 | data_analysis | ⚠️ WARN | spec | 5.0s | none | requires_plot: no valid PlotResult (none produced) |
| L1-3B.1 | data_analysis | ✅ PASS | spec | 4.2s | 1 plot | — |
| L1-3B.2 | data_analysis | ✅ PASS | spec | 3.6s | 1 plot | — |
| L1-3B.3 | data_analysis | ✅ PASS | spec | 3.0s | 1 plot | — |
| L1-3B.4 | data_analysis | ✅ PASS | spec | 3.0s | 1 plot | — |
| L1-3B.5 | data_analysis | ✅ PASS | spec | 8.0s | 1 plot | — |
| L1-3B.6 | data_analysis | ✅ PASS | spec | 3.7s | 1 plot | — |
| L1-3C.1 | data_analysis | ✅ PASS | spec | 10.4s | 1 plot | — |
| L1-3C.2 | data_analysis | ✅ PASS | spec | 3.6s | 1 plot | — |
| L1-3C.3 | data_analysis | ✅ PASS | spec | 3.4s | 1 plot | — |
| L1-3C.4 | data_analysis | ✅ PASS | spec | 4.5s | 1 plot | — |
| L1-3D.1 | data_analysis | ✅ PASS | spec | 16.1s | 1 plot | — |
| L1-3D.2 | data_analysis | ✅ PASS | spec | 15.4s | 1 plot | — |
| L1-3D.3 | data_analysis | ✅ PASS | spec | 6.8s | 1 plot | — |
| L1-3D.4 | data_analysis | ✅ PASS | spec | 6.0s | 1 plot | — |
| L1-3D.5 | data_analysis | ✅ PASS | spec | 16.3s | 1 plot | — |
| L1-4.1 | data_analysis | ✅ PASS | spec | 11.5s | 1 table | — |
| L1-4.2 | data_analysis | ✅ PASS | spec | 6.7s | 1 table | — |
| L1-5.1 | data_analysis | ✅ PASS | spec | 5.8s | 1 table | — |
| L1-5.2 | data_analysis | ✅ PASS | spec | 8.6s | 1 table | — |
| L1-5.3 | data_analysis | ✅ PASS | spec | 5.3s | 1 table | — |
| L1-6A.1 | data_analysis | ✅ PASS | spec | 2.7s | 1 plot | — |
| L1-6A.2 | data_analysis | ✅ PASS | spec | 3.1s | 1 plot | — |
| L1-6A.3 | data_analysis | ✅ PASS | spec | 3.4s | 1 plot | — |
| L1-6B.1 | data_analysis | ✅ PASS | spec | 3.9s | 1 plot | — |
| L1-6B.2 | data_analysis | ✅ PASS | spec | 4.0s | 1 plot | — |
| L1-6B.3 | data_analysis | ✅ PASS | spec | 4.3s | 1 plot | — |

---

## Top Failure Reasons

_No failures recorded._

---

## Detailed Results

### ✅ L1-1.1

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.31s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §1.1 — UMAP colored by cell type, no condition split.

**Response (truncated):** Here's the UMAP scatter plot showing 20,377 cells colored by their cell type. This visualization helps you explore the distribution and relationships among different cell types in your dataset (umap_plot).

---

### ✅ L1-1.2

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.47s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §1.2 — UMAP split by condition.

**Response (truncated):** Here's the UMAP scatter plot showing 20,377 cells colored by cell type, split by the original identity of the samples (you wrote 'condition', which was corrected to 'orig.ident'). This visualization helps you explore the distribution of cell types across different conditions.

---

### ✅ L1-1.3

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 2.91s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §1.3 — UMAP subset to specific cell types.

**Response (truncated):** Here's the UMAP scatter plot showing 3,794 cells, colored by cell type and subsetted to include only the categories 'CM.1' and 'Fib.1' (using `umap_plot`). I used 'celltype' for the subset key (you wrote 'cell type', case-corrected).

---

### ✅ L1-1.4

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 4.09s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §1.4 — UMAP subset + split by condition.

**Response (truncated):** Here's the UMAP scatter plot showing 3,794 cells colored by cell type, with the data split by the original identity of the samples. The plot is subsetted to include only the cell types 'CM.1' and 'Fib.1' (I used 'celltype' instead of 'cell type' for consistency). This visualization was created using the `umap_plot` tool.

---

### ✅ L1-2A.1

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.06s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §2A.1 — qc_violin_plot, single multi-panel figure.

**Response (truncated):** Here's the violin plot showing the quality control metrics (nCount_RNA, nFeature_RNA, and percent.mt) across different cell types for your dataset, which includes a total of 20,377 cells. This visualization was created using the `qc_violin_plot` tool.

---

### ✅ L1-2A.2

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 4.89s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §2A.2 — qc_violin_plot with split_by; metric × condition grid.

**Response (truncated):** Here's the QC violin plot showing three metrics (nCount_RNA, nFeature_RNA, and percent.mt) across cell types, split by the original identity of the conditions. This visualization includes data from 20,377 cells, allowing for a comparison of quality control metrics across the two conditions (qc_violin_plot). Note that I used 'orig.ident' for the split (you wrote 'condition', which I corrected).

---

### ✅ L1-2B.1

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 2.15s
- **Artifacts:** 1 table
- **Tool called:** True
- **Notes:** Yalu §2B.1 — qc_summary_table with automated flagging.

**Response (truncated):** Here's the QC summary for your dataset, which includes 20,377 cells across three metrics. All metrics are within typical ranges, indicating no issues with elevated mitochondrial fraction or low gene detection (qc_summary_table).

---

### ✅ L1-3A.1

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 2.42s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §3A.1 — single-gene UMAP overlay.

**Response (truncated):** Here's the feature plot of Tnnt2 expression on UMAP for 20,377 cells, showing how this gene is distributed across the dataset (feature_plot).

---

### ✅ L1-3A.2

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.35s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §3A.2 — single-gene + split by condition.

**Response (truncated):** Here's the feature plot of Tnnt2 expression on UMAP for 20,377 cells, split by the original identity of the samples (you wrote 'condition', which was adjusted to 'orig.ident'). This visualization highlights the distribution of Tnnt2 across different conditions.

---

### ✅ L1-3A.3

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 4.06s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §3A.3 — single-gene + cell-type subset.

**Response (truncated):** Here's the feature plot of Tnnt2 expression on UMAP, specifically for the CM.1 cell type, showing data from 1,412 cells. I used the `feature_plot` tool, and note that I canonicalized "cell type" to "celltype" for consistency.

---

### ⚠️ L1-3A.4

- **Category:** data_analysis
- **Expected status:** WARN
- **Actual status:** WARN
- **Duration:** 5.04s
- **Artifacts:** none
- **Tool called:** True
- **Notes:** Yalu §3A.4 — single-gene + subset + split. **Yalu spec gap (WARN):** §3A.3 / §3A.4 prompts drop "on the UMAP" keyword present in §3A.1 / §3A.2, making them linguistically indistinguishable from §3B.4 / §3B.5 (both share "in [Cell Type] only [, split by condition]" shape). LLM extractor cannot reliably disambiguate without the keyword. See local/from_yalu/scenario_gaps_2.md for proposed Yalu Layer 1 wording fix. Re-enable to PASS when Yalu spec adds "on the UMAP" to §3A.3 / §3A.4.


**Failures:**

- requires_plot: no valid PlotResult (none produced)

**Response (truncated):** The user asked about ['Tnnt2', 'CM.1'], but the tool output doesn't reference any of these entities. The tool may have ignored the parameters or produced output too terse to confirm they were acted on.

---

### ✅ L1-3B.1

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 4.19s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §3B.1 — violin, single gene, x-axis = cell type.

**Response (truncated):** Here's the violin plot of Tnnt2 expression across different cell types for 20,377 cells. I used the `violin_plot` tool, and note that I corrected "cell type" to "celltype" for consistency.

---

### ✅ L1-3B.2

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.58s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §3B.2 — violin split by condition.

**Response (truncated):** Here's the violin plot of Tnnt2 expression across different cell types, split by the original condition identifiers. I used the `violin_plot` tool for this analysis.

---

### ✅ L1-3B.3

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 2.97s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §3B.3 — violin, x-axis = condition.

**Response (truncated):** Here's the violin plot of Tnnt2 expression across different conditions (orig.ident) for 20,377 cells. I used the `violin_plot` tool, and please note that the condition you specified was canonicalized to match the dataset's identifiers.

---

### ✅ L1-3B.4

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 2.95s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §3B.4 — violin, subset to one cell type, x-axis = condition.

**Response (truncated):** Here's the violin plot of Tnnt2 expression across different conditions (orig.ident) for the subset of cells identified as CM.1. I used the `violin_plot` tool to visualize the distribution for 1,412 cells.

---

### ✅ L1-3B.5

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 8.05s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §3B.5 — violin, subset multi-celltype + split.

**Response (truncated):** Here's the violin plot of Tnnt2 expression across the specified cell types (CM.1, Fib.1, EC.1), split by condition (orig.ident). I used 'cell type' (you wrote 'celltype', case-corrected) for the grouping. This plot represents the distribution of Tnnt2 expression for 5,269 cells.

---

### ✅ L1-3B.6

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.75s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §3B.6 — violin, subset to one condition, x-axis = cell type.

**Response (truncated):** Here's the violin plot of Tnnt2 expression across different cell types, specifically for the Control condition. I used the `violin_plot` tool to visualize the distribution for 11,104 cells. Note that I adjusted the grouping to 'celltype' (you wrote 'cell type', case-corrected) and used 'orig.ident' for the condition (you wrote 'condition', canonicalized).

---

### ✅ L1-3C.1

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 10.37s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §3C.1 — multi-gene dot plot.

**Response (truncated):** Here's the dot plot showing the expression of Tnnt2, Pecam1, and Postn across different cell types. I used the `dot_plot` tool, and the cell types are represented on the y-axis with dot size indicating the percentage of cells expressing each gene and dot color reflecting the mean expression level. Note that I corrected "cell type" to "celltype" for consistency.

---

### ✅ L1-3C.2

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.56s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §3C.2 — dot plot split by condition.

**Response (truncated):** Here's the dot plot showing the expression of Tnnt2, Pecam1, and Postn across different cell types, split by the original identity of the samples. I used the `dot_plot` tool for this visualization.

---

### ✅ L1-3C.3

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.38s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §3C.3 — dot plot subset cell types.

**Response (truncated):** Here's the dot plot showing the expression of Tnnt2, Pecam1, and Postn across the cell types CM.1 and Fib.1. I used the `dot_plot` tool, and the data is subsetted to these specific cell types.

---

### ✅ L1-3C.4

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 4.51s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §3C.4 — dot plot subset + split.

**Response (truncated):** Here's the dot plot showing the expression of Tnnt2, Pecam1, and Postn across the specified cell types, split by the original identity of the conditions. The analysis focused on the cell types CM.1 and Fib.1 (dot_plot).

---

### ✅ L1-3D.1

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 16.08s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §3D.1 — multi-gene heatmap.

**Response (truncated):** Here's the heatmap showing the expression of Tnnt2, Pecam1, and Postn across different cell types. I used the `heatmap_plot` tool to visualize this data.

---

### ✅ L1-3D.2

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 15.41s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §3D.2 — heatmap split by condition.

**Response (truncated):** Here's the heatmap showing the expression of Tnnt2, Pecam1, and Postn across different cell types, split by the original identity of the conditions. I used the `heatmap_plot` tool for this visualization.

---

### ✅ L1-3D.3

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 6.78s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §3D.3 — heatmap subset.

**Response (truncated):** Here's the heatmap showing the expression of Tnnt2, Pecam1, and Postn across the selected cell types, CM.1 and Fib.1. I used the `heatmap_plot` tool to visualize this data, and the reference condition is anchored to blue.

---

### ✅ L1-3D.4

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 6.04s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §3D.4 — heatmap subset + split.

**Response (truncated):** Here's the heatmap showing the expression of Tnnt2, Pecam1, and Postn across cell types, split by the original identity of the samples. The analysis is focused on the cell types CM.1 and Fib.1 (heatmap_plot).

---

### ✅ L1-3D.5

- **Category:** data_analysis
- **Expected status:** WARN
- **Actual status:** PASS
- **Duration:** 16.31s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §3D.5 — spec-gap WARN. Yalu requests "individual cell level" rendering but sc.pl.heatmap aggregates by groupby by default. Per scenario_gaps.md, either Yalu spec needs cell-level rendering mode explicit, or §3D.5 collapses to §3D.1 with different groupby. Until resolved, the response is a regular aggregated heatmap with condition as the group axis.


**Response (truncated):** Here's the heatmap showing the expression of Tnnt2, Pecam1, and Postn across different conditions (orig.ident). I used the `heatmap_plot` tool, which includes a dendrogram for cluster similarity. Note that the grouping was adjusted from 'condition' to 'orig.ident' for this analysis.

---

### ✅ L1-4.1

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 11.46s
- **Artifacts:** 1 table
- **Tool called:** True
- **Notes:** Yalu §4.1 — find_all_markers; multi-group survey.

**Response (truncated):** I identified marker genes across 26 cell types, revealing a total of 46,617 significant markers (padj<0.05, logFC>0.0) using the `find_all_markers` tool. Each cell type has a list of genes that are specifically expressed compared to all other cells, which you can explore in the output provided.

---

### ✅ L1-4.2

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 6.65s
- **Artifacts:** 1 table
- **Tool called:** True
- **Notes:** Yalu §4.2 — find_markers; single-celltype one-vs-rest.

**Response (truncated):** I identified marker genes specifically expressed in the CM.1 cell type compared to all other cells, resulting in 4,747 significant markers (padj<0.05, logFC>0.0). The top 10 markers include Gm56685 (log2FC=26.37) and Gm57283 (log2FC=13.91), among others, using the `find_markers` tool.

---

### ✅ L1-5.1

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 5.83s
- **Artifacts:** 1 table
- **Tool called:** True
- **Notes:** Yalu §5.1 — run_de pairwise cell type vs cell type.

**Response (truncated):** I found the differential expression between CM.1 and Fib.1 based on cell type. A total of 7,067 genes were significant (padj<0.05), with 3,957 genes upregulated in CM.1 and 3,110 genes upregulated in Fib.1. The top genes up in CM.1 include Rbm20 and Fgf13, while Mrc2 and Gpr39 are among the top upregulated genes in Fib.1. Please note that these results are based on the Wilcoxon rank-sum test, whic

---

### ✅ L1-5.2

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 8.59s
- **Artifacts:** 1 table
- **Tool called:** True
- **Notes:** Yalu §5.2 — run_de pairwise condition vs condition.

**Response (truncated):** I performed pairwise differential expression analysis comparing the Control and Mutant groups based on the 'orig.ident' column. A total of 3,829 genes were found to be significant (adjusted p-value < 0.05), with 640 genes upregulated in Control and 3,189 genes upregulated in Mutant. The top upregulated genes in Control include Gm42743 (log2FC=7.67) and Il1rapl1 (log2FC=5.94), while in Mutant, H2-E

---

### ✅ L1-5.3

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 5.31s
- **Artifacts:** 1 table
- **Tool called:** True
- **Notes:** Yalu §5.3 — run_de pairwise condition vs condition + cell-type subset.

**Response (truncated):** I found the differentially expressed genes between the Control and Mutant groups within the CM.1 cell type. A total of 1,423 genes were significant (padj<0.05), with 335 genes upregulated in Control and 1,088 in Mutant. The top upregulated genes in Control include Afap1l1 (log2FC=28.61) and Gm56906 (log2FC=26.35), while in Mutant, Gm35533 (log2FC=-26.91) and Hsd17b7 (log2FC=-26.15) were the most p

---

### ✅ L1-6A.1

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 2.73s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §6A.1 — composition_barplot mode='count'.

**Response (truncated):** Here's the bar plot showing the cell type composition across the two conditions, displaying absolute cell numbers for each of the 26 cell types. I used the `composition_barplot` tool to create this visualization.

---

### ✅ L1-6A.2

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.10s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §6A.2 — composition_barplot mode='proportion_by_celltype'.

**Response (truncated):** Here's the bar plot showing the proportion of cell types across the two conditions in your dataset. I used the `composition_barplot` tool to visualize the composition of 26 cell types. Note that I canonicalized 'cell type' to 'celltype' and 'condition' to 'orig.ident' for the analysis.

---

### ✅ L1-6A.3

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.44s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §6A.3 — composition_barplot mode='proportion_by_condition'.

**Response (truncated):** Here's the bar plot showing the cell type composition across the two conditions, presented as proportions within each condition. I used the `composition_barplot` tool, which visualizes the distribution of 26 cell types across the specified conditions. Note that I canonicalized 'cell type' to 'celltype' and 'condition' to 'orig.ident' for clarity.

---

### ✅ L1-6B.1

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.87s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §6B.1 — composition mode='count' + celltypes subset.

**Response (truncated):** Here's the cell type composition bar plot showing the absolute counts of three selected cell types (CM.1, Fib.1, EC.1) across two conditions. I used the `composition_barplot` tool, with the cell type column canonicalized from 'cell type' to 'celltype' and the condition column from 'condition' to 'orig.ident'.

---

### ✅ L1-6B.2

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.99s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §6B.2 — composition mode='proportion_by_celltype' + subset.

**Response (truncated):** Here's the bar plot showing the cell type composition across conditions, displaying the proportions of the selected cell types (CM.1, Fib.1, EC.1) within each condition. I used the `composition_barplot` tool, with the cell type column canonicalized from 'cell type' to 'celltype' and the condition column from 'condition' to 'orig.ident'.

---

### ✅ L1-6B.3

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 4.25s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Yalu §6B.3 — composition mode='proportion_by_condition' + subset.

**Response (truncated):** Here's the bar plot showing the cell type composition across conditions, displaying the proportions of the selected cell types (CM.1, Fib.1, EC.1) within each condition. I used the `composition_barplot` tool, with the cell type column canonicalized from 'cell type' to 'celltype' and the condition column from 'condition' to 'orig.ident'.

---
