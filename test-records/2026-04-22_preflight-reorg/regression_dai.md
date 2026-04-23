# nvwa-mvp Regression Report

**Generated:** 2026-04-22 17:03:55 UTC  
**Dataset:** `/Users/yxchen/Documents/GitHub/nvwa-mvp/local/data/lung_combined_slim_v3.h5ad`  
**Model:** `gpt-4o-mini`  
**Duration:** 228.8s  

---

## Summary

| Result | Count |
|--------|-------|
| ✅ Pass  | 10 |
| ❌ Fail  | 8 |
| ⚠️ Warn  | 0 |
| 💥 Error | 0 |
| ⏭ Skip  | 0 |
| **Total** | **18** |

---

## Per-case Results

| case_id | category | status | duration | artifacts | failure reason |
|---------|----------|--------|----------|-----------|----------------|
| DAI-01 | data_analysis | ✅ PASS | 29.4s | none | — |
| DAI-02 | data_analysis | ✅ PASS | 3.1s | 1 plot | — |
| DAI-03 | data_analysis | ✅ PASS | 2.4s | 1 plot | — |
| DAI-04 | data_analysis | ✅ PASS | 5.2s | 4 plots | — |
| DAI-05 | data_analysis | ✅ PASS | 2.8s | 1 plot | — |
| DAI-07 | data_analysis | ❌ FAIL | 13.4s | none | must_contain: 'Fabp4' not in response; requires_plot: no valid PlotResult (none produced) |
| DAI-08 | data_analysis | ❌ FAIL | 2.0s | none | requires_plot: no valid PlotResult (none produced) |
| DAI-09 | data_analysis | ✅ PASS | 5.9s | 1 plot | — |
| DAI-10 | data_analysis | ✅ PASS | 18.2s | 1 plot | — |
| DAI-11 | data_analysis | ✅ PASS | 14.2s | 1 plot | — |
| DAI-12 | data_analysis | ✅ PASS | 65.4s | 1 plot | — |
| DAI-13 | data_analysis | ❌ FAIL | 4.1s | none | must_contain: 'heatmap' not in response; requires_plot: no valid PlotResult (none produced) |
| DAI-14 | data_analysis | ❌ FAIL | 1.5s | none | must_contain: 'heatmap' not in response; requires_plot: no valid PlotResult (none produced) |
| DAI-15 | data_analysis | ❌ FAIL | 6.3s | 1 plot | must_not_contain: forbidden string 'panels' found; must_not_contain: forbidden string 'each cell type' found |
| DAI-16 | out_of_scope | ✅ PASS | 2.0s | none | — |
| DAI-17 | data_analysis | ❌ FAIL | 5.8s | none | requires_plot: no valid PlotResult (none produced) |
| DAI-18 | data_analysis | ❌ FAIL | 25.0s | none | must_contain: 'volcano' not in response; requires_plot: no valid PlotResult (none produced) |
| DAI-19 | data_analysis | ❌ FAIL | 21.5s | none | must_contain: 'EC' not in response; must_contain: 'volcano' not in response |

---

## Top Failure Reasons

| Category | Count |
|----------|-------|
| `missing_artifact` | 7 |
| `other` | 1 |

---

## Detailed Results

### ✅ DAI-01

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 29.42s
- **Artifacts:** none
- **Notes:** Paraphrase of DA-01 (dataset overview). Lowercase, no punctuation.


**Response (truncated):** This dataset contains single-cell RNA sequencing data with the following details:  - **Total Cells**: 15,719 - **Total Genes**: 32,285 - **Observation Keys**:    - `celltype`   - `ident`   - `nCount_RNA`   - `nFeature_RNA`   - `orig.ident`   - `percent.mt`   - `seurat_clusters`    - **Sample Genes** (first 50 of 32,285):   - 0610005C13Rik   - 0610006L08Rik   - 0610009B22Rik   - 0610009E02Rik   - 0

---

### ✅ DAI-02

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.06s
- **Artifacts:** 1 plot
- **Notes:** Paraphrase of DA-02. No color-by specified — agent may default to clusters or cell_type; both acceptable.


**Response (truncated):** The UMAP plot has been generated, displaying the clusters identified in your dataset. This visualization provides insights into the spatial distribution of the different cell populations based on the clustering results.   Would you like to explore any specific clusters or perform further analyses?

---

### ✅ DAI-03

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 2.45s
- **Artifacts:** 1 plot
- **Notes:** UMAP colored by the cell_type obs column.


**Response (truncated):** The UMAP plot has been generated, displaying the cell type information. You can now visualize how different cell types are distributed across the UMAP space. If you have any further questions or need additional analyses, feel free to ask!

---

### ✅ DAI-04

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 5.22s
- **Artifacts:** 4 plots
- **Notes:** Two-gene expression, no plot type specified. Per C1 rule, feature + violin for each gene is the target; artifact_count not asserted so the test is tolerant of the agent batching differently.


**Response (truncated):** I've generated the gene expression plots for Fabp4 and Fabp5:  1. **Fabp4**:    - **Feature Plot**: Shows the spatial expression of Fabp4 on the UMAP.    - **Violin Plot**: Displays the distribution of Fabp4 expression across different seurat_clusters.  2. **Fabp5**:    - **Feature Plot**: Illustrates the spatial expression of Fabp5 on the UMAP.    - **Violin Plot**: Shows the distribution of Fabp

---

### ✅ DAI-05

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 2.78s
- **Artifacts:** 1 plot
- **Notes:** Tests Layer 2a parameter resolution: when user says "cell types", agent must group by the cell_type obs column, not default to seurat_clusters. Rewritten from a context-dependent follow-up in Dai's session so the prompt is self-contained.


**Response (truncated):** The violin plot of Fabp4 expression across different cell types has been generated successfully. This visualization allows you to see the distribution of Fabp4 expression levels within each cell type.   If you have any further questions or need additional analyses, feel free to ask!

---

### ❌ DAI-07

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 13.39s
- **Artifacts:** none
- **Notes:** Pairwise condition comparison. WT and CKO are condition labels in obs (not cluster names). Agent should subset to WT/CKO and plot — violin or DE both acceptable. Known-failing in Dai's session on condition resolution (see DAI-18/19).


**Failures:**

- must_contain: 'Fabp4' not in response
- requires_plot: no valid PlotResult (none produced)

**Response (truncated):** The pairwise differential expression analysis between CKO and WT has been completed. Here are the key results:  - **Total genes analyzed**: 32,285 - **Significant genes (adjusted p-value < 0.05)**: 1,929  ### Interpretation: - A positive log2 fold change (log2FC) indicates higher expression in CKO. - A negative log2FC indicates higher expression in WT.  If you would like to see the detailed result

---

### ❌ DAI-08

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 2.01s
- **Artifacts:** none
- **Notes:** Tests that WT and CKO are rendered on a single combined violin, not as separate panels. Rewritten from a context-dependent follow-up in Dai's session (where DAI-07 had produced split violins) so the "single combined plot" constraint is explicit in the prompt.


**Failures:**

- requires_plot: no valid PlotResult (none produced)

**Response (truncated):** To create a combined violin plot comparing Fabp4 expression between WT and CKO, I need to first confirm that "WT" and "CKO" are valid entries in the condition column. Let me check the available conditions in the dataset.   I'll proceed with that now.

---

### ✅ DAI-09

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 5.87s
- **Artifacts:** 1 plot
- **Notes:** "Different samples" is ambiguous — could mean orig.ident, condition, or another obs column. Agent should resolve or ask. Note: user typed FABP4 (uppercase) while the gene in the dataset is Fabp4 (mouse).


**Response (truncated):** The violin plot of Fabp4 expression across different samples (orig.ident) has been successfully generated. This visualization allows you to see the distribution of Fabp4 expression levels across the various samples.   If you have any further questions or need additional analysis, feel free to ask!

---

### ✅ DAI-10

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 18.22s
- **Artifacts:** 1 plot
- **Notes:** Three-condition comparison. All three labels must be resolved in the condition obs column. Paraphrase variant of DAI-07.


**Response (truncated):** The violin plot of Fabp4 expression across the conditions WT, CKO, and TKO has been generated successfully.   ### Insight: The plot illustrates the distribution of Fabp4 expression levels in each condition, allowing for a visual comparison of expression patterns. You can observe how the expression varies among the three groups.  Would you like to explore further or analyze something else?

---

### ✅ DAI-11

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 14.24s
- **Artifacts:** 1 plot
- **Notes:** Generic "by condition" — agent should auto-detect condition-like columns and split by all levels. Third phrasing of the same underlying request (DAI-07 / DAI-10 / DAI-11) — tests paraphrase robustness.


**Response (truncated):** The violin plot of Fabp4 expression across different conditions (orig.ident) has been generated successfully. This visualization allows you to see the distribution of Fabp4 expression levels in each condition.   If you have any further questions or need additional analysis, feel free to ask!

---

### ✅ DAI-12

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 65.43s
- **Artifacts:** 1 plot
- **Notes:** Markers of the EC (endothelial) cell type on a heatmap. Agent must resolve EC as a cell_type label and compute markers.


**Response (truncated):** The heatmap of unique markers for all clusters has been successfully generated. This visualization displays the expression levels of various genes across the different clusters in your dataset.  If you have any specific questions about the heatmap or if you would like to explore further analyses, please let me know!

---

### ❌ DAI-13

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 4.10s
- **Artifacts:** none
- **Notes:** Top-N markers variant of DAI-12. "how about" prefix is conversational; agent must still execute.


**Failures:**

- must_contain: 'heatmap' not in response
- requires_plot: no valid PlotResult (none produced)

**Response (truncated):** It seems that "EC" is not a recognized cluster in the dataset. The available clusters are indexed from 0 to 25.   Could you please specify the numeric cluster ID or provide the name of the cluster you are referring to? This will help me find the top markers for that specific cluster.

---

### ❌ DAI-14

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 1.53s
- **Artifacts:** none
- **Notes:** Per-cell-type top-3 markers. Dataset-generic request with no hard-coded genes.


**Failures:**

- must_contain: 'heatmap' not in response
- requires_plot: no valid PlotResult (none produced)

**Response (truncated):** To identify the unique top 3 markers for all cell types, I will first run the differential expression analysis to find the distinguishing genes for each cell type. Then, I will retrieve the top markers.  Let's proceed with that.

---

### ❌ DAI-15

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 6.34s
- **Artifacts:** 1 plot
- **Notes:** Subset UMAP to a single cell type (EC). Agent should subset, not highlight — "only show" implies subset. Known-failing in Dai's session (produced split panels instead of subset); negatives target the split-panel narration that previously slipped through.


**Failures:**

- must_not_contain: forbidden string 'panels' found
- must_not_contain: forbidden string 'each cell type' found

**Response (truncated):** The UMAP plot has been generated, displaying endothelial cells (EC) along with other cell types. Each cell type is represented in separate panels for clear comparison.   If you have any further questions or need additional analyses, feel free to ask!

---

### ✅ DAI-16

- **Category:** out_of_scope
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 2.04s
- **Artifacts:** none
- **Notes:** Jitter plot is not a supported tool. Agent must decline without crashing, or redirect to violin/strip.


**Response (truncated):** It seems you are looking for a jitter plot to visualize gene expression across all cell types. However, the available tools do not include a specific jitter plot function.   Instead, I can create a dot plot or a violin plot to show the distribution of gene expression across different cell types. Would you like to proceed with one of these options? If so, please specify the gene or genes you are in

---

### ❌ DAI-17

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 5.84s
- **Artifacts:** none
- **Notes:** Volcano plot for EC — requires a DE result first. Agent should run DE (EC vs rest, or prompt for comparator) then plot.


**Failures:**

- requires_plot: no valid PlotResult (none produced)

**Response (truncated):** It seems that the group "EC" is not found in the available clusters. The available clusters are: 1, 7, 8, 0, 5, 17, 2, 3, 23, 12, 14, 19, 13, 11, 16, 9, 4, 10, 6, 22, 18, 15, 24, 21, 20, and 25.  Could you please specify which cluster you would like to analyze for the volcano plot?

---

### ❌ DAI-18

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 25.01s
- **Artifacts:** none
- **Notes:** Subset to EC, pairwise DE between WT and CKO, volcano plot. Failed in Dai's session: condition labels WT/CKO not resolved against the condition obs column. Exact prompt that surfaced the bug.


**Failures:**

- must_contain: 'volcano' not in response
- requires_plot: no valid PlotResult (none produced)

**Response (truncated):** It appears that the group identifier "CKO" could not be resolved either. This suggests that the conditions might not be directly available in the expected format for the differential expression analysis.  Since the conditions are listed under the `orig.ident` column, I will filter the dataset to include only the endothelial cells (EC) and then perform the differential expression analysis between t

---

### ❌ DAI-19

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 21.53s
- **Artifacts:** none
- **Notes:** Paraphrase of DAI-18 with "condition" added as a disambiguator. Same bug class — useful to see whether the extra word flips the result.


**Failures:**

- must_contain: 'EC' not in response
- must_contain: 'volcano' not in response
- requires_plot: no valid PlotResult (none produced)

---
