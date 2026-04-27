# nvwa-mvp Regression Report

**Generated:** 2026-04-25 01:18:54 UTC  
**Dataset:** `/Users/yxchen/Documents/GitHub/nvwa-mvp/local/data/lung_combined_slim_v3.h5ad`  
**Model:** `gpt-4o-mini`  
**Duration:** 215.2s  

---

## Summary

| Result | Count |
|--------|-------|
| ✅ Pass  | 11 |
| ❌ Fail  | 7 |
| ⚠️ Warn  | 0 |
| 💥 Error | 0 |
| ⏭ Skip  | 0 |
| **Total** | **18** |

---

## Per-case Results

| case_id | category | status | duration | artifacts | failure reason |
|---------|----------|--------|----------|-----------|----------------|
| DAI-01 | data_analysis | ✅ PASS | 6.5s | none | — |
| DAI-02 | data_analysis | ✅ PASS | 4.8s | 1 plot | — |
| DAI-03 | data_analysis | ✅ PASS | 3.9s | 1 plot | — |
| DAI-04 | data_analysis | ✅ PASS | 6.0s | 4 plots | — |
| DAI-05 | data_analysis | ✅ PASS | 5.8s | 1 plot | — |
| DAI-07 | data_analysis | ❌ FAIL | 25.7s | none | must_contain: 'Fabp4' not in response; requires_plot: no valid PlotResult (none produced) |
| DAI-08 | data_analysis | ❌ FAIL | 3.1s | none | requires_plot: no valid PlotResult (none produced) |
| DAI-09 | data_analysis | ✅ PASS | 8.2s | 1 plot | — |
| DAI-10 | data_analysis | ✅ PASS | 16.1s | 1 plot | — |
| DAI-11 | data_analysis | ✅ PASS | 5.0s | 1 plot | — |
| DAI-12 | data_analysis | ✅ PASS | 41.0s | 1 plot | — |
| DAI-13 | data_analysis | ❌ FAIL | 5.1s | none | must_contain: 'heatmap' not in response; requires_plot: no valid PlotResult (none produced) |
| DAI-14 | data_analysis | ❌ FAIL | 15.7s | none | must_contain: 'heatmap' not in response; requires_plot: no valid PlotResult (none produced) |
| DAI-15 | data_analysis | ❌ FAIL | 7.5s | 1 plot | must_not_contain: forbidden string 'panels' found; must_not_contain: forbidden string 'each cell type' found |
| DAI-16 | out_of_scope | ✅ PASS | 4.7s | none | — |
| DAI-17 | data_analysis | ❌ FAIL | 5.0s | none | must_contain: 'volcano' not in response; requires_plot: no valid PlotResult (none produced) |
| DAI-18 | data_analysis | ✅ PASS | 26.2s | 1 plot | — |
| DAI-19 | data_analysis | ❌ FAIL | 24.4s | none | must_contain: 'volcano' not in response; requires_plot: no valid PlotResult (none produced) |

---

## Top Failure Reasons

| Category | Count |
|----------|-------|
| `missing_artifact` | 6 |
| `other` | 1 |

---

## Detailed Results

### ✅ DAI-01

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 6.51s
- **Artifacts:** none
- **Notes:** Paraphrase of DA-01 (dataset overview). Lowercase, no punctuation.


**Response (truncated):** This dataset contains single-cell RNA sequencing data with the following details:  - **Total Cells**: 15,719 - **Total Genes**: 32,285 - **Observation Keys**:    - `celltype`: Annotation of cell types   - `ident`: Identifier for cells   - `nCount_RNA`: Total RNA counts per cell   - `nFeature_RNA`: Number of genes detected per cell   - `orig.ident`: Original identifier for the sample   - `percent.m

---

### ✅ DAI-02

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 4.75s
- **Artifacts:** 1 plot
- **Notes:** Paraphrase of DA-02. No color-by specified — agent may default to clusters or cell_type; both acceptable.


**Response (truncated):** The UMAP plot has been generated, displaying the clusters identified in your dataset. This visualization provides insights into the spatial distribution of the different cell populations based on the clustering results.   Would you like to explore any specific clusters or perform further analyses?

---

### ✅ DAI-03

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.91s
- **Artifacts:** 1 plot
- **Notes:** UMAP colored by the cell_type obs column.


**Response (truncated):** The UMAP plot has been generated, displaying the cell type information. You can now visualize the distribution of different cell types across the UMAP space. If you have any further requests or need additional analysis, feel free to ask!

---

### ✅ DAI-04

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 6.04s
- **Artifacts:** 4 plots
- **Notes:** Two-gene expression, no plot type specified. Per C1 rule, feature + violin for each gene is the target; artifact_count not asserted so the test is tolerant of the agent batching differently.


**Response (truncated):** I've generated the gene expression plots for Fabp4 and Fabp5:  1. **Fabp4**:    - **Feature Plot**: Shows the spatial expression of Fabp4 on the UMAP.    - **Violin Plot**: Displays the distribution of Fabp4 expression across different seurat_clusters.  2. **Fabp5**:    - **Feature Plot**: Illustrates the spatial expression of Fabp5 on the UMAP.    - **Violin Plot**: Shows the distribution of Fabp

---

### ✅ DAI-05

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 5.79s
- **Artifacts:** 1 plot
- **Notes:** Tests Layer 2a parameter resolution: when user says "cell types", agent must group by the cell_type obs column, not default to seurat_clusters. Rewritten from a context-dependent follow-up in Dai's session so the prompt is self-contained.


**Response (truncated):** The violin plot of Fabp4 expression across different cell types has been generated successfully. This visualization allows you to see the distribution of Fabp4 expression levels within each cell type.   If you have any further questions or need additional analyses, feel free to ask!

---

### ❌ DAI-07

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 25.66s
- **Artifacts:** none
- **Notes:** Pairwise condition comparison. WT and CKO are condition labels in obs (not cluster names). Agent should subset to WT/CKO and plot — violin or DE both acceptable. Known-failing in Dai's session on condition resolution (see DAI-18/19).


**Failures:**

- must_contain: 'Fabp4' not in response
- requires_plot: no valid PlotResult (none produced)

---

### ❌ DAI-08

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 3.06s
- **Artifacts:** none
- **Notes:** Tests that WT and CKO are rendered on a single combined violin, not as separate panels. Rewritten from a context-dependent follow-up in Dai's session (where DAI-07 had produced split violins) so the "single combined plot" constraint is explicit in the prompt.


**Failures:**

- requires_plot: no valid PlotResult (none produced)

**Response (truncated):** To create a combined violin plot comparing Fabp4 expression between WT and CKO, I need to confirm the specific cell type or condition you want to focus on.   Could you please specify the cell type or confirm if you want to analyze all cells?

---

### ✅ DAI-09

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 8.19s
- **Artifacts:** 1 plot
- **Notes:** "Different samples" is ambiguous — could mean orig.ident, condition, or another obs column. Agent should resolve or ask. Note: user typed FABP4 (uppercase) while the gene in the dataset is Fabp4 (mouse).


**Response (truncated):** The violin plot of Fabp4 expression across different samples (orig.ident) has been generated successfully. This visualization allows you to see the distribution of Fabp4 expression levels across the various samples.   If you have any further questions or need additional analysis, feel free to ask!

---

### ✅ DAI-10

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 16.08s
- **Artifacts:** 1 plot
- **Notes:** Three-condition comparison. All three labels must be resolved in the condition obs column. Paraphrase variant of DAI-07.


**Response (truncated):** The violin plot of Fabp4 expression across the conditions (WT, CKO, and TKO) has been successfully generated.   ### Insight: The plot illustrates the distribution of Fabp4 expression levels in each condition, allowing for a visual comparison of expression patterns.   Would you like to explore further or analyze any specific aspects of this data?

---

### ✅ DAI-11

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 5.01s
- **Artifacts:** 1 plot
- **Notes:** Generic "by condition" — agent should auto-detect condition-like columns and split by all levels. Third phrasing of the same underlying request (DAI-07 / DAI-10 / DAI-11) — tests paraphrase robustness.


**Response (truncated):** The violin plot of Fabp4 expression across different conditions (orig.ident) has been generated successfully.   This visualization allows you to see the distribution of Fabp4 expression levels in each condition. If you have any further questions or need additional analysis, feel free to ask!

---

### ✅ DAI-12

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 40.98s
- **Artifacts:** 1 plot
- **Notes:** Markers of the EC (endothelial) cell type on a heatmap. Agent must resolve EC as a cell_type label and compute markers.


**Response (truncated):** The heatmap of unique markers for endothelial cells (EC) has been successfully generated. This heatmap displays the expression levels of various marker genes across the different clusters in your dataset.  ### Key Insights: - The heatmap provides a visual representation of how these markers are expressed in each cluster, allowing you to identify which clusters are enriched for endothelial cell mar

---

### ❌ DAI-13

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 5.15s
- **Artifacts:** none
- **Notes:** Top-N markers variant of DAI-12. "how about" prefix is conversational; agent must still execute.


**Failures:**

- must_contain: 'heatmap' not in response
- requires_plot: no valid PlotResult (none produced)

**Response (truncated):** It seems that "EC" is not recognized as a valid cluster in the dataset. The available clusters are indexed from 0 to 25.   Could you please specify the numeric cluster ID or provide the exact name of the cluster you are referring to? This will help me retrieve the top markers for the correct cluster.

---

### ❌ DAI-14

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 15.72s
- **Artifacts:** none
- **Notes:** Per-cell-type top-3 markers. Dataset-generic request with no hard-coded genes.


**Failures:**

- must_contain: 'heatmap' not in response
- requires_plot: no valid PlotResult (none produced)

**Response (truncated):** Here are the unique top 3 marker genes for each of the 26 cell types:  - **Cluster 0**: Chil3, Ctsd, Ccl6 - **Cluster 1**: Ramp2, Ly6a, Ier3 - **Cluster 2**: Rps24, Rps16, Rpl12 - **Cluster 3**: Trbc2, Cd3g, Rpl19 - **Cluster 4**: Gsn, Mgp, Serping1 - **Cluster 5**: Gpx3, Inmt, Mgp - **Cluster 6**: Sftpd, Sftpb, Sftpa1 - **Cluster 7**: Cxcl12, Sparcl1, Ly6c1 - **Cluster 8**: Cd79a, Igkc, Ighm - **

---

### ❌ DAI-15

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 7.51s
- **Artifacts:** 1 plot
- **Notes:** Subset UMAP to a single cell type (EC). Agent should subset, not highlight — "only show" implies subset. Known-failing in Dai's session (produced split panels instead of subset); negatives target the split-panel narration that previously slipped through.


**Failures:**

- must_not_contain: forbidden string 'panels' found
- must_not_contain: forbidden string 'each cell type' found

**Response (truncated):** The UMAP plot has been successfully generated, displaying endothelial cells (EC) along with other cell types. Each cell type is represented in separate panels for clear comparison.   If you have any further questions or need additional analyses, feel free to ask!

---

### ✅ DAI-16

- **Category:** out_of_scope
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 4.71s
- **Artifacts:** none
- **Notes:** Jitter plot is not a supported tool. Agent must decline without crashing, or redirect to violin/strip.


**Response (truncated):** This is not supported in this version. I've recorded your request for the development team to consider in the next iteration.

---

### ❌ DAI-17

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 5.01s
- **Artifacts:** none
- **Notes:** Volcano plot for EC — requires a DE result first. Agent should run DE (EC vs rest, or prompt for comparator) then plot.


**Failures:**

- must_contain: 'volcano' not in response
- requires_plot: no valid PlotResult (none produced)

**Response (truncated):** This is not supported in this version. I've recorded your request for the development team to consider in the next iteration.

---

### ✅ DAI-18

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 26.20s
- **Artifacts:** 1 plot
- **Notes:** Subset to EC, pairwise DE between WT and CKO, volcano plot. Failed in Dai's session: condition labels WT/CKO not resolved against the condition obs column. Exact prompt that surfaced the bug.


**Response (truncated):** The differential expression analysis between WT and CKO for endothelial cells (EC) has been successfully completed, revealing 1,929 significant genes (adjusted p-value < 0.05).  I have generated a volcano plot to visualize the results of this comparison. The plot highlights genes with significant expression differences, where positive log2 fold changes indicate higher expression in WT and negative

---

### ❌ DAI-19

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 24.40s
- **Artifacts:** none
- **Notes:** Paraphrase of DAI-18 with "condition" added as a disambiguator. Same bug class — useful to see whether the extra word flips the result.


**Failures:**

- must_contain: 'volcano' not in response
- requires_plot: no valid PlotResult (none produced)

**Response (truncated):** The cluster mapping does not provide specific names for the clusters, so I will need to identify which cluster(s) correspond to endothelial cells (EC) based on the previous metadata inspection.  Since we know that endothelial cells are present in the dataset, I will perform the differential expression analysis for all cells, but I will focus on the results for the EC cell type afterward.  Let's pr

---
