# nvwa-mvp Regression Report

**Generated:** 2026-04-22 22:23:35 UTC  
**Dataset:** `/Users/yxchen/Documents/GitHub/nvwa-mvp/local/data/lung_combined_slim_v3.h5ad`  
**Model:** `gpt-4o-mini`  
**Duration:** 247.0s  

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
| DAI-01 | data_analysis | ✅ PASS | 11.9s | none | — |
| DAI-02 | data_analysis | ✅ PASS | 3.8s | 1 plot | — |
| DAI-03 | data_analysis | ✅ PASS | 2.7s | 1 plot | — |
| DAI-04 | data_analysis | ✅ PASS | 7.2s | 4 plots | — |
| DAI-05 | data_analysis | ✅ PASS | 12.6s | 1 plot | — |
| DAI-07 | data_analysis | ❌ FAIL | 3.8s | none | requires_plot: no valid PlotResult (none produced) |
| DAI-08 | data_analysis | ✅ PASS | 9.7s | 1 plot | — |
| DAI-09 | data_analysis | ✅ PASS | 6.5s | 1 plot | — |
| DAI-10 | data_analysis | ✅ PASS | 17.2s | 1 plot | — |
| DAI-11 | data_analysis | ✅ PASS | 10.0s | 1 plot | — |
| DAI-12 | data_analysis | ❌ FAIL | 35.1s | 1 plot | must_contain: 'heatmap' not in response |
| DAI-13 | data_analysis | ❌ FAIL | 5.7s | none | must_contain: 'heatmap' not in response; requires_plot: no valid PlotResult (none produced) |
| DAI-14 | data_analysis | ❌ FAIL | 23.5s | none | requires_plot: no valid PlotResult (none produced) |
| DAI-15 | data_analysis | ❌ FAIL | 6.5s | 1 plot | must_not_contain: forbidden string 'panels' found; must_not_contain: forbidden string 'each cell type' found |
| DAI-16 | out_of_scope | ✅ PASS | 2.7s | none | — |
| DAI-17 | data_analysis | ❌ FAIL | 1.7s | none | requires_plot: no valid PlotResult (none produced) |
| DAI-18 | data_analysis | ❌ FAIL | 67.2s | none | must_contain: 'volcano' not in response; requires_plot: no valid PlotResult (none produced) |
| DAI-19 | data_analysis | ✅ PASS | 18.5s | 1 plot | — |

---

## Top Failure Reasons

| Category | Count |
|----------|-------|
| `missing_artifact` | 5 |
| `assertion_failed` | 1 |
| `other` | 1 |

---

## Detailed Results

### ✅ DAI-01

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 11.85s
- **Artifacts:** none
- **Notes:** Paraphrase of DA-01 (dataset overview). Lowercase, no punctuation.


**Response (truncated):** This dataset contains single-cell RNA sequencing data with the following details:  - **Total Cells**: 15,719 - **Total Genes**: 32,285 - **Observation Keys**:    - `celltype`   - `ident`   - `nCount_RNA` (total RNA counts)   - `nFeature_RNA` (number of genes detected)   - `orig.ident` (original sample identifier)   - `percent.mt` (percentage of mitochondrial genes)   - `seurat_clusters` (cluster a

---

### ✅ DAI-02

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.81s
- **Artifacts:** 1 plot
- **Notes:** Paraphrase of DA-02. No color-by specified — agent may default to clusters or cell_type; both acceptable.


**Response (truncated):** The UMAP plot has been generated, displaying the clusters identified in your dataset. This visualization provides insights into the spatial distribution of the different cell populations based on the clustering results.   Would you like to explore any specific clusters or perform further analyses?

---

### ✅ DAI-03

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 2.71s
- **Artifacts:** 1 plot
- **Notes:** UMAP colored by the cell_type obs column.


**Response (truncated):** The UMAP plot has been generated, displaying the cell type information. You can now visualize the distribution of different cell types across the UMAP space. If you have any further requests or need additional analysis, feel free to ask!

---

### ✅ DAI-04

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 7.17s
- **Artifacts:** 4 plots
- **Notes:** Two-gene expression, no plot type specified. Per C1 rule, feature + violin for each gene is the target; artifact_count not asserted so the test is tolerant of the agent batching differently.


**Response (truncated):** I've generated the gene expression plots for Fabp4 and Fabp5:  1. **Fabp4**:    - **Feature Plot**: This shows the spatial expression of Fabp4 on the UMAP.    - **Violin Plot**: This illustrates the distribution of Fabp4 expression across different seurat_clusters.  2. **Fabp5**:    - **Feature Plot**: This shows the spatial expression of Fabp5 on the UMAP.    - **Violin Plot**: This illustrates t

---

### ✅ DAI-05

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 12.55s
- **Artifacts:** 1 plot
- **Notes:** Tests Layer 2a parameter resolution: when user says "cell types", agent must group by the cell_type obs column, not default to seurat_clusters. Rewritten from a context-dependent follow-up in Dai's session so the prompt is self-contained.


**Response (truncated):** The violin plot of Fabp4 expression across different cell types has been generated successfully. This visualization allows you to see the distribution of Fabp4 expression levels within each cell type.   If you have any further questions or need additional analyses, feel free to ask!

---

### ❌ DAI-07

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 3.83s
- **Artifacts:** none
- **Notes:** Pairwise condition comparison. WT and CKO are condition labels in obs (not cluster names). Agent should subset to WT/CKO and plot — violin or DE both acceptable. Known-failing in Dai's session on condition resolution (see DAI-18/19).


**Failures:**

- requires_plot: no valid PlotResult (none produced)

**Response (truncated):** To compare the gene expression of Fabp4 between WT and CKO, I need to confirm the following:  1. The specific cell type (if any) you want to focus on for this comparison. 2. Ensure that both WT and CKO conditions exist in the dataset.  Please provide the specific cell type if you have one in mind, or let me know if you want to compare across all cell types.

---

### ✅ DAI-08

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 9.74s
- **Artifacts:** 1 plot
- **Notes:** Tests that WT and CKO are rendered on a single combined violin, not as separate panels. Rewritten from a context-dependent follow-up in Dai's session (where DAI-07 had produced split violins) so the "single combined plot" constraint is explicit in the prompt.


**Response (truncated):** The combined violin plot for Fabp4 expression in WT shows that there are 4,572 cells with a mean expression of approximately 0.193 and a median expression of 0.000.   However, since there are no CKO cells present in the plot, we cannot directly compare the expression levels between WT and CKO.   If you would like to explore further or analyze a different aspect, please let me know!

---

### ✅ DAI-09

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 6.54s
- **Artifacts:** 1 plot
- **Notes:** "Different samples" is ambiguous — could mean orig.ident, condition, or another obs column. Agent should resolve or ask. Note: user typed FABP4 (uppercase) while the gene in the dataset is Fabp4 (mouse).


**Response (truncated):** The violin plot of Fabp4 expression across different samples (orig.ident) has been generated successfully. This visualization allows you to see the distribution of Fabp4 expression levels across the various samples.   If you have any further questions or need additional analysis, feel free to ask!

---

### ✅ DAI-10

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 17.18s
- **Artifacts:** 1 plot
- **Notes:** Three-condition comparison. All three labels must be resolved in the condition obs column. Paraphrase variant of DAI-07.


**Response (truncated):** The violin plot of Fabp4 expression across the conditions (WT, CKO, and TKO) has been successfully generated.   ### Insight: The plot illustrates the distribution of Fabp4 expression levels in each condition, allowing for a visual comparison of expression patterns across WT, CKO, and TKO groups.   Would you like to explore further or analyze any specific aspects of this data?

---

### ✅ DAI-11

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 10.02s
- **Artifacts:** 1 plot
- **Notes:** Generic "by condition" — agent should auto-detect condition-like columns and split by all levels. Third phrasing of the same underlying request (DAI-07 / DAI-10 / DAI-11) — tests paraphrase robustness.


**Response (truncated):** The violin plot of Fabp4 expression across different conditions (orig.ident) has been generated successfully.   This visualization allows you to see the distribution of Fabp4 expression levels in each condition. If you have any further questions or need additional analyses, feel free to ask!

---

### ❌ DAI-12

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 35.06s
- **Artifacts:** 1 plot
- **Notes:** Markers of the EC (endothelial) cell type on a heatmap. Agent must resolve EC as a cell_type label and compute markers.


**Failures:**

- must_contain: 'heatmap' not in response

---

### ❌ DAI-13

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 5.71s
- **Artifacts:** none
- **Notes:** Top-N markers variant of DAI-12. "how about" prefix is conversational; agent must still execute.


**Failures:**

- must_contain: 'heatmap' not in response
- requires_plot: no valid PlotResult (none produced)

**Response (truncated):** It seems that "EC" is not recognized as a valid cluster in the dataset. The available clusters are indexed from 0 to 25.   Could you please specify the numeric cluster ID or provide the exact name of the cluster you are referring to? This will help me retrieve the top markers for that specific cluster.

---

### ❌ DAI-14

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 23.53s
- **Artifacts:** none
- **Notes:** Per-cell-type top-3 markers. Dataset-generic request with no hard-coded genes.


**Failures:**

- requires_plot: no valid PlotResult (none produced)

**Response (truncated):** The top 3 unique marker genes for each of the 26 cell types are as follows:  1. **Cluster 0**: Chil3, Ctsd, Ccl6 2. **Cluster 1**: Ramp2, Ly6a, Ier3 3. **Cluster 2**: Rps24, Rps16, Rpl12 4. **Cluster 3**: Trbc2, Cd3g, Rpl19 5. **Cluster 4**: Gsn, Mgp, Serping1 6. **Cluster 5**: Gpx3, Inmt, Mgp 7. **Cluster 6**: Sftpd, Sftpb, Sftpa1 8. **Cluster 7**: Cxcl12, Sparcl1, Ly6c1 9. **Cluster 8**: Cd79a, 

---

### ❌ DAI-15

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 6.46s
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
- **Duration:** 2.72s
- **Artifacts:** none
- **Notes:** Jitter plot is not a supported tool. Agent must decline without crashing, or redirect to violin/strip.


**Response (truncated):** It seems you are looking for a jitter plot to visualize gene expression across all cell types. However, I currently do not have a specific tool for generating jitter plots.   Instead, I can create a dot plot or a violin plot to show the distribution of gene expression across different cell types. Would you like to proceed with one of these options? If so, please specify the gene or genes you are i

---

### ❌ DAI-17

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 1.65s
- **Artifacts:** none
- **Notes:** Volcano plot for EC — requires a DE result first. Agent should run DE (EC vs rest, or prompt for comparator) then plot.


**Failures:**

- requires_plot: no valid PlotResult (none produced)

**Response (truncated):** To generate a volcano plot, I need to first confirm the context of "EC." Are you referring to a specific cluster or cell type? Please clarify so I can proceed accurately.

---

### ❌ DAI-18

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 67.23s
- **Artifacts:** none
- **Notes:** Subset to EC, pairwise DE between WT and CKO, volcano plot. Failed in Dai's session: condition labels WT/CKO not resolved against the condition obs column. Exact prompt that surfaced the bug.


**Failures:**

- must_contain: 'volcano' not in response
- requires_plot: no valid PlotResult (none produced)

**Response (truncated):** It appears that the group identifiers "CKO" and "WT" are not directly recognized in the dataset for the differential expression analysis. This suggests that the conditions may not be directly linked to the identifiers used in the analysis.  To resolve this, I will filter the dataset to include only the endothelial cells (EC) and then perform the differential expression analysis between the WT and 

---

### ✅ DAI-19

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 18.46s
- **Artifacts:** 1 plot
- **Notes:** Paraphrase of DAI-18 with "condition" added as a disambiguator. Same bug class — useful to see whether the extra word flips the result.


**Response (truncated):** The differential expression analysis between the WT and CKO conditions for endothelial cells (EC) has been successfully completed. A total of 1,929 genes were found to be significantly differentially expressed (adjusted p-value < 0.05).  The volcano plot visualizing the results has been generated. This plot highlights the genes with significant changes in expression between the two conditions, wit

---
