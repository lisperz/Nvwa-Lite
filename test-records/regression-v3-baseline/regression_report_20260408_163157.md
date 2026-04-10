# nvwa-mvp Regression Report

**Generated:** 2026-04-08 23:34:44 UTC  
**Dataset:** `/Users/yxchen/Downloads/pbmc_test.h5ad`  
**Model:** `gpt-4o-mini`  
**Duration:** 167.1s  

---

## Summary

| Result | Count |
|--------|-------|
| ✅ Pass  | 20 |
| ❌ Fail  | 4 |
| ⚠️ Warn  | 1 |
| 💥 Error | 0 |
| ⏭ Skip  | 0 |
| **Total** | **25** |

---

## Per-case Results

| case_id | category | status | duration | artifacts | failure reason |
|---------|----------|--------|----------|-----------|----------------|
| DA-01 | data_analysis | ✅ PASS | 5.1s | none | — |
| DA-02 | data_analysis | ✅ PASS | 4.3s | 1 plot | — |
| DA-03 | data_analysis | ✅ PASS | 2.9s | 1 plot | — |
| DA-04 | data_analysis | ✅ PASS | 4.8s | 1 plot | — |
| DA-05 | data_analysis | ✅ PASS | 3.4s | none | — |
| DA-06 | data_analysis | ✅ PASS | 2.0s | 1 table | — |
| DA-07 | data_analysis | ✅ PASS | 3.0s | 2 plots | — |
| DA-08 | data_analysis | ✅ PASS | 3.4s | 1 plot | — |
| DA-09 | data_analysis | ❌ FAIL | 5.8s | none | must_contain: 'not found' not in response |
| DA-10 | data_analysis | ❌ FAIL | 3.3s | none | must_contain: 'not exist' not in response |
| DA-11 | data_analysis | ✅ PASS | 1.2s | none | — |
| DA-12 | data_analysis | ✅ PASS | 1.6s | none | — |
| DA-13 | data_analysis | ✅ PASS | 9.5s | none | — |
| DA-14 | data_analysis | ❌ FAIL | 15.1s | none | must_contain: 'not found' not in response |
| DA-15 | data_analysis | ✅ PASS | 6.6s | none | — |
| DA-16 | data_analysis | ✅ PASS | 15.2s | none | — |
| DA-17 | data_analysis | ✅ PASS | 12.1s | none | — |
| INT-01 | interpretation | ❌ FAIL | 8.4s | 1 table | must_contain: 'general' not in response |
| ART-01 | artifact_interaction | ✅ PASS | 3.7s | none | — |
| OOS-01 | out_of_scope | ✅ PASS | 6.9s | none | — |
| OOS-02 | out_of_scope | ✅ PASS | 2.8s | none | — |
| RS-01 | resource_scale | ⚠️ WARN | 8.5s | none | must_contain: 'dot plot' not in response; requires_plot: no valid PlotResult (none produced) |
| RS-02 | resource_scale | ✅ PASS | 9.5s | none | — |
| RS-03 | resource_scale | ✅ PASS | 14.8s | 1 plot | — |
| RS-04 | resource_scale | ✅ PASS | 12.8s | 1 plot | — |

---

## Top Failure Reasons

| Category | Count |
|----------|-------|
| `gene_not_found` | 2 |
| `assertion_failed` | 2 |

---

## Detailed Results

### ✅ DA-01

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 5.09s
- **Artifacts:** none

**Response (truncated):** The dataset contains a total of **2,638 cells** and **1,838 processed genes**. It has undergone normalization, PCA computation, UMAP computation, and clustering (using the "louvain" key). Additionally, differential expression results are available for analysis.   If you need further insights or visualizations, feel free to ask!

---

### ✅ DA-02

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 4.30s
- **Artifacts:** 1 plot
- **Notes:** Requires X_umap in obsm and a leiden/louvain column in obs. Fails immediately (pre-flight) if X_umap absent — use a preprocessed h5ad.


**Response (truncated):** The UMAP plot has been successfully generated, colored by the Leiden cluster (louvain). You can now visualize the distribution of your clusters in the dataset. If you need further analysis or modifications, just let me know!

---

### ✅ DA-03

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 2.87s
- **Artifacts:** 1 plot
- **Notes:** User explicitly named a plot type — only 1 artifact expected. A feature plot alongside it is an error (unsolicited output). Per C1 rule: generate only the specified plot type when user names one.


**Response (truncated):** The violin plot for CD3D expression across all clusters (louvain groups) has been successfully generated. This visualization provides insights into the distribution of CD3D expression levels within each cluster.   If you have any further questions or need additional analyses, feel free to ask!

---

### ✅ DA-04

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 4.84s
- **Artifacts:** 1 plot
- **Notes:** Dataset-generic prompt — no hardcoded gene names. Agent discovers marker genes from the dataset and plots them. Works on any preprocessed dataset with clustering results.


**Response (truncated):** The dot plot has been successfully generated, displaying the expression of the top 3 marker genes for each cluster across the Louvain groups.   This visualization allows you to see the expression levels and the fraction of cells expressing each gene within the different clusters. If you have any further analyses or visualizations in mind, feel free to let me know!

---

### ✅ DA-05

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.35s
- **Artifacts:** none
- **Notes:** Requires a leiden/seurat clustering column in obs. Agent may ask for clarification (marker vs pairwise) — that is acceptable.


**Response (truncated):** The differential expression analysis has been completed using the Wilcoxon method across all clusters. A total of 8 groups were analyzed, including CD4 T cells, CD14+ Monocytes, B cells, CD8 T cells, and NK cells. The top 20 marker genes for each group have been stored.  Would you like to view the top markers or download the results?

---

### ✅ DA-06

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 2.04s
- **Artifacts:** 1 table
- **Notes:** QC stats returned as TableResult (CSV + display). must_contain checks agent response text only — mean/median live in the table, not the text.


**Response (truncated):** The QC metrics summary for the dataset with 2,638 cells includes three detected metrics. You can view the detailed statistics in the generated table.

---

### ✅ DA-07

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.04s
- **Artifacts:** 2 plots
- **Notes:** Generic gene expression request with no plot type specified — both a feature plot (spatial) and violin plot (by cell type) are expected. Per C1 rule: when user does not specify a plot type, generate both. Paraphrase of a supported request; previously miscategorized as out_of_scope.


**Response (truncated):** I've generated two visualizations for the expression levels of CD3D:  1. **Feature Plot**: This UMAP plot shows the spatial distribution of CD3D expression across the different cell populations. 2. **Violin Plot**: This plot illustrates the distribution of CD3D expression levels across the various clusters identified by the Louvain algorithm.  These visualizations provide a comprehensive view of h

---

### ✅ DA-08

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.42s
- **Artifacts:** 1 plot
- **Notes:** Paraphrase of a UMAP-by-cluster request. Agent must map this to umap_plot and succeed. Requires preprocessed dataset (X_umap present). Previously miscategorized as out_of_scope.


**Response (truncated):** The UMAP plot has been generated, showing the cells colored by their cluster assignment (louvain). This visualization provides a clear view of how the cells are grouped in the low-dimensional embedding space.   If you have any further analyses or visualizations in mind, feel free to let me know!

---

### ❌ DA-09

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** FAIL
- **Duration:** 5.78s
- **Artifacts:** none
- **Notes:** Gene validation: FAKEGENE999 does not exist. Agent must report the gene is missing (and ideally suggest alternatives).


**Failures:**

- must_contain: 'not found' not in response

**Response (truncated):** The gene "FAKEGENE999" is not present in the dataset. Here are some of the available genes:  - 7SK-2 - A1BG - A2M - AAAS - ABCA1 - ABCB1 - ABCD1 - AATF  Please let me know if you would like to visualize the expression of one of these genes or any other gene from the dataset!

---

### ❌ DA-10

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** FAIL
- **Duration:** 3.27s
- **Artifacts:** none
- **Notes:** obs-key validation: nonexistent_batch_col not in adata.obs. Agent must list available keys and refuse to plot.


**Failures:**

- must_contain: 'not exist' not in response

**Response (truncated):** The available categorical metadata column in the dataset is "louvain," which contains 8 unique values representing different cell types. Would you like to color the UMAP by this "louvain" column instead?

---

### ✅ DA-11

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 1.23s
- **Artifacts:** none
- **Notes:** obs-field mistaken as gene: leiden is a cluster label, not a gene. Explicit phrasing forces the LLM to look up "leiden" as a gene name so the gene validator fires. Without the constraint the LLM auto-recovers by picking real marker genes.


**Response (truncated):** The gene "leiden" is not found in the dataset. Please check the gene name and provide an alternative gene for the violin plot.

---

### ✅ DA-12

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 1.64s
- **Artifacts:** none
- **Notes:** Color-key validation on scatter plot. Previous runs showed agent silently substituting louvain for the bad key. "Do not fall back" + must_not_contain together catch silent substitution.


**Response (truncated):** The requested color key "imaginary_condition_col" is not available in the dataset. Therefore, I cannot generate the scatter plot as specified.   If you have another color key in mind or would like to explore a different visualization, please let me know!

---

### ✅ DA-13

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 9.52s
- **Artifacts:** none
- **Notes:** Invalid cluster ID: cluster 999 does not exist. Agent must not silently succeed — we check it did not produce DE results rather than asserting specific error phrasing (which varies by run).


**Response (truncated):** The available clusters in the dataset are as follows:  - Cluster 0: B cells - Cluster 1: CD14+ Monocytes - Cluster 2: CD4 T cells - Cluster 3: CD8 T cells - Cluster 4: Dendritic cells - Cluster 5: FCGR3A+ Monocytes - Cluster 6: Megakaryocytes - Cluster 7: NK cells  It appears that "cluster 999" does not exist. Please specify a valid cluster from the list above for the differential expression analy

---

### ❌ DA-14

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** FAIL
- **Duration:** 15.05s
- **Artifacts:** none
- **Notes:** Ribosomal RNA percentage is not a standard QC metric in most h5ad files. Agent must report the column is missing and list available QC columns. Must not produce a plot or silently substitute another metric. Source: Yalu NC-A1.


**Failures:**

- must_contain: 'not found' not in response

**Response (truncated):** It appears that the dataset does not support wildcard searches for gene names like "RPL*" or "RPS*". Instead, I can check for specific ribosomal protein genes that are commonly present in datasets.  Here are some common ribosomal protein genes: - RPL4 - RPL5 - RPL7 - RPL8 - RPS2 - RPS3 - RPS6 - RPS7  Would you like me to plot the expression levels for these specific ribosomal protein genes?

---

### ✅ DA-15

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 6.63s
- **Artifacts:** none
- **Notes:** CD3 is a protein complex encoding multiple genes (CD3E, CD3D, CD3G). Agent should recognize the ambiguity and ask which gene to visualize. Must not silently pick one gene and plot. Hard to assert specific clarification wording — must_not_contain guards against silent success. Source: Yalu NC-C2.


**Response (truncated):** It appears that the gene "CD3" is not found in the dataset. Here are some similar gene names you might consider: CHD3, CDV3, CDK3, CD93, CD83.   Please let me know if you'd like to explore one of these alternatives or if you have another gene in mind!

---

### ✅ DA-16

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 15.18s
- **Artifacts:** none
- **Notes:** Ghost_Cell_Type does not exist in cell_type annotation. Agent must validate the cell type before running and report which input failed. Must not loop or produce results using the wrong cell type. Source: Yalu NC-D1.


**Response (truncated):** It appears that "Ghost_Cell_Type" is not listed among the available cell types in the dataset. The existing cell types are:  - CD4 T cells: 1,144 cells (43.4%) - CD14+ Monocytes: 480 cells (18.2%) - B cells: 342 cells (13.0%) - CD8 T cells: 316 cells (12.0%) - NK cells: 154 cells (5.8%) - FCGR3A+ Monocytes: 150 cells (5.7%) - Dendritic cells: 37 cells (1.4%) - Megakaryocytes: 15 cells (0.6%)  Plea

---

### ✅ DA-17

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 12.08s
- **Artifacts:** none
- **Notes:** FakeCondition does not exist in the dataset. All three inputs (cell type, condition A, condition B) must be validated independently. Agent must report which specific input failed. Must not run DE on partial valid inputs. Source: Yalu NC-D2.


**Response (truncated):** It seems that the dataset does not contain a column for conditions or samples, which is necessary to compare "Early cardiomyocyte" between "FakeCondition" and "Control." The only available metadata column is "louvain," which contains cluster information.  To proceed, please confirm the following: 1. **Is "Early cardiomyocyte" a specific cluster in the dataset?** If so, please provide the correspon

---

### ❌ INT-01

- **Category:** interpretation
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 8.40s
- **Artifacts:** 1 table
- **Notes:** Biological interpretation question — no tool can assess dataset health. Agent must answer using general biological knowledge AND include a disclaimer that this is not an analysis of the user's specific data. PASS example: "This is based on general biological knowledge, not an analysis of your dataset." FAIL: confident verdict with no disclaimer.


**Failures:**

- must_contain: 'general' not in response

**Response (truncated):** The dataset has been assessed for quality control metrics. Here’s a summary of the findings for the 2,638 cells:  - **Total Counts (UMI depth)**: Indicates the overall sequencing depth. - **Number of Genes Expressed**: Reflects the complexity of the cells. - **Mitochondrial Percentage**: Provides insight into cell health, as high mitochondrial content can indicate stress or dying cells.  If you wo

---

### ✅ ART-01

- **Category:** artifact_interaction
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.66s
- **Artifacts:** none
- **Notes:** Covers the general class of artifact interaction prompts (modify color, resize, export prior result, etc.). No memory system exists — "the last plot" is an unresolvable pointer in a fresh session. Correct behavior: agent asks the user which plot or what analysis they are referring to. PASS = agent clarifies without claiming false success. FAIL = agent claims to have resized something, or flatly refuses with no guidance. Promote to multi-case suite when session memory is implemented.


**Response (truncated):** It seems there is no previous plot in this session. Please let me know what type of plot you would like to create, and I can assist you with that!

---

### ✅ OOS-01

- **Category:** out_of_scope
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 6.89s
- **Artifacts:** none
- **Notes:** Dataset integration / batch correction is out of scope. Agent should politely decline and redirect to supported analyses.


**Response (truncated):** I currently don't have the capability to directly integrate datasets or perform batch correction using Harmony. However, I can guide you through the process if you have access to the necessary tools and libraries in your environment.  Here’s a general workflow you can follow to integrate your scRNA-seq datasets and perform batch correction using Harmony:  1. **Load the Datasets**: Load both scRNA-

---

### ✅ OOS-02

- **Category:** out_of_scope
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 2.84s
- **Artifacts:** none
- **Notes:** Trajectory analysis and RNA velocity are out of scope. Agent should decline without crashing. On a raw dataset it may redirect to preprocessing rather than explicitly saying "not supported" — so we only check it did NOT perform the analysis, not the specific refusal wording.


**Response (truncated):** It seems that the current dataset does not support PAGA trajectory analysis or RNA velocity visualization directly. The available tools and analyses in this environment are focused on differential expression, clustering, and standard visualizations like UMAP and violin plots.  If you have specific questions about the dataset or would like to explore other analyses, please let me know!

---

### ⚠️ RS-01

- **Category:** resource_scale
- **Expected status:** WARN
- **Actual status:** WARN
- **Duration:** 8.55s
- **Artifacts:** none
- **Notes:** 10-gene panel smoke-test. Partial match (some genes not found) is acceptable for WARN.


**Failures:**

- must_contain: 'dot plot' not in response
- requires_plot: no valid PlotResult (none produced)
- expected_artifact_count: expected 1, got 0 (0 plot(s), 0 table(s))

**Response (truncated):** The gene "IL2" was not found in the dataset. Here are some alternatives you might consider: ILF2, IL32, IL27, IL24, PPIL2.   Would you like to proceed with the dotplot using the other genes (CD3D, MS4A1, NKG7, LYZ, CD8A, CD4, FOXP3, TNF, IFNG) or would you like to replace "IL2" with one of the suggested alternatives?

---

### ✅ RS-02

- **Category:** resource_scale
- **Expected status:** WARN
- **Actual status:** PASS
- **Duration:** 9.52s
- **Artifacts:** none
- **Notes:** Cluster-level aggregation. Should complete within time limits even for datasets with >10k cells.


**Response (truncated):** Here is the distribution of cells in each Leiden cluster:  - **CD4 T cells**: 1,144 cells (43.4%) - **CD14+ Monocytes**: 480 cells (18.2%) - **B cells**: 342 cells (13.0%) - **CD8 T cells**: 316 cells (12.0%) - **NK cells**: 154 cells (5.8%) - **FCGR3A+ Monocytes**: 150 cells (5.7%) - **Dendritic cells**: 37 cells (1.4%) - **Megakaryocytes**: 15 cells (0.6%)  If you need further analysis or visual

---

### ✅ RS-03

- **Category:** resource_scale
- **Expected status:** WARN
- **Actual status:** PASS
- **Duration:** 14.84s
- **Artifacts:** 1 plot
- **Notes:** Multi-step compound query. Agent may use multiple tool calls in one turn. WARN because not all steps may complete for every dataset.


**Response (truncated):** The dataset contains a total of **2,638 cells**.  ### Top 5 Marker Genes for Each Cluster: - **CD4 T cells**: LDHB, CD3D, RPS27, RPS25, RPS12 - **CD14+ Monocytes**: S100A9, S100A8, LYZ, FCN1, FTL - **B cells**: CD79A, CD74, CD79B, HLA-DRA, MS4A1 - **CD8 T cells**: CCL5, NKG7, CST7, GZMA, CTSW - **NK cells**: NKG7, GZMB, PRF1, GNLY, CTSW - **FCGR3A+ Monocytes**: LST1, FCER1G, AIF1, FTH1, IFITM2 - *

---

### ✅ RS-04

- **Category:** resource_scale
- **Expected status:** WARN
- **Actual status:** PASS
- **Duration:** 12.76s
- **Artifacts:** 1 plot
- **Notes:** Tests the heatmap tool, which has no other coverage in the suite. If rank_genes_groups is pre-computed the agent plots directly; otherwise it runs DE first then heatmap. WARN because dataset-dependent.


**Response (truncated):** The heatmap displaying the top 5 marker genes for each cluster has been successfully generated. It includes the following genes across the identified clusters:  - **CD4 T cells**: LDHB, CD3D, RPS27, RPS25, RPS12 - **CD14+ Monocytes**: S100A9, S100A8, LYZ, FCN1, FTL - **B cells**: CD79A, CD74, CD79B, HLA-DRA, MS4A1 - **CD8 T cells**: CCL5, NKG7, CST7, GZMA, CTSW - **NK cells**: NKG7, GZMB, PRF1, GN

---
