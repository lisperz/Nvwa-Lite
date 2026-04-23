# nvwa-mvp Regression Report

**Generated:** 2026-04-22 17:16:07 UTC  
**Dataset:** `/Users/yxchen/Documents/GitHub/nvwa-mvp/local/data/pbmc_test.h5ad`  
**Model:** `gpt-4o-mini`  
**Duration:** 228.1s  

---

## Summary

| Result | Count |
|--------|-------|
| ✅ Pass  | 21 |
| ❌ Fail  | 4 |
| ⚠️ Warn  | 1 |
| 💥 Error | 0 |
| ⏭ Skip  | 0 |
| **Total** | **26** |

---

## Per-case Results

| case_id | category | status | duration | artifacts | failure reason |
|---------|----------|--------|----------|-----------|----------------|
| DA-01 | data_analysis | ✅ PASS | 13.6s | none | — |
| DA-02 | data_analysis | ✅ PASS | 15.7s | 1 plot | — |
| DA-03 | data_analysis | ✅ PASS | 14.3s | 1 plot | — |
| DA-04 | data_analysis | ✅ PASS | 18.9s | 1 plot | — |
| DA-05 | data_analysis | ✅ PASS | 3.9s | none | — |
| DA-06 | data_analysis | ✅ PASS | 2.3s | 1 table | — |
| DA-07 | data_analysis | ❌ FAIL | 3.3s | 1 plot | expected_artifact_count: expected 2, got 1 (1 plot(s), 0 table(s)) |
| DA-08 | data_analysis | ✅ PASS | 3.4s | 1 plot | — |
| DA-09 | data_analysis | ❌ FAIL | 2.6s | none | must_contain: 'not present' not in response |
| DA-10 | data_analysis | ✅ PASS | 3.8s | none | — |
| DA-11 | data_analysis | ✅ PASS | 3.3s | none | — |
| DA-12 | data_analysis | ✅ PASS | 1.6s | none | — |
| DA-13 | data_analysis | ✅ PASS | 11.1s | none | — |
| DA-14 | data_analysis | ✅ PASS | 11.8s | none | — |
| DA-15 | data_analysis | ✅ PASS | 4.6s | none | — |
| DA-16 | data_analysis | ✅ PASS | 13.2s | none | — |
| DA-17 | data_analysis | ❌ FAIL | 47.5s | none | must_contain: 'FakeCondition' not in response |
| INT-01 | interpretation | ❌ FAIL | 5.1s | 1 table | expected_artifact_count: expected 0, got 1 (0 plot(s), 1 table(s)) |
| INT-02 | interpretation | ✅ PASS | 5.1s | none | — |
| ART-01 | artifact_interaction | ✅ PASS | 1.4s | none | — |
| OOS-01 | out_of_scope | ✅ PASS | 5.4s | none | — |
| OOS-02 | out_of_scope | ✅ PASS | 2.0s | none | — |
| RS-01 | resource_scale | ⚠️ WARN | 4.0s | none | must_contain: 'dot plot' not in response; requires_plot: no valid PlotResult (none produced) |
| RS-02 | resource_scale | ✅ PASS | 5.5s | none | — |
| RS-03 | resource_scale | ✅ PASS | 14.4s | 1 plot | — |
| RS-04 | resource_scale | ✅ PASS | 10.1s | 1 plot | — |

---

## Top Failure Reasons

| Category | Count |
|----------|-------|
| `other` | 2 |
| `gene_not_found` | 1 |
| `assertion_failed` | 1 |

---

## Detailed Results

### ✅ DA-01

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 13.55s
- **Artifacts:** none

**Response (truncated):** The dataset contains a total of **2,638 cells** and **1,838 processed genes**. It has undergone normalization, PCA computation, UMAP computation, and clustering (using the "louvain" key). Additionally, differential expression results are available for analysis.   If you need further insights or visualizations, feel free to ask!

---

### ✅ DA-02

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 15.71s
- **Artifacts:** 1 plot
- **Notes:** Requires X_umap in obsm and a leiden/louvain column in obs. Fails immediately (pre-flight) if X_umap absent — use a preprocessed h5ad.


**Response (truncated):** The UMAP plot has been generated, colored by the Leiden cluster (louvain). You can now visualize the distribution of your clusters in the dataset. If you have any further analyses or visualizations in mind, feel free to let me know!

---

### ✅ DA-03

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 14.28s
- **Artifacts:** 1 plot
- **Notes:** User explicitly named a plot type — only 1 artifact expected. A feature plot alongside it is an error (unsolicited output). Per C1 rule: generate only the specified plot type when user names one.


**Response (truncated):** The violin plot of CD3D expression across all clusters (louvain groups) has been successfully generated. This visualization provides insights into the distribution of CD3D expression levels within each cluster.   Would you like to explore further or analyze another gene?

---

### ✅ DA-04

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 18.86s
- **Artifacts:** 1 plot
- **Notes:** Dataset-generic prompt — no hardcoded gene names. Agent discovers marker genes from the dataset and plots them. Works on any preprocessed dataset with clustering results.


**Response (truncated):** The dot plot has been successfully generated, displaying the expression of the top 3 marker genes for each cluster across the Louvain groups. This visualization provides insights into the distinguishing gene expression patterns associated with each cluster.   If you have any further analyses or visualizations in mind, feel free to let me know!

---

### ✅ DA-05

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.87s
- **Artifacts:** none
- **Notes:** Requires a leiden/seurat clustering column in obs. Agent may ask for clarification (marker vs pairwise) — that is acceptable.


**Response (truncated):** The differential expression analysis has been completed using the Wilcoxon method. It identified the top 20 marker genes for each of the 8 clusters, including CD4 T cells, CD14+ Monocytes, B cells, CD8 T cells, and NK cells.   Would you like to see the top markers for a specific cluster or download the results?

---

### ✅ DA-06

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 2.29s
- **Artifacts:** 1 table
- **Notes:** QC stats returned as TableResult (CSV + display). must_contain checks agent response text only — mean/median live in the table, not the text.


**Response (truncated):** The QC metrics summary for the dataset with 2,638 cells includes three detected metrics. You can view the detailed statistics in the generated table.

---

### ❌ DA-07

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 3.34s
- **Artifacts:** 1 plot
- **Notes:** Generic gene expression request with no plot type specified — both a feature plot (spatial) and violin plot (by cell type) are expected. Per C1 rule: when user does not specify a plot type, generate both. Paraphrase of a supported request; previously miscategorized as out_of_scope.


**Failures:**

- expected_artifact_count: expected 2, got 1 (1 plot(s), 0 table(s))

**Response (truncated):** The violin plot of CD3D expression across different cell populations (louvain clusters) has been generated successfully. This visualization allows you to see the distribution of CD3D expression levels within each cluster.   If you have any further questions or need additional analyses, feel free to ask!

---

### ✅ DA-08

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.36s
- **Artifacts:** 1 plot
- **Notes:** Paraphrase of a UMAP-by-cluster request. Agent must map this to umap_plot and succeed. Requires preprocessed dataset (X_umap present). Previously miscategorized as out_of_scope.


**Response (truncated):** The UMAP plot has been generated, showing the cells colored by their cluster assignment (louvain). This visualization provides insight into how the cells are grouped in the low-dimensional embedding space.   If you have any further analyses or visualizations in mind, feel free to let me know!

---

### ❌ DA-09

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** FAIL
- **Duration:** 2.59s
- **Artifacts:** none
- **Notes:** Gene validation: FAKEGENE999 does not exist. Agent must report the gene is missing (and ideally suggest alternatives). Agent uses "not present" not "not found".


**Failures:**

- must_contain: 'not present' not in response

**Response (truncated):** The gene "FAKEGENE999" was not found in the dataset, and there are no close matches available. Please check the gene name for any typos or provide a different gene to visualize.

---

### ✅ DA-10

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 3.76s
- **Artifacts:** none
- **Notes:** obs-key validation: nonexistent_batch_col not in adata.obs. Agent lists available columns and offers alternatives — does not use a fixed phrase like "not exist". must_not_contain guards against silent success instead.


**Response (truncated):** The available categorical metadata column in the dataset is "louvain," which contains 8 unique values representing different cell types.   Would you like to color the UMAP by the "louvain" column instead?

---

### ✅ DA-11

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 3.33s
- **Artifacts:** none
- **Notes:** obs-field mistaken as gene: leiden is a cluster label, not a gene. Explicit phrasing forces the LLM to look up "leiden" as a gene name so the gene validator fires. Without the constraint the LLM auto-recovers by picking real marker genes.


**Response (truncated):** The gene name "leiden" is not found in the dataset. Please check the spelling or provide a different gene name for the violin plot.

---

### ✅ DA-12

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 1.58s
- **Artifacts:** none
- **Notes:** Color-key validation on scatter plot. Previous runs showed agent silently substituting louvain for the bad key. "Do not fall back" + must_not_contain together catch silent substitution.


**Response (truncated):** The requested color key "imaginary_condition_col" does not exist in the dataset. Therefore, I cannot generate the scatter plot as specified.   If you have another color key in mind or would like to explore a different visualization, please let me know!

---

### ✅ DA-13

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 11.06s
- **Artifacts:** none
- **Notes:** Invalid cluster ID: cluster 999 does not exist. Agent must not silently succeed — we check it did not produce DE results rather than asserting specific error phrasing (which varies by run).


**Response (truncated):** The available clusters in the dataset are as follows:  - Cluster 0: B cells - Cluster 1: CD14+ Monocytes - Cluster 2: CD4 T cells - Cluster 3: CD8 T cells - Cluster 4: Dendritic cells - Cluster 5: FCGR3A+ Monocytes - Cluster 6: Megakaryocytes - Cluster 7: NK cells  It appears that "cluster 999" does not exist. Please specify a valid cluster from the list above, and I can proceed with the different

---

### ✅ DA-14

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 11.78s
- **Artifacts:** none
- **Notes:** Ribosomal RNA percentage (pct_counts_rb) is not a standard QC metric in most h5ad files. Prompt is explicit about QC metric to avoid the agent interpreting it as a gene expression query (RPL/RPS genes). Agent must report the column is missing. Must not produce a plot. Source: Yalu NC-A1.


**Response (truncated):** The dataset contains a categorical metadata column for "louvain," which represents the clustering of cells into different types. However, it does not have a specific key for ribosomal RNA percentage.  If you would like, I can calculate the mitochondrial percentage instead, which is often used as a quality control metric. Would you like me to proceed with that?

---

### ✅ DA-15

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 4.60s
- **Artifacts:** none
- **Notes:** CD3 is a protein complex encoding multiple genes (CD3E, CD3D, CD3G). Agent should recognize the ambiguity and ask which gene to visualize. Must not silently pick one gene and plot. Hard to assert specific clarification wording — must_not_contain guards against silent success. Source: Yalu NC-C2.


**Response (truncated):** It seems that the gene "CD3" was not found in the dataset. Here are some similar gene names that you might be referring to: CHD3, CDV3, CDK3, CD93, CD83.   Please confirm if any of these are the gene you intended to analyze, or provide the correct gene name.

---

### ✅ DA-16

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 13.21s
- **Artifacts:** none
- **Notes:** Ghost_Cell_Type does not exist in cell_type annotation. Agent must validate the cell type before running and report which input failed. Must not loop or produce results using the wrong cell type. Source: Yalu NC-D1.


**Response (truncated):** It appears that I am unable to retrieve the conditions from the dataset. The only categorical metadata column available is "louvain," which contains the cell type information.  Since "Ghost_Cell_Type" is not a recognized cell type in the dataset, and I cannot confirm the conditions "disease" and "normal," I cannot proceed with the differential expression analysis.  Please provide a valid cell type

---

### ❌ DA-17

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** FAIL
- **Duration:** 47.52s
- **Artifacts:** none
- **Notes:** FakeCondition does not exist in the dataset. All three inputs (cell type, condition A, condition B) must be validated independently. Agent must report which specific input failed. Must not run DE on partial valid inputs. Source: Yalu NC-D2.


**Failures:**

- must_contain: 'FakeCondition' not in response

**Response (truncated):** It seems that the dataset does not contain a specific cell type labeled as "Early cardiomyocyte." The available cell types in the dataset are:  - CD4 T cells: 1,144 cells (43.4%) - CD14+ Monocytes: 480 cells (18.2%) - B cells: 342 cells (13.0%) - CD8 T cells: 316 cells (12.0%) - NK cells: 154 cells (5.8%) - FCGR3A+ Monocytes: 150 cells (5.7%) - Dendritic cells: 37 cells (1.4%) - Megakaryocytes: 15

---

### ❌ INT-01

- **Category:** interpretation
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 5.11s
- **Artifacts:** 1 table
- **Notes:** KNOWN FAIL — deferred pending T-013 (router wired into core.py). Root cause: query is genuinely ambiguous — "Is this a healthy dataset?" can legitimately be read as a QC request OR a general biology question. Model maps it to QC/Overview intent and calls summarize_qc_metrics before reaching the Interpretation rule. Two prompt patches attempted (Interpretation entry in Fast Track, carve-out in PROBE-FIRST rule) — neither held. Correct fix is a code-level gate in core.py that routes Layer 2b queries away from tool execution. Revisit after T-013.


**Failures:**

- expected_artifact_count: expected 0, got 1 (0 plot(s), 1 table(s))

**Response (truncated):** The quality control (QC) metrics summary for your dataset has been generated successfully. It includes statistics for 2,638 cells across three detected metrics.   To assess the health of the dataset, please review the following metrics: - Total counts (UMI depth) - Number of genes expressed - Mitochondrial percentage  If you have specific thresholds or criteria in mind for determining dataset heal

---

### ✅ INT-02

- **Category:** interpretation
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 5.12s
- **Artifacts:** none
- **Notes:** Unambiguous general biology question — no reasonable interpretation calls for running a tool. Unlike INT-01, this cannot be misconstrued as a data analysis request. Agent must answer from biological knowledge and include a disclaimer. PASS example: "This is based on general biological knowledge, not an analysis of your specific dataset." FAIL = model called a tool or answered without a disclaimer. Diagnostic value: if INT-02 passes and INT-01 fails, INT-01 is confirmed as an ambiguity edge case requiring T-013 (router), not a fundamental behavior gap.


**Response (truncated):** A high mitochondrial percentage in scRNA-seq data typically indicates that the cells may be experiencing stress or are in a state of apoptosis (programmed cell death). This is because mitochondrial genes are often upregulated in response to cellular stress, and a high proportion of mitochondrial RNA can suggest that the cell is not functioning optimally. In quality control assessments, cells with 

---

### ✅ ART-01

- **Category:** artifact_interaction
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 1.44s
- **Artifacts:** none
- **Notes:** Covers the general class of artifact interaction prompts (modify color, resize, export prior result, etc.). No memory system exists — "the last plot" is an unresolvable pointer in a fresh session. Correct behavior: agent asks the user which plot or what analysis they are referring to. PASS = agent clarifies without claiming false success. FAIL = agent claims to have resized something, or flatly refuses with no guidance. Promote to multi-case suite when session memory is implemented.


**Response (truncated):** There is currently no previous plot in this session. Please let me know what type of plot you would like to create, and I can assist you with that!

---

### ✅ OOS-01

- **Category:** out_of_scope
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 5.45s
- **Artifacts:** none
- **Notes:** Dataset integration / batch correction is out of scope. Agent should politely decline and redirect to supported analyses.


**Response (truncated):** I currently don't have the capability to directly integrate datasets or perform batch correction using Harmony. However, I can guide you through the process if you have access to the necessary tools and libraries in your environment.  Here’s a general workflow you can follow to integrate your scRNA-seq datasets and perform batch correction using Harmony:  1. **Load the Datasets**: Load both scRNA-

---

### ✅ OOS-02

- **Category:** out_of_scope
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 2.02s
- **Artifacts:** none
- **Notes:** Trajectory analysis and RNA velocity are out of scope. Agent should decline without crashing. On a raw dataset it may redirect to preprocessing rather than explicitly saying "not supported" — so we only check it did NOT perform the analysis, not the specific refusal wording.


**Response (truncated):** It seems that the current tools available do not support PAGA trajectory analysis or RNA velocity visualization directly. However, I can assist you with other analyses or visualizations related to your dataset.  Would you like to explore differential expression, visualize clusters, or analyze gene expression instead? Please let me know how you'd like to proceed!

---

### ⚠️ RS-01

- **Category:** resource_scale
- **Expected status:** WARN
- **Actual status:** WARN
- **Duration:** 3.99s
- **Artifacts:** none
- **Notes:** 10-gene panel smoke-test. Partial match (some genes not found) is acceptable for WARN.


**Failures:**

- must_contain: 'dot plot' not in response
- requires_plot: no valid PlotResult (none produced)
- expected_artifact_count: expected 1, got 0 (0 plot(s), 0 table(s))

**Response (truncated):** The gene "IL2" was not found in the dataset. Here are some similar gene suggestions: ILF2, IL32, IL27, IL24, PPIL2.   Would you like to proceed with the dotplot using one of these alternatives, or would you prefer to exclude "IL2" from the plot?

---

### ✅ RS-02

- **Category:** resource_scale
- **Expected status:** WARN
- **Actual status:** PASS
- **Duration:** 5.53s
- **Artifacts:** none
- **Notes:** Cluster-level aggregation. Should complete within time limits even for datasets with >10k cells.


**Response (truncated):** Here is the distribution of cells in each Leiden cluster:  - **CD4 T cells**: 1,144 cells (43.4%) - **CD14+ Monocytes**: 480 cells (18.2%) - **B cells**: 342 cells (13.0%) - **CD8 T cells**: 316 cells (12.0%) - **NK cells**: 154 cells (5.8%) - **FCGR3A+ Monocytes**: 150 cells (5.7%) - **Dendritic cells**: 37 cells (1.4%) - **Megakaryocytes**: 15 cells (0.6%)  If you need further analysis or visual

---

### ✅ RS-03

- **Category:** resource_scale
- **Expected status:** WARN
- **Actual status:** PASS
- **Duration:** 14.36s
- **Artifacts:** 1 plot
- **Notes:** Multi-step compound query. Agent may use multiple tool calls in one turn. WARN because not all steps may complete for every dataset.


**Response (truncated):** The dataset contains a total of **2,638 cells**.  I have identified the top 5 marker genes for each cluster:  - **CD4 T cells**: LDHB, CD3D, RPS27, RPS25, RPS12 - **CD14+ Monocytes**: S100A9, S100A8, LYZ, FCN1, FTL - **B cells**: CD79A, CD74, CD79B, HLA-DRA, MS4A1 - **CD8 T cells**: CCL5, NKG7, CST7, GZMA, CTSW - **NK cells**: NKG7, GZMB, PRF1, GNLY, CTSW - **FCGR3A+ Monocytes**: LST1, FCER1G, AIF

---

### ✅ RS-04

- **Category:** resource_scale
- **Expected status:** WARN
- **Actual status:** PASS
- **Duration:** 10.11s
- **Artifacts:** 1 plot
- **Notes:** Tests the heatmap tool, which has no other coverage in the suite. If rank_genes_groups is pre-computed the agent plots directly; otherwise it runs DE first then heatmap. WARN because dataset-dependent.


**Response (truncated):** The heatmap displaying the top 5 marker genes for each cluster has been successfully generated. It includes the following genes across the identified clusters:  - **CD4 T cells**: LDHB, CD3D, RPS27, RPS25, RPS12 - **CD14+ Monocytes**: S100A9, S100A8, LYZ, FCN1, FTL - **B cells**: CD79A, CD74, CD79B, HLA-DRA, MS4A1 - **CD8 T cells**: CCL5, NKG7, CST7, GZMA, CTSW - **NK cells**: NKG7, GZMB, PRF1, GN

---
