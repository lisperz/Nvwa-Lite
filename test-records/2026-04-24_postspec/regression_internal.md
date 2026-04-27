# nvwa-mvp Regression Report

**Generated:** 2026-04-25 01:05:48 UTC  
**Dataset:** `/Users/yxchen/Documents/GitHub/nvwa-mvp/local/data/pbmc_test.h5ad`  
**Model:** `gpt-4o-mini`  
**Duration:** 207.9s  

---

## Summary

| Result | Count |
|--------|-------|
| ✅ Pass  | 25 |
| ❌ Fail  | 2 |
| ⚠️ Warn  | 1 |
| 💥 Error | 0 |
| ⏭ Skip  | 0 |
| **Total** | **28** |

---

## Per-case Results

| case_id | category | status | duration | artifacts | failure reason |
|---------|----------|--------|----------|-----------|----------------|
| DA-01 | data_analysis | ✅ PASS | 2.7s | none | — |
| DA-02 | data_analysis | ✅ PASS | 4.2s | 1 plot | — |
| DA-03 | data_analysis | ✅ PASS | 5.9s | 1 plot | — |
| DA-04 | data_analysis | ✅ PASS | 7.0s | 1 plot | — |
| DA-05 | data_analysis | ✅ PASS | 5.7s | none | — |
| DA-06 | data_analysis | ✅ PASS | 3.4s | 1 table | — |
| DA-07 | data_analysis | ✅ PASS | 12.2s | 2 plots | — |
| DA-08 | data_analysis | ✅ PASS | 4.0s | 1 plot | — |
| DA-09 | data_analysis | ✅ PASS | 3.4s | none | — |
| DA-10 | data_analysis | ✅ PASS | 5.5s | none | — |
| DA-11 | data_analysis | ✅ PASS | 15.1s | none | — |
| DA-12 | data_analysis | ❌ FAIL | 2.6s | none | must_contain: 'imaginary_condition_col' not in response |
| DA-13 | data_analysis | ✅ PASS | 15.0s | none | — |
| DA-14 | data_analysis | ✅ PASS | 9.8s | none | — |
| DA-15 | data_analysis | ✅ PASS | 3.8s | none | — |
| DA-16 | data_analysis | ✅ PASS | 27.6s | none | — |
| DA-17 | data_analysis | ✅ PASS | 8.7s | none | — |
| INT-01 | interpretation | ❌ FAIL | 2.4s | 1 table | expected_artifact_count: expected 0, got 1 (0 plot(s), 1 table(s)) |
| INT-02 | interpretation | ✅ PASS | 8.8s | none | — |
| ART-01 | artifact_interaction | ✅ PASS | 1.5s | none | — |
| OOS-01 | out_of_scope | ✅ PASS | 2.4s | none | — |
| OOS-02 | out_of_scope | ✅ PASS | 2.0s | none | — |
| RS-01 | resource_scale | ⚠️ WARN | 7.1s | none | requires_plot: no valid PlotResult (none produced); expected_artifact_count: expected 1, got 0 (0 plot(s), 0 table(s)) |
| RS-02 | resource_scale | ✅ PASS | 7.3s | none | — |
| RS-03 | resource_scale | ✅ PASS | 13.4s | 1 plot | — |
| RS-04 | resource_scale | ✅ PASS | 14.3s | 1 plot | — |
| SP-01 | data_analysis | ✅ PASS | 3.9s | none | — |
| SP-02 | interpretation | ✅ PASS | 8.0s | none | — |

---

## Top Failure Reasons

| Category | Count |
|----------|-------|
| `assertion_failed` | 1 |
| `other` | 1 |

---

## Detailed Results

### ✅ DA-01

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 2.75s
- **Artifacts:** none

**Response (truncated):** I ran the `dataset_overview` tool to provide a summary of your single-cell RNA-seq dataset. The dataset contains 2,638 cells and 1,838 genes, and it has been identified as originating from human samples.

---

### ✅ DA-02

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 4.21s
- **Artifacts:** 1 plot
- **Notes:** Requires X_umap in obsm and a leiden/louvain column in obs. Fails immediately (pre-flight) if X_umap absent — use a preprocessed h5ad.


**Response (truncated):** The UMAP plot has been generated, colored by the Leiden cluster (louvain). This visualization allows you to see the distribution of cells across different clusters. If you have any further analyses or visualizations in mind, feel free to let me know!

---

### ✅ DA-03

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 5.91s
- **Artifacts:** 1 plot
- **Notes:** User explicitly named a plot type — only 1 artifact expected. A feature plot alongside it is an error (unsolicited output). Per C1 rule: generate only the specified plot type when user names one.


**Response (truncated):** The violin plot of CD3D expression across all clusters (louvain groups) has been successfully generated. This visualization provides insights into the distribution of CD3D expression levels within each cluster.   Would you like to explore further or analyze another gene?

---

### ✅ DA-04

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 7.01s
- **Artifacts:** 1 plot
- **Notes:** Dataset-generic prompt — no hardcoded gene names. Agent discovers marker genes from the dataset and plots them. Works on any preprocessed dataset with clustering results.


**Response (truncated):** The dot plot has been successfully generated, displaying the expression of the top 3 marker genes for each cluster across the louvain groups. This visualization provides insights into the distinct gene expression profiles associated with each cluster.   If you have any further analyses or visualizations in mind, feel free to let me know!

---

### ✅ DA-05

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 5.67s
- **Artifacts:** none
- **Notes:** Requires a leiden/seurat clustering column in obs. Agent may ask for clarification (marker vs pairwise) — that is acceptable.


**Response (truncated):** The differential expression analysis has been completed using the Wilcoxon method for all clusters. A total of 8 groups were analyzed, including CD4 T cells, CD14+ Monocytes, B cells, CD8 T cells, and NK cells. The top 20 marker genes for each group have been identified and stored for further exploration.   Would you like to see the top markers for a specific cluster or download the results?

---

### ✅ DA-06

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.38s
- **Artifacts:** 1 table
- **Notes:** QC stats returned as TableResult (CSV + display). must_contain checks agent response text only — mean/median live in the table, not the text.


**Response (truncated):** The QC metrics summary for the dataset with 2,638 cells includes three detected metrics. You can view the detailed statistics in the generated table.

---

### ✅ DA-07

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 12.22s
- **Artifacts:** 2 plots
- **Notes:** Generic gene expression request with no plot type specified — both a feature plot (spatial) and violin plot (by cell type) are expected. Per C1 rule: when user does not specify a plot type, generate both. Paraphrase of a supported request; previously miscategorized as out_of_scope.


**Response (truncated):** The expression levels of CD3D have been visualized in two ways:  1. **Violin Plot**: This plot shows the distribution of CD3D expression across different cell populations (louvain clusters). It provides insights into the expression levels within each cluster.  2. **Feature Plot**: This UMAP plot displays the spatial expression of CD3D across the entire dataset, highlighting where CD3D is expressed

---

### ✅ DA-08

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.98s
- **Artifacts:** 1 plot
- **Notes:** Paraphrase of a UMAP-by-cluster request. Agent must map this to umap_plot and succeed. Requires preprocessed dataset (X_umap present). Previously miscategorized as out_of_scope.


**Response (truncated):** The UMAP plot has been generated, showing the cells grouped in the low-dimensional embedding space, colored by their cluster assignment (louvain). This visualization allows you to see the distribution and separation of different clusters within the dataset.   If you have any further analyses or visualizations in mind, feel free to let me know!

---

### ✅ DA-09

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 3.44s
- **Artifacts:** none
- **Notes:** Gene validation: FAKEGENE999 does not exist. Agent must report the gene is missing (and ideally suggest alternatives). Updated 2026-04-24: agent consistently says "not found" (matches gene_lookup wording); previous "not present" assertion was a stale artifact of an earlier prompt rev.


**Response (truncated):** The gene "FAKEGENE999" was not found in the dataset, and there were no close matches available. Please check the gene name for any typos or provide a different gene to visualize.

---

### ✅ DA-10

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 5.54s
- **Artifacts:** none
- **Notes:** obs-key validation: nonexistent_batch_col not in adata.obs. Agent lists available columns and offers alternatives — does not use a fixed phrase like "not exist". must_not_contain guards against silent success instead.


**Response (truncated):** The available categorical metadata column in the dataset is "louvain," which contains 8 unique values representing different cell types. Would you like to color the UMAP by this "louvain" column instead?

---

### ✅ DA-11

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 15.07s
- **Artifacts:** none
- **Notes:** obs-field mistaken as gene: leiden is a cluster label, not a gene. Explicit phrasing forces the LLM to look up "leiden" as a gene name so the gene validator fires. Without the constraint the LLM auto-recovers by picking real marker genes.


**Response (truncated):** The gene name "leiden" is not found in the dataset, and there are no close matches available. Please check the gene name and provide a different one if possible.

---

### ❌ DA-12

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** FAIL
- **Duration:** 2.55s
- **Artifacts:** none
- **Notes:** Color-key validation on scatter plot. Previous runs showed agent silently substituting louvain for the bad key. "Do not fall back" + must_not_contain together catch silent substitution. KNOWN FAIL (2026-04-24) — deferred pending prompt cleanup. Rule 6 (OUT-OF-SCOPE DECLINE) over-triggers because the prompt's "imaginary" lexical cue + "do not fall back" directive align with the OOS attractor more strongly than PROBE-FIRST (Rule 2). Structural fix requires either eliminating Rule 6's verbatim reply template (collapses attractor) or moving PROBE-FIRST validation into code so it runs deterministically before the LLM ever sees the prompt. Tracked separately from T-038.


**Failures:**

- must_contain: 'imaginary_condition_col' not in response

**Response (truncated):** This is not supported in this version. I've recorded your request for the development team to consider in the next iteration.

---

### ✅ DA-13

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 14.96s
- **Artifacts:** none
- **Notes:** Invalid cluster ID: cluster 999 does not exist. Agent must not silently succeed — we check it did not produce DE results rather than asserting specific error phrasing (which varies by run).


**Response (truncated):** The available clusters in the dataset are as follows:  - Cluster 0: B cells - Cluster 1: CD14+ Monocytes - Cluster 2: CD4 T cells - Cluster 3: CD8 T cells - Cluster 4: Dendritic cells - Cluster 5: FCGR3A+ Monocytes - Cluster 6: Megakaryocytes - Cluster 7: NK cells  It appears that "cluster 999" does not exist. Please specify a valid cluster from the list above, and I can proceed with the different

---

### ✅ DA-14

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 9.83s
- **Artifacts:** none
- **Notes:** Ribosomal RNA percentage (pct_counts_rb) is not a standard QC metric in most h5ad files. Prompt is explicit about QC metric to avoid the agent interpreting it as a gene expression query (RPL/RPS genes). Agent must report the column is missing. Must not produce a plot. Source: Yalu NC-A1.


**Response (truncated):** The dataset does not currently have a calculated ribosomal RNA percentage. To obtain this information, I can calculate the percentage of counts from ribosomal genes for each cell. Would you like me to proceed with that calculation?

---

### ✅ DA-15

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 3.82s
- **Artifacts:** none
- **Notes:** CD3 is a protein complex encoding multiple genes (CD3E, CD3D, CD3G). Agent should recognize the ambiguity and ask which gene to visualize. Must not silently pick one gene and plot. Hard to assert specific clarification wording — must_not_contain guards against silent success. Source: Yalu NC-C2.


**Response (truncated):** It seems that the gene "CD3" was not found in the dataset. Here are some similar gene names that might be relevant: CHD3, CDV3, CDK3, CD93, CD83.   Please confirm if any of these are the gene you intended to analyze, or provide the correct gene name.

---

### ✅ DA-16

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 27.55s
- **Artifacts:** none
- **Notes:** Ghost_Cell_Type does not exist in cell_type annotation. Agent must validate the cell type before running and report which input failed. Must not loop or produce results using the wrong cell type. Source: Yalu NC-D1.


**Response (truncated):** It appears that I am unable to retrieve the condition metadata from the dataset. Currently, I only have access to the "louvain" column, which contains the cell type information.  Since "Ghost_Cell_Type" is not listed among the available cell types, and I cannot confirm the conditions "disease" and "normal," I cannot proceed with the differential expression analysis.  If you have specific condition

---

### ✅ DA-17

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 8.72s
- **Artifacts:** none
- **Notes:** FakeCondition does not exist in the dataset. All three inputs (cell type, condition A, condition B) must be validated independently. Agent must report which specific input failed. Must not run DE on partial valid inputs. Source: Yalu NC-D2.


**Response (truncated):** It seems that the dataset does not contain a column for conditions or samples, which is necessary to compare "Early cardiomyocyte" in "FakeCondition" vs "Control." The only available metadata column is "louvain," which contains cluster information.  To proceed, please confirm the following: 1. **Is "Early cardiomyocyte" a specific cluster in your dataset?** If so, please provide the corresponding 

---

### ❌ INT-01

- **Category:** interpretation
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 2.41s
- **Artifacts:** 1 table
- **Notes:** KNOWN FAIL — deferred pending T-013 (router wired into core.py). Root cause: query is genuinely ambiguous — "Is this a healthy dataset?" can legitimately be read as a QC request OR a general biology question. Model maps it to QC/Overview intent and calls summarize_qc_metrics before reaching the Interpretation rule. Two prompt patches attempted (Interpretation entry in Fast Track, carve-out in PROBE-FIRST rule) — neither held. Correct fix is a code-level gate in core.py that routes Layer 2b queries away from tool execution. Revisit after T-013.


**Failures:**

- expected_artifact_count: expected 0, got 1 (0 plot(s), 1 table(s))

**Response (truncated):** The quality control (QC) metrics summary for your dataset has been generated successfully. It includes statistics for 2,638 cells across three detected metrics.   Would you like to see the detailed QC metrics?

---

### ✅ INT-02

- **Category:** interpretation
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 8.79s
- **Artifacts:** none
- **Notes:** Unambiguous general biology question — no reasonable interpretation calls for running a tool. Unlike INT-01, this cannot be misconstrued as a data analysis request. Agent must answer from biological knowledge and include a disclaimer. PASS example: "This is based on general biological knowledge, not an analysis of your specific dataset." FAIL = model called a tool or answered without a disclaimer. Diagnostic value: if INT-02 passes and INT-01 fails, INT-01 is confirmed as an ambiguity edge case requiring T-013 (router), not a fundamental behavior gap.


**Response (truncated):** A high mitochondrial percentage in scRNA-seq data typically indicates that the cells may be experiencing stress or undergoing apoptosis (programmed cell death). This is because mitochondrial genes are often upregulated in response to cellular stress, and a high proportion of reads mapping to mitochondrial genes can suggest that the cell is not healthy. In quality control assessments, cells with ex

---

### ✅ ART-01

- **Category:** artifact_interaction
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 1.55s
- **Artifacts:** none
- **Notes:** Covers the general class of artifact interaction prompts (modify color, resize, export prior result, etc.). No memory system exists — "the last plot" is an unresolvable pointer in a fresh session. Correct behavior: agent asks the user which plot or what analysis they are referring to. PASS = agent clarifies without claiming false success. FAIL = agent claims to have resized something, or flatly refuses with no guidance. Promote to multi-case suite when session memory is implemented.


**Response (truncated):** There is currently no previous plot in this session. Please let me know what you'd like to visualize, and I can create a new plot for you!

---

### ✅ OOS-01

- **Category:** out_of_scope
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 2.42s
- **Artifacts:** none
- **Notes:** Dataset integration / batch correction is out of scope. Agent declines via Rule 6 (OUT-OF-SCOPE DECLINE PROTOCOL) with the standard wording "not supported in this version". must_not_contain still guards against fake-success substitution.


**Response (truncated):** This is not supported in this version. I've recorded your request for the development team to consider in the next iteration.

---

### ✅ OOS-02

- **Category:** out_of_scope
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 1.96s
- **Artifacts:** none
- **Notes:** Trajectory analysis and RNA velocity are out of scope. Agent should decline without crashing. On a raw dataset it may redirect to preprocessing rather than explicitly saying "not supported" — so we only check it did NOT perform the analysis, not the specific refusal wording.


**Response (truncated):** This is not supported in this version. I've recorded your request for the development team to consider in the next iteration.

---

### ⚠️ RS-01

- **Category:** resource_scale
- **Expected status:** WARN
- **Actual status:** WARN
- **Duration:** 7.09s
- **Artifacts:** none
- **Notes:** 10-gene panel smoke-test. Partial match (some genes not found) is acceptable for WARN.


**Failures:**

- requires_plot: no valid PlotResult (none produced)
- expected_artifact_count: expected 1, got 0 (0 plot(s), 0 table(s))

**Response (truncated):** The gene "IL2" was not found in the dataset. Here are some similar gene suggestions: ILF2, IL32, IL27, IL24, PPIL2.   Please confirm if you would like to proceed with one of these alternatives or if you want to exclude "IL2" from the dot plot.

---

### ✅ RS-02

- **Category:** resource_scale
- **Expected status:** WARN
- **Actual status:** PASS
- **Duration:** 7.27s
- **Artifacts:** none
- **Notes:** Cluster-level aggregation. Should complete within time limits even for datasets with >10k cells.


**Response (truncated):** Here is the distribution of cells in each Leiden cluster:  - **CD4 T cells**: 1,144 cells (43.4%) - **CD14+ Monocytes**: 480 cells (18.2%) - **B cells**: 342 cells (13.0%) - **CD8 T cells**: 316 cells (12.0%) - **NK cells**: 154 cells (5.8%) - **FCGR3A+ Monocytes**: 150 cells (5.7%) - **Dendritic cells**: 37 cells (1.4%) - **Megakaryocytes**: 15 cells (0.6%)  If you need further analysis or visual

---

### ✅ RS-03

- **Category:** resource_scale
- **Expected status:** WARN
- **Actual status:** PASS
- **Duration:** 13.38s
- **Artifacts:** 1 plot
- **Notes:** Multi-step compound query. Agent may use multiple tool calls in one turn. WARN because not all steps may complete for every dataset.


**Response (truncated):** The dataset contains a total of **2,638 cells**.  I have identified the top 5 marker genes for each cluster as follows:  - **CD4 T cells**: LDHB, CD3D, RPS27, RPS25, RPS12 - **CD14+ Monocytes**: S100A9, S100A8, LYZ, FCN1, FTL - **B cells**: CD79A, CD74, CD79B, HLA-DRA, MS4A1 - **CD8 T cells**: CCL5, NKG7, CST7, GZMA, CTSW - **NK cells**: NKG7, GZMB, PRF1, GNLY, CTSW - **FCGR3A+ Monocytes**: LST1, 

---

### ✅ RS-04

- **Category:** resource_scale
- **Expected status:** WARN
- **Actual status:** PASS
- **Duration:** 14.27s
- **Artifacts:** 1 plot
- **Notes:** Tests the heatmap tool, which has no other coverage in the suite. If rank_genes_groups is pre-computed the agent plots directly; otherwise it runs DE first then heatmap. WARN because dataset-dependent.


**Response (truncated):** The heatmap of the top 5 marker genes for each cluster has been successfully generated. It displays the expression levels of the following genes across the different clusters:  - **CD4 T cells**: LDHB, CD3D, RPS27, RPS25, RPS12 - **CD14+ Monocytes**: S100A9, S100A8, LYZ, FCN1, FTL - **B cells**: CD79A, CD74, CD79B, HLA-DRA, MS4A1 - **CD8 T cells**: CCL5, NKG7, CST7, GZMA, CTSW - **NK cells**: NKG7

---

### ✅ SP-01

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.89s
- **Artifacts:** none
- **Notes:** Spec-pipeline happy path. Router matches "data overview" → layer=2a, task_type=dataset_info. Extractor then picks dataset_overview from REGISTRY (only @register tool today; dataset_info is legacy and not in REGISTRY), resolver pass-through, validator ok, dispatch reads adata.uns["nvwa_meta"] populated by detect_species + classify_columns in scripts/run_tests.py:load_adata, responder narrates, gatekeeper passes. End-to-end new-path validation. Assertions are dataset-agnostic since content varies across the three target datasets. Prompt chosen to differ from DA-01 so router matches a different keyword and this case is independent.


**Response (truncated):** I ran the `dataset_overview` tool to provide a summary of your single-cell RNA-seq dataset. The dataset contains 2,638 cells and 1,838 genes, and it has been identified as belonging to the human species based on the gene-symbol case.

---

### ✅ SP-02

- **Category:** interpretation
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 7.99s
- **Artifacts:** none
- **Notes:** 2b plain-LLM path. Router matches "explain" (Tier-1 2b pattern) and returns layer=2b before Tier-2 keyword scan runs. Core.py routes to _run_plain_llm → ChatOpenAI.invoke without bind_tools → text response, no tool call. output_guard remains as a safety net for fabricated artifact claims (checked via must_not_contain phrases). Prompt intentionally distinct from INT-01/02 so SP-02 behavior doesn't couple to future INT prompt tuning. Avoid the word "embedding" here — it matches the umap_plot keyword set and would flip router to 2a.


**Response (truncated):** In single-cell RNA sequencing (scRNA-seq) data, similar cells tend to form tight clusters due to several key factors:  1. **Gene Expression Profiles**: Cells that share similar functions or developmental stages often exhibit comparable gene expression patterns. This similarity in gene expression leads to clustering in high-dimensional space, as cells with similar profiles are closer together.  2. 

---
