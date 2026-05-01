# nvwa-mvp Regression Report

**Generated:** 2026-05-01 17:33:07 UTC  
**Dataset:** `/Users/yxchen/GitHub/nvwa-mvp/local/data/pbmc_test.h5ad`  
**Model:** `gpt-4o-mini`  
**Duration:** 119.7s  

---

## Summary

| Result | Count |
|--------|-------|
| ✅ Pass  | 24 |
| ❌ Fail  | 3 |
| ⚠️ Warn  | 1 |
| 💥 Error | 0 |
| ⏭ Skip  | 0 |
| **Total** | **28** |

---

## Per-case Results

| case_id | category | status | routed_via | duration | artifacts | failure reason |
|---------|----------|--------|------------|----------|-----------|----------------|
| DA-01 | data_analysis | ✅ PASS | spec | 4.9s | none | — |
| DA-02 | data_analysis | ✅ PASS | spec | 3.5s | 1 plot | — |
| DA-03 | data_analysis | ✅ PASS | spec | 3.6s | 1 plot | — |
| DA-04 | data_analysis | ✅ PASS | legacy | 14.8s | 1 plot | — |
| DA-05 | data_analysis | ✅ PASS | legacy | 5.1s | none | — |
| DA-06 | data_analysis | ✅ PASS | spec | 3.5s | 1 table | — |
| DA-07 | data_analysis | ❌ FAIL | spec | 1.3s | none | requires_plot: no valid PlotResult (none produced); expected_artifact_count: expected 2, got 0 (0 plot(s), 0 table(s)) |
| DA-08 | data_analysis | ✅ PASS | spec | 4.3s | 1 plot | — |
| DA-09 | data_analysis | ✅ PASS | spec | 1.0s | none | — |
| DA-10 | data_analysis | ✅ PASS | legacy | 4.9s | none | — |
| DA-11 | data_analysis | ✅ PASS | spec | 1.0s | none | — |
| DA-12 | data_analysis | ❌ FAIL | legacy | 5.0s | none | must_contain: 'imaginary_condition_col' not in response |
| DA-13 | data_analysis | ✅ PASS | spec | 1.5s | none | — |
| DA-14 | data_analysis | ✅ PASS | spec | 3.0s | none | — |
| DA-15 | data_analysis | ✅ PASS | spec | 1.1s | none | — |
| DA-16 | data_analysis | ✅ PASS | spec | 1.6s | none | — |
| DA-17 | data_analysis | ✅ PASS | spec | 1.5s | none | — |
| INT-01 | interpretation | ❌ FAIL | spec | 2.9s | 1 table | expected_artifact_count: expected 0, got 1 (0 plot(s), 1 table(s)) |
| INT-02 | interpretation | ✅ PASS | legacy | 3.5s | none | — |
| ART-01 | artifact_interaction | ✅ PASS | legacy | 2.5s | none | — |
| OOS-01 | out_of_scope | ✅ PASS | legacy | 2.2s | none | — |
| OOS-02 | out_of_scope | ✅ PASS | legacy | 2.6s | none | — |
| RS-01 | resource_scale | ⚠️ WARN | spec | 2.1s | none | must_contain: 'dot plot' not in response; requires_plot: no valid PlotResult (none produced) |
| RS-02 | resource_scale | ✅ PASS | legacy | 5.6s | none | — |
| RS-03 | resource_scale | ✅ PASS | legacy | 15.4s | 1 plot | — |
| RS-04 | resource_scale | ✅ PASS | legacy | 12.8s | 1 plot | — |
| SP-01 | data_analysis | ✅ PASS | spec | 2.4s | none | — |
| SP-02 | interpretation | ✅ PASS | plain_llm | 5.7s | none | — |

---

## Top Failure Reasons

| Category | Count |
|----------|-------|
| `missing_artifact` | 1 |
| `assertion_failed` | 1 |
| `other` | 1 |

---

## Detailed Results

### ✅ DA-01

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 4.91s
- **Artifacts:** none
- **Tool called:** True

**Response (truncated):** Here's a summary of your dataset: it contains 2,638 cells and 1,838 genes, and the species identified is human. This overview provides a high-level look at the data structure (dataset_overview).

---

### ✅ DA-02

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.50s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Requires X_umap in obsm and a leiden/louvain column in obs. Fails immediately (pre-flight) if X_umap absent — use a preprocessed h5ad.


**Response (truncated):** Here's the UMAP scatter plot rendered for 2,638 cells, colored by the Louvain clustering results (I used UMAP to visualize the data). Note that I adjusted your input from 'leiden' to 'louvain' for consistency.

---

### ✅ DA-03

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.63s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** User explicitly named a plot type — only 1 artifact expected. A feature plot alongside it is an error (unsolicited output). Per C1 rule: generate only the specified plot type when user names one.
Gene changed from CD3D to MS4A1 on 2026-04-29 (violin_plot spec-pipeline migration). CD3D exists in pbmc_test's adata.raw but NOT in adata.var_names; the legacy violin_plot @tool silently fell back to adata.raw, masking the mismatch. The new gene_lookup primitive (src/domain/resolver/gene_lookup.py) restricts to adata.var_names per the C5 architectural lock — raw fallback is a hidden intent translation (G4) that gives the user data with different normalization without their knowledge. MS4A1 is a B-cell marker present in adata.var_names, preserving the test's intent (violin plot of A gene across clusters works) without depending on the raw-fallback anti-pattern.


**Response (truncated):** Here's the violin plot of MS4A1 expression across the louvain clusters for 2,638 cells (violin_plot). I used louvain (you wrote 'cluster', case-corrected) as the grouping variable.

---

### ✅ DA-04

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 14.82s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Dataset-generic prompt — no hardcoded gene names. Agent discovers marker genes from the dataset and plots them. Works on any preprocessed dataset with clustering results.


**Response (truncated):** The dot plot has been successfully generated, displaying the expression of the top 3 marker genes for each cluster across the louvain groups.   This visualization allows you to see the expression levels and the fraction of cells expressing each gene within the different clusters. If you have any further analyses or visualizations in mind, feel free to let me know!

---

### ✅ DA-05

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 5.12s
- **Artifacts:** none
- **Tool called:** True
- **Notes:** Requires a leiden/seurat clustering column in obs. Agent may ask for clarification (marker vs pairwise) — that is acceptable.


**Response (truncated):** The differential expression analysis has been completed using the Wilcoxon method across all clusters. A total of 8 groups were analyzed, including CD4 T cells, CD14+ Monocytes, B cells, CD8 T cells, and NK cells. The top 20 marker genes for each group have been identified and stored.   Would you like to see the top markers for a specific cluster or download the results?

---

### ✅ DA-06

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.48s
- **Artifacts:** 1 table
- **Tool called:** True
- **Notes:** QC stats returned as TableResult (CSV + display). must_contain checks agent response text only — mean/median live in the table, not the text.


**Response (truncated):** Here's the QC summary for your dataset, which includes 2,638 cells across three metrics. All metrics are within typical ranges, indicating no issues with elevated mitochondrial fraction or low gene detection (qc_summary_table).

---

### ❌ DA-07

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 1.34s
- **Artifacts:** none
- **Tool called:** False
- **Notes:** Generic gene expression request with no plot type specified — both a feature plot (spatial) and violin plot (by cell type) are expected. Per C1 rule: when user does not specify a plot type, generate both. Paraphrase of a supported request; previously miscategorized as out_of_scope.


**Failures:**

- requires_plot: no valid PlotResult (none produced)
- expected_artifact_count: expected 2, got 0 (0 plot(s), 0 table(s))

**Response (truncated):** I don't see a column named 'cell type' in the dataset. Available columns: n_genes, percent_mito, n_counts. I don't see 'CD3D' as a gene in the dataset. Did you mean: CDV3, CD93, CD37?

---

### ✅ DA-08

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 4.34s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Paraphrase of a UMAP-by-cluster request. Agent must map this to umap_plot and succeed. Requires preprocessed dataset (X_umap present). Previously miscategorized as out_of_scope.


**Response (truncated):** Here's the UMAP scatter plot showing 2,638 cells colored by their Louvain clusters. This visualization helps you explore the relationships and distributions of different cell types within your dataset (umap_plot).

---

### ✅ DA-09

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 1.03s
- **Artifacts:** none
- **Tool called:** False
- **Notes:** Gene validation: FAKEGENE999 does not exist. Agent must report the gene is missing. Asserts the user's literal input is echoed (proving the agent processed it without silent substitution) — more robust than coupling to specific failure-message wording. Updated 2026-04-30: switched from "not found" to "FAKEGENE999" echo (prior wording assertion kept rotting as renderer wording improved; behavior-level echo is stable across wording changes).


**Response (truncated):** I don't see 'FAKEGENE999' as a gene in the dataset.

---

### ✅ DA-10

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 4.91s
- **Artifacts:** none
- **Tool called:** True
- **Notes:** obs-key validation: nonexistent_batch_col not in adata.obs. Agent lists available columns and offers alternatives — does not use a fixed phrase like "not exist". must_not_contain guards against silent success instead.


**Response (truncated):** The available categorical metadata column in the dataset is "louvain," which contains 8 unique values representing different cell types.   Would you like to color the UMAP by the "louvain" column instead?

---

### ✅ DA-11

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 1.04s
- **Artifacts:** none
- **Tool called:** False
- **Notes:** obs-field mistaken as gene: leiden is a cluster label, not a gene. Explicit phrasing forces the LLM to look up "leiden" as a gene name so the gene validator fires. Without the constraint the LLM auto-recovers by picking real marker genes. Updated 2026-04-30: switched must_contain from "not found" to "leiden" echo (behavior-level echo is stable across wording changes; previous wording assertion was tied to a specific renderer phrasing that improved).


**Response (truncated):** I don't see 'leiden' as a gene in the dataset.

---

### ❌ DA-12

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** FAIL
- **Duration:** 4.99s
- **Artifacts:** none
- **Tool called:** False
- **Notes:** Color-key validation on scatter plot. Previous runs showed agent silently substituting louvain for the bad key. "Do not fall back" + must_not_contain together catch silent substitution. KNOWN FAIL (2026-04-24) — deferred pending prompt cleanup. Rule 6 (OUT-OF-SCOPE DECLINE) over-triggers because the prompt's "imaginary" lexical cue + "do not fall back" directive align with the OOS attractor more strongly than PROBE-FIRST (Rule 2). Structural fix requires either eliminating Rule 6's verbatim reply template (collapses attractor) or moving PROBE-FIRST validation into code so it runs deterministically before the LLM ever sees the prompt. Tracked separately from T-038.


**Failures:**

- must_contain: 'imaginary_condition_col' not in response

**Response (truncated):** This is not supported in this version. I've recorded your request for the development team to consider in the next iteration.

---

### ✅ DA-13

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 1.46s
- **Artifacts:** none
- **Tool called:** False
- **Notes:** Invalid cluster ID: cluster 999 does not exist. Agent must not silently succeed — we check it did not produce DE results rather than asserting specific error phrasing (which varies by run).


**Response (truncated):** I don't see '999' as a condition in the dataset. Available conditions: CD4 T cells, B cells, CD14+ Monocytes. I don't see 'all other cells' as a condition in the dataset. Available conditions: CD4 T cells, B cells, CD14+ Monocytes.

---

### ✅ DA-14

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 2.96s
- **Artifacts:** none
- **Tool called:** True
- **Notes:** Ribosomal RNA percentage (pct_counts_rb) is not a standard QC metric in most h5ad files. Prompt is explicit about QC metric to avoid the agent interpreting it as a gene expression query (RPL/RPS genes). Agent must report the column is missing. Must not produce a plot. Source: Yalu NC-A1.


**Response (truncated):** It seems there was an issue with your request. The quality control metric 'pct_counts_rb' you specified is not available in the dataset. The available metrics are 'n_counts', 'n_genes', and 'percent_mito'. Please check your input and try again with one of these metrics.

---

### ✅ DA-15

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 1.13s
- **Artifacts:** none
- **Tool called:** False
- **Notes:** CD3 is a protein complex encoding multiple genes (CD3E, CD3D, CD3G). Agent should recognize the ambiguity and ask which gene to visualize. Must not silently pick one gene and plot. Hard to assert specific clarification wording — must_not_contain guards against silent success. Source: Yalu NC-C2.


**Response (truncated):** I don't see 'CD3' as a gene in the dataset. Did you mean: CDV3, CD93, CD37?

---

### ✅ DA-16

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 1.59s
- **Artifacts:** none
- **Tool called:** False
- **Notes:** Ghost_Cell_Type does not exist in cell_type annotation. Agent must validate the cell type before running and report which input failed. Must not loop or produce results using the wrong cell type. Source: Yalu NC-D1.


**Response (truncated):** I don't see a column named 'condition' in the dataset. Available columns: n_genes, percent_mito, n_counts. I don't see a column named 'cell type' in the dataset. Available columns: n_genes, percent_mito, n_counts. I don't see 'disease' as a condition in the dataset. I don't see 'normal' as a condition in the dataset. I don't see 'Ghost_Cell_Type' as a cell type in the dataset.

---

### ✅ DA-17

- **Category:** data_analysis
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 1.49s
- **Artifacts:** none
- **Tool called:** False
- **Notes:** FakeCondition does not exist in the dataset. All three inputs (cell type, condition A, condition B) must be validated independently. Agent must report which specific input failed. Must not run DE on partial valid inputs. Source: Yalu NC-D2.


**Response (truncated):** I don't see a column named 'condition' in the dataset. Available columns: n_genes, percent_mito, n_counts. I don't see a column named 'cell type' in the dataset. Available columns: n_genes, percent_mito, n_counts. I don't see 'FakeCondition' as a condition in the dataset. I don't see 'Control' as a condition in the dataset. I don't see 'Early cardiomyocyte' as a cell type in the dataset.

---

### ❌ INT-01

- **Category:** interpretation
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 2.85s
- **Artifacts:** 1 table
- **Tool called:** True
- **Notes:** KNOWN FAIL — deferred pending T-013 (router wired into core.py). Root cause: query is genuinely ambiguous — "Is this a healthy dataset?" can legitimately be read as a QC request OR a general biology question. Model maps it to QC/Overview intent and calls summarize_qc_metrics before reaching the Interpretation rule. Two prompt patches attempted (Interpretation entry in Fast Track, carve-out in PROBE-FIRST rule) — neither held. Correct fix is a code-level gate in core.py that routes Layer 2b queries away from tool execution. Revisit after T-013.


**Failures:**

- expected_artifact_count: expected 0, got 1 (0 plot(s), 1 table(s))

**Response (truncated):** Here's the QC summary for your dataset, which includes 2,638 cells across three metrics. All metrics are within typical ranges, indicating no issues with elevated mitochondrial fraction or low gene detection (qc_summary_table).

---

### ✅ INT-02

- **Category:** interpretation
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.48s
- **Artifacts:** none
- **Tool called:** False
- **Notes:** Unambiguous general biology question — no reasonable interpretation calls for running a tool. Unlike INT-01, this cannot be misconstrued as a data analysis request. Agent must answer from biological knowledge and include a disclaimer. PASS example: "This is based on general biological knowledge, not an analysis of your specific dataset." FAIL = model called a tool or answered without a disclaimer. Diagnostic value: if INT-02 passes and INT-01 fails, INT-01 is confirmed as an ambiguity edge case requiring T-013 (router), not a fundamental behavior gap.


**Response (truncated):** A high mitochondrial percentage in scRNA-seq data typically indicates that the cells may be experiencing stress or are in a state of apoptosis (programmed cell death). This is because mitochondrial genes are often upregulated in response to cellular stress, and a high proportion of mitochondrial RNA can suggest that the cell is not functioning optimally. In quality control assessments, cells with 

---

### ✅ ART-01

- **Category:** artifact_interaction
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 2.55s
- **Artifacts:** none
- **Tool called:** False
- **Notes:** Covers the general class of artifact interaction prompts (modify color, resize, export prior result, etc.). No memory system exists — "the last plot" is an unresolvable pointer in a fresh session. Correct behavior: agent asks the user which plot or what analysis they are referring to. PASS = agent clarifies without claiming false success. FAIL = agent claims to have resized something, or flatly refuses with no guidance. Promote to multi-case suite when session memory is implemented.


**Response (truncated):** There is currently no previous plot in this session. Please let me know what type of plot you would like to create, and I can assist you with that!

---

### ✅ OOS-01

- **Category:** out_of_scope
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 2.16s
- **Artifacts:** none
- **Tool called:** False
- **Notes:** Dataset integration / batch correction is out of scope. Agent declines via Rule 6 (OUT-OF-SCOPE DECLINE PROTOCOL) with the standard wording "not supported in this version". must_not_contain still guards against fake-success substitution.


**Response (truncated):** This is not supported in this version. I've recorded your request for the development team to consider in the next iteration.

---

### ✅ OOS-02

- **Category:** out_of_scope
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 2.58s
- **Artifacts:** none
- **Tool called:** False
- **Notes:** Trajectory analysis and RNA velocity are out of scope. Agent should decline without crashing. On a raw dataset it may redirect to preprocessing rather than explicitly saying "not supported" — so we only check it did NOT perform the analysis, not the specific refusal wording.


**Response (truncated):** This is not supported in this version. I've recorded your request for the development team to consider in the next iteration.

---

### ⚠️ RS-01

- **Category:** resource_scale
- **Expected status:** WARN
- **Actual status:** WARN
- **Duration:** 2.14s
- **Artifacts:** none
- **Tool called:** False
- **Notes:** 10-gene panel smoke-test. Partial match (some genes not found) is acceptable for WARN.


**Failures:**

- must_contain: 'dot plot' not in response
- requires_plot: no valid PlotResult (none produced)
- expected_artifact_count: expected 1, got 0 (0 plot(s), 0 table(s))

**Response (truncated):** I don't see 'CD3D' as a gene in the dataset. Did you mean: CDV3, CD93, CD37? I don't see 'LYZ' as a gene in the dataset. I don't see 'CD8A' as a gene in the dataset. Did you mean: CDA, CD82, CD48? I don't see 'CD4' as a gene in the dataset. Did you mean: CD48, CD47, CDC40? I don't see 'FOXP3' as a gene in the dataset. Did you mean: SOX13, FOPNL, FOLR3? I don't see 'IL2' as a gene in the dataset. D

---

### ✅ RS-02

- **Category:** resource_scale
- **Expected status:** WARN
- **Actual status:** PASS
- **Duration:** 5.59s
- **Artifacts:** none
- **Tool called:** True
- **Notes:** Cluster-level aggregation. Should complete within time limits even for datasets with >10k cells.


**Response (truncated):** Here is the distribution of cells in each Leiden cluster:  - **CD4 T cells**: 1,144 cells (43.4%) - **CD14+ Monocytes**: 480 cells (18.2%) - **B cells**: 342 cells (13.0%) - **CD8 T cells**: 316 cells (12.0%) - **NK cells**: 154 cells (5.8%) - **FCGR3A+ Monocytes**: 150 cells (5.7%) - **Dendritic cells**: 37 cells (1.4%) - **Megakaryocytes**: 15 cells (0.6%)  If you need further analysis or visual

---

### ✅ RS-03

- **Category:** resource_scale
- **Expected status:** WARN
- **Actual status:** PASS
- **Duration:** 15.41s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Multi-step compound query. Agent may use multiple tool calls in one turn. WARN because not all steps may complete for every dataset.


**Response (truncated):** The dataset contains a total of **2,638 cells**.  I have identified the top 5 marker genes for each cluster:  - **CD4 T cells**: LDHB, CD3D, RPS27, RPS25, RPS12 - **CD14+ Monocytes**: S100A9, S100A8, LYZ, FCN1, FTL - **B cells**: CD79A, CD74, CD79B, HLA-DRA, MS4A1 - **CD8 T cells**: CCL5, NKG7, CST7, GZMA, CTSW - **NK cells**: NKG7, GZMB, PRF1, GNLY, CTSW - **FCGR3A+ Monocytes**: LST1, FCER1G, AIF

---

### ✅ RS-04

- **Category:** resource_scale
- **Expected status:** WARN
- **Actual status:** PASS
- **Duration:** 12.78s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Tests the heatmap tool, which has no other coverage in the suite. If rank_genes_groups is pre-computed the agent plots directly; otherwise it runs DE first then heatmap. WARN because dataset-dependent.


**Response (truncated):** The heatmap of the top 5 marker genes for each cluster has been successfully generated. It displays the expression levels of the following genes across the different clusters:  - **CD4 T cells**: LDHB, CD3D, RPS27, RPS25, RPS12 - **CD14+ Monocytes**: S100A9, S100A8, LYZ, FCN1, FTL - **B cells**: CD79A, CD74, CD79B, HLA-DRA, MS4A1 - **CD8 T cells**: CCL5, NKG7, CST7, GZMA, CTSW - **NK cells**: NKG7

---

### ✅ SP-01

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 2.38s
- **Artifacts:** none
- **Tool called:** True
- **Notes:** Spec-pipeline happy path. Router matches "data overview" → layer=2a, task_type=dataset_info. Extractor then picks dataset_overview from REGISTRY (only @register tool today; dataset_info is legacy and not in REGISTRY), resolver pass-through, validator ok, dispatch reads adata.uns["nvwa_meta"] populated by detect_species + classify_columns in scripts/run_tests.py:load_adata, responder narrates, gatekeeper passes. End-to-end new-path validation. Assertions are dataset-agnostic since content varies across the three target datasets. Prompt chosen to differ from DA-01 so router matches a different keyword and this case is independent.


**Response (truncated):** Here's a summary of your dataset: it contains 2,638 cells and 1,838 genes, with the species identified as human. This overview was generated using the dataset_overview tool.

---

### ✅ SP-02

- **Category:** interpretation
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 5.70s
- **Artifacts:** none
- **Tool called:** False
- **Notes:** 2b plain-LLM path. Router matches "explain" (Tier-1 2b pattern) and returns layer=2b before Tier-2 keyword scan runs. Core.py routes to _run_plain_llm → ChatOpenAI.invoke without bind_tools → text response, no tool call. output_guard remains as a safety net for fabricated artifact claims (checked via must_not_contain phrases). Prompt intentionally distinct from INT-01/02 so SP-02 behavior doesn't couple to future INT prompt tuning. Avoid the word "embedding" here — it matches the umap_plot keyword set and would flip router to 2a.


**Response (truncated):** In single-cell RNA sequencing (scRNA-seq) data, similar cells tend to form tight clusters due to several key factors:  1. **Gene Expression Profiles**: Cells that share similar functions or developmental stages often exhibit comparable gene expression patterns. This similarity in gene expression leads to clustering in high-dimensional space, as cells with similar profiles are closer together.  2. 

---
