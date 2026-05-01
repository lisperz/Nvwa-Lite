# nvwa-mvp Regression Report

**Generated:** 2026-05-01 03:51:08 UTC  
**Dataset:** `/Users/yxchen/GitHub/nvwa-mvp/local/data/lung_combined_slim_v3.h5ad`  
**Model:** `gpt-4o-mini`  
**Duration:** 81.4s  

---

## Summary

| Result | Count |
|--------|-------|
| ✅ Pass  | 14 |
| ❌ Fail  | 3 |
| ⚠️ Warn  | 1 |
| 💥 Error | 0 |
| ⏭ Skip  | 0 |
| **Total** | **18** |

---

## Per-case Results

| case_id | category | status | routed_via | duration | artifacts | failure reason |
|---------|----------|--------|------------|----------|-----------|----------------|
| DAI-01 | data_analysis | ✅ PASS | spec | 3.7s | none | — |
| DAI-02 | data_analysis | ✅ PASS | spec | 2.5s | 1 plot | — |
| DAI-03 | data_analysis | ✅ PASS | spec | 2.6s | 1 plot | — |
| DAI-04 | data_analysis | ✅ PASS | spec | 2.2s | 1 plot | — |
| DAI-05 | data_analysis | ✅ PASS | spec | 3.2s | 1 plot | — |
| DAI-07 | data_analysis | ✅ PASS | spec | 3.4s | 1 plot | — |
| DAI-08 | data_analysis | ✅ PASS | spec | 3.3s | 1 plot | — |
| DAI-09 | data_analysis | ✅ PASS | spec | 3.3s | 1 plot | — |
| DAI-10 | data_analysis | ✅ PASS | spec | 3.6s | 1 plot | — |
| DAI-11 | data_analysis | ✅ PASS | spec | 3.7s | 1 plot | — |
| DAI-12 | data_analysis | ❌ FAIL | spec | 7.1s | 1 table | must_contain: 'heatmap' not in response; requires_plot: no valid PlotResult (none produced) |
| DAI-13 | data_analysis | ❌ FAIL | spec | 7.0s | 1 table | must_contain: 'heatmap' not in response; requires_plot: no valid PlotResult (none produced) |
| DAI-14 | data_analysis | ❌ FAIL | spec | 9.7s | 1 table | must_contain: 'heatmap' not in response; requires_plot: no valid PlotResult (none produced) |
| DAI-15 | data_analysis | ✅ PASS | spec | 3.8s | 1 plot | — |
| DAI-16 | out_of_scope | ✅ PASS | legacy | 7.9s | none | — |
| DAI-17 | data_analysis | ⚠️ WARN | spec | 3.3s | none | must_contain: 'volcano' not in response; requires_plot: no valid PlotResult (none produced) |
| DAI-18 | data_analysis | ✅ PASS | spec | 5.1s | 1 table | — |
| DAI-19 | data_analysis | ✅ PASS | spec | 5.5s | 1 table | — |

---

## Top Failure Reasons

| Category | Count |
|----------|-------|
| `missing_artifact` | 3 |

---

## Detailed Results

### ✅ DAI-01

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.71s
- **Artifacts:** none
- **Tool called:** True
- **Notes:** Paraphrase of DA-01 (dataset overview). Lowercase, no punctuation.


**Response (truncated):** Here's a summary of your dataset: it contains 15,719 cells and 32,285 genes, with 19 annotated cell types including EC, B cell, T cell, and Fibroblast. I noticed that the 'ident' column has 26 unique values, and I'm unsure if this represents an experimental condition, sample ID, or batch label. Additionally, the gene symbols are in Title-case, which is common to either mouse or rat; could you clar

---

### ✅ DAI-02

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 2.48s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Paraphrase of DA-02. No color-by specified — agent may default to clusters or cell_type; both acceptable.


**Response (truncated):** Here's the UMAP scatter plot showing 15,719 cells colored by their cell type. This visualization helps you explore the distribution and relationships of different cell types within your dataset (umap_plot).

---

### ✅ DAI-03

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 2.63s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** UMAP colored by the cell_type obs column.


**Response (truncated):** Here's the UMAP scatter plot colored by cell type, showing a visualization of 15,719 cells (umap_plot). I used 'celltype' for the color coding (you wrote 'cell type', case-corrected).

---

### ✅ DAI-04

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 2.18s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Multi-gene expression visualization. Per Yalu §3C.1: "Multiple genes · All cell types · No condition split" → ONE dot_plot.
Updated 2026-04-29: previously the test accepted whatever artifact count the legacy LangChain path produced (typically 4 — feature + violin per gene = 2×2). New expectation: ONE dot_plot covering both genes side-by-side, per Yalu spec. Test will FAIL until dot_plot migrates to the spec pipeline (currently legacy `@tool`, falls back to LangChain which still applies the legacy "feature + violin per gene" rule from prompts.py:74). Fix lands when the dot_plot PR ships.


**Response (truncated):** Here's the dot plot showing the expression of Fabp4 and Fabp5 across different cell types. I used the `dot_plot` tool to visualize this data.

---

### ✅ DAI-05

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.25s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Tests Layer 2a parameter resolution: when user says "cell types", agent must group by the cell_type obs column, not default to seurat_clusters. Rewritten from a context-dependent follow-up in Dai's session so the prompt is self-contained.


**Response (truncated):** Here's the feature plot of Fabp4 on UMAP for 15,719 cells, showing the expression levels across the dataset. This visualization was created using the `feature_plot` tool.

---

### ✅ DAI-07

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.38s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Pairwise condition comparison. WT and CKO are condition labels in obs (not cluster names). Agent should subset to WT/CKO and plot — violin or DE both acceptable. Known-failing in Dai's session on condition resolution (see DAI-18/19).


**Response (truncated):** Here's the violin plot of Fabp4 expression across conditions for 9,640 cells, specifically subsetted to the 'WT' and 'CKO' groups. I used the `violin_plot` tool, and please note that the grouping was adjusted from 'condition' to 'orig.ident' for this analysis.

---

### ✅ DAI-08

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.29s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Tests that WT and CKO are rendered on a single combined violin, not as separate panels. Rewritten from a context-dependent follow-up in Dai's session (where DAI-07 had produced split violins) so the "single combined plot" constraint is explicit in the prompt.


**Response (truncated):** Here's the violin plot of Fabp4 expression across conditions for 9,640 cells, specifically subsetted to the 'WT' and 'CKO' groups. I used `violin_plot` for this visualization, and please note that the grouping was adjusted from 'condition' to 'orig.ident' for consistency.

---

### ✅ DAI-09

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.26s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** "Different samples" is ambiguous — could mean orig.ident, condition, or another obs column. Agent should resolve or ask. Note: user typed FABP4 (uppercase) while the gene in the dataset is Fabp4 (mouse).


**Response (truncated):** Here's the violin plot of Fabp4 expression across different conditions (orig.ident) for 15,719 cells. I used the `violin_plot` tool, and please note that the condition you specified was adjusted to the original identifier for clarity.

---

### ✅ DAI-10

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.57s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Three-condition comparison. All three labels must be resolved in the condition obs column. Paraphrase variant of DAI-07.


**Response (truncated):** Here's the violin plot of Fabp4 expression across different conditions (orig.ident) for 15,719 cells, specifically subsetted to the conditions WT, CKO, and TKO. I used the `violin_plot` tool.

---

### ✅ DAI-11

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.65s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Generic "by condition" — agent should auto-detect condition-like columns and split by all levels. Third phrasing of the same underlying request (DAI-07 / DAI-10 / DAI-11) — tests paraphrase robustness.


**Response (truncated):** Here's the violin plot of Fabp4 expression across different conditions (orig.ident) for 15,719 cells. I used the `violin_plot` tool, and please note that the condition you specified was canonicalized to 'orig.ident'.

---

### ❌ DAI-12

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 7.06s
- **Artifacts:** 1 table
- **Tool called:** True
- **Notes:** Markers of the EC (endothelial) cell type on a heatmap. Agent must resolve EC as a cell_type label and compute markers.


**Failures:**

- must_contain: 'heatmap' not in response
- requires_plot: no valid PlotResult (none produced)

**Response (truncated):** I identified marker genes specifically expressed in the endothelial cell (EC) type compared to all other cells. A total of 1849 significant markers were found (padj<0.05, logFC>0.0) using the one-vs-rest method. The top 10 markers include Ly6c1, Ramp2, and Ly6a, among others (find_markers).

---

### ❌ DAI-13

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 7.02s
- **Artifacts:** 1 table
- **Tool called:** True
- **Notes:** Top-N markers variant of DAI-12. "how about" prefix is conversational; agent must still execute.


**Failures:**

- must_contain: 'heatmap' not in response
- requires_plot: no valid PlotResult (none produced)

**Response (truncated):** I identified marker genes specifically expressed in the endothelial cell (EC) type compared to all other cells. A total of 1849 significant markers were found (padj<0.05, logFC>0.0) using the one-vs-rest method. The top 10 markers include Ly6c1, Ramp2, and Ly6a, with log2 fold changes indicating strong expression in ECs (find_markers).

---

### ❌ DAI-14

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 9.73s
- **Artifacts:** 1 table
- **Tool called:** True
- **Notes:** Per-cell-type top-3 markers. Dataset-generic request with no hard-coded genes.


**Failures:**

- must_contain: 'heatmap' not in response
- requires_plot: no valid PlotResult (none produced)

**Response (truncated):** Here are the marker genes identified across 19 cell types from your dataset, with a total of 35,360 significant markers found (padj<0.05, logFC>0.0). For example, the top markers for AT1 are Ager, Rtkn2, and Scnn1g, while B cells show Cd79a, Ms4a1, and Ighd. This analysis was performed using the `find_all_markers` tool.

---

### ✅ DAI-15

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.75s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Subset UMAP to a single cell type (EC). Agent should subset, not highlight — "only show" implies subset. Known-failing in Dai's session (produced split panels instead of subset); negatives target the split-panel narration that previously slipped through.


**Response (truncated):** Here's the UMAP scatter plot showing the distribution of 4,818 cells, specifically subsetted to the endothelial cell (EC) type. I used the `umap_plot` tool, which automatically colored the cells by their cell type.

---

### ✅ DAI-16

- **Category:** out_of_scope
- **Expected status:** FAIL
- **Actual status:** PASS
- **Duration:** 7.93s
- **Artifacts:** none
- **Tool called:** False
- **Notes:** Jitter plot is not a supported tool. Agent must decline without crashing, or redirect to violin/strip.


**Response (truncated):** This is not supported in this version. I've recorded your request for the development team to consider in the next iteration.

---

### ⚠️ DAI-17

- **Category:** data_analysis
- **Expected status:** WARN
- **Actual status:** WARN
- **Duration:** 3.29s
- **Artifacts:** none
- **Tool called:** True
- **Notes:** DEFERRED 2026-04-29 (DE migration). Set to WARN until two preconditions land: (1) find_markers (one-vs-rest DE for EC) ships in src/domain/de.py via the §4 PR; (2) multi-turn test infrastructure so turn 1 = find_markers, turn 2 = volcano follow-up. Original expectation was one-vs-rest DE for EC (via legacy LangChain path) then volcano. Failure here is expected and won't block the suite.


**Failures:**

- must_contain: 'volcano' not in response
- requires_plot: no valid PlotResult (none produced)

**Response (truncated):** The tool returned an error (output starts with 'Error:') but the response doesn't acknowledge it. The user would be misled about what happened.

---

### ✅ DAI-18

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 5.15s
- **Artifacts:** 1 table
- **Tool called:** True
- **Notes:** Single-turn DE-only prompt (volcano follow-up split off as DAI-18b when multi-turn infra lands). Tests run_de end-to-end: subset to EC celltype + pairwise DE WT vs CKO. Original 2026-04 prompt asked for DE+volcano in one turn; new spec-pipeline architecture is single-tool per turn, so the volcano half moves to a follow-up. The agent's response should end with a volcano offer per Path D (src/domain/de.py:_VOLCANO_OFFER) — that part is qualitative, not in must_contain.


**Response (truncated):** I found the differential expression between the WT and CKO conditions for the EC cell type. A total of 1,751 genes were significant (padj<0.05), with 516 genes upregulated in WT and 1,235 genes upregulated in CKO. The top upregulated genes in WT include Sntb1 (log2FC=5.62) and Gm50250 (log2FC=4.99), while in CKO, Srarp (log2FC=-27.45) and Rbp7 (log2FC=-7.45) were the most significantly downregulat

---

### ✅ DAI-19

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 5.49s
- **Artifacts:** 1 table
- **Tool called:** True
- **Notes:** Paraphrase of DAI-18 with "condition" added as a disambiguator — checks that the extra word still routes to run_de with subset to EC and groupby=condition. Volcano follow-up split off (see DAI-18 note).


**Response (truncated):** I performed pairwise differential expression analysis comparing the WT and CKO conditions within the EC cell type. The analysis revealed 1,751 significant genes (padj<0.05), with 516 genes upregulated in WT and 1,235 genes upregulated in CKO. Notably, the top upregulated genes in WT include Sntb1 and Gm50250, while in CKO, Srarp and Rbp7 are among the most significantly downregulated. This analysi

---
