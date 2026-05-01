# nvwa-mvp Regression Report

**Generated:** 2026-05-01 17:36:53 UTC  
**Dataset:** `/Users/yxchen/GitHub/nvwa-mvp/local/data/GSE223414_slim.h5ad`  
**Model:** `gpt-4o-mini`  
**Duration:** 133.5s  

---

## Summary

| Result | Count |
|--------|-------|
| ✅ Pass  | 16 |
| ❌ Fail  | 4 |
| ⚠️ Warn  | 0 |
| 💥 Error | 0 |
| ⏭ Skip  | 0 |
| **Total** | **20** |

---

## Per-case Results

| case_id | category | status | routed_via | duration | artifacts | failure reason |
|---------|----------|--------|------------|----------|-----------|----------------|
| ZH-01 | data_analysis | ✅ PASS | legacy | 6.9s | none | — |
| ZH-02 | data_analysis | ✅ PASS | legacy | 6.3s | 1 plot | — |
| ZH-03 | data_analysis | ✅ PASS | spec | 4.6s | 1 plot | — |
| ZH-04 | data_analysis | ✅ PASS | spec | 4.0s | 1 plot | — |
| ZH-05 | data_analysis | ✅ PASS | spec | 3.7s | 1 plot | — |
| ZH-06 | data_analysis | ❌ FAIL | spec | 1.4s | none | must_contain: 'DDX21' not in response; requires_plot: no valid PlotResult (none produced) |
| ZH-07 | data_analysis | ✅ PASS | spec | 3.1s | 1 plot | — |
| ZH-08 | data_analysis | ✅ PASS | spec | 3.5s | 1 plot | — |
| ZH-09 | data_analysis | ✅ PASS | spec | 4.0s | 1 plot | — |
| ZH-10 | data_analysis | ✅ PASS | legacy | 2.4s | none | — |
| ZH-11 | data_analysis | ✅ PASS | spec | 3.6s | 1 plot | — |
| ZH-12 | data_analysis | ✅ PASS | spec | 4.6s | 1 plot | — |
| ZH-13 | data_analysis | ✅ PASS | spec | 3.3s | 1 plot | — |
| ZH-14 | data_analysis | ✅ PASS | spec | 3.0s | 1 plot | — |
| ZH-15 | data_analysis | ✅ PASS | spec | 3.2s | 1 plot | — |
| ZH-16 | data_analysis | ✅ PASS | spec | 4.5s | 1 plot | — |
| ZH-17 | data_analysis | ✅ PASS | spec | 4.2s | 1 plot | — |
| ZH-18 | out_of_scope | ❌ FAIL | legacy | 6.4s | 3 plots | expected_artifact_count: expected 0, got 3 (3 plot(s), 0 table(s)) |
| ZH-19 | out_of_scope | ❌ FAIL | legacy | 51.4s | 3 plots | must_not_contain: forbidden string 'successfully generated' found; expected_artifact_count: expected 0, got 3 (3 plot(s), 0 table(s)) |
| ZH-20 | data_analysis | ❌ FAIL | spec | 4.0s | 1 plot | must_contain: 'D10' not in response; must_contain: 'group' not in response |

---

## Top Failure Reasons

| Category | Count |
|----------|-------|
| `out_of_scope` | 2 |
| `missing_artifact` | 1 |
| `assertion_failed` | 1 |

---

## Detailed Results

### ✅ ZH-01

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 6.89s
- **Artifacts:** none
- **Tool called:** False
- **Notes:** No gene specified — feature plot needs a gene. Agent should ask for a gene; must not silently pick one.


**Response (truncated):** To generate a feature plot, I need to know which specific gene you would like to visualize. Please provide the gene name you are interested in.

---

### ✅ ZH-02

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 6.31s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Paraphrase of DA-02 (no "the"). Should succeed identically.


**Response (truncated):** The UMAP plot has been generated, displaying the clusters identified in your dataset. This visualization provides insights into the spatial distribution of the different cell populations based on the clustering results.   Would you like to explore any specific clusters or perform further analyses?

---

### ✅ ZH-03

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 4.62s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** UMAP colored by cell_type obs column. Natural-language variant.


**Response (truncated):** Here's the UMAP scatter plot rendered for 45,460 cells, colored by cell type (you wrote 'cell type', which I corrected to 'cell_type'). This visualization helps you explore the distribution of different cell types in your dataset.

---

### ✅ ZH-04

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.99s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** UMAP split (one panel per condition level). Agent should detect a condition-like obs column.


**Response (truncated):** Here's the UMAP scatter plot showing 45,460 cells colored by cell type and split by the original identity of the samples (you specified 'condition', which I adjusted to 'orig.ident'). This visualization helps you explore the distribution of cell types across different conditions (umap_plot).

---

### ✅ ZH-05

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.70s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Feature plot of a single gene, split by condition (multi-panel).


**Response (truncated):** Here's the feature plot of DDX21 on UMAP for 45,460 cells, split by the original identity of the samples. I used 'orig.ident' (you wrote 'condition', which was adjusted for clarity) to categorize the data. This visualization highlights the expression of DDX21 across different conditions.

---

### ❌ ZH-06

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 1.39s
- **Artifacts:** none
- **Tool called:** False
- **Notes:** Ambiguous: "by cell type" on a feature plot could mean split-panel or just the single feature plot on a UMAP where cell types are visible. Agent should interpret reasonably.


**Failures:**

- must_contain: 'DDX21' not in response
- requires_plot: no valid PlotResult (none produced)

**Response (truncated):** I need to know the cell type for this analysis.

---

### ✅ ZH-07

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.14s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Single feature plot, no subsetting. "in all cells" is explicit scope confirmation.


**Response (truncated):** Here's the feature plot of DDX21 expression on UMAP for 45,460 cells, showing how this gene is distributed across the dataset (feature_plot).

---

### ✅ ZH-08

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.54s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Direct analogue of DA-03 on Zhang's dataset.


**Response (truncated):** Here's the violin plot of DDX21 expression across the different seurat_clusters for 45,460 cells. I used the `violin_plot` tool, and note that I canonicalized your input by changing 'cluster' to 'seurat_clusters'.

---

### ✅ ZH-09

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 4.02s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Same request as ZH-08 but grouped by cell_type instead of leiden.


**Response (truncated):** Here's the violin plot of DDX21 expression across different cell types for 45,460 cells. I used the `violin_plot` tool, and note that I canonicalized "cell type" to "cell_type" for consistency.

---

### ✅ ZH-10

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 2.35s
- **Artifacts:** none
- **Tool called:** False
- **Notes:** Volcano plot for the SHF cell-type subset. Single-gene phrasing is unusual for volcano (volcano shows all genes). Two behaviors acceptable: (1) agent asks for clarification (2b path), or (2) agent resolves SHF as a cell-type subset, runs DE on SHF vs rest, and plots the volcano with DDX21 highlighted. Silent substitution (e.g. scatter fallback, or false-success claim on the wrong plot) is the bug we guard against.


**Response (truncated):** This is not supported in this version. I've recorded your request for the development team to consider in the next iteration.

---

### ✅ ZH-11

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.60s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Imperative "label the cell types" — UMAP colored/annotated by cell_type. Grammatically odd but intent is clear. Silent-substitution guard: must not label by seurat_clusters (the cluster-id column) when the user asked for cell-type labels.


**Response (truncated):** Here's the UMAP scatter plot showing 45,460 cells colored by their cell type. The visualization provides insights into the distribution and clustering of different cell types within your dataset (umap_plot).

---

### ✅ ZH-12

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 4.64s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Subset to "SHF progenitor" cell type, feature plot of DDX21 split by condition. Cell-type label resolution (exact match vs partial) is the risk.


**Response (truncated):** Here's the feature plot of DDX21 expression on a UMAP for 2,873 SHF progenitor cells, split by the original identity of the samples. I used `feature_plot` for this visualization, with the cell type canonicalized from 'cell type' to 'cell_type' and the split by condition canonicalized to 'orig.ident'.

---

### ✅ ZH-13

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.28s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Violin for DDX21, subset to SHF, grouped by condition. User named the plot type — per C1 rule (see DA-03), only the violin is expected; a bonus feature plot is an unsolicited-output bug.


**Response (truncated):** Here's the violin plot of DDX21 expression across different conditions, specifically for the SHF progenitor cell type. I used the `violin_plot` tool to visualize the distribution for 2,873 cells.

---

### ✅ ZH-14

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 2.97s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Same shape as ZH-13 with a different gene — paraphrase robustness check. expected_artifact_count=1 per C1 rule; bonus plot is a scope-creep bug.


**Response (truncated):** Here's the violin plot of YAP1 expression across different conditions, specifically for the SHF progenitor cell type. I used the `violin_plot` tool to visualize the distribution for 2,873 cells.

---

### ✅ ZH-15

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.19s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** No cell-type subset; just condition grouping across all cells.


**Response (truncated):** Here's the violin plot of MKI67 expression across different conditions (orig.ident) for 45,460 cells. I used the `violin_plot` tool, and please note that the grouping was adjusted from 'condition' to 'orig.ident' for consistency.

---

### ✅ ZH-16

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 4.48s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Feature plot of MKI67 split by condition.


**Response (truncated):** Here's the feature plot of MKI67 on UMAP for 45,460 cells, split by the original identity of the samples (you specified 'condition', which was adjusted to 'orig.ident'). This visualization highlights the expression of MKI67 across different conditions.

---

### ✅ ZH-17

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 4.15s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Typo-variant of ZH-16 ("conditions" plural). Agent should treat identically.


**Response (truncated):** Here's the feature plot of MKI67 on UMAP for 45,460 cells, split by the original identity of the samples. I used `feature_plot` and adjusted the split to 'orig.ident' based on your input of 'condition'.

---

### ❌ ZH-18

- **Category:** out_of_scope
- **Expected status:** FAIL
- **Actual status:** FAIL
- **Duration:** 6.45s
- **Artifacts:** 3 plots
- **Tool called:** True
- **Notes:** Gene-gene correlation plot is not a standard tool in the current suite. If a scatter tool exists and handles this, promote to data_analysis. Verify against tool inventory before accepting as a legitimate FAIL. expected_artifact_count=0 is the strong guard; must_not_contain kept as a belt-and-suspenders phrase check.


**Failures:**

- expected_artifact_count: expected 0, got 3 (3 plot(s), 0 table(s))

**Response (truncated):** I've generated the following plots for you:  1. **Scatter Plot**: This plot shows the correlation between DDX21 and YAP1 expression levels across all cells. 2. **Violin Plot for DDX21**: This plot displays the distribution of DDX21 expression across different clusters. 3. **Violin Plot for YAP1**: This plot shows the distribution of YAP1 expression across the same clusters.  These visualizations w

---

### ❌ ZH-19

- **Category:** out_of_scope
- **Expected status:** FAIL
- **Actual status:** FAIL
- **Duration:** 51.43s
- **Artifacts:** 3 plots
- **Tool called:** True
- **Notes:** Same as ZH-18, subset to SHF. expected_artifact_count=0 is the strong guard; must_not_contain kept as a belt-and-suspenders phrase check.


**Failures:**

- must_not_contain: forbidden string 'successfully generated' found
- expected_artifact_count: expected 0, got 3 (3 plot(s), 0 table(s))

**Response (truncated):** The correlation plot between DDX21 and YAP1 in the SHF progenitor cells has been successfully generated.   ### Insights: - The scatter plot illustrates the relationship between the expression levels of DDX21 and YAP1, allowing you to assess any potential correlation between these two genes within the SHF progenitor population.  If you have any further analyses or visualizations in mind, feel free 

---

### ❌ ZH-20

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 4.04s
- **Artifacts:** 1 plot
- **Tool called:** True
- **Notes:** Subset to condition == D10, colored by "group" (another obs column). Two obs lookups in one prompt — disambiguation test. must_contain forces the agent to reference both D10 (subset) and group (color key); drops of either dimension in the response indicate silent substitution (e.g. the agent splits by condition instead of subsetting, or colors by something other than group).


**Failures:**

- must_contain: 'D10' not in response
- must_contain: 'group' not in response

**Response (truncated):** Here's the UMAP scatter plot showing 45,460 cells colored by cell type and split by the original identity of the samples (you wrote 'condition', which was corrected to 'orig.ident'). This visualization helps you explore the distribution of cell types across different conditions.

---
