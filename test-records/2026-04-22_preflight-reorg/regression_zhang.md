# nvwa-mvp Regression Report

**Generated:** 2026-04-22 17:07:20 UTC  
**Dataset:** `/Users/yxchen/Documents/GitHub/nvwa-mvp/local/data/GSE223414_slim.h5ad`  
**Model:** `gpt-4o-mini`  
**Duration:** 182.1s  

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

| case_id | category | status | duration | artifacts | failure reason |
|---------|----------|--------|----------|-----------|----------------|
| ZH-01 | data_analysis | ✅ PASS | 3.4s | none | — |
| ZH-02 | data_analysis | ✅ PASS | 5.9s | 1 plot | — |
| ZH-03 | data_analysis | ✅ PASS | 9.5s | 1 plot | — |
| ZH-04 | data_analysis | ✅ PASS | 9.1s | 1 plot | — |
| ZH-05 | data_analysis | ✅ PASS | 7.0s | 1 plot | — |
| ZH-06 | data_analysis | ✅ PASS | 2.7s | 1 plot | — |
| ZH-07 | data_analysis | ❌ FAIL | 4.3s | 2 plots | expected_artifact_count: expected 1, got 2 (2 plot(s), 0 table(s)) |
| ZH-08 | data_analysis | ✅ PASS | 7.4s | 1 plot | — |
| ZH-09 | data_analysis | ✅ PASS | 6.6s | 1 plot | — |
| ZH-10 | data_analysis | ✅ PASS | 8.7s | none | — |
| ZH-11 | data_analysis | ❌ FAIL | 6.2s | 1 plot | must_not_contain: forbidden string 'seurat_clusters' found |
| ZH-12 | data_analysis | ✅ PASS | 13.3s | 1 plot | — |
| ZH-13 | data_analysis | ✅ PASS | 17.8s | 1 plot | — |
| ZH-14 | data_analysis | ✅ PASS | 15.5s | 1 plot | — |
| ZH-15 | data_analysis | ✅ PASS | 5.9s | 1 plot | — |
| ZH-16 | data_analysis | ✅ PASS | 9.1s | 1 plot | — |
| ZH-17 | data_analysis | ✅ PASS | 8.8s | 1 plot | — |
| ZH-18 | out_of_scope | ❌ FAIL | 10.1s | 3 plots | must_not_contain: forbidden string 'successfully generated' found; expected_artifact_count: expected 0, got 3 (3 plot(s), 0 table(s)) |
| ZH-19 | out_of_scope | ❌ FAIL | 17.9s | 1 plot | expected_artifact_count: expected 0, got 1 (1 plot(s), 0 table(s)) |
| ZH-20 | data_analysis | ✅ PASS | 8.5s | 1 plot | — |

---

## Top Failure Reasons

| Category | Count |
|----------|-------|
| `other` | 2 |
| `out_of_scope` | 2 |

---

## Detailed Results

### ✅ ZH-01

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 3.44s
- **Artifacts:** none
- **Notes:** No gene specified — feature plot needs a gene. Agent should ask for a gene; must not silently pick one.


**Response (truncated):** To create a feature plot, I need to know which specific gene you would like to visualize. Please provide the gene name you are interested in.

---

### ✅ ZH-02

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 5.87s
- **Artifacts:** 1 plot
- **Notes:** Paraphrase of DA-02 (no "the"). Should succeed identically.


**Response (truncated):** The UMAP plot has been generated, displaying the clusters identified in your dataset. This visualization provides insights into the spatial distribution of the different cell populations based on the clustering results.   Would you like to explore any specific clusters or perform further analyses?

---

### ✅ ZH-03

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 9.53s
- **Artifacts:** 1 plot
- **Notes:** UMAP colored by cell_type obs column. Natural-language variant.


**Response (truncated):** The UMAP plot has been generated, displaying the cell types colored by their respective clusters (seurat_clusters). This visualization helps to identify the distribution and relationships among different cell types in your dataset.   Would you like to explore any specific clusters or perform further analysis?

---

### ✅ ZH-04

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 9.05s
- **Artifacts:** 1 plot
- **Notes:** UMAP split (one panel per condition level). Agent should detect a condition-like obs column.


**Response (truncated):** The UMAP plot has been successfully generated, displaying the clusters colored by `seurat_clusters` and split by `orig.ident` into 8 panels. This visualization allows you to observe the distribution of clusters across different conditions.   If you have any further analyses or visualizations in mind, feel free to let me know!

---

### ✅ ZH-05

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 6.99s
- **Artifacts:** 1 plot
- **Notes:** Feature plot of a single gene, split by condition (multi-panel).


**Response (truncated):** The feature plot of DDX21 expression on UMAP has been successfully generated, split by condition (orig.ident). The plot displays 8 panels, each representing a different condition, using a shared color scale for direct comparison.   If you have any further questions or need additional analyses, feel free to ask!

---

### ✅ ZH-06

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 2.69s
- **Artifacts:** 1 plot
- **Notes:** Ambiguous: "by cell type" on a feature plot could mean split-panel or just the single feature plot on a UMAP where cell types are visible. Agent should interpret reasonably.


**Response (truncated):** The feature plot of DDX21 expression on UMAP has been generated successfully. This visualization highlights the spatial distribution of DDX21 across the different cell types in your dataset.   Would you like to explore any specific aspects of this plot or analyze further?

---

### ❌ ZH-07

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 4.33s
- **Artifacts:** 2 plots
- **Notes:** Single feature plot, no subsetting. "in all cells" is explicit scope confirmation.


**Failures:**

- expected_artifact_count: expected 1, got 2 (2 plot(s), 0 table(s))

**Response (truncated):** The feature plot of DDX21 expression on UMAP has been generated, providing a spatial view of its expression across all cells. Additionally, a violin plot showing DDX21 expression distribution across different seurat_clusters has also been created.  These visualizations will help you understand the expression patterns of DDX21 in your dataset. Would you like to explore any specific aspects of these

---

### ✅ ZH-08

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 7.36s
- **Artifacts:** 1 plot
- **Notes:** Direct analogue of DA-03 on Zhang's dataset.


**Response (truncated):** The violin plot of DDX21 expression across all clusters has been successfully generated. This visualization allows you to see the distribution of DDX21 expression levels within each cluster.   If you have any further questions or need additional analyses, feel free to ask!

---

### ✅ ZH-09

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 6.55s
- **Artifacts:** 1 plot
- **Notes:** Same request as ZH-08 but grouped by cell_type instead of leiden.


**Response (truncated):** The violin plot of DDX21 expression across all cell types has been successfully generated. This visualization provides insights into the distribution of DDX21 expression levels within each cell type.   Would you like to explore further or analyze another gene?

---

### ✅ ZH-10

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 8.70s
- **Artifacts:** none
- **Notes:** Volcano plot for the SHF cell-type subset. Single-gene phrasing is unusual for volcano (volcano shows all genes). Two behaviors acceptable: (1) agent asks for clarification (2b path), or (2) agent resolves SHF as a cell-type subset, runs DE on SHF vs rest, and plots the volcano with DDX21 highlighted. Silent substitution (e.g. scatter fallback, or false-success claim on the wrong plot) is the bug we guard against.


**Response (truncated):** The gene **DDX21** is found in the dataset. However, it seems that the cluster mapping does not include a specific entry for "SHF."   Could you please confirm if "SHF" corresponds to a specific numeric cluster ID or if it is a broader term that might match multiple clusters? If you have a specific numeric cluster ID for SHF, please provide that, and I can proceed with generating the volcano plot.

---

### ❌ ZH-11

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** FAIL
- **Duration:** 6.22s
- **Artifacts:** 1 plot
- **Notes:** Imperative "label the cell types" — UMAP colored/annotated by cell_type. Grammatically odd but intent is clear. Silent-substitution guard: must not label by seurat_clusters (the cluster-id column) when the user asked for cell-type labels.


**Failures:**

- must_not_contain: forbidden string 'seurat_clusters' found

**Response (truncated):** The UMAP plot has been generated, with cell types labeled according to the `seurat_clusters`. This visualization allows you to see the distribution of different cell types in your dataset. If you have any further analyses or visualizations in mind, feel free to let me know!

---

### ✅ ZH-12

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 13.28s
- **Artifacts:** 1 plot
- **Notes:** Subset to "SHF progenitor" cell type, feature plot of DDX21 split by condition. Cell-type label resolution (exact match vs partial) is the risk.


**Response (truncated):** The feature plot of DDX21 expression in SHF progenitor cells (2,873 cells) has been generated, split by condition (orig.ident).   **Note**: Six of the eight conditions have no SHF progenitor cells and are not shown in the plot: Control-D14, Control-D30, Control-D5, PA-IVS-1v-D14, PA-IVS-1v-D30, and PA-IVS-1v-D5.  If you have any further questions or need additional analyses, feel free to ask!

---

### ✅ ZH-13

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 17.82s
- **Artifacts:** 1 plot
- **Notes:** Violin for DDX21, subset to SHF, grouped by condition. User named the plot type — per C1 rule (see DA-03), only the violin is expected; a bonus feature plot is an unsolicited-output bug.


**Response (truncated):** The violin plot of DDX21 expression in SHF progenitor cells, grouped by condition, has been successfully generated.   ### Summary of DDX21 Expression in SHF Progenitor: - **Control-D10**: 2,697 cells, mean expression = 0.484, median expression = 0.000 - **PA-IVS-1v-D10**: 176 cells, mean expression = 0.416, median expression = 0.000  ### Insights: - The mean expression of DDX21 is slightly higher 

---

### ✅ ZH-14

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 15.49s
- **Artifacts:** 1 plot
- **Notes:** Same shape as ZH-13 with a different gene — paraphrase robustness check. expected_artifact_count=1 per C1 rule; bonus plot is a scope-creep bug.


**Response (truncated):** The violin plot of YAP1 expression in SHF progenitor cells, grouped by condition, has been successfully generated.   ### Summary of YAP1 Expression in SHF Progenitor: - **Control-D10**:    - Number of cells: 2,697   - Mean expression: 0.287   - Median expression: 0.000 - **PA-IVS-1v-D10**:    - Number of cells: 176   - Mean expression: 0.379   - Median expression: 0.000  ### Insights: - The mean e

---

### ✅ ZH-15

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 5.93s
- **Artifacts:** 1 plot
- **Notes:** No cell-type subset; just condition grouping across all cells.


**Response (truncated):** The violin plot of MKI67 expression across different conditions (orig.ident) has been successfully generated. This visualization allows you to see the distribution of MKI67 expression levels in all cells grouped by condition.   If you have any further questions or need additional analyses, feel free to ask!

---

### ✅ ZH-16

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 9.06s
- **Artifacts:** 1 plot
- **Notes:** Feature plot of MKI67 split by condition.


**Response (truncated):** The feature plot of MKI67 expression on UMAP has been successfully generated, split by condition (orig.ident) into 8 panels, all using a shared color scale.   This visualization allows you to observe the spatial distribution of MKI67 expression across different conditions. If you have any further analyses or visualizations in mind, feel free to let me know!

---

### ✅ ZH-17

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 8.80s
- **Artifacts:** 1 plot
- **Notes:** Typo-variant of ZH-16 ("conditions" plural). Agent should treat identically.


**Response (truncated):** The feature plot of MKI67 expression on UMAP has been successfully generated, split by conditions (orig.ident) into 8 panels, all using a shared color scale.   This visualization allows you to observe the expression levels of MKI67 across different conditions. If you have any further analyses or visualizations in mind, feel free to let me know!

---

### ❌ ZH-18

- **Category:** out_of_scope
- **Expected status:** FAIL
- **Actual status:** FAIL
- **Duration:** 10.14s
- **Artifacts:** 3 plots
- **Notes:** Gene-gene correlation plot is not a standard tool in the current suite. If a scatter tool exists and handles this, promote to data_analysis. Verify against tool inventory before accepting as a legitimate FAIL. expected_artifact_count=0 is the strong guard; must_not_contain kept as a belt-and-suspenders phrase check.


**Failures:**

- must_not_contain: forbidden string 'successfully generated' found
- expected_artifact_count: expected 0, got 3 (3 plot(s), 0 table(s))

**Response (truncated):** The correlation analysis between DDX21 and YAP1 has been successfully generated. Here are the results:  1. **Scatter Plot**: This plot shows the correlation between DDX21 and YAP1 expression levels across the dataset.     2. **Violin Plot for DDX21**: This plot illustrates the distribution of DDX21 expression across different seurat_clusters.  3. **Violin Plot for YAP1**: This plot displays the di

---

### ❌ ZH-19

- **Category:** out_of_scope
- **Expected status:** FAIL
- **Actual status:** FAIL
- **Duration:** 17.92s
- **Artifacts:** 1 plot
- **Notes:** Same as ZH-18, subset to SHF. expected_artifact_count=0 is the strong guard; must_not_contain kept as a belt-and-suspenders phrase check.


**Failures:**

- expected_artifact_count: expected 0, got 1 (1 plot(s), 0 table(s))

**Response (truncated):** The scatter plot showing the correlation between DDX21 and YAP1 has been generated successfully. This visualization allows you to observe the relationship between the expression levels of these two genes. If you have any further questions or need additional analyses, feel free to ask!

---

### ✅ ZH-20

- **Category:** data_analysis
- **Expected status:** PASS
- **Actual status:** PASS
- **Duration:** 8.46s
- **Artifacts:** 1 plot
- **Notes:** Subset to condition == D10, colored by "group" (another obs column). Two obs lookups in one prompt — disambiguation test. must_contain forces the agent to reference both D10 (subset) and group (color key); drops of either dimension in the response indicate silent substitution (e.g. the agent splits by condition instead of subsetting, or colors by something other than group).


**Response (truncated):** The UMAP plot has been generated, displaying the conditions split by "D10" across different groups. Each panel represents a distinct condition, allowing for a clear visual comparison of the cell distributions.   If you have any further analyses or visualizations in mind, feel free to let me know!

---
