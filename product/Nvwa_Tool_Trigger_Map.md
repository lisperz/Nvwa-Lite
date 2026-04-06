# Nvwa Bio — Tool Trigger Map & Scenario Guide

*For: Yuxin (Layer 2a keyword coverage) + Chen (tool implementation scope)*

From: Yalu (CEO) · Date: April 4, 2026

---

## Purpose

This document defines every user scenario the agent must handle, what tools it should call, what keywords trigger each scenario, and what Layer 1 tests are required. It is the single source of truth for:

- **Yuxin** — use the keyword lists to verify Layer 2a test coverage. If a keyword pattern is not in the test suite, add it.
- **Chen** — use the tool names and Layer 1 test descriptions to verify every tool is implemented and tested on the Temple dataset.

*Note: Scenarios marked [tool missing?] need Chen to confirm whether the tool exists. Scenarios marked [flywheel] must save to the database from Day 1.*

---

## A — Quality Control

### A1 — General QC Overview

**Keywords:** QC · quality control · data quality · check my data · overview · how does my data look · is my data good

**Tools:** `violin_plot(nCount_RNA)` + `violin_plot(nFeature_RNA)` + `violin_plot(pct_counts_mt)`

**Behavior:** Always generate all 3 QC plots together. Never show just one unless the user specifically requests it.

**Layer 1 test:** All 3 violin tools run without error on Temple dataset. Each produces a valid image file.

---

### A2 — Dead Cells / Mitochondrial Concern

**Keywords:** dead cells · dying cells · mitochondria · mito · mt percent · cell stress · apoptosis · stressed cells

**Tools:** `violin_plot(pct_counts_mt)`

**Behavior:** Prioritize mito plot. Agent must interpret: median, IQR, and whether a filtering threshold was applied. Typical healthy cutoff: <20%.

**Layer 1 test:** `violin_plot` handles `pct_counts_mt` correctly even if column name varies (`percent.mito` vs `pct_counts_mt` vs `pct_mito`).

---

### A3 — Sequencing Depth / Low Quality Cells / Doublets

**Keywords:** sequencing depth · low quality · too few genes · doublets · nCount · nFeature · UMI · high counts

**Tools:** `violin_plot(nCount_RNA)` + `violin_plot(nFeature_RNA)` + `scatter_plot(nCount vs nFeature)` [tool missing?]

**Behavior:** For doublet suspicion: scatter plot of nCount vs nFeature is essential. High nCount + high nFeature = likely doublet. Confirm `scatter_plot` tool exists.

**Layer 1 test:** `scatter_plot` handles two continuous obs variables. Output is a valid image.

---

### A4 — Sample / Condition Balance

**Keywords:** balanced · distribution across samples · how many cells per condition · sample size · batch · imbalanced

**Tools:** `cell_count_table(groupby=condition)` + bar chart [tool missing?]

**Behavior:** Return both a numeric table and a bar chart. Agent should flag if any condition has >3x more cells than another.

**Layer 1 test:** Cell count aggregation returns consistent numbers across repeated runs. This was the state contamination bug in CEO test (45,460 vs 43,915).

---

## B — Clustering & UMAP

### B1 — Basic UMAP Colored by Cell Type or Cluster

**Keywords:** UMAP · show me the clusters · cell types · visualize cells · embedding · what does the data look like

**Tools:** `umap_plot(color_by=cell_type)` or `umap_plot(color_by=cluster)`

**Behavior:** Default to `cell_type` if annotation exists in h5ad. If user says "clusters" explicitly, use `leiden`/`seurat_clusters`. Always include legend in output.

**Layer 1 test:** `umap_plot` produces image with legend visible. Legend was broken in CEO test — verify legend is present in output file.

---

### B2 — UMAP Split by Condition + Colored by Cell Type

**Keywords:** split by condition · compare conditions · side by side · per sample · disease vs normal · split_by

**Tools:** `umap_plot(split_by=condition, color_by=cell_type)`

**Behavior:** CRITICAL: `split_by` and `color_by` are independent parameters. When user says "split by X and color by Y", both must be set simultaneously. Agent must NOT default to coloring by the split variable. This was the 3-round failure in CEO test.

**Layer 1 test:** `umap_plot` correctly handles `split_by` and `color_by` as independent parameters at the same time.

---

### B3 — UMAP Colored by Continuous Metadata

**Keywords:** color by number of genes · color by total counts · color by [metadata column] · recolor by

**Tools:** `umap_plot(color_by=obs_column)`

**Behavior:** Agent must probe available obs columns before plotting. If requested column does not exist, report what columns are available — do not substitute silently.

**Layer 1 test:** `umap_plot` accepts any obs column as `color_by` parameter. Returns error if column not found.

---

### B4 — Cell Composition per Condition

**Keywords:** composition · proportion · how many cells per cell type · cell type distribution · abundance · what fraction

**Tools:** `composition_barplot(groupby=condition, color_by=cell_type)` [tool missing?]

**Behavior:** Stacked bar chart showing proportion of each cell type per condition. Very common PI question. Confirm this tool exists — it may not be implemented yet.

**Layer 1 test:** Composition tool aggregates correctly by two variables. Output chart has readable legend.

---

### B5 — Re-run Clustering at Different Resolution

**Keywords:** resolution · re-cluster · clustering resolution · refine clusters · more clusters · fewer clusters · resolution 1.2

**Tools:** `run_clustering(resolution=X)` + `umap_plot()`

**Behavior:** Agent must confirm resolution parameter was applied. Output UMAP must reflect the new clustering, not reuse the old image.

**Layer 1 test:** `run_clustering` with different resolution values produces different cluster counts.

---

## C — Gene Expression Visualization

### C1 — Single Gene Expression Across Cell Types

**Keywords:** expression of [GENE] · where is [GENE] expressed · show me [GENE] · which cell type expresses · highest expression

**Tools:** `feature_plot(gene)` + `violin_plot(gene, groupby=cell_type)`

**Behavior:** If the user specifies a plot type (violin or feature), generate only that plot. If the user asks generically about gene expression without specifying a plot type, generate both feature plot (spatial) AND violin plot (by cell type). CRITICAL: when user asks "which cell type expresses highest" — answer MUST group by `cell_type`, never by cluster ID. This was the MKI67 bug in CEO test.

**Layer 1 test:** `violin_plot(gene, groupby=cell_type)` returns results with cell type names (not cluster IDs) as group labels.

---

### C2 — Gene Expression Split by Condition

**Keywords:** split by condition · disease vs normal · compare expression across conditions · [condition A] vs [condition B] · PA-IVS vs Control

**Tools:** `violin_plot(gene, groupby=condition)` or `feature_plot(gene, split_by=condition)`

**Behavior:** If user asks expression "across conditions", use `groupby=condition`. If they want spatial split on UMAP, use `split_by` on `feature_plot`. These are two different outputs.

**Layer 1 test:** Both `groupby` and `split_by` work independently on expression plots without interfering with each other.

---

### C3 — Gene Expression per Cell Type per Condition

**Keywords:** per cell type per condition · across all conditions and cell types · breakdown by both · split violin · grouped dot plot

**Tools:** `violin_plot(gene, groupby=cell_type, split=condition)` or `dot_plot(gene, groupby=cell_type, split=condition)`

**Behavior:** When dataset has many conditions (>4), dot plot is often cleaner than violin. Agent should suggest dot plot if conditions >4. Both x-axis (cell type) and color/split (condition) must be set correctly.

**Layer 1 test:** `violin_plot` and `dot_plot` both handle two grouping variables simultaneously without error.

---

### C4 — Gene Name Lookup and Synonym Handling

**Keywords:** (triggered when gene is not found, or user provides protein name / mouse gene name)

**Tools:** `gene_name_lookup(query)` [tool missing?]

**Behavior:** Three cases require gene name resolution: (1) User provides protein name instead of gene name (e.g. "CD3 protein" instead of CD3E). (2) User provides mouse gene name instead of human (e.g. Cd3e instead of CD3E). (3) User provides a gene that does not exist in the dataset. Agent must: look up the correct gene name, confirm with user before plotting, never silently substitute.

**Layer 1 test:** `gene_name_lookup` returns correct human gene symbol for common protein names and mouse orthologs. Returns clear "not found" message for non-existent genes.

---

### C5 — Non-existent Gene

**Keywords:** (triggered automatically when gene is not found in `adata.var` after lookup)

**Tools:** `gene_name_lookup(query)` — if still not found, return error

**Behavior:** If gene is not found after lookup: (1) Report clearly that the gene was not found. (2) Suggest similar gene names if any exist (e.g. TNSF11 → TNFSF11A). (3) Never plot an empty or zero-expression result without flagging it. Agent must not silently auto-substitute.

**Layer 1 test:** Agent returns error with suggestion when gene not found. Does not produce empty plots silently.

---

### C6 — Marker Gene Dot Plot — Full Pipeline

**Keywords:** plot top N marker genes from each cluster · dot plot marker genes · top markers per cell type · top 3 per cluster · top 5 per cell type

**Tools:** `de_all_clusters()` → `extract_top_n_per_celltype(n)` → `dot_plot(genes, groupby=cell_type)`

**Behavior:** CRITICAL — this is a 3-step pipeline, not a single tool call. Step 1: run Mode 1 DE to get ranked gene list per cell type. Step 2: extract top N genes per cell type separately (not top N globally). Step 3: merge all per-cell-type top genes into one list, then `dot_plot`. The most common error is taking top N from the global list — this was the dot plot bug in CEO test (3 global genes instead of top 3 per cell type).

**Layer 1 test:** Pipeline produces one set of top N genes per cell type. For 8 cell types with top 3, output should have up to 24 unique genes (may be less if shared). Verify gene count in dot plot x-axis.

---

### C7 — Pathway / Gene Family Expression

**Keywords:** interferon genes · HLA genes · ribosomal genes · pathway · gene set · genes related to [function]

**Tools:** `dot_plot(gene_list, groupby=cell_type)` or `heatmap(gene_list, groupby=cell_type)`

**Behavior:** Agent should recognize common gene family keywords and construct the gene list automatically (e.g. "interferon genes" → IFNG, IFNGR1, IFNAR1, ISG15...). For unknown pathways, ask user to provide the gene list.

**Layer 1 test:** `dot_plot` and `heatmap` handle a list of genes as input. Output shows all genes in list.

---

## D — Differential Expression Analysis

### DE Mode Definitions

There are four distinct DE modes. The agent must correctly identify which mode to use based on the user's request. Using the wrong mode is a silent error — the result looks valid but answers the wrong question.

| Mode | What it compares | Required inputs | Typical trigger phrases |
|---|---|---|---|
| **Mode 1** | Each cell type vs all others (one-vs-rest) | None — runs on full dataset | *find marker genes / top genes per cluster / what defines each cell type* |
| **Mode 2** | Within one cell type, condition A vs condition B | Cell type + condition A + condition B | *DE in [cell type] between disease and normal / what changes in cardiomyocytes* |
| **Mode 3a** | All cells, condition A vs condition B | Condition A + condition B | *disease vs normal overall / export full DEG table / upregulated in disease* |
| **Mode 3b** | Cell type A vs cell type B | Cell type A + cell type B | *B cells vs Monocytes / compare [cell type A] and [cell type B]* |

---

### D1 — Mode 1 — Find Marker Genes for All Cell Types

**Keywords:** marker genes · what defines each cell type · top genes per cluster · find markers · characterize clusters · what makes each cell type unique

**Tools:** `de_all_clusters()`

**Behavior:** One-vs-rest comparison. No condition filtering. Returns ranked gene list per cell type. Typically followed by `dot_plot` or `heatmap` (see C6 pipeline).

**Layer 1 test:** `de_all_clusters()` completes without error. Returns table with columns: cluster/cell_type, gene, score, log2fc, pval, pval_adj.

---

### D2 — Mode 2 — DE Within One Cell Type Between Conditions

**Keywords:** in [cell type] between · within cardiomyocytes · DE in [cell type] · what changes in [cell type] between disease and normal · differentially expressed in [cell type]

**Tools:** `de_within_celltype(cell_type, condition_A, condition_B)`

**Behavior:** Requires 3 inputs: cell type + condition A + condition B. Agent must validate all 3 exist before running. If any input not found, report error immediately — do not loop or guess. This was the complete failure in CEO test (Early cardiomyocyte, PA-IVS vs Control).

**Layer 1 test:** `de_within_celltype` correctly subsets by `cell_type` first, then by condition. Test the exact failing case: "Early cardiomyocyte" in PA-IVS vs Control on Temple dataset. Must return result or clear error — not loop.

---

### D3a — Mode 3a — DE Between Conditions Across All Cells

**Keywords:** disease vs normal overall · export full DEG table · all cells · overall difference between conditions · compare [condition A] vs [condition B] · upregulated in disease

**Tools:** `de_pairwise_conditions(condition_A, condition_B)`

**Behavior:** No cell type subset. Uses all cells. Requires condition A + condition B only. When user does not specify a cell type, default to this mode. Output is typically large — offer CSV export immediately.

**Layer 1 test:** `de_pairwise_conditions` runs on full dataset without cell type filter. Returns valid result table. Export to CSV works.

---

### D3b — Mode 3b — DE Between Two Cell Types

**Keywords:** compare [cell type A] vs [cell type B] · difference between monocytes and T cells · B cells vs Monocytes · what genes distinguish [cell type A] from [cell type B]

**Tools:** `de_pairwise_celltypes(cell_type_A, cell_type_B)`

**Behavior:** Subsets to two specified cell types, compares directly. No condition involved. Distinguish from Mode 2: Mode 2 = same cell type, two conditions. Mode 3b = two cell types, no condition filter.

**Layer 1 test:** `de_pairwise_celltypes` correctly subsets to two cell types and returns valid DE result.

---

### D4 — Volcano Plot

**Keywords:** volcano plot · significant genes · fold change · p-value · upregulated · downregulated · show significance

**Tools:** `volcano_plot(de_result)`

**Behavior:** Always requires a prior DE result as input. Cannot run standalone. If DE has not been run yet, agent should run the appropriate DE mode first, then generate volcano plot. Label top significant genes on the plot.

**Layer 1 test:** `volcano_plot` handles edge cases: all genes non-significant, single significant gene, very large dataset. Also test dataset formats that caused silent crash in pre-launch testing.

---

## E — Export & Flywheel

### E1 — Export DE Results to CSV

**Keywords:** export · download · save · CSV · table · full list · all genes · give me the data

**Tools:** `export_csv(de_result)`

**Behavior:** Export requires a prior DE result. If multiple DE analyses were run in the session, agent must ask which result to export. Output file must be non-empty and have correct columns.

**Layer 1 test:** CSV export produces valid file with columns: gene, log2fc, pval, pval_adj. File is non-empty. Filename includes the analysis type and date.

---

### E2 — Save PI Annotation Correction [flywheel]

**Keywords:** that's wrong · this should be · rename · update annotation · re-annotate · change cell type · correct this

**Tools:** `update_annotation(cluster_id, new_label)` + `save_correction_to_db()`

**Behavior:** FLYWHEEL REQUIREMENT — Day 1. When PI corrects a cell type label, the correction must be saved to the corrections table in the database. Agent confirms the update and logs it. Subsequent UMAP must reflect the corrected label.

**Layer 1 test:** Correction saves to DB correctly. Subsequent UMAP plot reflects the updated label without requiring a page reload.

---

### E3 — Save Next Experiment Plan [flywheel]

**Keywords:** next steps · what should I do next · follow-up experiment · I want to validate · we should · plan to · I'm thinking of

**Tools:** `save_research_intent_to_db()`

**Behavior:** FLYWHEEL REQUIREMENT — Day 1. When PI mentions future experiments or next steps, agent should offer to log it. This feeds the `research_intents` table. Agent should confirm what was saved.

**Layer 1 test:** Research intent saves to DB. Entry is retrievable in a subsequent session.

---

## Appendix — Known Bugs from CEO Test (April 1, 2026)

These bugs must be fixed before first customer delivery.

| Bug | Description | Scenario |
|---|---|---|
| **State contamination** | Same query returned 45,460 vs 43,915 cells on two runs. Results from one test contaminating the next. | A4 |
| **Pairwise DE loop** | `de_within_celltype` for "Early cardiomyocyte" PA-IVS vs Control looped 4 times with no result. Cell type filter was broken. | D2 |
| **Dot plot top N bug** | Asked for top 3 per cell type; returned 3 global genes (TGFBI, S100A11, IGFBP7) instead. | C6 |
| **MKI67 cluster vs cell type** | Asked which cell type expresses MKI67 highest; returned cluster IDs twice with contradictory values. Layer 1 failure (tool returned wrong column); root cause is Layer 2b (agent didn't understand user was asking by cell type, not cluster). | C1 |
| **UMAP legend removed** | Asked to add legend; agent internally removed labels but claimed success. | B1 |
| **UMAP split parameters** | `split_by` and `color_by` treated as the same parameter; required 3 correction rounds. Layer 2a symptom; root cause is Layer 2b (agent misunderstood split\_by vs color\_by as the same intent). | B2 |
