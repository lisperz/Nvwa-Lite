# Nvwa Bio — Tool Trigger Map & Scenario Guide
*For: Yuxin (Layer 2a keyword coverage) + Chen (tool implementation scope)*
*From: Yalu (CEO) · Last updated: April 8, 2026*

---

## Purpose

This document defines every user scenario the agent must handle, what tools it should call, what keywords trigger each scenario, and what Layer 1 tests are required. It is the single source of truth for:

- **Yuxin** — use the keyword lists to verify Layer 2a test coverage. If a keyword pattern is not in the test suite, add it.
- **Chen** — use the tool names and Layer 1 test descriptions to verify every tool is implemented and tested on the Temple dataset.

> Scenarios marked **[defer V1.5]** are intentionally out of scope for V1.0. Do not implement.
> Scenarios marked **[flywheel]** must save to the database from Day 1.

---

## A — Quality Control

### A1 — General QC overview

| | |
|---|---|
| **Keywords** | QC · quality control · data quality · check my data · overview · how does my data look · is my data good |
| **Tools** | `violin_plot(nCount_RNA)` + `violin_plot(nFeature_RNA)` + `violin_plot(pct_counts_mt)` |
| **Behavior** | Always generate all 3 QC plots together. Never show just one unless the user specifically requests it. |
| **Layer 1 test** | All 3 violin tools run without error on Temple dataset. Each produces a valid image file. |

### A2 — Dead cells / mitochondrial concern

| | |
|---|---|
| **Keywords** | dead cells · dying cells · mitochondria · mito · mt percent · cell stress · apoptosis · stressed cells |
| **Tools** | `violin_plot(pct_counts_mt)` |
| **Behavior** | Prioritize mito plot. Agent must interpret: median, IQR, and whether a filtering threshold was applied. Typical healthy cutoff: <20%. |
| **Layer 1 test** | violin_plot handles pct_counts_mt correctly even if column name varies (percent.mito vs pct_counts_mt vs pct_mito). |

### A3 — Sequencing depth / low quality cells / doublets

| | |
|---|---|
| **Keywords** | sequencing depth · low quality · too few genes · doublets · nCount · nFeature · UMI · high counts |
| **Tools** | `violin_plot(nCount_RNA)` + `violin_plot(nFeature_RNA)` |
| **Behavior** | Generate both violin plots for nCount and nFeature. For doublet suspicion, agent should note that a scatter plot of nCount vs nFeature would be informative but is not available in V1.0. |
| **Layer 1 test** | Both violin tools run without error. Output is valid image. |
| **Note** | scatter_plot(nCount vs nFeature) — **[defer V1.5]** |

### A4 — Sample / condition balance

| | |
|---|---|
| **Keywords** | balanced · distribution across samples · how many cells per condition · sample size · batch · imbalanced |
| **Tools** | `cell_count_table(groupby=condition)` |
| **Behavior** | Return a numeric table showing cell counts per condition. Agent should flag if any condition has >3x more cells than another. A bar chart is not required — the table is sufficient for this use case. |
| **Layer 1 test** | Cell count aggregation returns consistent numbers across repeated runs. This was the state contamination bug in CEO test (45,460 vs 43,915). |
| **Note** | Bar chart visualization — **[defer V1.5]** |

---

## B — Clustering & UMAP

### B1 — Basic UMAP colored by cell type or cluster

| | |
|---|---|
| **Keywords** | UMAP · show me the clusters · cell types · visualize cells · embedding · what does the data look like |
| **Tools** | `umap_plot(color_by=cell_type)` or `umap_plot(color_by=cluster)` |
| **Behavior** | Default to cell_type if annotation exists in h5ad. If user says 'clusters' explicitly, use leiden/seurat_clusters. Always include legend in output. |
| **Layer 1 test** | umap_plot produces image with legend visible. Legend was broken in CEO test — verify legend is present in output file. |

### B2 — UMAP split by condition + colored by cell type

| | |
|---|---|
| **Keywords** | split by condition · compare conditions · side by side · per sample · disease vs normal · split_by |
| **Tools** | `umap_plot(split_by=condition, color_by=cell_type)` |
| **Behavior** | CRITICAL: split_by and color_by are independent parameters. When user says 'split by X and color by Y', both must be set simultaneously. Agent must NOT default to coloring by the split variable. This was the 3-round failure in CEO test. |
| **Layer 1 test** | umap_plot correctly handles split_by and color_by as independent parameters at the same time. |

### B3 — UMAP colored by continuous metadata (e.g. nGenes, total counts)

| | |
|---|---|
| **Keywords** | color by number of genes · color by total counts · color by [metadata column] · recolor by |
| **Tools** | `umap_plot(color_by=obs_column)` |
| **Behavior** | Agent must probe available obs columns before plotting. If requested column does not exist, report what columns are available — do not substitute silently. |
| **Layer 1 test** | umap_plot accepts any obs column as color_by parameter. Returns error if column not found. |

### B4 — Cell composition per condition

| | |
|---|---|
| **Keywords** | composition · proportion · how many cells per cell type · cell type distribution · abundance · what fraction |
| **Tools** | `cell_count_table(groupby=condition, by=cell_type)` |
| **Behavior** | Return a numeric table showing cell counts per cell type per condition. Agent should flag dominant or absent cell types per condition. |
| **Layer 1 test** | cell_count_table aggregates correctly by two variables. Returns consistent results across runs. |
| **Note** | Stacked composition barplot — **[defer V1.5]** |

### B5 — Re-run clustering at different resolution

| | |
|---|---|
| **Keywords** | resolution · re-cluster · clustering resolution · refine clusters · more clusters · fewer clusters · resolution 1.2 |
| **Tools** | `run_clustering(resolution=X)` + `umap_plot()` |
| **Behavior** | Agent must confirm resolution parameter was applied. Output UMAP must reflect the new clustering, not reuse the old image. |
| **Layer 1 test** | run_clustering with different resolution values produces different cluster counts. |

---

## C — Gene Expression Visualization

### C1 — Single gene expression across cell types

| | |
|---|---|
| **Keywords** | expression of [GENE] · where is [GENE] expressed · show me [GENE] · which cell type expresses · highest expression |
| **Tools** | `feature_plot(gene)` + `violin_plot(gene, groupby=cell_type)` |
| **Behavior** | Plot behavior depends on how the user asks. If the user mentions a specific plot type (e.g. "show me a violin plot" or "feature plot for GENE"), generate only that plot. If the user asks generically about gene expression without specifying a plot type (e.g. "show me expression of GENE" / "where is GENE expressed"), generate both feature plot and violin plot. CRITICAL: when user asks "which cell type expresses highest" — answer MUST group by cell_type, never by cluster ID. This was the MKI67 bug in CEO test. |
| **Layer 1 test** | violin_plot(gene, groupby=cell_type) returns results with cell type names (not cluster IDs) as group labels. |

### C2 — Gene expression split by condition

| | |
|---|---|
| **Keywords** | split by condition · disease vs normal · compare expression across conditions · [condition A] vs [condition B] · PA-IVS vs Control |
| **Tools** | `violin_plot(gene, groupby=condition)` or `feature_plot(gene, split_by=condition)` |
| **Behavior** | If user asks expression 'across conditions', use groupby=condition. If they want spatial split on UMAP, use split_by on feature_plot. These are two different outputs. |
| **Layer 1 test** | Both groupby and split_by work independently on expression plots without interfering with each other. |

### C3 — Gene expression per cell type per condition (many conditions)

| | |
|---|---|
| **Keywords** | per cell type per condition · across all conditions and cell types · breakdown by both · split violin · grouped dot plot |
| **Tools** | `violin_plot(gene, groupby=cell_type, split=condition)` or `dot_plot(gene, groupby=cell_type, split=condition)` |
| **Behavior** | When dataset has many conditions (>4), dot plot is often cleaner than violin. Agent should suggest dot plot if conditions >4. Both x-axis (cell type) and color/split (condition) must be set correctly. |
| **Layer 1 test** | violin_plot and dot_plot both handle two grouping variables simultaneously without error. |

### C4 — Gene name lookup and synonym handling

| | |
|---|---|
| **Keywords** | (triggered when gene is not found, or user provides protein name / mouse gene name) |
| **Tools** | `gene_name_lookup(query)` |
| **Behavior** | Three cases require gene name resolution: (1) User provides protein name instead of gene name (e.g. 'CD3 protein' instead of CD3E). (2) User provides mouse gene name instead of human (e.g. Cd3e instead of CD3E). (3) User provides a gene that does not exist in the dataset. Agent must: look up the correct gene name, confirm with user before plotting, never silently substitute. |
| **Layer 1 test** | gene_name_lookup returns correct human gene symbol for common protein names and mouse orthologs. Returns clear 'not found' message for non-existent genes. |

### C5 — Non-existent gene

| | |
|---|---|
| **Keywords** | (triggered automatically when gene is not found in adata.var after lookup) |
| **Tools** | `gene_name_lookup(query)` — if still not found, return error |
| **Behavior** | If gene is not found after lookup: (1) Report clearly that the gene was not found. (2) Suggest similar gene names if any exist (e.g. TNSF11 → TNFSF11A). (3) Never plot an empty or zero-expression result without flagging it. Agent must not silently auto-substitute. |
| **Layer 1 test** | Agent returns error with suggestion when gene not found. Does not produce empty plots silently. |

### C6 — Marker gene dot plot — full pipeline

| | |
|---|---|
| **Keywords** | plot top N marker genes from each cluster · dot plot marker genes · top markers per cell type · top 3 per cluster · top 5 per cell type |
| **Tools** | `de_all_clusters()` → `extract_top_n_per_celltype(n)` → `dot_plot(genes, groupby=cell_type)` |
| **Behavior** | CRITICAL — this is a 3-step pipeline, not a single tool call. Step 1: run Mode 1 DE to get ranked gene list per cell type. Step 2: extract top N genes per cell type separately (not top N globally). Step 3: merge all per-cell-type top genes into one list, then dot_plot. The most common error is taking top N from the global list — this was the dot plot bug in CEO test (3 global genes instead of top 3 per cell type). |
| **Layer 1 test** | Pipeline produces one set of top N genes per cell type. For 8 cell types with top 3, output should have up to 24 unique genes (may be less if shared). Verify gene count in dot plot x-axis. |

### C7 — Pathway / gene family expression

| | |
|---|---|
| **Keywords** | interferon genes · HLA genes · ribosomal genes · pathway · gene set · genes related to [function] |
| **Tools** | `dot_plot(gene_list, groupby=cell_type)` or `heatmap(gene_list, groupby=cell_type)` |
| **Behavior** | Agent should recognize common gene family keywords and construct the gene list automatically (e.g. 'interferon genes' → IFNG, IFNGR1, IFNAR1, ISG15...). For unknown pathways, ask user to provide the gene list. |
| **Layer 1 test** | dot_plot and heatmap handle a list of genes as input. Output shows all genes in list. |

---

## D — Differential Expression Analysis

### DE mode definitions

There are four distinct DE modes. The agent must correctly identify which mode to use based on the user's request. Using the wrong mode is a silent error — the result looks valid but answers the wrong question.

| Mode | What it compares | Required inputs | Typical trigger phrases |
|------|-----------------|-----------------|------------------------|
| **Mode 1** | Each cell type vs all others (one-vs-rest) | None — runs on full dataset | *find marker genes / top genes per cluster / what defines each cell type* |
| **Mode 2** | Within one cell type, condition A vs condition B | Cell type + condition A + condition B | *DE in [cell type] between disease and normal / what changes in cardiomyocytes* |
| **Mode 3a** | All cells, condition A vs condition B | Condition A + condition B | *disease vs normal overall / export full DEG table / upregulated in disease* |
| **Mode 3b** | Cell type A vs cell type B | Cell type A + cell type B | *B cells vs Monocytes / compare [cell type A] and [cell type B]* |

### D1 — Mode 1 — Find marker genes for all cell types

| | |
|---|---|
| **Keywords** | marker genes · what defines each cell type · top genes per cluster · find markers · characterize clusters · what makes each cell type unique |
| **Tools** | `de_all_clusters()` |
| **Behavior** | One-vs-rest comparison. No condition filtering. Returns ranked gene list per cell type. Typically followed by dot_plot or heatmap (see C6 pipeline). |
| **Layer 1 test** | de_all_clusters() completes without error. Returns table with columns: cluster/cell_type, gene, score, log2fc, pval, pval_adj. |

### D2 — Mode 2 — DE within one cell type between conditions

| | |
|---|---|
| **Keywords** | in [cell type] between · within cardiomyocytes · DE in [cell type] · what changes in [cell type] between disease and normal · differentially expressed in [cell type] |
| **Tools** | `de_within_celltype(cell_type, condition_A, condition_B)` |
| **Behavior** | Requires 3 inputs: cell type + condition A + condition B. Agent must validate all 3 exist before running. If any input not found, report error immediately — do not loop or guess. This was the complete failure in CEO test (Early cardiomyocyte, PA-IVS vs Control). |
| **Layer 1 test** | de_within_celltype correctly subsets by cell_type first, then by condition. Test the exact failing case: 'Early cardiomyocyte' in PA-IVS vs Control on Temple dataset. Must return result or clear error — not loop. |

### D3a — Mode 3a — DE between conditions across all cells

| | |
|---|---|
| **Keywords** | disease vs normal overall · export full DEG table · all cells · overall difference between conditions · compare [condition A] vs [condition B] · upregulated in disease |
| **Tools** | `de_pairwise_conditions(condition_A, condition_B)` |
| **Behavior** | No cell type subset. Uses all cells. Requires condition A + condition B only. When user does not specify a cell type, default to this mode. Output is typically large — offer CSV export immediately. |
| **Layer 1 test** | de_pairwise_conditions runs on full dataset without cell type filter. Returns valid result table. Export to CSV works. |

### D3b — Mode 3b — DE between two cell types

| | |
|---|---|
| **Keywords** | compare [cell type A] vs [cell type B] · difference between monocytes and T cells · B cells vs Monocytes · what genes distinguish [cell type A] from [cell type B] |
| **Tools** | `de_pairwise_celltypes(cell_type_A, cell_type_B)` |
| **Behavior** | Subsets to two specified cell types, compares directly. No condition involved. Distinguish from Mode 2: Mode 2 = same cell type, two conditions. Mode 3b = two cell types, no condition filter. |
| **Layer 1 test** | de_pairwise_celltypes correctly subsets to two cell types and returns valid DE result. |

### D4 — Volcano plot

| | |
|---|---|
| **Keywords** | volcano plot · significant genes · fold change · p-value · upregulated · downregulated · show significance |
| **Tools** | `volcano_plot(de_result)` |
| **Behavior** | Always requires a prior DE result as input. Cannot run standalone. If DE has not been run yet, agent should run the appropriate DE mode first, then generate volcano plot. Label top significant genes on the plot. |
| **Layer 1 test** | volcano_plot handles edge cases: all genes non-significant, single significant gene, very large dataset. Also test dataset formats that caused silent crash in pre-launch testing. |

---

## E — Export & Flywheel

### E1 — Export DE results to CSV

| | |
|---|---|
| **Keywords** | export · download · save · CSV · table · full list · all genes · give me the data |
| **Tools** | `export_csv(de_result)` |
| **Behavior** | Export requires a prior DE result. If multiple DE analyses were run in the session, agent must ask which result to export. Output file must be non-empty and have correct columns. |
| **Layer 1 test** | CSV export produces valid file with columns: gene, log2fc, pval, pval_adj. File is non-empty. Filename includes the analysis type and date. |

### E2 — PI annotation correction [internal only — defer V1.5]

| | |
|---|---|
| **Keywords** | that's wrong · this should be · rename · update annotation · re-annotate · change cell type · correct this |
| **Tools** | `save_correction_to_db()` |
| **Behavior** | INTERNAL ONLY — V1.0 assumes users upload data with accurate, finalized cell type annotations. If a PI mentions a correction, Nvwa can acknowledge it and silently log it to the corrections table for future model improvement. No user-facing confirmation required. update_annotation() is not implemented in V1.0 — UMAP labels will not update in session. |
| **Layer 1 test** | Correction saves to DB without error. No user-facing update required. |
| **Note** | User-facing annotation correction (with live UMAP update) — **[defer V1.5]** |

### E3 — Research intent logging [internal only — defer V1.5]

| | |
|---|---|
| **Keywords** | next steps · what should I do next · follow-up experiment · I want to validate · we should · plan to · I'm thinking of |
| **Tools** | `save_research_intent_to_db()` |
| **Behavior** | INTERNAL ONLY — If a PI mentions future experiments or next steps, Nvwa can silently log the intent to the research_intents table for future model improvement. No user-facing confirmation required. Users are not informed this is being logged. |
| **Layer 1 test** | Research intent saves to DB without error. |
| **Note** | User-facing research intent tracking (with cross-session memory) — **[defer V1.5]** |

---

## Appendix — Known bugs from CEO test (April 1, 2026)

These bugs must be fixed before first customer delivery.

| Bug | Description | Layer | Scenario |
|-----|-------------|-------|----------|
| **State contamination** | Same query returned 45,460 vs 43,915 cells on two runs. Results from one test contaminating the next. | Layer 1 | A4 |
| **Pairwise DE loop** | de_within_celltype for 'Early cardiomyocyte' PA-IVS vs Control looped 4 times with no result. Cell type filter was broken. | Layer 1 | D2 |
| **Dot plot top N bug** | Asked for top 3 per cell type; returned 3 global genes (TGFBI, S100A11, IGFBP7) instead. | Layer 1 / 2a | C6 |
| **MKI67 cluster vs cell type** | Asked which cell type expresses MKI67 highest; returned cluster IDs twice with contradictory values. Root cause: Layer 2a (tool returned wrong dimension) + Layer 2b (agent did not understand user was asking by cell type, not cluster). | Layer 2a + 2b | C1 |
| **UMAP legend removed** | Asked to add legend; agent internally removed labels but claimed success. | Layer 2a | B1 |
| **UMAP split parameters** | split_by and color_by treated as the same parameter; required 3 correction rounds. Root cause: Layer 2b (agent misunderstood split_by vs color_by as the same intent) — not a Layer 2a artifact failure. | Layer 2b | B2 |

---

## Deferred to V1.5

| Scenario | Tool | Reason |
|----------|------|--------|
| A3 — Doublet scatter plot | `scatter_plot(nCount vs nFeature)` | Not critical for V1.0 launch |
| A4 — Condition balance bar chart | Bar chart from cell_count_table | Number table sufficient for V1.0 |
| B4 — Cell composition barplot | `composition_barplot` | Stacked chart deferred; table sufficient for V1.0 |
| E2 — User-facing annotation correction | `update_annotation()` | V1.0 assumes finalized annotations; no live UMAP update needed |
| E3 — User-facing research intent tracking | Cross-session memory | Requires cross-session memory feature not yet built |
