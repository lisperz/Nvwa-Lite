# Yalu Workflow Capture — Responses
*April 2026*

---

## Section 1 — Module logic

### 1A. `_apply_subset(subset_key, subset_value)`

**Typical `subset_key` values:**
- `celltype` — most common, one or multiple values
- `condition` — varies by dataset: could be `genotype`, `age`, `timepoint`, or `sample`
- The exact column name depends on the dataset — agent should probe metadata on ingestion

**Normalization of `subset_value` strings:**
- No silent auto-correction. Agent screens the target obs column, finds the closest match, and asks the user to confirm before proceeding.
- Example: *"I didn't find 'cd8t' in the celltype column, but I found 'CD8 T'. Is that what you meant?"*

**Subset on 2+ keys:**
- Sequential subset preferred over simultaneous — one key at a time for clarity.
- Scanpy supports simultaneous subset with `&` operator, but sequential is cleaner for agent logic.
- Example: `subset_data(adata, CELLTYPE_COL, 'CD8 T')` → `subset_data(adata_sub, CONDITION_COL, 'CTRL')`

**When `subset_value` doesn't match:**
- Same as normalization — agent probes the column, surfaces the closest match, and confirms with user before running.

---

### 1B. `_resolve_identifier` — 4 modes

**`gene` — ID system:**

Standard input is HGNC gene symbol (e.g. `CD8A`, `FOXP3`). Edge cases to handle:

- **Protein name → gene symbol:** Agent should convert (e.g. "CD8 protein" → `CD8A`)
- **Cross-species symbol:** If dataset is mouse but user queries human symbol, agent should convert (e.g. `CD8A` → `Cd8a`)
- **Typo / partial recall:** Agent screens all feature names, finds closest match, confirms with user
- **Ensembl ID dataset:** Convert to gene symbol on data ingestion — not required for current MVP, defer to post-MVP
- **Note:** Chen has done some work on gene name resolution in the current version — coordinate with him directly

**`cluster` — integer IDs vs biological labels:**

- Seurat datasets: cluster column typically named `seurat_clusters` or `xxx_snn_resxxx`
- Scanpy datasets: typically `leiden` or `louvain`
- Cell type annotation column: typically named `celltype`, `cell_type`, `annotation`, or `seurat_annotations`

**Key design decision:** Users never care about integer cluster IDs — they always want biological labels. On data ingestion, agent must:
1. Screen all obs columns
2. Identify the column most likely to contain cell type annotations (semantic match on column name + content)
3. Set that as the default `CELLTYPE_COL`
4. If uncertain, ask the user to confirm

Priority column names to match: `celltype`, `cell_type`, `annotation`, `cell_annotation`, `seurat_annotations`
Columns to avoid as default: `leiden`, `louvain`, `seurat_clusters`, `xxx_snn_res`

This directly addresses Dr. Dai's feedback that the agent kept plotting by cluster ID instead of cell type.

**`obs_category` — fuzzy vs exact match:**
- Same as 1A — no silent auto-correction. Agent finds closest match and confirms with user.

**`de_group` — A vs B definition:**
- Covered in detail in the scenario documentation (`Nvwa_All_Features_Layer3.md`). Three main patterns:
  1. Cell type vs cell type (e.g. CD8 T vs CD4 Naive T)
  2. Condition vs condition across all cells (e.g. CTRL vs STIM)
  3. Condition vs condition within a single cell type (e.g. CD8 T_CTRL vs CD8 T_STIM)

---

### 1C. `_run_de` — 3 flavors

**`survey` (`find_all_markers`):**
- Method: `wilcoxon`
- Only positive markers: `only_pos=True`
- Filter: `pvals_adj < 0.05`, `logfoldchanges > 0`
- Default: top 10 genes per cell type, ranked by logfoldchanges
- Full results saved to CSV

**`one_vs_rest` (`find_markers`):**
- Method: `wilcoxon`
- Single cell type vs all other cells as reference
- Filter: `pvals_adj < 0.05`, `logfoldchanges > 0`
- No pre-subset by default. If user requests subset before finding markers (e.g. "find markers for CD8 T cells in CTRL only"), this becomes a workflow: `subset_data` → `find_markers`

**`pairwise` (`run_de`):**
- Method: `wilcoxon`
- Used when comparing two specific groups: cell type vs cell type, condition vs condition, or condition vs condition within a single cell type
- Filter: `pvals_adj < 0.05`
- logfc threshold: `0.25` (Seurat default)
- Returns both upregulated and downregulated genes with summary statistics
- Statistical disclaimer always shown: results are single-cell level testing, not pseudo-bulk. For publication, consider DESeq2.

---

## Section 2 — Top 7 tools

All Section 2 content is already documented and validated in the shared folder `Nvwa_Agent_Redesign_April2026`:

- **Code + parameter overrides:** `Nvwa_Complete_ToolRegistry.md` — complete code for every tool with all non-default parameters annotated
- **Good output images + notebook cells:** `Nvwa_Validation_Notebook_pbmc_ifnb.ipynb` — all scenarios validated on pbmc_ifnb dataset with output figures

---

## Section 3 — Sanity check per analysis type

This section is difficult to answer comprehensively because quality issues in scRNA-seq analysis are highly context-dependent — there is no fixed checklist that covers all failure modes. The most reliable sanity check is biological prior knowledge: does the result make sense given what we know about the biology?

For the current MVP, we recommend focusing on getting the tools and workflows right first. Quality check heuristics can be added iteratively as we encounter real failure cases from beta users.

If specific sanity check logic is needed for a particular analysis type, the best approach is to discuss with Yalu case by case as real examples come up.

---

## Section 4 — Per scenario notebooks

All scenario notebook cells are in `Nvwa_Validation_Notebook_pbmc_ifnb.ipynb`.

Input state, parameter overrides, and reasoning are documented in:
- `Nvwa_All_Features_Layer3.md` — workflow spec per scenario
- `Nvwa_Complete_ToolRegistry.md` — tool-level parameter details

---
