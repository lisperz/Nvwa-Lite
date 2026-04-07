# Fix Plan: Condition-Split Feature Plot + Cell-Type Subset Analysis

## Problem Statement

When a user asks: *"Plot TNNT2 expression split by condition. Is there a difference between PA-IVS and normal cardiomyocytes?"*

The agent should:
1. Generate a **condition-split feature plot** (one UMAP panel per condition, colored by TNNT2)
2. **Subset to cardiomyocytes only** for the comparison
3. Compare TNNT2 expression in cardiomyocytes between PA-IVS vs Control conditions
4. Produce a conclusion that answers the specific restricted question

The agent currently:
1. Generates a global (unsplit) feature plot
2. Generates a violin plot across all conditions (not restricted to cardiomyocytes)
3. Gives a broad conclusion about all cell types

## Root Cause Analysis

### Gap 1: `feature_plot` tool lacks `split_by` support

`plot_umap()` in `executor.py` already supports `split_by` with full grid logic (lines 102–168). But `plot_feature()` is a bare-bones wrapper (lines 265–286) that only accepts `(adata, gene)` — no `split_by` parameter.

The `feature_plot` **tool** in `tools.py` only exposes `gene` as its argument. When the agent receives "plot TNNT2 split by condition", it has no way to pass `split_by` through `feature_plot`. It would need to use `umap_plot(color_by="TNNT2", split_by="orig.ident")` instead, but the prompt doesn't guide it to do so.

### Gap 2: No cell-type subsetting capability

No tool exists that can filter `adata` to a specific cell type before analysis. When the user asks about "cardiomyocytes only", the agent has no mechanism to:
- Temporarily subset the data to cardiomyocytes
- Run violin plot / expression analysis on just that subset
- Compare conditions within that subset

### Gap 3: No prompt guidance for compound queries

The system prompt lacks guidance for queries that combine:
- A gene expression visualization request
- A split-by-condition constraint
- A cell-type subset constraint
- A comparison question

The existing intent mapping handles simple single-dimension queries well, but doesn't address compound queries like *"[gene] expression split by [condition], comparing [group A] vs [group B] in [cell type]"*.

## Repair Plan

### Change 1: Add `split_by` to `feature_plot` tool

**Files:** `src/agent/tools.py`

Add `split_by` parameter to the `feature_plot` tool definition so the agent can request condition-split feature plots directly:

```python
@tool
def feature_plot(gene: str, split_by: str = "") -> str:
    """Generate a feature plot showing gene expression on UMAP (viridis colormap).

    Args:
        gene: The gene name to visualize.
        split_by: Optional observation key to split the plot into separate panels
                  (e.g. 'orig.ident' to show one panel per condition/sample).
    """
```

**Implementation:** Reuse `plot_umap()` internally since it already has `split_by` grid logic, passing `cmap="viridis"` to get the continuous gene expression colormap. OR add `split_by` support to `plot_feature()` in `executor.py`.

**Preferred approach:** Modify `plot_feature()` in `executor.py` to accept `split_by`, reusing the same grid-panel logic from `plot_umap()`. This keeps the feature plot's viridis colormap semantics while adding split capability.

### Change 2: Add `subset_cells` tool for cell-type-restricted analysis

**Files:** `src/agent/tools.py` (new tool), `src/analysis/subset.py` (new file, optional)

Add a tool that temporarily subsets the dataset and runs a violin plot or expression comparison on just that subset:

```python
@tool
def subset_violin_plot(
    genes: str,
    subset_key: str,
    subset_value: str,
    groupby: str = "",
) -> str:
    """Generate a violin plot for a SUBSET of cells filtered by a metadata column.

    Use this when the user wants to compare gene expression WITHIN a specific
    cell type across conditions, or within a specific condition across cell types.

    Examples:
    - "Show TNNT2 in cardiomyocytes across conditions"
      → subset_key="cell_type", subset_value="Cardiomyocytes", groupby="orig.ident"
    - "Compare CD3E between PA-IVS and Control in T cells"
      → subset_key="cell_type", subset_value="T cells", groupby="orig.ident"

    Args:
        genes: Gene name(s) to plot. Comma-separated for multiple.
        subset_key: The metadata column to filter on (e.g., 'cell_type').
        subset_value: The value to filter for (e.g., 'Cardiomyocytes').
        groupby: The observation key to group the subset by (e.g., 'orig.ident').
    """
```

**Implementation:**
1. Validate `subset_key` and `subset_value` exist in `adata.obs`
2. Create subset: `adata_sub = adata[adata.obs[subset_key] == subset_value]`
3. If `groupby` is empty, auto-detect clustering key
4. Call `plot_violin(adata_sub, genes=gene_list, groupby=groupby)`
5. Annotate the result with the subset information
6. Support fuzzy matching on `subset_value` (case-insensitive, partial match) — crucial since users may say "cardiomyocytes" while the data has "Cardiomyocyte" or "CM"

### Change 3: Add `subset_feature_plot` tool for cell-type-restricted spatial view

**Files:** `src/agent/tools.py` (new tool)

For showing spatial gene expression restricted to a cell type, with optional condition split:

```python
@tool
def subset_feature_plot(
    gene: str,
    subset_key: str,
    subset_value: str,
    split_by: str = "",
) -> str:
    """Generate a feature plot (UMAP colored by gene expression) for a SUBSET of cells.

    Use this when the user wants to see gene expression on UMAP for just one
    cell type or one condition, optionally split by another variable.

    Args:
        gene: Gene name to visualize.
        subset_key: The metadata column to filter on (e.g., 'cell_type').
        subset_value: The value to filter for (e.g., 'Cardiomyocytes').
        split_by: Optional observation key to split the subset into panels.
    """
```

### Change 4: Update system prompt with compound query intent mapping

**File:** `src/agent/prompts.py`

Add a new section after the existing "GENE EXPRESSION ANALYSIS INTENT MAPPING":

```
## COMPOUND GENE EXPRESSION QUERIES (CRITICAL)

When user asks about gene expression with BOTH a visualization constraint AND a biological comparison:

### Pattern: "[gene] expression split by [condition]"
- Map "split by condition" → use `split_by=<condition_column>` parameter
- For feature plots: `feature_plot(gene="TNNT2", split_by="orig.ident")`
- For UMAP plots: `umap_plot(color_by="TNNT2", split_by="orig.ident")`
- ALWAYS split when user says "split by", "per condition", "for each sample", "across conditions"

### Pattern: "Is there a difference in [gene] between [groupA] and [groupB] in [cell type]?"
- This requires CELL-TYPE-RESTRICTED analysis
- Step 1: Identify the cell type mentioned (e.g., "cardiomyocytes")
- Step 2: Identify the comparison groups (e.g., "PA-IVS" vs "normal/Control")
- Step 3: Use `subset_violin_plot(genes="TNNT2", subset_key="cell_type", subset_value="Cardiomyocytes", groupby="orig.ident")` to compare within that cell type only
- Step 4: Optionally use `subset_feature_plot(gene="TNNT2", subset_key="cell_type", subset_value="Cardiomyocytes", split_by="orig.ident")` for spatial context
- Step 5: Provide conclusion specifically about that cell type, NOT all cells

### Pattern: Combined — "[gene] split by [condition], difference between [groupA] and [groupB] [cell type]"
Full workflow:
1. `feature_plot(gene="TNNT2", split_by="orig.ident")` — global spatial view split by condition
2. `subset_violin_plot(genes="TNNT2", subset_key="cell_type", subset_value="Cardiomyocytes", groupby="orig.ident")` — quantitative comparison restricted to cardiomyocytes
3. Conclusion must state: "In [cell_type] specifically, [gene] shows [higher/lower/similar] expression in [groupA] compared to [groupB]"

### CRITICAL: Constraint Propagation
When a user query contains multiple constraints, you MUST carry ALL of them:
- Visualization constraint ("split by condition") → applies to plot tool parameters
- Subset constraint ("in cardiomyocytes") → applies to which cells are analyzed
- Comparison constraint ("PA-IVS vs normal") → applies to what groups are compared
- Gene constraint ("TNNT2") → applies to which gene is measured

NEVER drop a constraint. If the user says "split by condition", the plot MUST be split.
If the user says "in cardiomyocytes", the analysis MUST be restricted to cardiomyocytes.
```

### Change 5: Update `plot_feature()` in `executor.py` to support `split_by`

**File:** `src/plotting/executor.py`

Add `split_by` support to `plot_feature()` by extracting the grid-panel logic from `plot_umap()` into a shared helper, or by directly adding the split logic:

```python
def plot_feature(
    adata: AnnData,
    gene: str,
    split_by: str | None = None,
) -> PlotResult:
```

**Implementation:** When `split_by` is provided:
1. Get unique groups from `adata.obs[split_by]`
2. Create subplot grid (same layout logic as `plot_umap`)
3. For each group, subset adata and call `sc.pl.umap(subset, color=gene, cmap="viridis", ax=ax, show=False)`
4. Use consistent `vmin`/`vmax` across all panels (compute globally first) so colors are comparable
5. Add a shared colorbar

**Critical detail:** `vmin`/`vmax` must be computed from the FULL dataset before splitting, so the color scale is consistent across panels. This is essential for visual comparison.

## File Change Summary

| File | Change Type | Description |
|---|---|---|
| `src/plotting/executor.py` | Modify | Add `split_by` + shared `vmin/vmax` to `plot_feature()` |
| `src/agent/tools.py` | Modify | Add `split_by` param to `feature_plot` tool; add `subset_violin_plot` and `subset_feature_plot` tools; register in `get_all_tools()` |
| `src/agent/prompts.py` | Modify | Add COMPOUND QUERY and CONSTRAINT PROPAGATION prompt sections |

**Optional new file (only if tools.py exceeds 300 lines too much):**
| `src/agent/subset_tools.py` | New | Extract `subset_violin_plot` and `subset_feature_plot` if tools.py gets too large |

## Tool Count Impact
Current: 28 tools → After: 30 tools (adding `subset_violin_plot` and `subset_feature_plot`)

## Design Decisions

### Why not just update the prompt to use existing tools?
The agent could theoretically be told "use `umap_plot(color_by='TNNT2', split_by='orig.ident')` instead of `feature_plot`" for split requests. But:
1. **Feature plot uses viridis colormap** specifically designed for continuous gene expression. `umap_plot` uses scanpy's default categorical palette when coloring by gene, which doesn't give the same visual clarity.
2. **No subsetting tool exists at all** — prompt guidance alone can't fix the cell-type restriction gap since there's no tool that filters cells.

### Why fuzzy matching on subset values?
Users will say "cardiomyocytes" but the data might have "Cardiomyocyte", "CM", "cardiomyocyte", etc. The subset tools should:
1. Try exact match first
2. Try case-insensitive match
3. Try substring match (e.g., "cardio" matches "Cardiomyocyte")
4. If no match, list available values and ask for clarification

### Why shared vmin/vmax for split feature plots?
Without consistent color scales, a low-expression panel could use the same bright yellow as a high-expression panel, making visual comparison meaningless. The reviewer's code hint uses `sc.pl.umap(obj, color="TNNT2")` per panel — we improve on this by ensuring comparable color ranges across panels.

## Expected Behavior After Fix

**User:** "Plot TNNT2 expression split by condition. Is there a difference between PA-IVS and normal cardiomyocytes?"

**Agent should:**
1. Call `feature_plot(gene="TNNT2", split_by="orig.ident")` → generates 8-panel feature plot (one per condition), viridis colormap, consistent color scale
2. Call `subset_violin_plot(genes="TNNT2", subset_key="cell_type", subset_value="Cardiomyocytes", groupby="orig.ident")` → violin plot showing TNNT2 expression in cardiomyocytes only, grouped by condition
3. Provide conclusion: "Looking specifically at cardiomyocytes, TNNT2 expression in PA-IVS conditions shows [higher/lower/similar] levels compared to Control conditions. The median expression in PA-IVS cardiomyocytes is X vs Y in Control cardiomyocytes."

## Testing Checklist

- [ ] "Plot TNNT2 expression split by condition" → feature plot with one panel per condition
- [ ] "Is there a difference between PA-IVS and normal cardiomyocytes?" → subset analysis restricted to cardiomyocytes
- [ ] "Show CD3E in T cells across conditions" → subset violin of CD3E in T cells grouped by condition
- [ ] "Compare MKI67 between disease and healthy in NK cells" → subset violin restricted to NK cells
- [ ] Split feature plot has consistent color scale across all panels
- [ ] Subset tool handles case-insensitive cell type matching
- [ ] Subset tool provides helpful error when cell type not found
- [ ] Conclusion specifically addresses the cell-type-restricted comparison
- [ ] Global (unsplit) feature_plot still works when split_by is not provided
