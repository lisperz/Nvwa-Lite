# Yalu Scenario Gaps — Discovered During Migration

*Running list of gaps in Yalu's Layer 1 / Layer 2 / Layer 3 docs surfaced during the L4 spec-pipeline tool migration. To take to Yalu for a doc update.*

---

## Missing scenarios (top-10 customer-used tools with no Yalu coverage)

- **`inspect_metadata`** (5 calls, audit rank #5) — no Yalu Layer 1 scenario. Real prompts that route to it: "How many cells per cluster?", "Show me the distribution of cells across conditions", "Are samples balanced?" (cell-count interpretation). Single-column distribution is a discovery primitive used before plotting; Yalu Layer 1 doesn't include it.
- **`get_cluster_mapping`** (4 calls, audit rank #7) — no Yalu Layer 1 scenario. Real intent: "what does cluster N correspond to in cell-type annotation?" — a discovery utility for datasets where clusters are numeric IDs and need cross-reference to named cell types.

## Mismatch — Yalu Layer 3 tool count vs migration map collapse decision

- **UMAP**: Yalu Layer 3 defines `umap_plot` and `umap_plot_split` as two separate tools with substantially different bodies (split version has unified-axis logic). Migration map collapses into one `umap_plot(split_by=...)`. Yalu Layer 3 needs an update to reflect "one tool per intent, internal branching for variants" if that's the locked architecture.
- **Heatmap**: same pattern — `heatmap_plot` + `heatmap_plot_split` in Layer 3, collapsed to one in implementation.
- **QC violin**: same — `qc_violin_plot` + `qc_violin_plot_split` in Layer 3.
- **Feature plot**: same — `feature_plot` + `feature_plot_split`.
- **Violin plot**: same — `violin_plot` + `violin_plot_split`.

## Mismatch — Yalu Layer 2 "Workflow" type vs server-side composition

- **§1.3, §1.4 (UMAP subset)**: Yalu marks as `Type: Workflow` with explicit step list (`subset_data` then `umap_plot`). Implementation collapses subset into `umap_plot(subset_key, subset_value)` via internal `_apply_subset` (server-side composition). Yalu's "Workflow" framing implies a multi-step LLM tool chain; implementation is single-step. Need to align Yalu Layer 2 to "Workflow becomes single-tool with subset params" if that's the locked architecture.
- **§3D.3, §3D.4 (Heatmap subset)**: same pattern.
- **§3A.3, §3A.4 (Feature plot subset)**: same pattern.
- **§3B.4, §3B.5 (Violin subset)**: same pattern.
- **§3C.3, §3C.4 (Dot plot subset)**: same pattern.

## Symbolic vs concrete column names

- Yalu Layer 2/3 use symbolic constants `CELLTYPE_COL`, `CONDITION_COL` for obs columns. Real datasets use varied names (`cell_type`, `celltype`, `annotation`, `label` for cell-type; `orig.ident`, `condition`, `sample`, `batch` for condition). Yalu's docs need a "column resolution standard" section explaining how implementations map the symbolic constants to actual column names per dataset (or accept that the migration handles it via resolver).

## Underspecified scenarios

- **§3D.5 "Multiple genes · Grouped by condition · Individual cell level"** — Yalu notes "Individual cell level, not aggregated", but `sc.pl.heatmap(groupby=CONDITION_COL)` aggregates by default. Either the spec should call out the rendering mode (cell-level vs aggregated) explicitly, or §3D.5 collapses to §3D.1 with a different groupby.
- **`reference` param** in Yalu split tools (heatmap_split, feature_split, violin_split) — no user prompt examples in Layer 1 show how the user signals the reference/control condition. Without a prompt convention, the LLM can't extract `reference`. Yalu needs example prompts: "compare A and B with A as control" → `reference="A"`.

## LLM-extraction conventions not specified

- **Gene name normalization** — Yalu Layer 3 uses bare symbols (`CD3D`, `MS4A1`). Real users may pass lowercase (`cd3d`), aliases, or partial matches (`CD3` for CD3 family). No spec on whether broad-term expansion (`CD3` → `CD3D, CD3E, CD3G`) should clarify or assume.
- **Cell-type broad terms** — Yalu §1.3 example "Show the UMAP for [Cell Type A] and [Cell Type B] only" uses placeholder names. No guidance on what happens when user says "cardiomyocytes" against a dataset with `Early cardiomyocyte, Ventricular cardiomyocyte` — clarification turn or assume all-matching? (Implementation chose clarification.)

## Subset_data tool exposure

- Yalu Layer 2 lists `subset_data` under "Shared Reusable Tools" but no user-facing Layer 1 scenario routes to it directly. With server-side composition, `subset_data` may not need to be a registered LLM-facing tool at all. Yalu should clarify whether `subset_data` is an LLM-callable workflow primitive or a server-side helper only.

---

## How to use this doc

When a new gap surfaces during migration, append a bullet under the relevant section. Send to Yalu in batches when ready for a Layer 1/2/3 doc-update conversation.
