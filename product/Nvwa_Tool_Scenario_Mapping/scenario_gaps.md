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

## Mismatch — Yalu Layer 3 output contract vs implementation

- **§6 composition: legacy CSV+plot → spec-pipeline image-only** — Legacy `composition_analysis` returned both a count CSV table AND a stacked bar chart. Yalu Layer 3 spec for `composition_barplot` shows image-only (matplotlib `plt.savefig`); the new spec-pipeline implementation matches Yalu spec → no CSV. Post-legacy purge, users lose CSV access for composition. Yalu hasn't stated whether CSV output is intentionally dropped or just unspecified. **Proposed fix:** Yalu confirm intent — image-only is final, OR add a §6 CSV scenario (e.g. "Show the cell counts table for cell type × condition" → `composition_barplot` returns CSV instead of image). (Surfaced 2026-04-29 via §6 PR shipping decision.)

## Prompt ambiguity / disambiguator keywords

- **§3A.3 / §3A.4 vs §3B.4 / §3B.5** — §3A.1 / §3A.2 disambiguate from §3B via "on the UMAP". §3A.3 / §3A.4 drop this keyword, leaving the differentiator as "split by" (§3A.4) vs "grouped by" (§3B.4) — too subtle for reliable LLM extraction across variance. Empirically the extractor mis-routes L1-3A.4 to `violin_plot`. Adding a §3A.4 few-shot caused over-generalization that broke §3B.4 / §3B.5 routing — closed feedback loop, no fixed point via few-shot patching. **Proposed fix:** add "on the UMAP" to §3A.3 / §3A.4 prompts to restore the keyword disambiguator already used in §3A.1 / §3A.2. (Surfaced 2026-04-29 via L1-3A.4 regression failure.)

## Symbolic vs concrete column names

- Yalu Layer 2/3 use symbolic constants `CELLTYPE_COL`, `CONDITION_COL` for obs columns. Real datasets use varied names (`cell_type`, `celltype`, `annotation`, `label` for cell-type; `orig.ident`, `condition`, `sample`, `batch` for condition). Yalu's docs need a "column resolution standard" section explaining how implementations map the symbolic constants to actual column names per dataset (or accept that the migration handles it via resolver).

## Underspecified scenarios

- **§3D.5 "Multiple genes · Grouped by condition · Individual cell level"** — Yalu notes "Individual cell level, not aggregated", but `sc.pl.heatmap(groupby=CONDITION_COL)` aggregates by default. Either the spec should call out the rendering mode (cell-level vs aggregated) explicitly, or §3D.5 collapses to §3D.1 with a different groupby. Re-confirmed 2026-04-29 by L1-3D.5 regression test (set to `expected_status: WARN`); current implementation falls back to aggregated rendering.
- **`reference` param** in Yalu split tools (heatmap_split, feature_split, violin_split) — no user prompt examples in Layer 1 show how the user signals the reference/control condition. Without a prompt convention, the LLM can't extract `reference`. Yalu needs example prompts: "compare A and B with A as control" → `reference="A"`.
- **Post-`find_all_markers` visualization workflow** — Yalu Layer 1 lists "top marker genes" as a valid input alternative in §3C.1, §3C.2, §3D.1, §3D.2 plus the §4.3 workflow itself, implying an "after `find_all_markers`, plot top markers as dot plot / heatmap" flow. But Yalu Layer 2 documents no explicit follow-up offer (compare to §5's "Would you like a volcano plot?" pattern). Open: does `find_all_markers` emit a "show as dot plot / heatmap?" offer line? Does the agent silently auto-run when user requests "top markers" plot without prior state? **Proposed fix:** Yalu §4 add a "Nvwa will offer to plot the top markers as a dot plot or heatmap" sentence — symmetric to §5's volcano-offer pattern. (Surfaced 2026-04-29 via §3C/§3D/§4 cross-reference.)
- **Plot-tool subsetting beyond cell-type — sample / condition / timepoint values** — Yalu §1.3, §1.4 (and §3A.3/§3A.4, §3B.4/§3B.5, §3D.3/§3D.4) frame subsetting strictly as "for [Cell Type]" / "in [Cell Type] only". Real customer prompts use sample / condition / timepoint values as filters — e.g. ZH-20: "show me UMAP of condition D10 by group" against the Zhang dataset where `D10` is a timepoint substring of `orig.ident` values like `Control-D10`, `PA-IVS-1v-D10`. The plot tools (`umap_plot`, `violin_plot`, `heatmap_plot`) declare `subset_value.field_type="cell_type"`, biasing the LLM extractor to drop non-cell-type values silently — even though the runtime resolver supports cross-field dispatch by sibling `subset_key`'s column role (added with violin_plot 2026-04-29). Three layers compound the gap: (1) extractor drops the value, (2) tool field_type is stale, (3) `_dispatch_subset_value`'s SAMPLE_ID branch does case-insensitive exact match — partial values like `D10` don't match compound `Control-D10` even when extraction is correct. **Proposed fix:** Yalu confirm whether sample / condition / timepoint subsetting is a real customer workflow. If yes, add explicit Layer 1 scenarios (e.g. "Show the UMAP for [Condition Value] only" / "Show expression of [Gene] in [Sample/Timepoint] only") AND clarify substring/semantic matching policy for compound values (does `D10` resolve to both `Control-D10` and `PA-IVS-1v-D10` — auto-pick, ambiguity-flag, or substring expansion?). If no, ZH-20 test reframes to use cell-type subsetting only. (Surfaced 2026-04-30 via ZH-20 regression + Bucket A investigation.)

## LLM-extraction conventions not specified

- **Gene name normalization** — Yalu Layer 3 uses bare symbols (`CD3D`, `MS4A1`). Real users may pass lowercase (`cd3d`), aliases, or partial matches (`CD3` for CD3 family). No spec on whether broad-term expansion (`CD3` → `CD3D, CD3E, CD3G`) should clarify or assume.
- **Cell-type broad terms** — Yalu §1.3 example "Show the UMAP for [Cell Type A] and [Cell Type B] only" uses placeholder names. No guidance on what happens when user says "cardiomyocytes" against a dataset with `Early cardiomyocyte, Ventricular cardiomyocyte` — clarification turn or assume all-matching? (Implementation chose clarification.)
- **Multi-gene "show expression" without plot-type keyword — default unspecified** — Real customer prompts often skip the plot-type keyword (e.g. DAI-04 "can you show me the gene expression for Fabp4 and Fabp5"). The extractor currently defaults to `dot_plot` (matching §3C's "efficient multi-gene overview" framing). Yalu Layer 1 doesn't articulate this default. **Proposed fix:** Yalu §3 (intro) or §3C add: "When asking to visualize multiple genes without specifying a plot type, Nvwa defaults to a dot plot. Specify 'heatmap' or 'violin' for alternatives." (Surfaced 2026-04-29 via DAI-04.)
- **Single-gene "compare expression of [gene] between [conds]" — viz vs DE compute** — Yalu §5 covers genome-wide DE between conditions (table output). But "compare the gene expression of [single gene] between [conds]" (DAI-07 / DAI-09 / DAI-10 / DAI-11 phrasings) linguistically resembles §5 while semantically being a §3B-shape visualization request. The extractor can mis-route to `run_de` (table, not single-gene plot). **Proposed fix:** Yalu §3B add an explicit example "compare the gene expression of [gene] between [conds]" → `violin_plot`, NOT `run_de`. Or add a global note distinguishing single-gene visual comparison from DE genome-wide compute. (Surfaced 2026-04-29 via DAI-07/09/10/11.)

## Subset_data tool exposure

- Yalu Layer 2 lists `subset_data` under "Shared Reusable Tools" but no user-facing Layer 1 scenario routes to it directly. With server-side composition, `subset_data` may not need to be a registered LLM-facing tool at all. Yalu should clarify whether `subset_data` is an LLM-callable workflow primitive or a server-side helper only.

---

## How to use this doc

When a new gap surfaces during migration, append a bullet under the relevant section. Send to Yalu in batches when ready for a Layer 1/2/3 doc-update conversation.
