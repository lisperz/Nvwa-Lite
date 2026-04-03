# Fix: Inconsistent Response for "Distribution of Cells Across Conditions"

## Problem

When a user asks "Show me the distribution of cells across different conditions and samples. Are they balanced?", the agent gives inconsistent responses:

1. **Bad response** (Image 1): Generates 2 QC violin plots (`nCount_RNA`, `nFeature_RNA` grouped by `orig.ident`) — interprets "distribution" as QC metric distributions
2. **Good response** (Image 2): Returns cell count breakdown per condition with percentages — correctly interprets "distribution" as cell composition

The root cause is **ambiguity in the word "distribution"** with no explicit intent mapping in the system prompt. The LLM sometimes interprets it as:
- **Cell count distribution** (correct): how many cells per condition/sample
- **Expression/QC metric distribution** (wrong): violin plots of numeric values grouped by condition

## Root Cause Analysis

1. **Missing intent mapping**: The `USER INTENT MAPPING` section in `prompts.py` has no entry for "distribution of cells across conditions/samples/batches"
2. **Correct tool exists but isn't mapped**: `inspect_metadata` already returns exactly what's needed — categorical column values with cell counts and percentages
3. **LLM defaults to visualization**: Without explicit guidance, the LLM gravitates toward visualization tools (violin plots) when it sees "distribution" + "conditions"

## Dataset Context

The test dataset (`GSE223414_slim.h5ad`) has these relevant `obs` columns:
- `orig.ident` (8 categories): Control-D5, Control-D10, Control-D14, Control-D30, PA-IVS-1v-D5, PA-IVS-1v-D10, PA-IVS-1v-D14, PA-IVS-1v-D30
- `cell_type` (12 categories): Cardiac mesoderm, Early cardiomyocyte, etc.
- `nCount_RNA`, `nFeature_RNA`, `percent.mito`: QC metrics

## Proposed Solution

### Approach: Add explicit intent mapping in the system prompt

Add a new entry to the `USER INTENT MAPPING` section in `src/agent/prompts.py` that maps "distribution/composition/breakdown of cells across conditions/samples/batches" queries to the `inspect_metadata` tool.

### Changes

**File: `src/agent/prompts.py`**

Add the following to the `USER INTENT MAPPING (MVP SPECIAL)` section (after the "What's in my data?" entry):

```
- **"Cell distribution" / "How many cells per condition?" / "Are samples balanced?" / "Cell composition across samples/conditions/batches"** -> Use `inspect_metadata()` to show the distribution of cells across categorical metadata columns (conditions, samples, batches, cell types). This returns cell counts and percentages for each category. Do NOT use violin plots for this — violin plots show QC metric distributions, not cell composition.
  - IMPORTANT: When a user says "distribution of cells across conditions/samples", they want CELL COUNTS per group, not QC metric violin plots.
  - Example: "Show me the distribution of cells across conditions" -> `inspect_metadata()`
  - Example: "Are my samples balanced?" -> `inspect_metadata()`
  - Example: "How many cells per sample?" -> `inspect_metadata()`
  - Example: "Cell composition by condition" -> `inspect_metadata()`
```

### Why this approach

1. **Minimal change**: Only modifies the system prompt, no tool code changes needed
2. **Explicit disambiguation**: Directly addresses the ambiguity by telling the LLM what "distribution of cells across conditions" means
3. **Negative instruction**: Explicitly says "Do NOT use violin plots for this" to prevent the bad interpretation
4. **Follows existing pattern**: The same approach is used for other ambiguous queries (e.g., DE analysis disambiguation)
5. **Low risk**: Prompt additions are additive — they won't break existing behavior for other queries

### Risk Assessment

- **Prompt fragility**: The system prompt is documented as fragile. However, this is a pure addition to the intent mapping section (not modifying existing prompt text), which is the safest type of prompt change.
- **Over-correction**: The negative instruction ("Do NOT use violin plots") is scoped specifically to "cell distribution" queries, so it won't affect legitimate violin plot requests.

### Testing

After the change, verify with these queries:
1. "Show me the distribution of cells across different conditions" → should use `inspect_metadata`
2. "Are the samples balanced?" → should use `inspect_metadata`
3. "Show a violin plot for nCount_RNA across conditions" → should still generate violin plots (not affected)
4. "Show QC metrics" → should still use `summarize_qc_metrics_tool` (not affected)
