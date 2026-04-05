# Implementation Roadmap - Boss's Requirements

## Overview
Implementing all 4 levels of testing requirements while maintaining CLAUDE.md principles (files <300 lines, no code smells, strong typing).

## Current File Sizes
- `src/agent/tools.py`: 302 lines ⚠️ (already at limit)
- `src/plotting/executor.py`: 261 lines (safe)
- `src/agent/prompts.py`: 248 lines (safe)
- `src/ui/app.py`: 228 lines (safe)

## Implementation Strategy

### Phase 1: Level 1 Completion (Multi-gene Violin Plot)
**Goal**: Support multiple genes in violin plot

**Changes:**
- `src/plotting/executor.py`: Add `plot_violin_multi()` function
- `src/agent/tools.py`: Update `violin_plot` tool to accept comma-separated genes

**Estimated lines**: +30 lines to executor.py (291 total - safe)

### Phase 2: Level 2 Completion (Scatter Plot)
**Goal**: Add gene-gene correlation scatter plot

**Changes:**
- `src/plotting/executor.py`: Add `plot_scatter()` function
- `src/agent/tools.py`: Add `scatter_plot` tool

**Estimated lines**: +40 lines to executor.py (331 total - ⚠️ EXCEEDS LIMIT)

**Solution**: Split `executor.py` into multiple files:
```
src/plotting/
  ├── executor.py (core functions: umap, violin, feature)
  ├── comparison.py (dotplot, heatmap, scatter)
  └── advanced.py (future: split views, highlighting)
```

### Phase 3: Level 3 Customization
**Goal**: Add plot customization (labels, split views, colors)

**Changes:**
- `src/plotting/executor.py`: Add optional parameters to `plot_umap()`
  - `show_labels: bool = False`
  - `show_legend: bool = True`
  - `split_by: str | None = None`
- `src/agent/tools.py`: Update `umap_plot` tool signature

**Estimated lines**: +50 lines (need to refactor first)

### Phase 4: Level 4 Reasoning Tools
**Goal**: Add calculation and annotation capabilities

**New file needed**: `src/analysis/calculations.py`
```python
def calculate_average_expression(adata, gene, groupby)
def find_highest_expressing_cluster(adata, gene, groupby)
def rename_cluster(adata, old_name, new_name, groupby)
```

**Changes:**
- `src/agent/tools.py`: Add 3 new tools (⚠️ will exceed 300 lines)

**Solution**: Split `tools.py`:
```
src/agent/
  ├── tools.py (core: dataset_info, preprocess, DE)
  ├── plot_tools.py (all plotting tools)
  └── analysis_tools.py (reasoning: calculate, highlight, rename)
```

## Refactoring Plan (To Avoid Code Smells)

### 1. Split `src/plotting/executor.py` (261 → 3 files)
```
executor.py (basic plots: ~150 lines)
  - plot_umap()
  - plot_violin()
  - plot_feature()
  - validation helpers

comparison.py (comparison plots: ~120 lines)
  - plot_dotplot()
  - plot_heatmap()
  - plot_scatter() [NEW]

advanced.py (advanced features: ~80 lines)
  - plot_umap_split() [NEW]
  - plot_umap_highlighted() [NEW]
```

### 2. Split `src/agent/tools.py` (302 → 3 files)
```
tools.py (core tools: ~120 lines)
  - dataset_info
  - check_data_status
  - preprocess_data
  - differential_expression
  - Global state management

plot_tools.py (plotting tools: ~150 lines)
  - umap_plot
  - violin_plot
  - dotplot
  - feature_plot
  - heatmap_plot
  - volcano_plot_tool
  - scatter_plot [NEW]

analysis_tools.py (reasoning tools: ~80 lines)
  - calculate_average_expression [NEW]
  - find_highest_expression [NEW]
  - highlight_cluster [NEW]
  - rename_cluster [NEW]
```

### 3. Create `src/analysis/calculations.py` (~100 lines)
```python
- calculate_cluster_averages()
- find_top_expressing_cluster()
- get_cluster_statistics()
- rename_cluster_labels()
```

## Implementation Order

### Step 1: Refactor existing code (no new features)
1. Split `executor.py` → `executor.py` + `comparison.py`
2. Split `tools.py` → `tools.py` + `plot_tools.py`
3. Update imports in `core.py` and `app.py`
4. Test that everything still works

### Step 2: Add Level 1 & 2 features
1. Enhance `plot_violin()` for multiple genes
2. Add `plot_scatter()` to `comparison.py`
3. Add `scatter_plot` tool to `plot_tools.py`
4. Update system prompt with new capabilities

### Step 3: Add Level 3 customization
1. Add parameters to `plot_umap()`: labels, legend, split_by
2. Update `umap_plot` tool to accept these parameters
3. Update system prompt with customization examples

### Step 4: Add Level 4 reasoning
1. Create `src/analysis/calculations.py`
2. Create `src/agent/analysis_tools.py`
3. Add 4 new reasoning tools
4. Update system prompt with reasoning capabilities

## Code Quality Checks

### CLAUDE.md Compliance
- ✅ All files <300 lines
- ✅ Each folder ≤8 files
- ✅ Strong typing (no unstructured dicts)
- ✅ No code smells:
  - No rigidity (modular design)
  - No redundancy (DRY principle)
  - No circular dependencies
  - No fragility (isolated changes)
  - No obscurity (clear naming)
  - No data clumps (use dataclasses)
  - No needless complexity (YAGNI)

### Testing Strategy
After each step:
1. Run local tests with sample data
2. Push to GitHub
3. Verify Streamlit deployment
4. Test with boss's prompts

## Estimated Timeline
- Step 1 (Refactor): 1 hour
- Step 2 (Level 1-2): 30 minutes
- Step 3 (Level 3): 45 minutes
- Step 4 (Level 4): 1 hour
- **Total**: ~3-4 hours of implementation

## Risk Mitigation
- Each step is independently testable
- Can rollback to previous commit if issues arise
- Refactoring first ensures clean foundation
- All changes maintain backward compatibility
