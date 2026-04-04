# Critical Bug Fix: Dynamic Marker Gene Detection

## Problem

The system was hardcoding immune cell markers (CD3E, MS4A1, NKG7, etc.) in the system prompt, causing two critical issues:

1. **Wrong cell type suggestions**: When users uploaded non-immune datasets (e.g., forebrain data), the agent would still reference immune markers that don't exist in the dataset
2. **Confusing user experience**: Users saw mentions of "CD3E", "MS4A1" etc. even though these genes weren't in their data

### Example Issue
- User uploaded: `/Users/zhuchen/Documents/Sub UMap of the forebrain.h5ad`
- Dataset contains: 49,637 cells, 32,031 genes (Ensembl IDs: ENSDARG...)
- Immune markers present: **NONE**
- System still mentioned: CD3E, MS4A1, CD79A, NKG7, etc.

## Solution

### 1. Dynamic Tissue Detection
The system now automatically detects tissue type by checking for marker genes:

```python
tissue_markers = {
    "immune": ["CD3E", "CD3D", "MS4A1", "CD79A", "NKG7", ...],
    "brain": ["SLC17A7", "GAD1", "GAD2", "AQP4", "MBP", ...],
    "liver": ["ALB", "AFP", "CYP3A4", ...],
    "kidney": ["NPHS1", "NPHS2", "SLC12A1", ...],
    "heart": ["MYH6", "MYH7", "TNNT2", ...],
    "lung": ["SFTPC", "SFTPB", "SCGB1A1", ...],
}
```

If ≥3 markers from a tissue type are found, the system identifies the tissue and provides relevant marker knowledge.

### 2. Ensembl ID Detection
The system detects when datasets use Ensembl gene IDs (e.g., ENSDARG00000000001) and adjusts its behavior:

```python
uses_ensembl = any(g.startswith("ENS") for g in sample_genes_list)
```

When Ensembl IDs are detected, the system tells the agent to rely on differential expression rather than assuming cell types.

### 3. Expanded Marker Knowledge
The system prompt now includes markers for multiple tissue types:

**Immune cells:**
- T cells, B cells, NK cells, Monocytes, Dendritic cells, Macrophages

**Brain cells:**
- Excitatory neurons, Inhibitory neurons, Astrocytes, Oligodendrocytes, Microglia, Dopaminergic neurons

**Other tissues:**
- Hepatocytes, Kidney, Heart, Lung

### 4. Conditional Interpretation
The agent is now instructed to:
- **Only use marker knowledge for genes that actually exist in the dataset**
- **Rely on differential expression when markers are absent**
- **Describe patterns without assuming cell types when appropriate**

## Changes Made

### File: `src/agent/prompts.py`

**Function `build_system_prompt()`:**
- Added tissue type detection logic
- Added Ensembl ID format detection
- Replaced hardcoded `available_markers` with dynamic `marker_context`
- Context now includes detected tissue type or appropriate fallback message

**System Prompt Template:**
- Replaced hardcoded immune marker list with `{marker_genes}` placeholder
- Added comprehensive marker knowledge for multiple tissue types
- Added **IMPORTANT** warning to only use markers that exist in dataset
- Updated interpretation guidelines to handle Ensembl IDs and unknown tissues

## Testing

### Before Fix:
```
User uploads forebrain.h5ad (no immune markers)
Agent: "Available marker genes: CD3E, MS4A1, NKG7..."
User: "Show me CD3E"
Agent: "CD3E not found. Did you mean..."
```

### After Fix:
```
User uploads forebrain.h5ad
Agent: "Dataset uses Ensembl gene IDs. Marker gene identification
       will require differential expression analysis."
User: "What cell types do I have?"
Agent: "Let me run differential expression to identify the marker
       genes that distinguish each cluster..."
```

## Impact

### Positive:
- ✅ Works with any tissue type (immune, brain, liver, kidney, heart, lung)
- ✅ Handles Ensembl IDs gracefully
- ✅ No more confusing references to non-existent genes
- ✅ Agent provides appropriate guidance based on dataset characteristics

### No Breaking Changes:
- Existing immune cell datasets still work perfectly
- All existing tools remain unchanged
- System prompt structure preserved

## Related to Boss's Requirements

This fix addresses the foundation needed for the boss's testing framework:
- **Level 1-2**: Agent can now handle diverse datasets
- **Level 4**: Better reasoning about cell types based on actual data

## Next Steps

See `discuss/boss-requirements-analysis.md` for:
1. Missing tools needed for Level 3-4 testing
2. Customization features (split.by, labels, colors)
3. Reasoning capabilities (find highest expression, rename clusters)
