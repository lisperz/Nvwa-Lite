# Bug Fix: Infinite Loop on Preprocessed Data

## Problem
When users uploaded already-preprocessed data (like `pbmc_test.h5ad`), the agent would:
1. Detect that data needs preprocessing (incorrectly)
2. Try to run preprocessing
3. Fail because data is already normalized
4. Loop back and try again infinitely

## Root Cause
The `pbmc_test.h5ad` file was already preprocessed:
- ✅ Has UMAP coordinates (`X_umap`)
- ✅ Has clustering (`louvain` column)
- ✅ Already log-normalized (values range -2.85 to 10.0, not raw counts)
- ✅ Has raw layer saved

The agent didn't check if preprocessing was already done before attempting to preprocess again.

## Solution

### 1. Added State Check in `preprocess_data` Tool
```python
# Check if data is already preprocessed
if _dataset_state is not None:
    if _dataset_state.has_umap and _dataset_state.has_clustering:
        return (
            f"Error: Data is already preprocessed!\n\n"
            f"Current state:\n{_dataset_state.summary()}\n\n"
            f"The dataset already has UMAP coordinates and clustering.\n"
            f"You can directly use visualization tools.\n"
            f"No need to preprocess again."
        )
```

### 2. Updated System Prompt
Added explicit guidance at the top of the decision tree:

**"If Data is ALREADY Preprocessed (has UMAP/clustering):"**
- DO NOT run preprocess_data again!
- Acknowledge data is ready
- Offer visualization options directly
- NEVER suggest preprocessing

### 3. Reordered Decision Tree
Changed from:
```
1. If NOT preprocessed → preprocess
2. If preprocessed → visualize
```

To:
```
1. FIRST: Check if ALREADY preprocessed → visualize directly
2. ONLY IF not preprocessed → offer to preprocess
```

## Testing

### Before Fix:
```
User: "Show me the data"
Agent: "Let me preprocess..." → ERROR → "Let me preprocess..." → ERROR → LOOP
```

### After Fix:
```
User: "Show me the data"
Agent: "Your data is already preprocessed! Would you like me to show a UMAP?"
User: "Yes"
Agent: [Shows UMAP successfully]
```

## Files Changed
- `src/agent/tools.py`: Added state check in `preprocess_data()`
- `src/agent/prompts.py`: Reordered and emphasized state checking

## Deployment
- Commit: `7d3f9bf`
- Pushed to GitHub: ✅
- Streamlit auto-deploy: ⏳ 2-5 minutes

## How to Test
1. Upload `pbmc_test.h5ad` (already preprocessed)
2. Ask: "Show me the data" or "Visualize this"
3. Agent should recognize it's preprocessed and show UMAP directly
4. No infinite loop!

## Related Files
- Already preprocessed datasets: `pbmc_test.h5ad`, any Seurat-exported files
- Raw count datasets: Still work as before, will be preprocessed correctly
