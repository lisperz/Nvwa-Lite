# Quick Reference: Fixed DE Intent Handling

## Three DE Analysis Modes

### Mode 1: All Clusters (One-vs-Rest for Each)
**User Intent**: "Full analysis report for all clusters"

**Agent Workflow**:
```
1. differential_expression() → runs one-vs-rest for ALL clusters
2. get_de_results_table() → exports ALL clusters
```

**Expected Output**:
- Message: "Differential expression results table generated for ALL clusters"
- Clusters: 8 (B cells, CD14+ Monocytes, CD4 T cells, ...)
- CSV contains all clusters

---

### Mode 2: Single Cluster (One-vs-Rest for One)
**User Intent**: "DEGs for Cluster 3" or "Report for cluster 4"

**Agent Workflow**:
```
1. get_cluster_degs(cluster="3") → resolves to CD8 T cells, runs one-vs-rest
2. get_de_results_table(target_cluster="3") → exports CD8 T cells ONLY
```

**Expected Output**:
- Message: "Marker gene analysis complete (one-vs-rest for CD8 T cells)"
- CSV contains CD8 T cells ONLY

---

### Mode 3: Pairwise Comparison
**User Intent**: "Compare Cluster 1 vs Cluster 5"

**Agent Workflow**:
```
1. get_cluster_mapping() → show mapping
2. compare_groups_de(group1="1", group2="5") → pairwise comparison
3. get_pairwise_de_table() → export pairwise results
```

**Expected Output**:
- Message: "Pairwise differential expression complete. Comparison: 1 vs 5"
- CSV contains pairwise comparison results

---

## Test Queries

### ✅ Should Work Correctly Now

1. **"Please prepare the full analysis report and the marker tables for download."**
   - Mode: 1 (All Clusters)
   - Expected: ALL 8 clusters exported

2. **"Please prepare the full analysis report for cluster 4 and the marker tables for download."**
   - Mode: 2 (Single Cluster)
   - Expected: Cluster 4 (Dendritic cells) ONLY

3. **"Please prepare the full analysis report for cluster 6 and the marker tables for download."**
   - Mode: 2 (Single Cluster)
   - Expected: Cluster 6 (Megakaryocytes) ONLY

4. **"Can I get a CSV file containing all the differentially expressed genes for Cluster 3?"**
   - Mode: 2 (Single Cluster)
   - Expected: Cluster 3 (CD8 T cells) one-vs-rest, NOT pairwise

5. **"Compare gene expression between Cluster 1 and Cluster 5"**
   - Mode: 3 (Pairwise)
   - Expected: CD14+ Monocytes vs FCGR3A+ Monocytes pairwise comparison

---

## Cluster Index Mapping (pbmc3k_processed)

```
0: B cells
1: CD14+ Monocytes
2: CD4 T cells
3: CD8 T cells
4: Dendritic cells
5: FCGR3A+ Monocytes
6: Megakaryocytes
7: NK cells
```

---

## Key Fixes

1. ✅ Added `get_cluster_degs()` tool for single-cluster one-vs-rest
2. ✅ Added `target_cluster` parameter to `get_de_results_table()`
3. ✅ Added `resolve_analysis_scope()` for centralized scope handling
4. ✅ Updated prompts with clear intent disambiguation rules
5. ✅ No more default fallback to B cells

---

## Validation

Run automated tests:
```bash
cd /Users/zhuchen/Downloads/Nvwa\ Bio\ Technical\ Challange/nvwa-lite
uv run python scripts/test_intent_disambiguation.py
```

All tests should pass with correct cluster resolution.
