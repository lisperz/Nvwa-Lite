# Tools vs Skills in Your Project

## Current Status: Tools Only (No Skills)

**Your project uses TOOLS, not SKILLS.**

### What You Have: 21 Tools

**Location:** `src/agent/tools.py`

**Your tools include:**
1. `dataset_info` - Get dataset metadata
2. `preprocess_data` - Run QC and normalization
3. `umap_plot` - Generate UMAP visualization
4. `violin_plot` - Create violin plots
5. `feature_plot` - Show gene expression on UMAP
6. `heatmap_plot` - Generate heatmap
7. `differential_expression` - Run DE for all clusters
8. `get_cluster_degs` - Run DE for one cluster (NEW)
9. `compare_groups_de` - Pairwise DE comparison
10. `get_top_markers` - Get top marker genes
11. `get_de_results_table` - Export DE results as CSV
12. `get_pairwise_de_table` - Export pairwise DE as CSV
13. `find_highest_expression` - Find cluster with highest gene expression
14. `rename_cluster` - Rename cluster labels
15. `get_cluster_mapping` - Get cluster ID to name mapping
16. `plot_scatter` - Create scatter plots
17. `summarize_qc_metrics_tool` - Summarize QC metrics
18. `summarize_obs_column` - Summarize any obs column
19. `inspect_metadata` - Inspect categorical metadata
20. `volcano_plot` - Create volcano plot
21. `query_cell_type_markers` - Query cell type markers (if available)

**How they work:**
```python
# Agent sees tools as function definitions
tools = get_all_tools()  # Returns list of LangChain tools
llm_with_tools = llm.bind_tools(tools)  # LLM can call these tools

# Agent decides which tool to call based on user query
# Example: "Show me a UMAP" → calls umap_plot()
```

## What Are Skills? (Claude Code Concept)

**Skills are NOT the same as tools.**

### Skills (Claude Code Feature)
- **Higher-level workflows** composed of multiple steps
- **Reusable patterns** for common tasks
- **User-invocable** via slash commands (e.g., `/commit`, `/review-pr`)
- **Defined in `.claude/` directory** (you don't have this)
- **Examples from Claude Code:**
  - `/commit` - Create git commit with smart message
  - `/review-pr` - Review pull request
  - `/simplify` - Refactor code for simplicity

### Tools (Your Current Implementation)
- **Single-purpose functions** that do one thing
- **Called by the agent** based on user intent
- **Defined in code** (`src/agent/tools.py`)
- **Examples from your project:**
  - `umap_plot()` - Just creates a UMAP
  - `differential_expression()` - Just runs DE analysis
  - `get_de_results_table()` - Just exports results

## Do You Need Skills?

**Short answer: No, not right now.**

### Why You Don't Need Skills

1. **Your agent already has workflows**
   - The system prompt defines multi-step workflows
   - Example: "Show markers" → `differential_expression()` → `get_top_markers()`
   - The agent orchestrates tool calls automatically

2. **Your tools are well-designed**
   - Each tool has a clear, single purpose
   - Tools compose well together
   - Agent can chain them based on user intent

3. **No repetitive manual workflows**
   - Users don't need to type complex command sequences
   - Agent handles the complexity via prompt instructions

4. **Skills add complexity**
   - Need to maintain `.claude/` directory
   - Need to define skill prompts and workflows
   - Need to document when to use each skill
   - Your current approach is simpler and works well

### When You MIGHT Want Skills

**Consider adding skills if:**

1. **Repetitive multi-step workflows emerge**
   - Example: "Full analysis report" that always does the same 10 steps
   - Could create `/full-report` skill

2. **Users want shortcuts**
   - Example: `/qc` to run all QC checks
   - Example: `/annotate` to run full cell type annotation workflow

3. **Complex domain-specific workflows**
   - Example: `/trajectory-analysis` for pseudotime analysis
   - Example: `/batch-correction` for batch effect removal

4. **Integration with external tools**
   - Example: `/export-to-cellxgene` to export to CellxGene format
   - Example: `/sync-to-s3` to backup results

## Your Current Architecture (Tools-Based)

```
User Query
    ↓
Agent (LLM + System Prompt)
    ↓
Intent Recognition
    ↓
Tool Selection
    ↓
Tool Execution (1 or more tools)
    ↓
Result Synthesis
    ↓
Response to User
```

**This works well because:**
- System prompt provides workflow guidance
- Agent is smart enough to chain tools
- Tools are atomic and composable
- No need for predefined skill workflows

## Comparison: Tools vs Skills

| Aspect | Tools (Your Current) | Skills (Claude Code) |
|--------|---------------------|---------------------|
| **Definition** | Python functions | Workflow definitions |
| **Invocation** | Agent decides | User types `/skill-name` |
| **Complexity** | Single-purpose | Multi-step workflows |
| **Location** | `src/agent/tools.py` | `.claude/skills/` |
| **Flexibility** | Agent composes | Fixed workflow |
| **Maintenance** | Code changes | Prompt changes |
| **Use Case** | Atomic operations | Common workflows |

## Example: How Your System Works Without Skills

**User:** "I want to see marker genes for all clusters"

**Agent reasoning:**
1. Parse intent: User wants marker gene analysis
2. Check system prompt: "Show markers" → `differential_expression` → `get_top_markers`
3. Call tools in sequence:
   - `differential_expression()` - Runs one-vs-rest for all clusters
   - `get_top_markers()` - Extracts top genes per cluster
4. Synthesize response with results

**No skill needed** - Agent handles the workflow automatically.

## Recommendation

**Keep your current tools-based approach.**

### Reasons:
1. ✓ It's working well
2. ✓ Simpler to maintain
3. ✓ More flexible (agent can adapt)
4. ✓ No additional infrastructure needed
5. ✓ Easier to debug

### Future consideration:
- If you find users repeatedly asking for the same complex multi-step workflows
- If you want to add user-invocable shortcuts (e.g., `/full-report`)
- If you need to standardize certain analysis pipelines

**Then** consider adding skills. But for now, your tools are sufficient.

## Summary

**What you have:** 21 well-designed tools that the agent orchestrates intelligently

**What you don't have:** Skills (predefined workflows invoked by slash commands)

**Do you need skills?** No, not for your current use case

**Your architecture is correct** - Tools are the right abstraction for your single-cell analysis agent.
