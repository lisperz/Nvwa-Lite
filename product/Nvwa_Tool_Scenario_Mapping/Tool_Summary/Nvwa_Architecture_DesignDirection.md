# Nvwa Agent Rebuild — Design Direction
*Yalu — April 2026 | For team discussion*

---

## Background

Our current agent has a core problem: **for the same feature, it frequently calls the wrong tool or executes the wrong logic when the use case is slightly more complex.** For example, feature plots split by condition render each condition at a different figure size; subset logic is hardcoded into specific tools (e.g. `sub_violin`) instead of being a reusable step that any visualization can use.

The goal of this rebuild is to fix the root cause — not to keep patching individual tools.

---

## Core Design Principle

**Start from scenarios, not from tools.**

The root cause of our current problems is that we defined tools first, then asked the agent to call them. But the tool boundaries were never aligned with how analysis actually works in practice, so the agent doesn't know how to handle complex scenarios.

The correct order is:

```
Scenario definition → Workflow design → Tool boundaries → Tool implementation
```

---

## Three-Layer Architecture

### Layer 1 — Scenarios
For each analysis feature, what are the real use cases a PI would encounter? Each scenario is defined by two dimensions:
- **What the user wants to see** (which gene, which cell type, whether to group by condition)
- **What data operations are required** (subset or not, x-axis grouping, split or not)

The scenario document serves two purposes:
1. User education material (teaching users how to phrase their requests)
2. Foundation for agent training data (to be developed later)

### Layer 2 — Workflows (pre-defined)
Each scenario maps to a pre-defined workflow. The agent's job is to **identify the scenario and call the corresponding workflow** — not to decide on the fly how to combine tools.

Workflows are composed of steps that can be reused across scenarios. For example:
- `subset_data` is a universal step that can precede any visualization
- `split_by_condition` subplot logic is handled once, centrally — not re-implemented in each individual tool

### Layer 3 — Tools
Tools are the smallest execution unit called by workflows. Tools execute; they do not make decisions. Parameters are passed in by the workflow; the tool itself does not determine the scenario.

**Key design decisions:**
- Subset is not a standalone tool — it is a parameter on every visualization tool
- Split-by-condition subplot logic is encapsulated once and reused, not duplicated per tool
- Tool boundaries are defined after scenarios and workflows are finalized — not locked in advance

---

## Current Progress

We have completed **Layer 1 scenario mapping for Gene Expression Visualization**, covering:

- Feature Plot: 4 scenarios
- Violin Plot: 6 scenarios
- Dot Plot: 4 scenarios
- Heatmap: 5 scenarios

19 scenarios total, each with a clear scenario description and a standard user prompt. Output as both Word and Markdown.

**Remaining features to map (roughly low to high complexity):**
1. UMAP visualization
2. QC analysis & visualization
3. Find markers
4. Cell type composition analysis
5. Differential expression analysis

---

## On Yuxin's Tool Implementation Form

Yuxin's form (parameters, thresholds, quality standards, notebook code) is extremely valuable — it is the core input for Layer 3.

**Recommended timing: fill the form after scenarios and workflows are finalized.**

The form assumes tool boundaries are already settled. If we fill it now, we risk investing effort in an architecture that hasn't been aligned with the new scenario map yet, and will need to redo parts of it.

Suggested sequence:
1. **Now:** Complete scenario mapping for all features (Layer 1)
2. **Alignment meeting:** Confirm scenario completeness, discuss workflow boundaries
3. **After:** Use Yuxin's form to capture implementation details per tool (Layer 3)

---

## Goals for Today's Meeting

1. Confirm that the scenario-first design direction is aligned across the team
2. Discuss the workflow layer — how granular should pre-defined workflows be?
3. Confirm when Yuxin's form should be filled in
4. Decide on priorities and pace for the next phase

