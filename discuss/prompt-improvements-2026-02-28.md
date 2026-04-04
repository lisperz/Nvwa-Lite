# System Prompt Improvements - 2026-02-28

## Problem Statement

The original system prompt was too passive and technical, making the platform difficult for non-experts to use. Users needed to know specific biological terms and methods to interact effectively with the agent.

## Key Issues Addressed

1. **Passive behavior**: Agent waited for explicit commands instead of guiding users
2. **Technical jargon**: Used terms like "Leiden clustering", "HVG selection" without explanation
3. **No biological interpretation**: Generated plots without explaining what they mean
4. **Poor workflow guidance**: Didn't help users understand the typical analysis pipeline
5. **Limited intent understanding**: Couldn't interpret vague requests like "show me my data"

## Changes Made

### 1. Added "Your Role" Section
- Reframed the agent as a "knowledgeable colleague" not just a tool executor
- Emphasized proactive guidance and plain language explanations
- Set expectation for biological interpretation of all results

### 2. Added "Common User Goals" Mapping
- Maps vague user requests to concrete analysis goals
- Examples:
  - "Show me my data" → Preprocess + UMAP visualization
  - "What cell types do I have?" → Clustering + marker gene identification
  - "Is gene X expressed?" → Feature plot + violin plot + interpretation

### 3. Added "Proactive Workflow Guidance"
- Decision tree for handling different data states
- Automatic preprocessing suggestions when data isn't ready
- Automatic follow-up actions (e.g., show UMAP after preprocessing)
- Clear explanations of what each step does in plain language

### 4. Added "Biological Interpretation" Requirements
- REQUIRED biological context after every plot
- Marker gene knowledge base for common cell types:
  - T cells: CD3E, CD3D, CD8A, CD4
  - B cells: MS4A1, CD79A, CD79B
  - NK cells: NKG7, GNLY, NCAM1
  - Monocytes: CD14, LYZ, CST3
  - Dendritic cells: FCER1A, CST3
  - Megakaryocytes: PPBP, PF4
- Guidelines for interpreting different plot types

### 5. Rewrote Interaction Rules
- Focus on understanding intent first
- Be proactive about preprocessing
- Explain in plain language (avoid jargon)
- Always provide biological context
- Guide lost users through the workflow
- Be conversational ("I", "you", "let's")

### 6. Added Example Interactions
- Shows how to handle common scenarios
- Demonstrates the conversational, helpful tone
- Illustrates proactive behavior and biological interpretation

### 7. Enhanced Processing State Feedback
- Added visual indicators (⚠️ for issues, ✓ for ready)
- Provides actionable guidance based on current state
- Examples:
  - "⚠️ DATA NOT PREPROCESSED - User will need preprocessing before visualization"
  - "✓ Ready for visualization. Suggest running differential_expression to identify cell types"

## Expected Impact

### For Non-Expert Users:
- Can ask vague questions like "show me my data" and get helpful guidance
- Don't need to know technical terms like "UMAP" or "Leiden clustering"
- Receive biological interpretation of results, not just plots
- Get guided through the typical analysis workflow

### For Expert Users:
- Still have full control with specific requests
- Get richer biological context for results
- Benefit from proactive suggestions for next steps

### Example User Experience Improvements:

**Before:**
```
User: "Show me my data"
Agent: "Error: UMAP not computed. Run preprocessing first."
```

**After:**
```
User: "Show me my data"
Agent: "I'll create a UMAP visualization showing your cells. First, I need to
preprocess the data to group similar cells together - this filters low-quality
cells, normalizes expression, and identifies clusters. This takes about 30 seconds.
Proceeding now..."
[Runs preprocessing]
[Shows UMAP]
"I've grouped your 2,700 cells into 8 clusters based on gene expression similarity.
Each color represents a different cell population. The clear separation suggests
you have distinct cell types in your dataset."
```

## File Modified

- `src/agent/prompts.py`
  - `SYSTEM_PROMPT_TEMPLATE`: Completely rewritten (11 lines → 152 lines)
  - `build_system_prompt()`: Enhanced with actionable state feedback

## Testing Recommendations

Test with these user queries to verify improvements:

1. **Vague requests**:
   - "What can I do?"
   - "Show me my data"
   - "Analyze this"

2. **Non-technical language**:
   - "What cell types do I have?"
   - "Find important genes"
   - "Compare the groups"

3. **Specific genes**:
   - "Is CD3E expressed?"
   - "Show me MS4A1"

4. **Lost users**:
   - "I don't know what to do"
   - "Help"

Expected behavior: Agent should understand intent, provide guidance, execute appropriate analyses, and explain results in plain language with biological context.

## Next Steps

Consider these additional improvements:

1. **Add more marker genes**: Expand the knowledge base with tissue-specific markers
2. **Add pathway information**: Help interpret gene sets (e.g., "high inflammatory markers")
3. **Add quality metrics interpretation**: Explain what good/bad QC metrics mean
4. **Add comparison capabilities**: Help users compare datasets or conditions
5. **Add export guidance**: Suggest how to use results in publications

## Compliance with CLAUDE.md

- File remains under 300 lines (200 lines total)
- Strong typing maintained
- No code smells introduced
- Follows clean architecture principles
