"""LLM-based spec extractor — converts a user prompt into a Spec via OpenAI JSON output.

Inline prompt (coupled to the Spec + tool catalog). prompts.py owns the
conversational/system prompt; this file owns the spec-extraction prompt only.

Provider portability: the only OpenAI-specific code is `_llm_call_json` below.
Schema build, prompt rendering, and Spec construction are provider-agnostic.
Switching providers means swapping the body of that one function.
"""

from __future__ import annotations

import json
import logging
import os

from openai import OpenAI

# Importing src.domain triggers @register side-effects from family modules so REGISTRY
# is populated before extractor reads it. Safe to keep even if app already imported it.
import src.domain  # noqa: F401

from src.core.spec import Spec
from src.core.registry import REGISTRY, get_tool_names

logger = logging.getLogger(__name__)


_DEFAULT_MODEL = "gpt-4o-mini"


_EXTRACTOR_PROMPT_TMPL = """You convert a user query about a single-cell RNA-seq dataset into a structured Spec.

Output a JSON object with EXACTLY these three keys:
- "scenario_id": a short string label for the task type. Use "plot", "compare", "explore", "qc", or "unknown".
- "tool_name": one of the tool names listed below, OR the literal string "none" if no listed tool fits the user's intent.
- "pre_canonical_params": an object mapping the chosen tool's parameter names to the user's values. Use the user's words verbatim (e.g. "gapdh" not "GAPDH"); a downstream resolver canonicalizes them. Use {{}} when tool_name is "none".

Available tools:
{tool_catalog}

Phrasing → param mapping (applies to any plot tool with subset_key / subset_value / groupby / split_by):
- "across <X>" / "by <X>" / "grouped by <X>" where X names a category dimension (cell types, conditions, clusters) → groupby = the matching obs column.
- "split by <X>" (the LITERAL word "split" must appear) → split_by = the matching obs column. Without "split", "by <X>" is groupby, NOT split_by.
- "in <V>" / "just in <V>" / "restricted to <V>" / "<V> only" where V is a SPECIFIC NAMED value (a cell-type name like "SHF" or "T cell"; a condition label like "WT" or "D10") → subset_key + subset_value: pick the column whose values include V (cell-type column for cell-type names; condition column for condition labels). Put the verbatim V into subset_value as a list.
- A short token like "SHF" or "EC" inside "in <X>" is almost always a cell-type abbreviation, NOT a gene — map to subset_value, not gene.

Multi-gene plot-tool selection (Yalu §3C / §3D):
- Multiple gene names + "show expression" / "show gene expression for X and Y" / "show X, Y, Z across cell types" — when no specific plot type is named → dot_plot (Yalu §3C is the efficient multi-gene overview default).
- Multiple gene names + "heatmap" / "show a heatmap" / clustering / dendrogram framing → heatmap_plot.
- Multiple gene names + "dot plot" / "dotplot" → dot_plot.
- Single gene "show expression" → feature_plot or violin_plot per single-gene phrasing rules above.

Examples:
  USER: "Show CD3D across cell types"
  OUT:  {{"scenario_id":"plot","tool_name":"violin_plot","pre_canonical_params":{{"gene":"CD3D","groupby":"cell type"}}}}

  USER: "Show MKI67 in SHF by condition"
  OUT:  {{"scenario_id":"plot","tool_name":"violin_plot","pre_canonical_params":{{"gene":"MKI67","subset_key":"cell type","subset_value":["SHF"],"groupby":"condition"}}}}
  (Note: "by condition" without "split" → groupby, not split_by.)

  USER: "Show CD3D in WT only across cell types"
  OUT:  {{"scenario_id":"plot","tool_name":"violin_plot","pre_canonical_params":{{"gene":"CD3D","subset_key":"condition","subset_value":["WT"],"groupby":"cell type"}}}}

  USER: "Show CD3D across cell types split by condition"
  OUT:  {{"scenario_id":"plot","tool_name":"violin_plot","pre_canonical_params":{{"gene":"CD3D","groupby":"cell type","split_by":"condition"}}}}
  (Note: "split by condition" — the word "split" present → split_by.)

  USER: "Show UMAP for cardiomyocytes only"
  OUT:  {{"scenario_id":"plot","tool_name":"umap_plot","pre_canonical_params":{{"subset_key":"cell type","subset_value":["cardiomyocyte"]}}}}

  USER: "Find DE genes between WT and CKO"
  OUT:  {{"scenario_id":"compare","tool_name":"run_de","pre_canonical_params":{{"groupby":"condition","group1":"WT","group2":"CKO"}}}}

  USER: "DE between T cell and B cell"
  OUT:  {{"scenario_id":"compare","tool_name":"run_de","pre_canonical_params":{{"groupby":"cell type","group1":"T cell","group2":"B cell"}}}}

  USER: "Run DE between WT and CKO in EC cells"
  OUT:  {{"scenario_id":"compare","tool_name":"run_de","pre_canonical_params":{{"groupby":"condition","group1":"WT","group2":"CKO","subset_key":"cell type","subset_value":["EC"]}}}}
  (Note: "in EC cells" = subset to EC celltype before pairwise DE; groupby stays the condition column.)

  USER: "generate the different gene expression of EC between WT and CKO"
  OUT:  {{"scenario_id":"compare","tool_name":"run_de","pre_canonical_params":{{"groupby":"condition","group1":"WT","group2":"CKO","subset_key":"cell type","subset_value":["EC"]}}}}
  (Note: "different gene expression of [celltype] between [cond_a] and [cond_b]" — celltype is the SUBSET focus, conditions are the comparison axis. Equivalent shape to "DE between conds in [celltype] cells".)

  USER: "can you compare the gene expression of Fabp4 between WT and CKO"
  OUT:  {{"scenario_id":"plot","tool_name":"violin_plot","pre_canonical_params":{{"gene":"Fabp4","groupby":"condition","subset_key":"condition","subset_value":["WT","CKO"]}}}}
  (Note: "compare gene expression of [gene] between/across [conds]" — SINGLE-gene visual comparison → violin_plot. NOT run_de — run_de finds DE genes globally; here the user names a specific gene to visualize. NOT find_markers either — gene is named, not a celltype.)

  USER: "Find marker genes for B cells"
  OUT:  {{"scenario_id":"explore","tool_name":"find_markers","pre_canonical_params":{{"groupby":"cell type","celltype":["B cell"]}}}}
  (Note: "marker genes for [single cell type]" → find_markers, NOT run_de. find_markers is one-vs-rest for the named cell type.)

  USER: "Find marker genes for all cell types"
  OUT:  {{"scenario_id":"explore","tool_name":"find_all_markers","pre_canonical_params":{{"groupby":"cell type"}}}}
  (Note: "for all cell types" / "for every cluster" → find_all_markers, the multi-group survey. No celltype param.)

  USER: "What genes define EC?"
  OUT:  {{"scenario_id":"explore","tool_name":"find_markers","pre_canonical_params":{{"groupby":"cell type","celltype":["EC"]}}}}
  (Note: "what genes define [X]" / "what makes [X] unique" / "markers for [X]" — same pattern as find_markers.)

  USER: "can you show me the gene expression for Fabp4 and Fabp5"
  OUT:  {{"scenario_id":"plot","tool_name":"dot_plot","pre_canonical_params":{{"genes":["Fabp4","Fabp5"]}}}}
  (Note: multi-gene "show expression" without specific plot type → dot_plot per Yalu §3C default. NOT feature/violin per gene; NOT heatmap unless user said "heatmap".)

  USER: "Show a dot plot of CD3D, MS4A1, NKG7 across cell types"
  OUT:  {{"scenario_id":"plot","tool_name":"dot_plot","pre_canonical_params":{{"genes":["CD3D","MS4A1","NKG7"],"groupby":"cell type"}}}}

  USER: "Show the cell number for each cell type across conditions"
  OUT:  {{"scenario_id":"plot","tool_name":"composition_barplot","pre_canonical_params":{{"celltype_col":"cell type","condition_col":"condition","mode":"count"}}}}
  (Note: "cell number" / "cell count" / "how many cells per X" → mode='count'. composition_barplot has TWO obs columns: celltype_col + condition_col, NOT groupby/split_by.)

  USER: "Show the proportion of each cell type across conditions"
  OUT:  {{"scenario_id":"plot","tool_name":"composition_barplot","pre_canonical_params":{{"celltype_col":"cell type","condition_col":"condition","mode":"proportion_by_celltype"}}}}
  (Note: "proportion of each cell type" → mode='proportion_by_celltype' (within-celltype normalization).)

  USER: "Show the cell type composition for each condition"
  OUT:  {{"scenario_id":"plot","tool_name":"composition_barplot","pre_canonical_params":{{"celltype_col":"cell type","condition_col":"condition","mode":"proportion_by_condition"}}}}
  (Note: "composition for each condition" / "composition per condition" → mode='proportion_by_condition' (stacked bars per condition).)

  USER: "Show the cell number for T cell and B cell across conditions"
  OUT:  {{"scenario_id":"plot","tool_name":"composition_barplot","pre_canonical_params":{{"celltype_col":"cell type","condition_col":"condition","mode":"count","celltypes":["T cell","B cell"]}}}}
  (Note: Yalu §6B subset variant — celltypes list filters before plotting.)

  USER: "Show me the QC metrics for all cells"
  OUT:  {{"scenario_id":"qc","tool_name":"qc_violin_plot","pre_canonical_params":{{}}}}
  (Note: Yalu §2A.1 — multi-panel QC metric violins. Auto-detects metrics; defaults groupby to cell-type column. Use qc_violin_plot — NOT violin_plot — whenever the user names QC quantities (`nCount_RNA` / `nFeature_RNA` / `pct_counts_mt` / "QC metrics" / "quality") instead of a gene.)

  USER: "Show the QC metrics across cell types split by condition"
  OUT:  {{"scenario_id":"qc","tool_name":"qc_violin_plot","pre_canonical_params":{{"groupby":"cell type","split_by":"condition"}}}}
  (Note: Yalu §2A.2 — same QC violin, gridded per condition. "split by condition" → split_by.)

  USER: "yes" / "yes please" / "go ahead" (after the assistant offered a volcano plot in the prior turn)
  OUT:  {{"scenario_id":"plot","tool_name":"volcano_plot","pre_canonical_params":{{}}}}
  (volcano_plot with no params plots the most recent run_de result from adata.uns.)

When to use tool_name="none":
- The user's intent cannot be fulfilled by any tool listed above.
- Examples: trajectory/velocity analysis, dataset integration, plot types not listed, unrelated visualization requests.
- Do NOT force-pick the closest-named tool. Returning "none" routes the request to the legacy fallback path; force-picking would mis-dispatch and confuse the user.

Do not invent parameter values not implied by the user query. If a required parameter is missing from the user query, leave it out of pre_canonical_params — the validator will surface it as needing input.

Output ONLY the JSON object, no commentary."""


def extract(
    user_prompt: str,
    chat_history: list[tuple[str, str]] | None = None,
    model: str = _DEFAULT_MODEL,
) -> Spec:
    """Extract a Spec from the user prompt via structured-output LLM call.

    The extractor populates spec.pre_canonical_params with the user's verbatim
    values; the resolver downstream builds spec.params (additive pattern).
    """
    system = _EXTRACTOR_PROMPT_TMPL.format(tool_catalog=_format_tool_catalog())
    response_dict = _llm_call_json(
        system=system,
        user=user_prompt,
        history=chat_history or [],
        model=model,
    )
    if "tool_name" not in response_dict:
        raise RuntimeError("Extractor: response missing 'tool_name' key.")
    return Spec(
        scenario_id=str(response_dict.get("scenario_id", "unknown")),
        tool_name=str(response_dict["tool_name"]),
        pre_canonical_params=dict(response_dict.get("pre_canonical_params", {})),
    )


def _format_tool_catalog() -> str:
    """Render registered tools + their LLM-facing params for the extractor prompt."""
    if not REGISTRY:
        return "(no tools registered)"
    lines: list[str] = []
    for entry in REGISTRY.values():
        lines.append(f"- {entry.name}: {entry.description}")
        for p in entry.params:
            req = "required" if p.required else f"default={p.default!r}"
            enum = f" [one of {p.enum}]" if p.enum else ""
            desc = f" — {p.description}" if p.description else ""
            lines.append(f"    {p.name}: {p.type}, {req}{enum}{desc}")
    return "\n".join(lines)


def _llm_call_json(
    system: str,
    user: str,
    history: list[tuple[str, str]],
    model: str,
) -> dict:
    """Single seam for the LLM provider. Returns parsed JSON dict.

    Switching to a different provider (Anthropic, local model, etc.) means
    swapping the body of this function only.
    """
    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        raise RuntimeError("Extractor: OPENAI_API_KEY not set in environment.")

    client = OpenAI(api_key=api_key)

    messages: list[dict] = [{"role": "system", "content": system}]
    for role, content in history:
        api_role = "assistant" if role == "assistant" else "user"
        messages.append({"role": api_role, "content": content})
    messages.append({"role": "user", "content": user})

    response = client.chat.completions.create(
        model=model,
        messages=messages,
        response_format={"type": "json_object"},
        temperature=0,
    )
    content = response.choices[0].message.content
    if not content:
        raise RuntimeError("Extractor: empty response from LLM.")
    return json.loads(content)
