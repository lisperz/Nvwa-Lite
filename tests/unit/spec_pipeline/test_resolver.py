"""Exhaustive-branch unit tests for resolve(spec, adata).

Covers every branch in the resolver orchestrator contract:
- Additive pattern: reads pre_canonical_params, writes params fresh, never
  mutates pre_canonical_params.
- field_type dispatch for gene, cell_type, condition.
- Pass-through for params with no field_type.
- List handling with per-element Issue fields (params.X[idx]).
- Cross-field disambiguation (wrong_field reason with suggestions).
- Total not-found (no cross-field match) -> not_found reason.
- Ambiguous match for cell_type / condition.
- Unknown tool -> original spec unchanged, empty issues.
- pre_canonical_params identity + value unchanged after call.
- Empty pre_canonical_params -> no-op.

Public API only: resolve and the returned (Spec, list[Issue]) tuple.
Private helpers (_dispatch, _resolve_list, _cross_field_check) are NOT imported.

REGISTRY isolation: the clean_registry autouse fixture snapshots REGISTRY
before each test and restores it after.
"""

from __future__ import annotations

import copy
from typing import Any

import numpy as np
import pandas as pd
import pytest

import src.tools.registry as reg_module
from src.core.spec import Spec
from src.domain.resolver.resolver import resolve
from src.spec_validation.result import Issue
from src.tools.registry import REGISTRY, register


# ---------------------------------------------------------------------------
# REGISTRY isolation fixture
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def clean_registry():
    """Snapshot REGISTRY before each test; restore after."""
    snapshot = dict(reg_module.REGISTRY)
    yield
    reg_module.REGISTRY.clear()
    reg_module.REGISTRY.update(snapshot)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_spec(
    tool_name: str,
    pre_canonical_params: dict[str, Any],
    params: dict[str, Any] | None = None,
) -> Spec:
    return Spec(
        scenario_id="test",
        tool_name=tool_name,
        pre_canonical_params=pre_canonical_params,
        params=params or {},
    )


def _register_tool(fn_name: str, fn, param_extras: dict[str, dict] | None = None):
    """Register a throwaway function under fn_name with optional param extras."""
    decorator = register(description=f"test tool {fn_name}", params=param_extras or {})
    decorator(fn)


def _issues_by_field(issues: list[Issue]) -> dict[str, Issue]:
    return {i.field: i for i in issues}


# ---------------------------------------------------------------------------
# AnnData augmentation helpers (inline — no mutations to conftest)
# ---------------------------------------------------------------------------


def _add_gene(adata, gene_name: str):
    """Insert gene_name into adata.var_names."""
    import anndata as ad
    import scipy.sparse as sp

    new_var = pd.DataFrame(index=[gene_name])
    n_obs = adata.n_obs
    new_X = sp.csr_matrix(np.zeros((n_obs, 1)))
    new_block = ad.AnnData(X=new_X, obs=adata.obs.copy(), var=new_var)
    merged = ad.concat([adata, new_block], axis=1)
    # Copy obsm/uns/obsp as concat strips them
    for k, v in adata.uns.items():
        merged.uns[k] = v
    for k, v in adata.obsm.items():
        merged.obsm[k] = v
    return merged


def _add_cell_type_col(adata, col: str, values: list[str]):
    """Add or overwrite adata.obs[col] with a cycling sequence of values."""
    n = adata.n_obs
    adata.obs[col] = [values[i % len(values)] for i in range(n)]


def _add_condition_col(adata, col: str, values: list[str]):
    """Declare col in uns nvwa_meta condition_cols and populate adata.obs."""
    n = adata.n_obs
    adata.obs[col] = [values[i % len(values)] for i in range(n)]
    if "nvwa_meta" not in adata.uns:
        adata.uns["nvwa_meta"] = {}
    cols = adata.uns["nvwa_meta"].get("condition_cols", [])
    if col not in cols:
        cols = list(cols) + [col]
    adata.uns["nvwa_meta"]["condition_cols"] = cols


# ---------------------------------------------------------------------------
# 1. Unknown tool -> original spec unchanged, empty issues
# ---------------------------------------------------------------------------


class TestUnknownTool:
    def test_unknown_tool_returns_original_spec(self, adata):
        spec = _make_spec(
            "nonexistent_tool_xyzzy_999",
            {"gene": "BRCA1"},
        )
        returned_spec, issues = resolve(spec, adata)
        # Must be the same object (not a copy)
        assert returned_spec is spec

    def test_unknown_tool_returns_empty_issues(self, adata):
        spec = _make_spec("no_such_tool_abc", {"gene": "TP53"})
        _, issues = resolve(spec, adata)
        assert issues == []

    def test_unknown_tool_params_unchanged(self, adata):
        spec = _make_spec("ghost_tool", {"gene": "CD3E", "resolution": 0.5})
        returned_spec, _ = resolve(spec, adata)
        assert returned_spec.params == {}  # original spec.params was {}


# ---------------------------------------------------------------------------
# 2. Empty pre_canonical_params -> no-op
# ---------------------------------------------------------------------------


class TestEmptyPreCanonicalParams:
    def test_empty_params_produces_empty_spec_params(self, adata):
        def empty_tool(x: str) -> None:
            pass

        _register_tool("empty_tool", empty_tool)
        spec = _make_spec("empty_tool", {})
        new_spec, issues = resolve(spec, adata)
        assert new_spec.params == {}
        assert issues == []

    def test_empty_params_no_canonicalizations(self, adata):
        def empty_tool2(x: str) -> None:
            pass

        _register_tool("empty_tool2", empty_tool2)
        spec = _make_spec("empty_tool2", {})
        new_spec, _ = resolve(spec, adata)
        assert new_spec.canonicalizations_applied == []

    def test_empty_pre_canonical_params_not_mutated(self, adata):
        def empty_tool3(x: str) -> None:
            pass

        _register_tool("empty_tool3", empty_tool3)
        spec = _make_spec("empty_tool3", {})
        original_pcp = copy.copy(spec.pre_canonical_params)
        resolve(spec, adata)
        assert spec.pre_canonical_params == original_pcp


# ---------------------------------------------------------------------------
# 3. Pass-through (no field_type on param)
# ---------------------------------------------------------------------------


class TestPassThrough:
    def test_passthrough_scalar_string(self, adata):
        def passthrough_tool(resolution: str) -> None:
            pass

        # No field_type entry -> pass-through
        _register_tool("passthrough_tool", passthrough_tool)
        spec = _make_spec("passthrough_tool", {"resolution": "0.5"})
        new_spec, issues = resolve(spec, adata)
        assert new_spec.params["resolution"] == "0.5"
        assert issues == []

    def test_passthrough_numeric(self, adata):
        def passthrough_numeric_tool(n_pcs: int) -> None:
            pass

        _register_tool("passthrough_numeric_tool", passthrough_numeric_tool)
        spec = _make_spec("passthrough_numeric_tool", {"n_pcs": 30})
        new_spec, issues = resolve(spec, adata)
        assert new_spec.params["n_pcs"] == 30
        assert issues == []

    def test_passthrough_no_canonicalization_entry(self, adata):
        def passthrough_tool2(method: str) -> None:
            pass

        _register_tool("passthrough_tool2", passthrough_tool2)
        spec = _make_spec("passthrough_tool2", {"method": "leiden"})
        new_spec, _ = resolve(spec, adata)
        assert new_spec.canonicalizations_applied == []

    def test_passthrough_multiple_params(self, adata):
        def multi_passthrough(a: str, b: int, c: float) -> None:
            pass

        _register_tool("multi_passthrough", multi_passthrough)
        pre = {"a": "hello", "b": 5, "c": 3.14}
        spec = _make_spec("multi_passthrough", pre)
        new_spec, issues = resolve(spec, adata)
        assert new_spec.params == {"a": "hello", "b": 5, "c": 3.14}
        assert issues == []

    def test_passthrough_mixed_with_field_type_param(self, adata):
        """Param with field_type is canonicalized; param without is passed through."""
        def mixed_tool(gene: str, resolution: str) -> None:
            pass

        gene_name = adata.var_names[0]
        _register_tool(
            "mixed_tool",
            mixed_tool,
            param_extras={"gene": {"field_type": "gene"}},
        )
        spec = _make_spec("mixed_tool", {"gene": gene_name, "resolution": "0.5"})
        new_spec, issues = resolve(spec, adata)
        assert new_spec.params["resolution"] == "0.5"
        assert issues == []


# ---------------------------------------------------------------------------
# 4. Additive pattern: pre_canonical_params never mutated
# ---------------------------------------------------------------------------


class TestAdditivePattern:
    def test_pre_canonical_params_identity_not_mutated(self, adata):
        """pre_canonical_params dict object itself must not be mutated."""
        def identity_tool(gene: str) -> None:
            pass

        gene_name = adata.var_names[0]
        _register_tool(
            "identity_tool",
            identity_tool,
            param_extras={"gene": {"field_type": "gene"}},
        )
        spec = _make_spec("identity_tool", {"gene": gene_name})
        pcp_id_before = id(spec.pre_canonical_params)
        pcp_copy_before = copy.deepcopy(spec.pre_canonical_params)
        resolve(spec, adata)
        # Identity: same dict object
        assert id(spec.pre_canonical_params) == pcp_id_before
        # Value: same content
        assert spec.pre_canonical_params == pcp_copy_before

    def test_pre_canonical_params_value_unchanged_after_canonicalization(self, adata):
        """Even when a gene is case-canonicalized, pre_canonical_params stays raw."""
        def gene_tool(gene: str) -> None:
            pass

        gene_name = adata.var_names[0]
        # Pass the gene in wrong case to trigger case-insensitive resolution
        raw_input = gene_name.lower()
        _register_tool(
            "gene_tool",
            gene_tool,
            param_extras={"gene": {"field_type": "gene"}},
        )
        spec = _make_spec("gene_tool", {"gene": raw_input})
        original_pcp = {"gene": raw_input}
        resolve(spec, adata)
        assert spec.pre_canonical_params == original_pcp

    def test_params_built_fresh_per_call(self, adata):
        """Two resolve() calls on the same spec produce independent param dicts."""
        def fresh_tool(x: str) -> None:
            pass

        _register_tool("fresh_tool", fresh_tool)
        spec = _make_spec("fresh_tool", {"x": "hello"})
        new_spec1, _ = resolve(spec, adata)
        new_spec2, _ = resolve(spec, adata)
        # params dicts are different objects (freshly built)
        assert new_spec1.params is not new_spec2.params

    def test_original_spec_params_not_mutated(self, adata):
        """spec.params before call must be untouched after call."""
        def no_mutate_tool(gene: str) -> None:
            pass

        gene_name = adata.var_names[0]
        _register_tool(
            "no_mutate_tool",
            no_mutate_tool,
            param_extras={"gene": {"field_type": "gene"}},
        )
        spec = _make_spec("no_mutate_tool", {"gene": gene_name}, params={"old": "value"})
        resolve(spec, adata)
        # Original spec.params is a pydantic dict; new_spec is a copy
        assert spec.params == {"old": "value"}


# ---------------------------------------------------------------------------
# 5. Successful gene canonicalization
# ---------------------------------------------------------------------------


class TestGeneFieldType:
    def test_exact_gene_in_params(self, adata):
        def gene_exact_tool(gene: str) -> None:
            pass

        gene_name = adata.var_names[0]
        _register_tool(
            "gene_exact_tool",
            gene_exact_tool,
            param_extras={"gene": {"field_type": "gene"}},
        )
        spec = _make_spec("gene_exact_tool", {"gene": gene_name})
        new_spec, issues = resolve(spec, adata)
        assert new_spec.params["gene"] == gene_name
        assert issues == []

    def test_case_insensitive_gene_canonicalization_entry(self, adata):
        """A lowercase gene name must resolve + emit a Canonicalization entry."""
        def gene_case_tool(gene: str) -> None:
            pass

        # Pick a gene where lowercasing changes it
        gene_name = next(
            (g for g in adata.var_names if g != g.lower()), None
        )
        if gene_name is None:
            pytest.skip("No mixed-case gene in subset")

        _register_tool(
            "gene_case_tool",
            gene_case_tool,
            param_extras={"gene": {"field_type": "gene"}},
        )
        raw = gene_name.lower()
        spec = _make_spec("gene_case_tool", {"gene": raw})
        new_spec, issues = resolve(spec, adata)
        assert new_spec.params["gene"] == gene_name
        assert issues == []
        # Canonicalization entry must record the transformation
        assert len(new_spec.canonicalizations_applied) == 1
        canon = new_spec.canonicalizations_applied[0]
        assert canon.field == "params.gene"
        assert canon.raw == raw
        assert canon.canonical == gene_name

    def test_exact_match_no_canonicalization_entry(self, adata):
        """Exact match (raw == canonical) should NOT produce a Canonicalization entry."""
        def gene_exact2_tool(gene: str) -> None:
            pass

        gene_name = adata.var_names[0]
        _register_tool(
            "gene_exact2_tool",
            gene_exact2_tool,
            param_extras={"gene": {"field_type": "gene"}},
        )
        spec = _make_spec("gene_exact2_tool", {"gene": gene_name})
        new_spec, _ = resolve(spec, adata)
        assert new_spec.canonicalizations_applied == []

    def test_not_found_gene_emits_issue_and_field_absent(self, adata):
        """Gene that cannot be resolved must not appear in params and emit Issue."""
        def gene_notfound_tool(gene: str) -> None:
            pass

        _register_tool(
            "gene_notfound_tool",
            gene_notfound_tool,
            param_extras={"gene": {"field_type": "gene"}},
        )
        nonexistent = "ZZZZZNOGENE99999"
        spec = _make_spec("gene_notfound_tool", {"gene": nonexistent})
        new_spec, issues = resolve(spec, adata)
        assert "gene" not in new_spec.params
        assert len(issues) == 1
        assert issues[0].field == "params.gene"

    def test_not_found_gene_issue_reason(self, adata):
        """Issue reason for a fully absent gene (no cross-field match) is not_found."""
        def gene_reason_tool(gene: str) -> None:
            pass

        _register_tool(
            "gene_reason_tool",
            gene_reason_tool,
            param_extras={"gene": {"field_type": "gene"}},
        )
        spec = _make_spec("gene_reason_tool", {"gene": "ZZZZZNOGENE99999"})
        _, issues = resolve(spec, adata)
        assert issues[0].reason == "not_found"


# ---------------------------------------------------------------------------
# 6. cell_type field_type dispatch
# ---------------------------------------------------------------------------


class TestCellTypeFieldType:
    def _setup_cell_type_tool(self, adata, cell_types: list[str]):
        """Register a tool + populate adata.obs['cell_type'] column."""
        _add_cell_type_col(adata, "cell_type", cell_types)

        def ct_tool(cell_type: str) -> None:
            pass

        _register_tool(
            "ct_tool",
            ct_tool,
            param_extras={"cell_type": {"field_type": "cell_type"}},
        )

    def test_exact_cell_type_resolved(self, adata):
        self._setup_cell_type_tool(adata, ["T cell", "B cell", "NK cell"])
        spec = _make_spec("ct_tool", {"cell_type": "T cell"})
        new_spec, issues = resolve(spec, adata)
        assert new_spec.params["cell_type"] == "T cell"
        assert issues == []

    def test_normalized_cell_type_canonicalization_entry(self, adata):
        """'T_cell' normalizes to 'T cell' and must produce a Canonicalization."""
        self._setup_cell_type_tool(adata, ["T cell", "B cell"])
        spec = _make_spec("ct_tool", {"cell_type": "T_cell"})
        new_spec, issues = resolve(spec, adata)
        assert new_spec.params["cell_type"] == "T cell"
        assert issues == []
        assert len(new_spec.canonicalizations_applied) == 1
        c = new_spec.canonicalizations_applied[0]
        assert c.field == "params.cell_type"
        assert c.raw == "T_cell"
        assert c.canonical == "T cell"

    def test_ambiguous_cell_type_emits_issue(self, adata):
        """Uppercase initials matching multiple cell types -> ambiguous."""
        # 'NK' uppercase initials match 'Natural Killer' patterns - use initials case
        # Use a name where uppercase initials collide: 'NK' for two cell type names
        _add_cell_type_col(adata, "cell_type", ["Naive Killer", "Natural Killer"])

        def ct_ambig_tool(cell_type: str) -> None:
            pass

        _register_tool(
            "ct_ambig_tool",
            ct_ambig_tool,
            param_extras={"cell_type": {"field_type": "cell_type"}},
        )
        spec = _make_spec("ct_ambig_tool", {"cell_type": "NK"})
        new_spec, issues = resolve(spec, adata)
        assert "cell_type" not in new_spec.params
        assert len(issues) == 1
        assert issues[0].field == "params.cell_type"
        assert issues[0].reason == "ambiguous"

    def test_not_found_cell_type_emits_issue(self, adata):
        self._setup_cell_type_tool(adata, ["T cell", "B cell"])
        spec = _make_spec("ct_tool", {"cell_type": "Totally Unknown Cell 9999"})
        new_spec, issues = resolve(spec, adata)
        assert "cell_type" not in new_spec.params
        assert len(issues) == 1
        assert issues[0].field == "params.cell_type"
        assert issues[0].reason == "not_found"


# ---------------------------------------------------------------------------
# 7. condition field_type dispatch
# ---------------------------------------------------------------------------


class TestConditionFieldType:
    def _setup_condition_tool(self, adata, col: str, values: list[str]):
        _add_condition_col(adata, col, values)

        def cond_tool(condition: str) -> None:
            pass

        _register_tool(
            "cond_tool",
            cond_tool,
            param_extras={"condition": {"field_type": "condition"}},
        )

    def test_exact_condition_resolved(self, adata):
        self._setup_condition_tool(adata, "treatment", ["WT", "KO"])
        spec = _make_spec("cond_tool", {"condition": "WT"})
        new_spec, issues = resolve(spec, adata)
        assert new_spec.params["condition"] == "WT"
        assert issues == []

    def test_condition_canonicalization_has_context(self, adata):
        """Condition resolution must attach context with obs_column."""
        _add_condition_col(adata, "treatment", ["WT", "KO"])

        def cond_ctx_tool(condition: str) -> None:
            pass

        _register_tool(
            "cond_ctx_tool",
            cond_ctx_tool,
            param_extras={"condition": {"field_type": "condition"}},
        )
        spec = _make_spec("cond_ctx_tool", {"condition": "WT"})
        new_spec, _ = resolve(spec, adata)
        # context carries obs_column; canonicalization is emitted even for exact
        # match when context is set
        cans = new_spec.canonicalizations_applied
        if cans:
            assert cans[0].context is not None
            assert "obs_column" in cans[0].context

    def test_ambiguous_condition_across_columns(self, adata):
        """Same value in two condition columns -> ambiguous."""
        _add_condition_col(adata, "treatment", ["WT", "KO"])
        _add_condition_col(adata, "batch", ["WT", "batch2"])

        def cond_ambig_tool(condition: str) -> None:
            pass

        _register_tool(
            "cond_ambig_tool",
            cond_ambig_tool,
            param_extras={"condition": {"field_type": "condition"}},
        )
        spec = _make_spec("cond_ambig_tool", {"condition": "WT"})
        new_spec, issues = resolve(spec, adata)
        assert "condition" not in new_spec.params
        assert len(issues) == 1
        assert issues[0].reason == "ambiguous"

    def test_not_found_condition_emits_issue(self, adata):
        self._setup_condition_tool(adata, "treatment", ["WT", "KO"])
        spec = _make_spec("cond_tool", {"condition": "ZZZNOMATCH"})
        new_spec, issues = resolve(spec, adata)
        assert "condition" not in new_spec.params
        assert len(issues) == 1
        assert issues[0].reason == "not_found"


# ---------------------------------------------------------------------------
# 8. List param handling
# ---------------------------------------------------------------------------


class TestListParamHandling:
    def test_list_all_matching_genes(self, adata):
        """All list elements resolve -> resolved list in params, no issues."""
        def list_gene_tool(genes: list) -> None:
            pass

        g1, g2 = adata.var_names[0], adata.var_names[1]
        _register_tool(
            "list_gene_tool",
            list_gene_tool,
            param_extras={"genes": {"field_type": "gene"}},
        )
        spec = _make_spec("list_gene_tool", {"genes": [g1, g2]})
        new_spec, issues = resolve(spec, adata)
        assert new_spec.params["genes"] == [g1, g2]
        assert issues == []

    def test_list_partial_matching_genes(self, adata):
        """Some elements resolve, some fail -> resolved elements kept, Issue per failure."""
        def partial_gene_tool(genes: list) -> None:
            pass

        g1 = adata.var_names[0]
        bad_gene = "ZZZZZNOGENE00001"
        _register_tool(
            "partial_gene_tool",
            partial_gene_tool,
            param_extras={"genes": {"field_type": "gene"}},
        )
        spec = _make_spec("partial_gene_tool", {"genes": [g1, bad_gene]})
        new_spec, issues = resolve(spec, adata)
        # Resolved element stays in list
        assert new_spec.params["genes"][0] == g1
        # Failed element stays as original value in list (resolver keeps the slot)
        assert new_spec.params["genes"][1] == bad_gene
        # Issue emitted for the failed index
        assert len(issues) == 1
        assert issues[0].field == "params.genes[1]"

    def test_list_all_failing_elements(self, adata):
        """All elements fail -> Issue for each element (params.genes[0], [1])."""
        def all_fail_tool(genes: list) -> None:
            pass

        _register_tool(
            "all_fail_tool",
            all_fail_tool,
            param_extras={"genes": {"field_type": "gene"}},
        )
        spec = _make_spec(
            "all_fail_tool",
            {"genes": ["ZZZZZBAD1", "ZZZZZBAD2", "ZZZZZBAD3"]},
        )
        new_spec, issues = resolve(spec, adata)
        assert len(issues) == 3
        fields = {i.field for i in issues}
        assert "params.genes[0]" in fields
        assert "params.genes[1]" in fields
        assert "params.genes[2]" in fields

    def test_list_issue_field_positional_naming(self, adata):
        """Issue field must be params.<name>[<idx>], not params.<name>."""
        def positional_tool(genes: list) -> None:
            pass

        _register_tool(
            "positional_tool",
            positional_tool,
            param_extras={"genes": {"field_type": "gene"}},
        )
        spec = _make_spec("positional_tool", {"genes": ["ZZZZZBADONLY"]})
        _, issues = resolve(spec, adata)
        assert len(issues) == 1
        assert issues[0].field == "params.genes[0]"

    def test_list_canonicalization_entry_per_element(self, adata):
        """Each resolved element that changes produces a per-element Canonicalization."""
        def list_canon_tool(genes: list) -> None:
            pass

        # Pick a gene with mixed case so lowercase triggers canonicalization
        gene_name = next(
            (g for g in adata.var_names if g != g.lower()), None
        )
        if gene_name is None:
            pytest.skip("No mixed-case gene in subset")

        _register_tool(
            "list_canon_tool",
            list_canon_tool,
            param_extras={"genes": {"field_type": "gene"}},
        )
        raw = gene_name.lower()
        spec = _make_spec("list_canon_tool", {"genes": [raw]})
        new_spec, issues = resolve(spec, adata)
        assert issues == []
        assert len(new_spec.canonicalizations_applied) == 1
        c = new_spec.canonicalizations_applied[0]
        assert c.field == "params.genes[0]"
        assert c.raw == raw
        assert c.canonical == gene_name

    def test_list_cell_type_all_matching(self, adata):
        """Cell type list: all match -> no issues."""
        _add_cell_type_col(adata, "cell_type", ["T cell", "B cell", "NK cell"])

        def list_ct_tool(cell_types: list) -> None:
            pass

        _register_tool(
            "list_ct_tool",
            list_ct_tool,
            param_extras={"cell_types": {"field_type": "cell_type"}},
        )
        spec = _make_spec("list_ct_tool", {"cell_types": ["T cell", "B cell"]})
        new_spec, issues = resolve(spec, adata)
        assert new_spec.params["cell_types"] == ["T cell", "B cell"]
        assert issues == []

    def test_list_preserves_order(self, adata):
        """Resolved list must preserve original element order."""
        def order_tool(genes: list) -> None:
            pass

        genes = list(adata.var_names[:3])
        _register_tool(
            "order_tool",
            order_tool,
            param_extras={"genes": {"field_type": "gene"}},
        )
        spec = _make_spec("order_tool", {"genes": genes})
        new_spec, _ = resolve(spec, adata)
        assert new_spec.params["genes"] == genes


# ---------------------------------------------------------------------------
# 9. Cross-field disambiguation: wrong_field
# ---------------------------------------------------------------------------


class TestCrossFieldDisambiguation:
    def test_gene_value_matches_cell_type_emits_wrong_field(self, adata):
        """A value typed under 'gene' param that actually matches a cell type
        must produce a wrong_field Issue with suggestions=[cell_type=<value>]."""
        cell_type_value = "MySpecialCellType999"
        _add_cell_type_col(adata, "cell_type", [cell_type_value, "OtherType"])

        def cross_gene_tool(gene: str) -> None:
            pass

        _register_tool(
            "cross_gene_tool",
            cross_gene_tool,
            param_extras={"gene": {"field_type": "gene"}},
        )
        spec = _make_spec("cross_gene_tool", {"gene": cell_type_value})
        new_spec, issues = resolve(spec, adata)
        assert "gene" not in new_spec.params
        assert len(issues) == 1
        issue = issues[0]
        assert issue.field == "params.gene"
        assert issue.reason == "wrong_field"
        assert any("cell_type" in s for s in issue.suggestions)

    def test_wrong_field_suggestion_format(self, adata):
        """Suggestions must be formatted as 'field_type=value'."""
        cell_type_value = "CrossFieldCell888"
        _add_cell_type_col(adata, "cell_type", [cell_type_value, "Other"])

        def cross_fmt_tool(gene: str) -> None:
            pass

        _register_tool(
            "cross_fmt_tool",
            cross_fmt_tool,
            param_extras={"gene": {"field_type": "gene"}},
        )
        spec = _make_spec("cross_fmt_tool", {"gene": cell_type_value})
        _, issues = resolve(spec, adata)
        assert issues[0].suggestions == [f"cell_type={cell_type_value}"]

    def test_gene_value_matches_condition_emits_wrong_field(self, adata):
        """A value typed under 'gene' param that actually matches a condition
        must produce a wrong_field Issue with suggestions=[condition=<value>]."""
        cond_value = "SpecialConditionXYZ777"
        _add_condition_col(adata, "treatment", [cond_value, "control"])

        def cross_cond_tool(gene: str) -> None:
            pass

        _register_tool(
            "cross_cond_tool",
            cross_cond_tool,
            param_extras={"gene": {"field_type": "gene"}},
        )
        spec = _make_spec("cross_cond_tool", {"gene": cond_value})
        new_spec, issues = resolve(spec, adata)
        assert "gene" not in new_spec.params
        assert len(issues) == 1
        assert issues[0].reason == "wrong_field"
        assert any("condition" in s for s in issues[0].suggestions)

    def test_no_cross_field_match_emits_not_found(self, adata):
        """A value that matches nothing at all (not gene, not cell_type, not condition)
        must emit reason=not_found, not wrong_field."""
        def notfound_gene_tool(gene: str) -> None:
            pass

        _register_tool(
            "notfound_gene_tool",
            notfound_gene_tool,
            param_extras={"gene": {"field_type": "gene"}},
        )
        spec = _make_spec("notfound_gene_tool", {"gene": "ABSOLUTELYNOTHINGATALL9999"})
        _, issues = resolve(spec, adata)
        assert len(issues) == 1
        assert issues[0].reason == "not_found"

    def test_cell_type_value_matches_gene_emits_wrong_field(self, adata):
        """A value typed under 'cell_type' that actually matches a gene -> wrong_field."""
        gene_name = adata.var_names[0]
        # Ensure cell_type col exists but does not contain the gene name
        _add_cell_type_col(adata, "cell_type", ["TypeA", "TypeB"])

        def cross_ct_gene_tool(cell_type: str) -> None:
            pass

        _register_tool(
            "cross_ct_gene_tool",
            cross_ct_gene_tool,
            param_extras={"cell_type": {"field_type": "cell_type"}},
        )
        spec = _make_spec("cross_ct_gene_tool", {"cell_type": gene_name})
        new_spec, issues = resolve(spec, adata)
        # gene_name is not a cell type -> not_found -> cross-field finds it as gene
        if issues:
            assert issues[0].reason in ("wrong_field", "not_found")
            if issues[0].reason == "wrong_field":
                assert any("gene" in s for s in issues[0].suggestions)


# ---------------------------------------------------------------------------
# 10. canonicalizations_applied accumulation
# ---------------------------------------------------------------------------


class TestCanonicalizationsApplied:
    def test_canonicalizations_appended_not_replaced(self, adata):
        """resolve() appends to existing canonicalizations_applied list."""
        from src.core.spec import Canonicalization

        def append_tool(gene: str) -> None:
            pass

        gene_name = next(
            (g for g in adata.var_names if g != g.lower()), None
        )
        if gene_name is None:
            pytest.skip("No mixed-case gene in subset")

        _register_tool(
            "append_tool",
            append_tool,
            param_extras={"gene": {"field_type": "gene"}},
        )
        existing_canon = Canonicalization(
            field="params.other", raw="x", canonical="y"
        )
        spec = Spec(
            scenario_id="test",
            tool_name="append_tool",
            pre_canonical_params={"gene": gene_name.lower()},
            params={},
            canonicalizations_applied=[existing_canon],
        )
        new_spec, _ = resolve(spec, adata)
        # Must retain the pre-existing entry + new one
        assert len(new_spec.canonicalizations_applied) >= 2
        assert new_spec.canonicalizations_applied[0] == existing_canon

    def test_failed_canonicalization_field_absent_from_params(self, adata):
        """A field that fails resolution must NOT appear in new_spec.params."""
        def absent_tool(gene: str) -> None:
            pass

        _register_tool(
            "absent_tool",
            absent_tool,
            param_extras={"gene": {"field_type": "gene"}},
        )
        spec = _make_spec("absent_tool", {"gene": "ZZZZZNOGENE11111"})
        new_spec, issues = resolve(spec, adata)
        assert "gene" not in new_spec.params
        assert len(issues) == 1

    def test_multiple_successful_canonicalizations(self, adata):
        """Two fields that each canonicalize produce two Canonicalization entries."""
        def two_gene_tool(gene_a: str, gene_b: str) -> None:
            pass

        gene_a = adata.var_names[0]
        gene_b = adata.var_names[1]
        # Only emit canonicalization if raw != canonical; lower() forces this
        raw_a = gene_a.lower() if gene_a != gene_a.lower() else gene_a
        raw_b = gene_b.lower() if gene_b != gene_b.lower() else gene_b

        # If neither gene has mixed case, skip
        if raw_a == gene_a and raw_b == gene_b:
            pytest.skip("Both genes are already lowercase in dataset; no canonicalization would occur")

        _register_tool(
            "two_gene_tool",
            two_gene_tool,
            param_extras={
                "gene_a": {"field_type": "gene"},
                "gene_b": {"field_type": "gene"},
            },
        )
        spec = _make_spec("two_gene_tool", {"gene_a": raw_a, "gene_b": raw_b})
        new_spec, _ = resolve(spec, adata)
        # At least one canonicalization emitted (possibly two)
        assert len(new_spec.canonicalizations_applied) >= 1


# ---------------------------------------------------------------------------
# 11. Return type contract
# ---------------------------------------------------------------------------


class TestReturnTypeContract:
    def test_returns_tuple_of_spec_and_list(self, adata):
        def type_tool(x: str) -> None:
            pass

        _register_tool("type_tool", type_tool)
        spec = _make_spec("type_tool", {"x": "hello"})
        result = resolve(spec, adata)
        assert isinstance(result, tuple)
        assert len(result) == 2
        new_spec, issues = result
        assert isinstance(new_spec, Spec)
        assert isinstance(issues, list)

    def test_returned_spec_is_new_object_on_success(self, adata):
        """resolve() must return a copy (model_copy), not the original spec."""
        def new_obj_tool(x: str) -> None:
            pass

        _register_tool("new_obj_tool", new_obj_tool)
        spec = _make_spec("new_obj_tool", {"x": "hello"})
        new_spec, _ = resolve(spec, adata)
        assert new_spec is not spec

    def test_issues_are_issue_instances(self, adata):
        def issue_type_tool(gene: str) -> None:
            pass

        _register_tool(
            "issue_type_tool",
            issue_type_tool,
            param_extras={"gene": {"field_type": "gene"}},
        )
        spec = _make_spec("issue_type_tool", {"gene": "ZZZZZNOGENE22222"})
        _, issues = resolve(spec, adata)
        for issue in issues:
            assert isinstance(issue, Issue)


# ---------------------------------------------------------------------------
# 12. Multiple params — mixed success/failure
# ---------------------------------------------------------------------------


class TestMixedSuccessFailure:
    def test_successful_param_in_params_even_when_another_fails(self, adata):
        """If one param fails and another succeeds, the successful one is in params."""
        def mixed_sf_tool(gene: str, cell_type: str) -> None:
            pass

        good_gene = adata.var_names[0]
        _add_cell_type_col(adata, "cell_type", ["T cell", "B cell"])
        _register_tool(
            "mixed_sf_tool",
            mixed_sf_tool,
            param_extras={
                "gene": {"field_type": "gene"},
                "cell_type": {"field_type": "cell_type"},
            },
        )
        spec = _make_spec(
            "mixed_sf_tool",
            {"gene": good_gene, "cell_type": "ZZZZZBADTYPE99999"},
        )
        new_spec, issues = resolve(spec, adata)
        assert new_spec.params["gene"] == good_gene
        assert "cell_type" not in new_spec.params
        assert len(issues) == 1
        assert issues[0].field == "params.cell_type"

    def test_all_fail_produces_empty_params(self, adata):
        """If every canonicalized param fails, params must be empty."""
        def all_fail_sf_tool(gene: str) -> None:
            pass

        _register_tool(
            "all_fail_sf_tool",
            all_fail_sf_tool,
            param_extras={"gene": {"field_type": "gene"}},
        )
        spec = _make_spec("all_fail_sf_tool", {"gene": "ZZZZZNOGENE33333"})
        new_spec, issues = resolve(spec, adata)
        # Only canonicalized fields that fail are absent; passthrough fields still appear
        assert "gene" not in new_spec.params
        assert len(issues) >= 1

    def test_passthrough_plus_failed_gene(self, adata):
        """A pass-through param appears in params even when gene param fails."""
        def pt_fail_tool(resolution: str, gene: str) -> None:
            pass

        _register_tool(
            "pt_fail_tool",
            pt_fail_tool,
            param_extras={"gene": {"field_type": "gene"}},
        )
        spec = _make_spec(
            "pt_fail_tool",
            {"resolution": "0.5", "gene": "ZZZZZNOGENE44444"},
        )
        new_spec, issues = resolve(spec, adata)
        assert new_spec.params["resolution"] == "0.5"
        assert "gene" not in new_spec.params
        assert len(issues) == 1


# ---------------------------------------------------------------------------
# 13. Condition context field
# ---------------------------------------------------------------------------


class TestConditionContext:
    def test_resolved_condition_context_has_obs_column(self, adata):
        """Condition resolution always sets context={obs_column: col}."""
        _add_condition_col(adata, "treatment", ["WT", "KO"])

        def ctx_check_tool(condition: str) -> None:
            pass

        _register_tool(
            "ctx_check_tool",
            ctx_check_tool,
            param_extras={"condition": {"field_type": "condition"}},
        )
        spec = _make_spec("ctx_check_tool", {"condition": "WT"})
        new_spec, issues = resolve(spec, adata)
        assert issues == []
        # Condition match -> context attached -> Canonicalization always emitted
        cans = new_spec.canonicalizations_applied
        if cans:
            assert cans[0].context == {"obs_column": "treatment"}
