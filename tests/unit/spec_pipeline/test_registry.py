"""Exhaustive-branch unit tests for src/tools/registry.py.

Covers every public function and every documented behavioral branch of the
internal _build_param_specs / _type_tag helpers, exercised only through the
public API (no direct import of private helpers).

REGISTRY isolation: the ``clean_registry`` autouse fixture snapshots and
restores the module-level REGISTRY dict so no test leaks state.
"""

from __future__ import annotations

import inspect
import sys
from typing import List, Optional

import pytest

import src.tools.registry as reg_module
from src.tools.registry import (
    REGISTRY,
    ParamSpec,
    ToolEntry,
    get_default_value,
    get_tool,
    get_tool_description_block,
    get_tool_names,
    register,
)


# ---------------------------------------------------------------------------
# Isolation fixture
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def clean_registry():
    """Save REGISTRY state before each test, restore it after."""
    snapshot = dict(reg_module.REGISTRY)
    yield
    reg_module.REGISTRY.clear()
    reg_module.REGISTRY.update(snapshot)


# ---------------------------------------------------------------------------
# Helpers — register a function and return its ToolEntry
# ---------------------------------------------------------------------------


def _reg(fn, *, description="test tool", kind="atomic", params=None) -> ToolEntry:
    """Register *fn* and return the resulting ToolEntry."""
    decorator = register(description=description, kind=kind, params=params or {})
    decorator(fn)
    return reg_module.REGISTRY[fn.__name__]


# ---------------------------------------------------------------------------
# @register — basic registration
# ---------------------------------------------------------------------------


class TestRegisterBasic:
    def test_entry_stored_under_function_name(self):
        def my_tool(x: str) -> None:
            pass

        _reg(my_tool)
        assert "my_tool" in reg_module.REGISTRY

    def test_entry_name_matches_function_name(self):
        def tool_alpha(x: str) -> None:
            pass

        entry = _reg(tool_alpha)
        assert entry.name == "tool_alpha"

    def test_entry_description_stored(self):
        def tool_beta(x: str) -> None:
            pass

        entry = _reg(tool_beta, description="beta description")
        assert entry.description == "beta description"

    def test_entry_callable_is_original_function(self):
        def tool_gamma(x: str) -> None:
            pass

        entry = _reg(tool_gamma)
        assert entry.callable is tool_gamma

    def test_register_returns_original_function(self):
        def tool_delta(x: str) -> None:
            pass

        result = register(description="d")(tool_delta)
        assert result is tool_delta

    def test_register_overwrites_on_duplicate_name(self):
        def dup_tool(x: str) -> None:
            pass

        _reg(dup_tool, description="first")

        def dup_tool(x: str) -> None:  # noqa: F811 — intentional redefinition
            pass

        entry = _reg(dup_tool, description="second")
        assert entry.description == "second"
        assert len([k for k in reg_module.REGISTRY if k == "dup_tool"]) == 1


# ---------------------------------------------------------------------------
# @register — kind field
# ---------------------------------------------------------------------------


class TestRegisterKind:
    def test_default_kind_is_atomic(self):
        def atomic_tool(x: str) -> None:
            pass

        entry = _reg(atomic_tool)
        assert entry.kind == "atomic"

    def test_kind_workflow(self):
        def workflow_tool(x: str) -> None:
            pass

        entry = _reg(workflow_tool, kind="workflow")
        assert entry.kind == "workflow"

    def test_kind_atomic_explicit(self):
        def another_atomic(x: str) -> None:
            pass

        entry = _reg(another_atomic, kind="atomic")
        assert entry.kind == "atomic"


# ---------------------------------------------------------------------------
# _build_param_specs — parameter filtering
# ---------------------------------------------------------------------------


class TestBuildParamSpecsFiltering:
    def test_self_is_skipped(self):
        class FakeClass:
            def my_method(self, x: str) -> None:
                pass

        entry = _reg(FakeClass.my_method)
        names = [p.name for p in entry.params]
        assert "self" not in names
        assert "x" in names

    def test_cls_is_skipped(self):
        class FakeClass:
            @classmethod
            def my_classmethod(cls, x: str) -> None:
                pass

        # classmethod wraps the function; access __func__ to avoid descriptor
        entry = _reg(FakeClass.my_classmethod.__func__)
        names = [p.name for p in entry.params]
        assert "cls" not in names

    def test_adata_is_skipped(self):
        def tool_with_adata(adata, x: str) -> None:
            pass

        entry = _reg(tool_with_adata)
        names = [p.name for p in entry.params]
        assert "adata" not in names
        assert "x" in names

    def test_var_positional_is_skipped(self):
        def tool_with_args(x: str, *args) -> None:
            pass

        entry = _reg(tool_with_args)
        names = [p.name for p in entry.params]
        assert "args" not in names
        assert "x" in names

    def test_var_keyword_is_skipped(self):
        def tool_with_kwargs(x: str, **kwargs) -> None:
            pass

        entry = _reg(tool_with_kwargs)
        names = [p.name for p in entry.params]
        assert "kwargs" not in names
        assert "x" in names

    def test_all_skip_targets_together(self):
        def tool_mixed(adata, x: str, *args, **kwargs) -> None:
            pass

        entry = _reg(tool_mixed)
        names = [p.name for p in entry.params]
        assert names == ["x"]

    def test_no_params_produces_empty_list(self):
        def tool_no_params() -> None:
            pass

        entry = _reg(tool_no_params)
        assert entry.params == []

    def test_only_adata_produces_empty_list(self):
        def tool_only_adata(adata) -> None:
            pass

        entry = _reg(tool_only_adata)
        assert entry.params == []


# ---------------------------------------------------------------------------
# _type_tag — type annotation variants (via public API)
# ---------------------------------------------------------------------------


class TestTypeTagViaPublicAPI:
    def test_str_annotation(self):
        def tool_str(x: str) -> None:
            pass

        entry = _reg(tool_str)
        assert entry.params[0].type == "str"

    def test_int_annotation(self):
        def tool_int(x: int) -> None:
            pass

        entry = _reg(tool_int)
        assert entry.params[0].type == "int"

    def test_float_annotation(self):
        def tool_float(x: float) -> None:
            pass

        entry = _reg(tool_float)
        assert entry.params[0].type == "float"

    def test_bool_annotation(self):
        def tool_bool(x: bool) -> None:
            pass

        entry = _reg(tool_bool)
        assert entry.params[0].type == "bool"

    def test_no_annotation_returns_any(self):
        def tool_no_ann(x) -> None:
            pass

        entry = _reg(tool_no_ann)
        assert entry.params[0].type == "Any"

    def test_optional_str_annotation(self):
        def tool_opt_str(x: Optional[str]) -> None:
            pass

        entry = _reg(tool_opt_str)
        # "typing." prefix should be stripped; result contains "Optional[str]"
        assert "Optional" in entry.params[0].type
        assert "typing." not in entry.params[0].type

    def test_optional_int_annotation(self):
        def tool_opt_int(x: Optional[int]) -> None:
            pass

        entry = _reg(tool_opt_int)
        assert "Optional" in entry.params[0].type
        assert "typing." not in entry.params[0].type

    def test_list_builtin_generic(self):
        def tool_list(x: list[str]) -> None:
            pass

        entry = _reg(tool_list)
        # builtin list[str] — no "typing." prefix; check no class wrapper remains
        assert "list" in entry.params[0].type
        assert "<class" not in entry.params[0].type

    def test_typing_list_generic(self):
        def tool_typing_list(x: List[str]) -> None:
            pass

        entry = _reg(tool_typing_list)
        assert "List" in entry.params[0].type
        assert "typing." not in entry.params[0].type

    def test_custom_class_annotation(self):
        class MyCustomClass:
            pass

        def tool_custom(x: MyCustomClass) -> None:
            pass

        entry = _reg(tool_custom)
        # _type_tag strips the "<class '...'>" wrapper
        assert "<class" not in entry.params[0].type
        assert "MyCustomClass" in entry.params[0].type

    def test_bare_list_annotation(self):
        def tool_bare_list(x: list) -> None:
            pass

        entry = _reg(tool_bare_list)
        assert entry.params[0].type == "list"


# ---------------------------------------------------------------------------
# _build_param_specs — required / default
# ---------------------------------------------------------------------------


class TestRequiredAndDefault:
    def test_no_default_is_required(self):
        def tool_required(x: str) -> None:
            pass

        entry = _reg(tool_required)
        p = entry.params[0]
        assert p.required is True
        assert p.default is None

    def test_with_default_is_not_required(self):
        def tool_optional_param(x: str = "hello") -> None:
            pass

        entry = _reg(tool_optional_param)
        p = entry.params[0]
        assert p.required is False
        assert p.default == "hello"

    def test_default_none_is_not_required(self):
        def tool_none_default(x: Optional[str] = None) -> None:
            pass

        entry = _reg(tool_none_default)
        p = entry.params[0]
        assert p.required is False
        assert p.default is None

    def test_default_int(self):
        def tool_int_default(n: int = 10) -> None:
            pass

        entry = _reg(tool_int_default)
        p = entry.params[0]
        assert p.required is False
        assert p.default == 10

    def test_default_bool_false(self):
        def tool_bool_default(flag: bool = False) -> None:
            pass

        entry = _reg(tool_bool_default)
        p = entry.params[0]
        assert p.required is False
        assert p.default is False

    def test_mixed_required_and_default(self):
        def tool_mixed(a: str, b: int = 5) -> None:
            pass

        entry = _reg(tool_mixed)
        params = {p.name: p for p in entry.params}
        assert params["a"].required is True
        assert params["b"].required is False
        assert params["b"].default == 5

    def test_param_order_preserved(self):
        def tool_ordered(first: str, second: int, third: float = 1.0) -> None:
            pass

        entry = _reg(tool_ordered)
        names = [p.name for p in entry.params]
        assert names == ["first", "second", "third"]


# ---------------------------------------------------------------------------
# _build_param_specs — extras dict (enum, field_type, description)
# ---------------------------------------------------------------------------


class TestExtrasDict:
    def test_enum_override(self):
        def tool_enum(resolution: str) -> None:
            pass

        extras = {"resolution": {"enum": ["low", "medium", "high"]}}
        entry = _reg(tool_enum, params=extras)
        p = entry.params[0]
        assert p.enum == ["low", "medium", "high"]

    def test_field_type_override(self):
        def tool_field_type(gene_name: str) -> None:
            pass

        extras = {"gene_name": {"field_type": "gene"}}
        entry = _reg(tool_field_type, params=extras)
        p = entry.params[0]
        assert p.field_type == "gene"

    def test_description_override(self):
        def tool_desc(cluster: str) -> None:
            pass

        extras = {"cluster": {"description": "cluster label to highlight"}}
        entry = _reg(tool_desc, params=extras)
        p = entry.params[0]
        assert p.description == "cluster label to highlight"

    def test_all_three_extras_together(self):
        def tool_full_extras(cell_type: str) -> None:
            pass

        extras = {
            "cell_type": {
                "enum": ["T cell", "B cell", "NK"],
                "field_type": "cell_type",
                "description": "target cell type",
            }
        }
        entry = _reg(tool_full_extras, params=extras)
        p = entry.params[0]
        assert p.enum == ["T cell", "B cell", "NK"]
        assert p.field_type == "cell_type"
        assert p.description == "target cell type"

    def test_extras_for_unknown_param_is_ignored(self):
        """Extras keyed on a param name not in the signature must not crash."""
        def tool_simple(x: str) -> None:
            pass

        extras = {"nonexistent": {"enum": ["a", "b"]}}
        # Should register without error; x gets no extras
        entry = _reg(tool_simple, params=extras)
        p = entry.params[0]
        assert p.name == "x"
        assert p.enum is None

    def test_missing_extras_key_uses_defaults(self):
        """Param with no extras entry gets None enum/field_type, empty description."""
        def tool_no_extras(x: str) -> None:
            pass

        entry = _reg(tool_no_extras)
        p = entry.params[0]
        assert p.enum is None
        assert p.field_type is None
        assert p.description == ""

    def test_partial_extras_only_enum(self):
        """Only 'enum' key in extras; field_type and description stay default."""
        def tool_partial(x: str) -> None:
            pass

        extras = {"x": {"enum": ["a", "b"]}}
        entry = _reg(tool_partial, params=extras)
        p = entry.params[0]
        assert p.enum == ["a", "b"]
        assert p.field_type is None
        assert p.description == ""


# ---------------------------------------------------------------------------
# get_tool
# ---------------------------------------------------------------------------


class TestGetTool:
    def test_hit_returns_tool_entry(self):
        def registered_fn(x: str) -> None:
            pass

        _reg(registered_fn)
        result = get_tool("registered_fn")
        assert isinstance(result, ToolEntry)
        assert result.name == "registered_fn"

    def test_miss_returns_none(self):
        result = get_tool("definitely_not_registered_xyzzy")
        assert result is None

    def test_hit_returns_correct_entry(self):
        def fn_a(a: str) -> None:
            pass

        def fn_b(b: int) -> None:
            pass

        _reg(fn_a, description="fn_a desc")
        _reg(fn_b, description="fn_b desc")
        assert get_tool("fn_a").description == "fn_a desc"
        assert get_tool("fn_b").description == "fn_b desc"


# ---------------------------------------------------------------------------
# get_tool_names
# ---------------------------------------------------------------------------


class TestGetToolNames:
    def test_empty_registry_returns_empty_list(self):
        reg_module.REGISTRY.clear()
        assert get_tool_names() == []

    def test_returns_registered_names(self):
        reg_module.REGISTRY.clear()

        def tool_x(x: str) -> None:
            pass

        def tool_y(y: int) -> None:
            pass

        _reg(tool_x)
        _reg(tool_y)
        names = get_tool_names()
        assert set(names) == {"tool_x", "tool_y"}

    def test_returns_list_type(self):
        assert isinstance(get_tool_names(), list)

    def test_count_matches_registry(self):
        reg_module.REGISTRY.clear()

        def fn1(x: str) -> None:
            pass

        def fn2(x: str) -> None:
            pass

        def fn3(x: str) -> None:
            pass

        _reg(fn1)
        _reg(fn2)
        _reg(fn3)
        assert len(get_tool_names()) == 3


# ---------------------------------------------------------------------------
# get_tool_description_block
# ---------------------------------------------------------------------------


class TestGetToolDescriptionBlock:
    def test_empty_registry_returns_empty_string(self):
        reg_module.REGISTRY.clear()
        assert get_tool_description_block() == ""

    def test_contains_tool_name_and_description(self):
        reg_module.REGISTRY.clear()

        def my_describe_tool(x: str) -> None:
            pass

        _reg(my_describe_tool, description="a helpful description")
        block = get_tool_description_block()
        assert "my_describe_tool" in block
        assert "a helpful description" in block

    def test_contains_param_name_and_type(self):
        reg_module.REGISTRY.clear()

        def tool_with_param(gene: str) -> None:
            pass

        _reg(tool_with_param)
        block = get_tool_description_block()
        assert "gene" in block
        assert "str" in block

    def test_required_param_shows_required(self):
        reg_module.REGISTRY.clear()

        def tool_req(gene: str) -> None:
            pass

        _reg(tool_req)
        block = get_tool_description_block()
        assert "required" in block

    def test_optional_param_shows_default(self):
        reg_module.REGISTRY.clear()

        def tool_opt(resolution: float = 0.5) -> None:
            pass

        _reg(tool_opt)
        block = get_tool_description_block()
        assert "default=" in block
        assert "0.5" in block

    def test_enum_appears_in_block(self):
        reg_module.REGISTRY.clear()

        def tool_enum_block(kind: str) -> None:
            pass

        extras = {"kind": {"enum": ["umap", "tsne", "pca"]}}
        _reg(tool_enum_block, params=extras)
        block = get_tool_description_block()
        assert "umap" in block

    def test_param_description_appears_in_block(self):
        reg_module.REGISTRY.clear()

        def tool_desc_block(gene: str) -> None:
            pass

        extras = {"gene": {"description": "gene symbol to plot"}}
        _reg(tool_desc_block, params=extras)
        block = get_tool_description_block()
        assert "gene symbol to plot" in block

    def test_no_param_tool_still_shows_tool_line(self):
        reg_module.REGISTRY.clear()

        def no_param_tool() -> None:
            pass

        _reg(no_param_tool, description="no params here")
        block = get_tool_description_block()
        assert "no_param_tool" in block
        assert "no params here" in block

    def test_multiple_tools_both_present(self):
        reg_module.REGISTRY.clear()

        def alpha_tool(x: str) -> None:
            pass

        def beta_tool(y: int) -> None:
            pass

        _reg(alpha_tool, description="alpha desc")
        _reg(beta_tool, description="beta desc")
        block = get_tool_description_block()
        assert "alpha_tool" in block
        assert "beta_tool" in block
        assert "alpha desc" in block
        assert "beta desc" in block

    def test_returns_string_type(self):
        assert isinstance(get_tool_description_block(), str)

    def test_no_enum_no_bracket_in_block(self):
        """When enum is None, the '[one of ...]' substring must not appear."""
        reg_module.REGISTRY.clear()

        def plain_tool(x: str) -> None:
            pass

        _reg(plain_tool)
        block = get_tool_description_block()
        assert "[one of" not in block

    def test_no_description_no_dash_in_block(self):
        """When param description is empty, the ' -- ' separator must not appear."""
        reg_module.REGISTRY.clear()

        def plain_tool2(x: str) -> None:
            pass

        _reg(plain_tool2)
        block = get_tool_description_block()
        # The em-dash separator only appears when description is non-empty
        assert " — " not in block


# ---------------------------------------------------------------------------
# get_default_value
# ---------------------------------------------------------------------------


class TestGetDefaultValue:
    def test_known_tool_known_param_with_default(self):
        def tool_dv(resolution: float = 0.8) -> None:
            pass

        _reg(tool_dv)
        assert get_default_value("tool_dv", "resolution") == 0.8

    def test_known_tool_required_param_returns_none(self):
        """Required params have default=None on ParamSpec (by construction)."""
        def tool_req_dv(gene: str) -> None:
            pass

        _reg(tool_req_dv)
        # default is None for required params — the function should return None
        assert get_default_value("tool_req_dv", "gene") is None

    def test_unknown_tool_returns_none(self):
        assert get_default_value("nonexistent_tool_xyzzy", "param") is None

    def test_unknown_param_returns_none(self):
        def tool_for_unknown_param(x: str) -> None:
            pass

        _reg(tool_for_unknown_param)
        assert get_default_value("tool_for_unknown_param", "nonexistent_param") is None

    def test_default_none_value_returned(self):
        def tool_none_dv(x: Optional[str] = None) -> None:
            pass

        _reg(tool_none_dv)
        result = get_default_value("tool_none_dv", "x")
        assert result is None

    def test_default_string_value(self):
        def tool_str_dv(method: str = "leiden") -> None:
            pass

        _reg(tool_str_dv)
        assert get_default_value("tool_str_dv", "method") == "leiden"

    def test_default_int_value(self):
        def tool_int_dv(n_neighbors: int = 15) -> None:
            pass

        _reg(tool_int_dv)
        assert get_default_value("tool_int_dv", "n_neighbors") == 15

    def test_default_bool_false(self):
        def tool_bool_dv(normalize: bool = False) -> None:
            pass

        _reg(tool_bool_dv)
        assert get_default_value("tool_bool_dv", "normalize") is False

    def test_multiple_params_correct_one_returned(self):
        def tool_multi(a: str = "alpha", b: int = 42) -> None:
            pass

        _reg(tool_multi)
        assert get_default_value("tool_multi", "a") == "alpha"
        assert get_default_value("tool_multi", "b") == 42


# ---------------------------------------------------------------------------
# Registry isolation check (meta-test)
# ---------------------------------------------------------------------------


class TestRegistryIsolation:
    def test_registrations_in_one_test_dont_bleed_into_next(self):
        """This test registers a tool; the next test confirms it's gone."""
        def canary_tool(x: str) -> None:
            pass

        _reg(canary_tool)
        assert "canary_tool" in reg_module.REGISTRY

    def test_canary_tool_not_present(self):
        """Canary registered in the previous test must be gone (fixture worked)."""
        assert "canary_tool" not in reg_module.REGISTRY


# ---------------------------------------------------------------------------
# Edge cases — func wrapper (hasattr fn, "func")
# ---------------------------------------------------------------------------


class TestFuncWrapperUnwrap:
    def test_functools_partial_like_wrapper(self):
        """register() unwraps fn.func when present (e.g. functools.partial).

        The ToolEntry should be keyed on the inner function's __name__.
        """
        import functools

        def inner_tool(x: str, y: int = 0) -> None:
            pass

        # partial has a .func attribute pointing to the wrapped callable
        wrapped = functools.partial(inner_tool, y=5)

        decorator = register(description="wrapped")
        decorator(wrapped)

        # Registry key is inner_tool.__name__ because raw = fn.func
        assert "inner_tool" in reg_module.REGISTRY
        entry = reg_module.REGISTRY["inner_tool"]
        assert entry.callable is inner_tool


# ---------------------------------------------------------------------------
# ParamSpec dataclass field defaults
# ---------------------------------------------------------------------------


class TestParamSpecDefaults:
    def test_paramspec_enum_defaults_none(self):
        p = ParamSpec(name="x", type="str", required=True)
        assert p.enum is None

    def test_paramspec_field_type_defaults_none(self):
        p = ParamSpec(name="x", type="str", required=True)
        assert p.field_type is None

    def test_paramspec_description_defaults_empty_string(self):
        p = ParamSpec(name="x", type="str", required=True)
        assert p.description == ""

    def test_paramspec_default_defaults_none(self):
        p = ParamSpec(name="x", type="str", required=True)
        assert p.default is None


# ---------------------------------------------------------------------------
# ToolEntry dataclass structure
# ---------------------------------------------------------------------------


class TestToolEntryStructure:
    def test_tool_entry_has_expected_fields(self):
        def entry_tool(x: str) -> None:
            pass

        entry = _reg(entry_tool)
        assert hasattr(entry, "name")
        assert hasattr(entry, "kind")
        assert hasattr(entry, "description")
        assert hasattr(entry, "params")
        assert hasattr(entry, "callable")

    def test_params_is_list(self):
        def list_tool(x: str) -> None:
            pass

        entry = _reg(list_tool)
        assert isinstance(entry.params, list)

    def test_params_elements_are_paramspec(self):
        def paramspec_tool(x: str) -> None:
            pass

        entry = _reg(paramspec_tool)
        for p in entry.params:
            assert isinstance(p, ParamSpec)
