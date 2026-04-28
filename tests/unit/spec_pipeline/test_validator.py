"""Exhaustive-branch unit tests for src/spec_validation/validator.py.

Public API only: validate(spec) -> tuple[Spec, ValidationResult].
_type_matches is NOT imported directly.

Every documented branch is covered:
  - unknown_tool -> needs_input with unknown_tool Issue + suggestions=get_tool_names()
  - all required params present + correct types -> ok, no issues
  - required param missing -> missing Issue
  - present param wrong type (shallow) -> wrong_type Issue
  - Optional[X] accepts None and X
  - list[str] accepts list instances (generic subscript stripped)
  - Any type tag always passes
  - unknown custom-class type tag -> pass (no false positive)
  - int accepts int; bool quirk (bool is int subclass)
  - float accepts int (ints are accepted where floats expected)
  - spec returned is the SAME object (identity preserved)
  - canonicalizations_applied and pre_canonical_params unchanged
  - multiple issues accumulate in issues list

Bug policy: xfail + report; do NOT edit src/ files.
"""

from __future__ import annotations

from typing import Any, Optional

import pytest

import src.core.registry as reg_module
from src.core.spec import Spec, Canonicalization
from src.spec_validation.result import Issue, ValidationResult
from src.spec_validation.validator import validate
from src.core.registry import REGISTRY, get_tool_names, register


# ---------------------------------------------------------------------------
# REGISTRY isolation fixture
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def clean_registry():
    """Snapshot REGISTRY before each test; restore after so registrations don't leak."""
    snapshot = dict(reg_module.REGISTRY)
    yield
    reg_module.REGISTRY.clear()
    reg_module.REGISTRY.update(snapshot)


# ---------------------------------------------------------------------------
# Throwaway tool registration helpers
# ---------------------------------------------------------------------------


def _register(fn, *, description="test tool", kind="atomic", params=None):
    """Register fn and return the ToolEntry."""
    decorator = register(description=description, kind=kind, params=params or {})
    decorator(fn)
    return reg_module.REGISTRY[fn.__name__]


def _make_spec(tool_name: str, params: dict | None = None, **kwargs) -> Spec:
    """Construct a minimal valid Spec for testing."""
    return Spec(
        scenario_id="test_scenario",
        tool_name=tool_name,
        params=params or {},
        **kwargs,
    )


# ---------------------------------------------------------------------------
# TestUnknownTool — unknown_tool issue branch
# ---------------------------------------------------------------------------


class TestUnknownTool:
    def test_unknown_tool_status_is_needs_input(self):
        spec = _make_spec("definitely_not_registered_xyzzy_abc")
        _spec, result = validate(spec)
        assert result.status == "needs_input"

    def test_unknown_tool_stage_is_validator(self):
        spec = _make_spec("nonexistent_tool_123")
        _spec, result = validate(spec)
        assert result.stage == "validator"

    def test_unknown_tool_has_one_issue(self):
        spec = _make_spec("ghost_tool")
        _spec, result = validate(spec)
        assert len(result.issues) == 1

    def test_unknown_tool_issue_field_is_tool_name(self):
        spec = _make_spec("ghost_tool")
        _spec, result = validate(spec)
        assert result.issues[0].field == "tool_name"

    def test_unknown_tool_issue_reason_is_unknown_tool(self):
        spec = _make_spec("ghost_tool")
        _spec, result = validate(spec)
        assert result.issues[0].reason == "unknown_tool"

    def test_unknown_tool_suggestions_equals_get_tool_names(self):
        """Suggestions list must equal the current registry's tool names."""
        # Register a throwaway tool so registry is non-empty
        def anchor_tool(x: str) -> None:
            pass

        _register(anchor_tool)

        spec = _make_spec("not_real_tool")
        _spec, result = validate(spec)
        assert set(result.issues[0].suggestions) == set(get_tool_names())

    def test_unknown_tool_suggestions_empty_when_registry_empty(self):
        reg_module.REGISTRY.clear()
        spec = _make_spec("any_tool")
        _spec, result = validate(spec)
        assert result.issues[0].suggestions == []

    def test_unknown_tool_returns_same_spec_object(self):
        spec = _make_spec("ghost_tool")
        returned_spec, _result = validate(spec)
        assert returned_spec is spec


# ---------------------------------------------------------------------------
# TestOkPath — no issues, status "ok"
# ---------------------------------------------------------------------------


class TestOkPath:
    def test_no_params_tool_ok(self):
        def no_param_tool() -> None:
            pass

        _register(no_param_tool)
        spec = _make_spec("no_param_tool")
        _spec, result = validate(spec)
        assert result.status == "ok"
        assert result.issues == []

    def test_all_required_params_present_and_correct_type_ok(self):
        def my_tool(gene: str, n: int) -> None:
            pass

        _register(my_tool)
        spec = _make_spec("my_tool", params={"gene": "BRCA1", "n": 10})
        _spec, result = validate(spec)
        assert result.status == "ok"
        assert result.issues == []

    def test_optional_param_absent_is_ok(self):
        """Missing optional param with default is NOT an issue — validator defers to dispatch."""
        def tool_with_optional(gene: str, resolution: float = 0.5) -> None:
            pass

        _register(tool_with_optional)
        spec = _make_spec("tool_with_optional", params={"gene": "BRCA1"})
        _spec, result = validate(spec)
        assert result.status == "ok"

    def test_stage_is_validator_on_ok(self):
        def simple_ok_tool(x: str) -> None:
            pass

        _register(simple_ok_tool)
        spec = _make_spec("simple_ok_tool", params={"x": "hello"})
        _spec, result = validate(spec)
        assert result.stage == "validator"

    def test_ok_returns_empty_issues_list(self):
        def ok_tool(x: str) -> None:
            pass

        _register(ok_tool)
        spec = _make_spec("ok_tool", params={"x": "hello"})
        _spec, result = validate(spec)
        assert result.issues == []


# ---------------------------------------------------------------------------
# TestMissingIssue — required param absent
# ---------------------------------------------------------------------------


class TestMissingIssue:
    def test_missing_required_param_status_needs_input(self):
        def req_tool(gene: str) -> None:
            pass

        _register(req_tool)
        spec = _make_spec("req_tool", params={})
        _spec, result = validate(spec)
        assert result.status == "needs_input"

    def test_missing_issue_reason_is_missing(self):
        def req_tool2(gene: str) -> None:
            pass

        _register(req_tool2)
        spec = _make_spec("req_tool2", params={})
        _spec, result = validate(spec)
        issues = [i for i in result.issues if i.reason == "missing"]
        assert len(issues) == 1

    def test_missing_issue_field_has_params_prefix(self):
        def req_tool3(cluster: str) -> None:
            pass

        _register(req_tool3)
        spec = _make_spec("req_tool3", params={})
        _spec, result = validate(spec)
        issue = result.issues[0]
        assert issue.field == "params.cluster"

    def test_two_required_params_both_missing_two_issues(self):
        def two_req_tool(gene: str, cluster: str) -> None:
            pass

        _register(two_req_tool)
        spec = _make_spec("two_req_tool", params={})
        _spec, result = validate(spec)
        missing_issues = [i for i in result.issues if i.reason == "missing"]
        assert len(missing_issues) == 2

    def test_missing_one_of_two_required_params_one_issue(self):
        def partial_tool(gene: str, cluster: str) -> None:
            pass

        _register(partial_tool)
        spec = _make_spec("partial_tool", params={"gene": "BRCA1"})
        _spec, result = validate(spec)
        missing_issues = [i for i in result.issues if i.reason == "missing"]
        assert len(missing_issues) == 1
        assert missing_issues[0].field == "params.cluster"

    def test_missing_issue_suggestions_is_empty_list(self):
        """missing issues carry no suggestions (field=value would be too specific)."""
        def miss_sug_tool(gene: str) -> None:
            pass

        _register(miss_sug_tool)
        spec = _make_spec("miss_sug_tool", params={})
        _spec, result = validate(spec)
        issue = result.issues[0]
        assert issue.suggestions == []

    def test_optional_param_absent_no_missing_issue(self):
        def opt_only_tool(resolution: float = 0.5) -> None:
            pass

        _register(opt_only_tool)
        spec = _make_spec("opt_only_tool", params={})
        _spec, result = validate(spec)
        missing_issues = [i for i in result.issues if i.reason == "missing"]
        assert missing_issues == []


# ---------------------------------------------------------------------------
# TestWrongTypeIssue — present param fails shallow isinstance
# ---------------------------------------------------------------------------


class TestWrongTypeIssue:
    def test_wrong_type_str_param_given_int_status_needs_input(self):
        def str_tool(name: str) -> None:
            pass

        _register(str_tool)
        spec = _make_spec("str_tool", params={"name": 42})
        _spec, result = validate(spec)
        assert result.status == "needs_input"

    def test_wrong_type_issue_reason_is_wrong_type(self):
        def wt_tool(name: str) -> None:
            pass

        _register(wt_tool)
        spec = _make_spec("wt_tool", params={"name": 42})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert len(wt_issues) == 1

    def test_wrong_type_issue_field_has_params_prefix(self):
        def wt_field_tool(cluster: str) -> None:
            pass

        _register(wt_field_tool)
        spec = _make_spec("wt_field_tool", params={"cluster": 99})
        _spec, result = validate(spec)
        issue = result.issues[0]
        assert issue.field == "params.cluster"

    def test_wrong_type_suggestions_contains_type_tag(self):
        """wrong_type suggestions = [p.type] (the registry type tag string)."""
        def wt_sug_tool(gene: str) -> None:
            pass

        _register(wt_sug_tool)
        spec = _make_spec("wt_sug_tool", params={"gene": 123})
        _spec, result = validate(spec)
        issue = result.issues[0]
        assert "str" in issue.suggestions

    def test_str_value_for_int_param_wrong_type(self):
        def int_tool(n: int) -> None:
            pass

        _register(int_tool)
        spec = _make_spec("int_tool", params={"n": "not_an_int"})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert len(wt_issues) == 1

    def test_str_value_for_list_param_wrong_type(self):
        def list_tool(items: list) -> None:
            pass

        _register(list_tool)
        spec = _make_spec("list_tool", params={"items": "not_a_list"})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert len(wt_issues) == 1

    def test_list_for_dict_param_wrong_type(self):
        def dict_tool(config: dict) -> None:
            pass

        _register(dict_tool)
        spec = _make_spec("dict_tool", params={"config": ["a", "b"]})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert len(wt_issues) == 1


# ---------------------------------------------------------------------------
# TestOptionalTypeTag — Optional[X] unwrapping
# ---------------------------------------------------------------------------


class TestOptionalTypeTag:
    def test_optional_str_accepts_none(self):
        def opt_str_tool(label: Optional[str] = None) -> None:
            pass

        _register(opt_str_tool)
        spec = _make_spec("opt_str_tool", params={"label": None})
        _spec, result = validate(spec)
        assert result.status == "ok"
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert wt_issues == []

    def test_optional_str_accepts_str(self):
        def opt_str_tool2(label: Optional[str] = None) -> None:
            pass

        _register(opt_str_tool2)
        spec = _make_spec("opt_str_tool2", params={"label": "hello"})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert wt_issues == []

    def test_optional_str_rejects_int(self):
        """Optional[str] unwraps to str; an int value must fail the str check."""
        def opt_str_tool3(label: Optional[str] = None) -> None:
            pass

        _register(opt_str_tool3)
        spec = _make_spec("opt_str_tool3", params={"label": 99})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert len(wt_issues) == 1

    def test_optional_int_accepts_none(self):
        def opt_int_tool(n: Optional[int] = None) -> None:
            pass

        _register(opt_int_tool)
        spec = _make_spec("opt_int_tool", params={"n": None})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert wt_issues == []

    def test_optional_int_accepts_int(self):
        def opt_int_tool2(n: Optional[int] = None) -> None:
            pass

        _register(opt_int_tool2)
        spec = _make_spec("opt_int_tool2", params={"n": 5})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert wt_issues == []


# ---------------------------------------------------------------------------
# TestGenericSubscriptStripping — list[str] -> list, dict[str,int] -> dict
# ---------------------------------------------------------------------------


class TestGenericSubscriptStripping:
    def test_list_str_annotated_param_accepts_list(self):
        """list[str] type tag strips to list; a list value must pass."""
        def list_str_tool(genes: list[str]) -> None:
            pass

        _register(list_str_tool)
        spec = _make_spec("list_str_tool", params={"genes": ["BRCA1", "TP53"]})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert wt_issues == []

    def test_list_str_annotated_param_rejects_str(self):
        """list[str] strips to list; a bare str must fail."""
        def list_str_tool2(genes: list[str]) -> None:
            pass

        _register(list_str_tool2)
        spec = _make_spec("list_str_tool2", params={"genes": "BRCA1"})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert len(wt_issues) == 1

    def test_list_int_annotated_param_accepts_list(self):
        def list_int_tool(counts: list[int]) -> None:
            pass

        _register(list_int_tool)
        spec = _make_spec("list_int_tool", params={"counts": [1, 2, 3]})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert wt_issues == []

    def test_bare_list_annotated_param_accepts_list(self):
        def bare_list_tool(items: list) -> None:
            pass

        _register(bare_list_tool)
        spec = _make_spec("bare_list_tool", params={"items": [1, 2, 3]})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert wt_issues == []


# ---------------------------------------------------------------------------
# TestAnyTypeTag — Any always passes
# ---------------------------------------------------------------------------


class TestAnyTypeTag:
    def test_any_accepts_str(self):
        def any_tool(value) -> None:
            pass

        _register(any_tool)
        spec = _make_spec("any_tool", params={"value": "hello"})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert wt_issues == []

    def test_any_accepts_int(self):
        def any_tool2(value) -> None:
            pass

        _register(any_tool2)
        spec = _make_spec("any_tool2", params={"value": 42})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert wt_issues == []

    def test_any_accepts_none(self):
        def any_tool3(value) -> None:
            pass

        _register(any_tool3)
        spec = _make_spec("any_tool3", params={"value": None})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert wt_issues == []

    def test_any_accepts_list(self):
        def any_tool4(value) -> None:
            pass

        _register(any_tool4)
        spec = _make_spec("any_tool4", params={"value": [1, 2, 3]})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert wt_issues == []

    def test_any_accepts_dict(self):
        def any_tool5(value) -> None:
            pass

        _register(any_tool5)
        spec = _make_spec("any_tool5", params={"value": {"key": "val"}})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert wt_issues == []


# ---------------------------------------------------------------------------
# TestUnknownCustomTypeTag — unknown tags default to pass (no false positive)
# ---------------------------------------------------------------------------


class TestUnknownCustomTypeTag:
    def test_custom_class_tag_does_not_produce_wrong_type(self):
        """A type tag the validator doesn't recognize must not false-positive."""
        class MySpecialClass:
            pass

        def custom_tag_tool(obj: MySpecialClass) -> None:
            pass

        _register(custom_tag_tool)
        # Pass a plain dict — validator must not know what MySpecialClass is
        spec = _make_spec("custom_tag_tool", params={"obj": {"key": "value"}})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert wt_issues == []

    def test_unknown_string_tag_does_not_produce_wrong_type(self):
        """If someone injects a weird type tag string, validator must still pass."""
        from src.core.registry import ParamSpec, ToolEntry
        import src.core.registry as reg_module

        def stub_fn() -> None:
            pass

        entry = ToolEntry(
            name="weird_tag_tool",
            kind="atomic",
            description="test",
            params=[ParamSpec(name="x", type="WeirdUnknownType", required=True)],
            callable=stub_fn,
        )
        reg_module.REGISTRY["weird_tag_tool"] = entry
        spec = _make_spec("weird_tag_tool", params={"x": "anything"})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert wt_issues == []


# ---------------------------------------------------------------------------
# TestBoolAndIntQuirks — bool is int subclass, float accepts int
# ---------------------------------------------------------------------------


class TestBoolAndIntQuirks:
    def test_int_param_accepts_int(self):
        def int_tool2(n: int) -> None:
            pass

        _register(int_tool2)
        spec = _make_spec("int_tool2", params={"n": 10})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert wt_issues == []

    def test_bool_param_accepts_true(self):
        """bool is the expected type; True passes isinstance(True, (bool,))."""
        def bool_tool(flag: bool) -> None:
            pass

        _register(bool_tool)
        spec = _make_spec("bool_tool", params={"flag": True})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert wt_issues == []

    def test_bool_param_accepts_false(self):
        def bool_tool2(flag: bool) -> None:
            pass

        _register(bool_tool2)
        spec = _make_spec("bool_tool2", params={"flag": False})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert wt_issues == []

    def test_int_param_accepts_bool_because_bool_is_int_subclass(self):
        """bool is an int subclass — isinstance(True, (int,)) is True.
        Validator source comment: 'bool is an int subclass; accept intentionally'.
        This is the documented quirk — NOT a bug.
        """
        def int_tool3(n: int) -> None:
            pass

        _register(int_tool3)
        spec = _make_spec("int_tool3", params={"n": True})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        # bool IS an int subclass; the validator accepts this intentionally
        assert wt_issues == []

    def test_bool_param_rejects_int_zero(self):
        """Bool param + int 0 is rejected as wrong_type.

        Asymmetric with the int-accepts-bool direction: the source comment on
        validator.py's simple_types map marks bool as accepted for int params
        ("bool is an int subclass"), but the reverse is NOT true — bool params
        use `(bool,)` isinstance, and `isinstance(0, (bool,))` is False.
        Callers passing 0/1 for a bool param will get a wrong_type issue.
        """
        def bool_only_tool(flag: bool) -> None:
            pass

        _register(bool_only_tool)
        spec = _make_spec("bool_only_tool", params={"flag": 0})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert len(wt_issues) == 1
        assert wt_issues[0].field == "params.flag"

    def test_float_param_accepts_float(self):
        def float_tool(resolution: float) -> None:
            pass

        _register(float_tool)
        spec = _make_spec("float_tool", params={"resolution": 0.5})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert wt_issues == []

    def test_float_param_accepts_int(self):
        """ints are accepted where floats are expected — documented in source."""
        def float_tool2(resolution: float) -> None:
            pass

        _register(float_tool2)
        spec = _make_spec("float_tool2", params={"resolution": 1})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert wt_issues == []

    def test_float_param_rejects_str(self):
        def float_tool3(resolution: float) -> None:
            pass

        _register(float_tool3)
        spec = _make_spec("float_tool3", params={"resolution": "high"})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert len(wt_issues) == 1


# ---------------------------------------------------------------------------
# TestSpecIdentityPreservation — validator never mutates spec
# ---------------------------------------------------------------------------


class TestSpecIdentityPreservation:
    def test_returned_spec_is_same_object_ok_path(self):
        def identity_ok_tool(x: str) -> None:
            pass

        _register(identity_ok_tool)
        spec = _make_spec("identity_ok_tool", params={"x": "hello"})
        returned_spec, _ = validate(spec)
        assert returned_spec is spec

    def test_returned_spec_is_same_object_needs_input_path(self):
        """Even when validation fails, the same spec object is returned."""
        def identity_err_tool(x: str) -> None:
            pass

        _register(identity_err_tool)
        spec = _make_spec("identity_err_tool", params={})  # missing required x
        returned_spec, _ = validate(spec)
        assert returned_spec is spec

    def test_returned_spec_is_same_object_unknown_tool(self):
        spec = _make_spec("ghost_tool_for_identity")
        returned_spec, _ = validate(spec)
        assert returned_spec is spec

    def test_spec_params_not_mutated(self):
        def no_mutate_tool(x: str) -> None:
            pass

        _register(no_mutate_tool)
        original_params = {"x": "hello"}
        spec = _make_spec("no_mutate_tool", params=dict(original_params))
        validate(spec)
        assert spec.params == original_params

    def test_spec_tool_name_not_mutated(self):
        def no_mutate_tn_tool(x: str) -> None:
            pass

        _register(no_mutate_tn_tool)
        spec = _make_spec("no_mutate_tn_tool", params={"x": "hi"})
        validate(spec)
        assert spec.tool_name == "no_mutate_tn_tool"

    def test_spec_scenario_id_not_mutated(self):
        def no_mutate_sid_tool(x: str) -> None:
            pass

        _register(no_mutate_sid_tool)
        spec = Spec(scenario_id="my_scenario", tool_name="no_mutate_sid_tool", params={"x": "hi"})
        validate(spec)
        assert spec.scenario_id == "my_scenario"

    def test_canonicalizations_applied_not_mutated(self):
        def canon_tool(gene: str) -> None:
            pass

        _register(canon_tool)
        canon = Canonicalization(field="gene", raw="brca1", canonical="BRCA1")
        spec = Spec(
            scenario_id="test",
            tool_name="canon_tool",
            params={"gene": "BRCA1"},
            canonicalizations_applied=[canon],
        )
        validate(spec)
        assert spec.canonicalizations_applied == [canon]

    def test_pre_canonical_params_not_mutated(self):
        def precanon_tool(gene: str) -> None:
            pass

        _register(precanon_tool)
        pre = {"gene": "brca1"}
        spec = Spec(
            scenario_id="test",
            tool_name="precanon_tool",
            params={"gene": "BRCA1"},
            pre_canonical_params=dict(pre),
        )
        validate(spec)
        assert spec.pre_canonical_params == pre


# ---------------------------------------------------------------------------
# TestMultipleIssuesAccumulate — both missing + wrong_type accumulate
# ---------------------------------------------------------------------------


class TestMultipleIssuesAccumulate:
    def test_two_missing_params_two_issues(self):
        def multi_miss_tool(a: str, b: int, c: float = 0.5) -> None:
            pass

        _register(multi_miss_tool)
        spec = _make_spec("multi_miss_tool", params={})
        _spec, result = validate(spec)
        missing_issues = [i for i in result.issues if i.reason == "missing"]
        assert len(missing_issues) == 2

    def test_missing_and_wrong_type_both_accumulate(self):
        """One param missing, one present but wrong type -> both issues reported."""
        def mix_issue_tool(gene: str, n: int) -> None:
            pass

        _register(mix_issue_tool)
        # n is present but wrong type; gene is absent
        spec = _make_spec("mix_issue_tool", params={"n": "not_an_int"})
        _spec, result = validate(spec)
        reasons = {i.reason for i in result.issues}
        assert "missing" in reasons
        assert "wrong_type" in reasons
        assert len(result.issues) == 2

    def test_two_wrong_type_params_two_issues(self):
        def two_wt_tool(a: str, b: int) -> None:
            pass

        _register(two_wt_tool)
        spec = _make_spec("two_wt_tool", params={"a": 99, "b": "not_int"})
        _spec, result = validate(spec)
        wt_issues = [i for i in result.issues if i.reason == "wrong_type"]
        assert len(wt_issues) == 2

    def test_multiple_issues_status_still_needs_input(self):
        def multi_status_tool(a: str, b: int) -> None:
            pass

        _register(multi_status_tool)
        spec = _make_spec("multi_status_tool", params={"a": 123, "b": "bad"})
        _spec, result = validate(spec)
        assert result.status == "needs_input"


# ---------------------------------------------------------------------------
# TestValidationResultStructure — result is a proper ValidationResult
# ---------------------------------------------------------------------------


class TestValidationResultStructure:
    def test_result_is_validation_result_instance(self):
        def struct_tool(x: str) -> None:
            pass

        _register(struct_tool)
        spec = _make_spec("struct_tool", params={"x": "hello"})
        _spec, result = validate(spec)
        assert isinstance(result, ValidationResult)

    def test_issues_is_list_of_issue_instances(self):
        def issue_list_tool(gene: str) -> None:
            pass

        _register(issue_list_tool)
        spec = _make_spec("issue_list_tool", params={})
        _spec, result = validate(spec)
        assert isinstance(result.issues, list)
        for issue in result.issues:
            assert isinstance(issue, Issue)

    def test_returns_tuple(self):
        def tuple_tool(x: str) -> None:
            pass

        _register(tuple_tool)
        spec = _make_spec("tuple_tool", params={"x": "hi"})
        output = validate(spec)
        assert isinstance(output, tuple)
        assert len(output) == 2

    def test_defaults_disclosed_never_emitted(self):
        """Validator must never return status='defaults_disclosed'."""
        def defaults_tool(x: str, y: int = 5) -> None:
            pass

        _register(defaults_tool)
        spec = _make_spec("defaults_tool", params={"x": "hello"})
        _spec, result = validate(spec)
        assert result.status != "defaults_disclosed"

    def test_stage_always_validator(self):
        """All results from validate() must have stage='validator'."""
        def stage_check_tool(x: str) -> None:
            pass

        _register(stage_check_tool)

        # ok path
        spec_ok = _make_spec("stage_check_tool", params={"x": "hi"})
        _, r1 = validate(spec_ok)
        assert r1.stage == "validator"

        # needs_input path
        spec_bad = _make_spec("stage_check_tool", params={})
        _, r2 = validate(spec_bad)
        assert r2.stage == "validator"

        # unknown tool path
        spec_unk = _make_spec("totally_unknown_tool_zzz")
        _, r3 = validate(spec_unk)
        assert r3.stage == "validator"


# ---------------------------------------------------------------------------
# TestCanonExtraFieldsUnchanged — canonicalizations_applied / pre_canonical_params
# ---------------------------------------------------------------------------


class TestCanonExtraFieldsUnchanged:
    def test_empty_canonicalizations_applied_stays_empty(self):
        def canon_empty_tool(x: str) -> None:
            pass

        _register(canon_empty_tool)
        spec = _make_spec("canon_empty_tool", params={"x": "val"})
        returned_spec, _ = validate(spec)
        assert returned_spec.canonicalizations_applied == []

    def test_non_empty_canonicalizations_applied_unchanged(self):
        def canon_nonempty_tool(gene: str) -> None:
            pass

        _register(canon_nonempty_tool)
        canon = Canonicalization(field="gene", raw="brca1", canonical="BRCA1")
        spec = Spec(
            scenario_id="test",
            tool_name="canon_nonempty_tool",
            params={"gene": "BRCA1"},
            canonicalizations_applied=[canon],
        )
        returned_spec, _ = validate(spec)
        assert len(returned_spec.canonicalizations_applied) == 1
        assert returned_spec.canonicalizations_applied[0].raw == "brca1"
        assert returned_spec.canonicalizations_applied[0].canonical == "BRCA1"

    def test_pre_canonical_params_empty_stays_empty(self):
        def pre_empty_tool(x: str) -> None:
            pass

        _register(pre_empty_tool)
        spec = _make_spec("pre_empty_tool", params={"x": "val"})
        returned_spec, _ = validate(spec)
        assert returned_spec.pre_canonical_params == {}

    def test_pre_canonical_params_non_empty_unchanged(self):
        def pre_nonempty_tool(gene: str) -> None:
            pass

        _register(pre_nonempty_tool)
        spec = Spec(
            scenario_id="test",
            tool_name="pre_nonempty_tool",
            params={"gene": "BRCA1"},
            pre_canonical_params={"gene": "brca1"},
        )
        returned_spec, _ = validate(spec)
        assert returned_spec.pre_canonical_params == {"gene": "brca1"}
