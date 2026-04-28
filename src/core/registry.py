"""Tool registry — self-registering catalog for the spec pipeline.

Lives in src/core/ alongside spec.py / results.py / types.py because it's pure
infrastructure (a dict of ToolEntry dataclasses + a decorator that populates it).
Bio-knowledge tool bodies in src/domain/<section>.py import @register from here.
"""

import inspect
from dataclasses import dataclass
from typing import Any, Callable, Literal, Optional


@dataclass
class ParamSpec:
    name: str
    type: str
    required: bool
    default: Any = None
    enum: Optional[list[Any]] = None
    field_type: Optional[str] = None
    description: str = ""


@dataclass
class ToolEntry:
    name: str
    kind: Literal["atomic", "workflow"]
    description: str
    params: list[ParamSpec]
    callable: Callable


REGISTRY: dict[str, ToolEntry] = {}


def register(
    *,
    description: str,
    kind: Literal["atomic", "workflow"] = "atomic",
    params: Optional[dict[str, dict]] = None,
) -> Callable:
    def decorator(fn: Callable) -> Callable:
        raw = fn.func if hasattr(fn, "func") else fn
        REGISTRY[raw.__name__] = ToolEntry(
            name=raw.__name__,
            kind=kind,
            description=description,
            params=_build_param_specs(raw, params or {}),
            callable=raw,
        )
        return fn
    return decorator


def _build_param_specs(fn: Callable, extras: dict[str, dict]) -> list[ParamSpec]:
    specs = []
    for name, param in inspect.signature(fn).parameters.items():
        if param.kind in (inspect.Parameter.VAR_POSITIONAL, inspect.Parameter.VAR_KEYWORD):
            continue
        # "adata" is injected by the dispatcher from session state; not an
        # extractor-facing param. Skipping it keeps tool catalog + JSON schema clean.
        if name in ("self", "cls", "adata"):
            continue
        extra = extras.get(name, {})
        specs.append(ParamSpec(
            name=name,
            type=_type_tag(param.annotation),
            required=param.default is inspect.Parameter.empty,
            default=None if param.default is inspect.Parameter.empty else param.default,
            enum=extra.get("enum"),
            field_type=extra.get("field_type"),
            description=extra.get("description", ""),
        ))
    return specs


def _type_tag(ann: Any) -> str:
    if ann is inspect.Parameter.empty:
        return "Any"
    tag = str(ann).replace("typing.", "")
    if tag.startswith("<class '") and tag.endswith("'>"):
        tag = tag[len("<class '"):-2]
    return tag


def get_tool(name: str) -> Optional[ToolEntry]:
    return REGISTRY.get(name)


def get_tool_names() -> list[str]:
    return list(REGISTRY.keys())


def get_tool_description_block() -> str:
    lines = []
    for entry in REGISTRY.values():
        lines.append(f"{entry.name}: {entry.description}")
        for p in entry.params:
            req = "required" if p.required else f"default={p.default!r}"
            enum = f" [one of {p.enum}]" if p.enum else ""
            desc = f" — {p.description}" if p.description else ""
            lines.append(f"  {p.name}: {p.type}, {req}{enum}{desc}")
    return "\n".join(lines)


def get_default_value(tool_name: str, param: str) -> Any:
    entry = REGISTRY.get(tool_name)
    if entry is None:
        return None
    for p in entry.params:
        if p.name == param:
            return p.default
    return None
