"""Spec and Canonicalization — LLM output contract for the prompt pipeline."""

from typing import Any, Optional

from pydantic import BaseModel


class Canonicalization(BaseModel):
    field: str
    raw: str
    canonical: str
    context: Optional[dict[str, Any]] = None


class Spec(BaseModel):
    scenario_id: str
    tool_name: str
    params: dict[str, Any] = {}
    pre_canonical_params: dict[str, Any] = {}
    canonicalizations_applied: list[Canonicalization] = []
