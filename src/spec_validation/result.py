"""ValidationResult and Issue — validator output contract."""

from typing import Literal, Optional

from pydantic import BaseModel


class Issue(BaseModel):
    field: str
    reason: Literal[
        "not_found",
        "ambiguous",
        "wrong_type",
        "missing",
        "wrong_field",
        "unknown_tool",
        "unknown_scenario",
    ]
    suggestions: list[str] = []


class ValidationResult(BaseModel):
    status: Literal["ok", "defaults_disclosed", "needs_input"]
    stage: Optional[Literal["resolver", "validator"]] = None
    issues: list[Issue] = []
