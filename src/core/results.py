"""Tool-output schema for the spec pipeline.

Discriminated union: TextResult (text-only) | ArtifactResult (image/csv + metadata).
Replaces the legacy str-prefix "Error:" status convention for @register tools;
legacy LangChain-loop tools in src/agent/tools.py still use the str convention
(_classify_tool_status in src/agent/core.py keeps handling that path).

Bytes carried on ArtifactResult are bridged into the legacy artifact channel
by src/agent/artifacts.consume_artifact_result, so UI consumption in
src/ui/app.py stays unchanged during migration.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

ToolStatus = Literal["success", "error", "warning"]
ArtifactKind = Literal["image", "csv"]


@dataclass(frozen=True)
class TextResult:
    """Tool returns this when output is text-only (no plot, no table)."""

    text: str
    status: ToolStatus = "success"
    error_message: str | None = None


@dataclass(frozen=True)
class ArtifactResult:
    """Tool returns this when output includes an image or csv artifact.

    Field population by artifact_kind:
      - "image": image_bytes required; code optional (snippet for UI display).
      - "csv":   csv_data required; display_df optional (markdown for UI); code optional.

    text is the responder/gatekeeper-facing narrative — appears in chat and is
    what the gatekeeper scans for honesty checks.
    """

    text: str
    artifact_kind: ArtifactKind
    tool_name: str
    params_used: dict
    entities_acted_on: list[str] = field(default_factory=list)
    image_bytes: bytes | None = None
    csv_data: str | None = None
    display_df: str | None = None
    code: str | None = None
    status: ToolStatus = "success"
    error_message: str | None = None


ToolResult = TextResult | ArtifactResult


class ToolExecutionError(Exception):
    """Raise from inside a @register tool body to signal a structured failure.

    The dispatcher (AgentRunner._dispatch_registered) catches this and wraps to
    a status='error' tool_output string. Use this instead of returning an
    "Error: ..." string — the latter is the legacy LangChain-loop convention
    and is not recognized by the new spec-pipeline contract.
    """

    def __init__(self, message: str, *, tool_name: str | None = None):
        super().__init__(message)
        self.tool_name = tool_name
