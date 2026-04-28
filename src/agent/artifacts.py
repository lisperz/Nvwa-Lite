"""New-path artifact buffer + ArtifactResult → legacy-shape bridge.

Owns the singleton list for artifacts produced by spec-pipeline @register tools.
Lives separate from the legacy src/agent/tools.py singletons so that file can be
deleted entirely when L4 migration completes — this module is the long-term home.

UI ([src/ui/app.py]) aggregates this module's getters with the legacy
[src/agent/tools.py] getters during migration. After legacy tools are deleted,
the legacy getters go away and this becomes the sole artifact channel.
"""

from __future__ import annotations

from src.core.results import ArtifactResult
from src.domain.plotting.executor import PlotResult, TableResult


_plot_results: list[PlotResult] = []
_table_results: list[TableResult] = []


def consume_artifact_result(result: ArtifactResult) -> None:
    """Translate a new-path ArtifactResult into the legacy singleton format.

    Called by AgentRunner._dispatch_registered after a @register tool returns
    an ArtifactResult. UI continues to read get_plot_results() / get_table_results()
    (aggregated across legacy + this module) without knowing which producer wrote.
    """
    if result.artifact_kind == "image":
        _plot_results.append(PlotResult(
            image=result.image_bytes or b"",
            code=result.code or "",
            message=result.text,
        ))
    elif result.artifact_kind == "csv":
        _table_results.append(TableResult(
            csv_data=result.csv_data or "",
            code=result.code or "",
            message=result.text,
            display_df=result.display_df or "",
        ))


def get_plot_results() -> list[PlotResult]:
    """Return all PlotResults produced by spec-pipeline tools this turn."""
    return list(_plot_results)


def get_table_results() -> list[TableResult]:
    """Return all TableResults produced by spec-pipeline tools this turn."""
    return list(_table_results)


def clear() -> None:
    """Clear both buffers between turns. Called from src/ui/app.py."""
    _plot_results.clear()
    _table_results.clear()
