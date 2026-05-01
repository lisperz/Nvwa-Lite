"""Re-export shim for legacy callers; deleted in Stage 3.

Canonical home is src/domain/qc_metrics.py since Stage 1 of the legacy
decoupling plan (local/strategy/legacy_decoupling_plan_2026-04-30.md).
This file exists so legacy src/agent/tools.py keeps working until Stage 3
deletes the legacy modules and this shim together.
"""

from src.domain.qc_metrics import (  # noqa: F401
    get_obs_column_statistics,
    resolve_qc_metric_column,
    summarize_qc_metrics,
)
