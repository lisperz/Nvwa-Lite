"""Upload-time classifier for adata.obs column roles.

Writes adata.uns["nvwa_meta"] with per-column role labels, derived condition_cols,
ambiguous_cols, and a user_resolutions dict (populated externally by UI/Redis).
Heuristic-only: no LLM. Columns the heuristic cannot confidently label are marked
AMBIGUOUS with candidates, surfaced to the user via dataset_overview().

Preserves any keys already present in adata.uns["nvwa_meta"] (e.g. species written
by species_detector.py running earlier in the load_dataset path) and any existing
user_resolutions. Overwrites only schema / condition_cols / ambiguous_cols.
"""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass
from enum import Enum
from typing import TYPE_CHECKING, Any, Optional

import pandas as pd

if TYPE_CHECKING:
    from anndata import AnnData

logger = logging.getLogger(__name__)


class Role(str, Enum):
    CONDITION = "condition"
    CELL_TYPE = "cell_type"
    CLUSTERING = "clustering"
    QC_METRIC = "qc_metric"
    SAMPLE_ID = "sample_id"
    BATCH = "batch"
    AMBIGUOUS = "ambiguous"
    OTHER = "other"


@dataclass
class SchemaEntry:
    role: Role
    confidence: str  # "high" | "low"
    n_unique: int
    values: Optional[list[str]] = None
    candidates: Optional[list[str]] = None


_KNOWN_NAMES: dict[str, Role] = {
    # Cell type
    "cell_type": Role.CELL_TYPE,
    "celltype": Role.CELL_TYPE,
    "cell_types": Role.CELL_TYPE,
    # Clustering
    "leiden": Role.CLUSTERING,
    "louvain": Role.CLUSTERING,
    "cluster": Role.CLUSTERING,
    "clusters": Role.CLUSTERING,
    "seurat_clusters": Role.CLUSTERING,
    # Batch
    "batch": Role.BATCH,
    # Sample id
    "sample_id": Role.SAMPLE_ID,
    "donor": Role.SAMPLE_ID,
    "donor_id": Role.SAMPLE_ID,
    "subject": Role.SAMPLE_ID,
    "subject_id": Role.SAMPLE_ID,
    "orig.ident": Role.SAMPLE_ID,
    # Condition
    "genotype": Role.CONDITION,
    "treatment": Role.CONDITION,
    "condition": Role.CONDITION,
    "timepoint": Role.CONDITION,
    "time_point": Role.CONDITION,
    "stimulus": Role.CONDITION,
    "cohort": Role.CONDITION,
    "sex": Role.CONDITION,
    "phase": Role.CONDITION,  # Seurat CellCycleScoring output
}

_AMBIGUOUS_HINTS: set[str] = {"sample", "group", "ident"}

_QC_PATTERNS: tuple[str, ...] = (
    "ncount_",
    "nfeature_",
    "percent.mt",
    "pct_counts_",
    "total_counts",
    "n_genes",
    "n_counts",
)

_CLUSTERING_PATTERNS: tuple[str, ...] = ("_snn_res.",)

_LOW_CARDINALITY_THRESHOLD: int = 20

_MAX_VALUES_PREVIEW: int = 20


def classify(adata: "AnnData") -> None:
    """Classify adata.obs columns; merge result into adata.uns["nvwa_meta"].

    Idempotent. Preserves existing nvwa_meta keys outside this module's scope
    (species, user_resolutions). Overwrites schema, condition_cols, ambiguous_cols.
    """
    schema: dict[str, dict] = {}
    for col in adata.obs.columns:
        entry = _classify_one(str(col), adata.obs[col])
        schema[str(col)] = _entry_to_dict(entry)

    condition_cols = [
        c for c, e in schema.items()
        if e["role"] == Role.CONDITION.value and e["confidence"] == "high"
    ]
    ambiguous_cols = [
        c for c, e in schema.items() if e["role"] == Role.AMBIGUOUS.value
    ]

    merged = _existing_meta(adata)
    merged["schema"] = schema
    merged["condition_cols"] = condition_cols
    merged["ambiguous_cols"] = ambiguous_cols
    merged.setdefault("user_resolutions", {})

    adata.uns["nvwa_meta"] = merged
    logger.info(
        "column_classifier: %d cols classified (%d condition, %d ambiguous)",
        len(schema),
        len(condition_cols),
        len(ambiguous_cols),
    )


def _classify_one(name: str, series: pd.Series) -> SchemaEntry:
    n_unique = int(series.nunique(dropna=True))
    name_lower = name.lower()

    # 1. Ambiguous name hints — always ask user
    if name_lower in _AMBIGUOUS_HINTS:
        return SchemaEntry(
            role=Role.AMBIGUOUS,
            confidence="low",
            n_unique=n_unique,
            values=_values_preview(series),
            candidates=[Role.CONDITION.value, Role.SAMPLE_ID.value, Role.BATCH.value],
        )

    # 2. Known exact name (case-insensitive)
    if name_lower in _KNOWN_NAMES:
        role = _KNOWN_NAMES[name_lower]
        return SchemaEntry(
            role=role,
            confidence="high",
            n_unique=n_unique,
            values=_values_preview(series) if n_unique <= _LOW_CARDINALITY_THRESHOLD else None,
        )

    # 3. QC pattern substring
    if any(p in name_lower for p in _QC_PATTERNS):
        return SchemaEntry(role=Role.QC_METRIC, confidence="high", n_unique=n_unique)

    # 4. Clustering pattern substring (e.g. Seurat RNA_snn_res.0.5)
    if any(p in name_lower for p in _CLUSTERING_PATTERNS):
        return SchemaEntry(role=Role.CLUSTERING, confidence="high", n_unique=n_unique)

    # 5. Numeric dtype — high cardinality is continuous QC, low is ambiguous
    if pd.api.types.is_numeric_dtype(series):
        if n_unique > _LOW_CARDINALITY_THRESHOLD:
            return SchemaEntry(role=Role.QC_METRIC, confidence="high", n_unique=n_unique)
        return SchemaEntry(
            role=Role.AMBIGUOUS,
            confidence="low",
            n_unique=n_unique,
            values=_values_preview(series),
            candidates=[Role.CONDITION.value, Role.CLUSTERING.value],
        )

    # 6. Low-cardinality categorical/string, unknown name → ask user
    if n_unique <= _LOW_CARDINALITY_THRESHOLD:
        return SchemaEntry(
            role=Role.AMBIGUOUS,
            confidence="low",
            n_unique=n_unique,
            values=_values_preview(series),
            candidates=[Role.CONDITION.value, Role.SAMPLE_ID.value, Role.BATCH.value],
        )

    # 7. Fallthrough — high-cardinality string, unknown name
    return SchemaEntry(role=Role.OTHER, confidence="low", n_unique=n_unique)


def _values_preview(series: pd.Series) -> list[str]:
    """Return up to _MAX_VALUES_PREVIEW unique string values from the series."""
    unique = [str(v) for v in series.dropna().unique()]
    return unique[:_MAX_VALUES_PREVIEW]


def _entry_to_dict(entry: SchemaEntry) -> dict[str, Any]:
    """Convert SchemaEntry to a plain dict with a plain-str role (uns-safe)."""
    d = asdict(entry)
    d["role"] = entry.role.value
    return d


def _existing_meta(adata: "AnnData") -> dict[str, Any]:
    """Return a mutable copy of adata.uns['nvwa_meta'] if it's a dict, else empty."""
    uns: Any = getattr(adata, "uns", {}) or {}
    meta = uns.get("nvwa_meta", {}) if hasattr(uns, "get") else {}
    return dict(meta) if isinstance(meta, dict) else {}
