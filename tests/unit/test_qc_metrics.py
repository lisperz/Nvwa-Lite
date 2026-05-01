"""Unit tests for QC metrics calculation and grouping."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from src.domain.qc_metrics import summarize_qc_metrics


@pytest.fixture
def adata_with_qc():
    """Create a synthetic AnnData with QC metrics and condition column."""
    n_cells = 100
    n_genes = 50

    X = np.random.poisson(5, size=(n_cells, n_genes))

    obs = pd.DataFrame({
        "total_counts": np.random.randint(1000, 5000, n_cells),
        "n_genes_by_counts": np.random.randint(500, 2000, n_cells),
        "pct_counts_mt": np.random.uniform(0, 10, n_cells),
        "condition": ["Control"] * 50 + ["Mutant"] * 50,
        "orig.ident": ["Sample_A"] * 30 + ["Sample_B"] * 30 + ["Sample_C"] * 40,
    })

    var = pd.DataFrame(index=[f"Gene_{i}" for i in range(n_genes)])

    return AnnData(X=X, obs=obs, var=var)


def test_summarize_qc_metrics_global(adata_with_qc):
    """Global QC summary returns one row per metric."""
    df = summarize_qc_metrics(adata_with_qc)

    assert len(df) == 3
    assert "metric" in df.columns
    assert set(df["metric"]) == {"total_counts", "n_genes_by_counts", "pct_counts_mt"}
    assert "n_cells" in df.columns
    assert "mean" in df.columns
    assert "median" in df.columns
    assert "group" not in df.columns


def test_summarize_qc_metrics_grouped_by_condition(adata_with_qc):
    """Grouped QC returns one row per (group × metric), with 'group' column first."""
    df = summarize_qc_metrics(adata_with_qc, groupby="condition")

    assert "group" in df.columns
    assert df.columns[0] == "group"
    assert df.columns[1] == "metric"

    groups = df["group"].unique().tolist()
    assert set(groups) == {"Control", "Mutant"}

    # 2 groups × 3 metrics = 6 rows
    assert len(df) == 6

    control_rows = df[df["group"] == "Control"]
    assert len(control_rows) == 3
    assert control_rows["n_cells"].iloc[0] == 50

    mutant_rows = df[df["group"] == "Mutant"]
    assert len(mutant_rows) == 3
    assert mutant_rows["n_cells"].iloc[0] == 50


def test_summarize_qc_metrics_grouped_three_groups(adata_with_qc):
    """Grouped QC with three groups returns 3 × 3 = 9 rows."""
    df = summarize_qc_metrics(adata_with_qc, groupby="orig.ident")

    assert "group" in df.columns
    assert len(df["group"].unique()) == 3
    assert len(df) == 9


def test_summarize_qc_metrics_grouped_invalid_column(adata_with_qc):
    """Invalid groupby column raises ValueError with helpful message."""
    with pytest.raises(ValueError, match="not found in adata.obs"):
        summarize_qc_metrics(adata_with_qc, groupby="nonexistent_column")


def test_summarize_qc_metrics_grouped_preserves_stats(adata_with_qc):
    """Per-group stats are computed only on cells in that group."""
    df = summarize_qc_metrics(adata_with_qc, groupby="condition")

    control_total = df[(df["group"] == "Control") & (df["metric"] == "total_counts")]
    assert len(control_total) == 1

    # Verify n_cells matches actual group size
    assert int(control_total["n_cells"].iloc[0]) == 50

    # Verify mean is computed only on Control cells
    expected_mean = float(
        adata_with_qc.obs.loc[adata_with_qc.obs["condition"] == "Control", "total_counts"].mean()
    )
    assert abs(float(control_total["mean"].iloc[0]) - expected_mean) < 1e-6
