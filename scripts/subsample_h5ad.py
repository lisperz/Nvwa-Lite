#!/usr/bin/env python3
"""
Subsample a large .h5ad file to a smaller one for fast regression testing.

Usage
-----
  uv run python scripts/subsample_h5ad.py --input /path/to/full.h5ad --output /path/to/subset.h5ad
  uv run python scripts/subsample_h5ad.py --input /path/to/full.h5ad --output /path/to/subset.h5ad --n-cells 2000
"""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Subsample an h5ad file for testing.")
    parser.add_argument("--input", required=True, help="Path to the source .h5ad file.")
    parser.add_argument("--output", required=True, help="Path for the subsampled .h5ad file.")
    parser.add_argument("--n-cells", type=int, default=2000, help="Number of cells to keep (default: 2000).")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42).")
    parser.add_argument("--drop-raw", action="store_true", help="Drop the .raw layer to reduce file size further.")
    args = parser.parse_args()

    import scanpy as sc

    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    print(f"Reading {input_path} ...")
    adata = sc.read_h5ad(input_path)
    print(f"  Original: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    n_cells = min(args.n_cells, adata.n_obs)
    sc.pp.subsample(adata, n_obs=n_cells, random_state=args.seed)
    print(f"  Subsampled to: {adata.n_obs:,} cells")

    # Keep only the embeddings needed for tests; drop bulky extras
    keep_obsm = {"X_umap", "X_pca"}
    for key in list(adata.obsm.keys()):
        if key not in keep_obsm:
            del adata.obsm[key]

    if args.drop_raw:
        adata.raw = None
        print("  Dropped .raw layer")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path, compression="gzip")

    size_mb = output_path.stat().st_size / 1024 / 1024
    print(f"  Written to: {output_path}  ({size_mb:.1f} MB)")
    print(adata)


if __name__ == "__main__":
    main()
