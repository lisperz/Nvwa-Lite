# The `src/analysis/` Folder Explained

The `analysis/` folder contains the **bioinformatics computation layer** — pure data processing logic with no UI or plotting concerns. It sits between the agent tools and the raw scanpy library, providing clean, validated wrappers that the agent can call safely.

---

## Background: Single-Cell RNA-seq Data Formats

### What is `.h5ad`?

`.h5ad` is the standard file format for single-cell RNA-seq data in Python. It stands for **HDF5 Annotated Data**.

**Structure:**

```
.h5ad file
├── X                    — Main data matrix (cells × genes), usually sparse
├── obs                  — Cell metadata (DataFrame: cell_id, cell_type, cluster, etc.)
├── var                  — Gene metadata (DataFrame: gene_name, highly_variable, etc.)
├── obsm                 — Multi-dimensional cell annotations (UMAP coords, PCA, etc.)
├── varm                 — Multi-dimensional gene annotations (PCA loadings, etc.)
├── uns                  — Unstructured metadata (DE results, parameters, etc.)
├── layers               — Alternative data matrices (raw counts, normalized, etc.)
└── raw                  — Backup of full gene set before subsetting
```

**Why HDF5?**

- **Efficient storage**: Single-cell datasets can have millions of cells × 20,000 genes. HDF5 compresses sparse matrices efficiently.
- **Partial loading**: You can read just the metadata without loading the full matrix into memory.
- **Cross-platform**: Works on Linux, macOS, Windows. Compatible with R (via `anndata2ri`).

### Other common single-cell file formats

| Format | Used by | Description |
|--------|---------|-------------|
| `.h5ad` | Python (scanpy, scvi-tools) | Standard Python format |
| `.rds` / `.rda` | R (Seurat) | R's native serialization format |
| `.loom` | Python/R (loompy, Seurat) | HDF5-based, designed for streaming large datasets |
| `.mtx` + `barcodes.tsv` + `genes.tsv` | 10x Genomics CellRanger | Sparse matrix in Matrix Market format + metadata files |
| `.h5` (10x format) | 10x Genomics | HDF5 with 10x-specific structure (not the same as `.h5ad`) |

**Conversion tools:**

- `scanpy.read_10x_mtx()` — reads 10x `.mtx` format into `.h5ad`
- `scanpy.read_loom()` — reads `.loom` into `.h5ad`
- `anndata2ri` — converts between Python `.h5ad` and R Seurat objects

---

## Real-World Use Cases

### Use Case 1: Cancer Research — Tumor Microenvironment Analysis

**Scenario:** A researcher has a tumor biopsy and wants to understand the immune cell composition.

**Workflow:**

1. **Upload `.h5ad`** — tumor sample with 50,000 cells
2. **Preprocess** — QC → normalize → cluster
3. **Identify cell types** — Run DE analysis, check marker genes:
   - Cluster 0: high `CD3E`, `CD3D` → T cells
   - Cluster 1: high `MS4A1`, `CD79A` → B cells
   - Cluster 2: high `LYZ`, `CD14` → Monocytes
   - Cluster 3: high `EPCAM`, `KRT18` → Tumor cells
4. **Visualize** — UMAP colored by cell type, violin plots for key markers
5. **Export** — Download plots for publication

**Why this matters:** Understanding which immune cells infiltrate the tumor helps predict treatment response (e.g., immunotherapy works better with high T cell infiltration).

---

### Use Case 2: Drug Discovery — Target Validation

**Scenario:** A pharma company wants to validate a drug target gene (`GENE_X`) in diseased vs. healthy tissue.

**Workflow:**

1. **Upload two datasets** — healthy lung (dataset A) and fibrotic lung (dataset B)
2. **Preprocess both** — normalize, cluster
3. **Compare expression** — Violin plot of `GENE_X` across conditions
4. **DE analysis** — Find genes co-expressed with `GENE_X` in disease
5. **Volcano plot** — Identify significantly upregulated genes in fibrotic cells

**Why this matters:** If `GENE_X` is highly expressed only in diseased cells, it's a good drug target (minimizes side effects on healthy tissue).

---

### Use Case 3: Developmental Biology — Cell Lineage Tracing

**Scenario:** A biologist studies how stem cells differentiate into neurons.

**Workflow:**

1. **Upload time-series data** — cells collected at day 0, 3, 7, 14
2. **Preprocess** — UMAP shows a trajectory from stem cells → progenitors → neurons
3. **DE analysis** — Find genes that turn on/off during differentiation
4. **Heatmap** — Show gene expression changes across pseudotime
5. **Feature plot** — Overlay key transcription factors on UMAP

**Why this matters:** Identifies the molecular switches that control cell fate decisions, useful for regenerative medicine.

---

### Use Case 4: Clinical Diagnostics — Disease Subtyping

**Scenario:** A hospital wants to classify leukemia patients into subtypes based on single-cell profiles.

**Workflow:**

1. **Upload patient sample** — blood cells from a leukemia patient
2. **Preprocess** — cluster cells
3. **Compare to reference** — Check if clusters match known leukemia subtypes (AML, CML, ALL)
4. **Marker validation** — Confirm subtype by checking expression of diagnostic markers
5. **Report** — Generate plots for clinical record

**Why this matters:** Different leukemia subtypes require different treatments. Single-cell analysis provides more precise diagnosis than bulk RNA-seq.

---

### Use Case 5: Quality Control — Dataset Validation Before Publication

**Scenario:** A researcher wants to verify their scRNA-seq data is high quality before submitting to a journal.

**Workflow:**

1. **Upload raw data** — check cell counts, gene counts
2. **QC metrics** — violin plots of `n_genes_by_counts`, `pct_counts_mt`
3. **Filter** — remove low-quality cells (< 200 genes, > 20% mitochondrial)
4. **Clustering** — verify expected cell types are present
5. **DE analysis** — confirm marker genes match literature

**Why this matters:** Journals require rigorous QC. Nvwa-Lite automates the standard checks and generates publication-ready figures.

---

## File 1: `preprocessing.py`

### What it does

Takes a **raw count matrix** (an `AnnData` object where each value is an integer read count) and runs the standard single-cell RNA-seq processing pipeline, returning a fully analysis-ready dataset.

### Why it exists

Raw scRNA-seq data cannot be visualized or clustered directly. It needs several transformation steps first. This file encapsulates the entire pipeline in one function so the agent can trigger it with a single tool call.

### The pipeline, step by step

```
Raw counts (cells × genes)
        ↓
1. QC metrics          — calculate % mitochondrial genes, total counts per cell
        ↓
2. Filter cells        — remove cells with too few genes (likely empty droplets)
   Filter genes        — remove genes expressed in too few cells (likely noise)
        ↓
3. Mito filter         — remove cells with high mitochondrial % (likely dying cells)
        ↓
4. Normalize           — scale each cell to 10,000 total counts (removes sequencing depth bias)
   Log1p transform     — apply log(x+1) to compress the dynamic range
        ↓
5. Save .raw           — store the log-normalized full gene set for DE analysis later
        ↓
6. Highly Variable Genes (HVG) — keep only the top 2,000 most informative genes
        ↓
7. Scale               — zero-mean, unit-variance per gene (required for PCA)
        ↓
8. PCA                 — reduce to ~50 principal components (dimensionality reduction)
        ↓
9. Neighbors graph     — build a k-nearest-neighbor graph in PCA space
   UMAP                — project to 2D for visualization
        ↓
10. Leiden clustering  — detect communities in the neighbor graph (cell type clusters)
        ↓
Processed AnnData (cells × HVGs, with UMAP coords + cluster labels)
```

### Key design detail: why `.raw` is saved after log1p

DE analysis (`rank_genes_groups`) needs the **log-normalized full gene set**, not the HVG-subset scaled data. So `.raw` is saved at step 5 — after normalization but before HVG subsetting. This is the standard scanpy convention.

### Return value

Returns a `(new_adata, PreprocessingResult)` tuple. The `PreprocessingResult` dataclass carries a human-readable summary message (cell counts before/after, steps performed) that the agent relays back to the user.

---

## File 2: `differential.py`

### What it does

Finds **marker genes** — genes that are significantly more (or less) expressed in one cluster compared to all others. This is how you identify what cell type each cluster represents.

### Why it exists

After clustering, you have groups of cells but no labels. Differential expression (DE) analysis answers: *"What genes define cluster 3?"* If cluster 3 highly expresses `MS4A1` and `CD79A`, it's likely B cells.

### Two functions

**`run_differential_expression(adata, groupby, method, n_genes)`**

Calls `sc.tl.rank_genes_groups()` — scanpy's built-in DE engine. Supports three statistical methods:

| Method | When to use |
|--------|-------------|
| `wilcoxon` | Default. Non-parametric, robust, works on most datasets |
| `t-test` | Faster, assumes roughly normal distribution |
| `logreg` | Logistic regression; good for very large datasets |

Results are stored back into `adata.uns["rank_genes_groups"]` in-place (scanpy convention).

**`get_de_dataframe(adata, group)`**

Extracts the stored DE results for a **specific cluster** into a tidy pandas DataFrame:

| Column | Meaning |
|--------|---------|
| `gene` | Gene name |
| `log2fc` | Log2 fold change vs. rest of cells |
| `pval` | Raw p-value |
| `pval_adj` | Adjusted p-value (corrected for multiple testing) |

This DataFrame is what `volcano.py` consumes to draw the volcano plot.

---

## How the two files connect to the rest of the system

```
User: "Preprocess the data"
        ↓
tools.py: preprocess_data tool
        ↓
analysis/preprocessing.py: run_preprocessing()
        ↓
Returns new AnnData → stored in session state

User: "Run differential expression"
        ↓
tools.py: differential_expression tool
        ↓
analysis/differential.py: run_differential_expression()
        ↓
Results stored in adata.uns["rank_genes_groups"]

User: "Show volcano plot for cluster 0"
        ↓
tools.py: volcano_plot_tool
        ↓
analysis/differential.py: get_de_dataframe()  ← extracts cluster 0 results
        ↓
plotting/volcano.py: plot_volcano()            ← draws the plot
```

---

## Summary

| File | Input | Output | Purpose |
|------|-------|--------|---------|
| `preprocessing.py` | Raw count AnnData | Processed AnnData + UMAP + clusters | Make data ready for visualization |
| `differential.py` | Processed AnnData + cluster key | Marker gene table per cluster | Identify what each cluster is |
