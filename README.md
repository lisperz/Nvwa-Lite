# Nvwa-Lite: Single-Cell RNA-seq Visualization Agent

A Dockerized web app that lets biologists explore scRNA-seq data through natural language. Powered by OpenAI, LangChain, Scanpy, and Streamlit.

## Quick Start

```bash
# 1. Clone the repo
git clone <repo-url> && cd nvwa-lite

# 2. Set your OpenAI API key
cp .env.example .env
# Edit .env: OPENAI_API_KEY=sk-...

# 3. Run with Docker
docker-compose up

# 4. Open browser at http://localhost:8501
```

## Local Development

```bash
uv sync
./scripts/start.sh   # start
./scripts/stop.sh    # stop
```

## Features

### Upload & Analyze Your Own Data
- Upload any `.h5ad` file via the sidebar
- Built-in PBMC3k demo dataset always available
- Dataset state tracking (knows what processing has been done)

### Preprocessing Pipeline
Ask: *"Preprocess the data"* — runs the full pipeline automatically:
- QC metrics → filter cells/genes → mito filter
- Normalize → log1p → highly variable genes
- PCA → UMAP → Leiden clustering

### 10 Analysis Tools

| Tool | What it does |
|------|-------------|
| `dataset_info` | Dataset metadata, cell/gene counts, obs keys |
| `check_data_status` | What preprocessing steps have been applied |
| `preprocess_data` | Full QC → UMAP → clustering pipeline |
| `differential_expression` | Marker genes per cluster (Wilcoxon/t-test) |
| `umap_plot` | UMAP colored by cluster or gene |
| `violin_plot` | Gene expression distribution across groups |
| `dotplot` | Multi-gene dot plot |
| `feature_plot` | Gene expression overlaid on UMAP (viridis) |
| `heatmap_plot` | Gene expression heatmap |
| `volcano_plot_tool` | DE volcano plot with non-overlapping labels |

### Download Plots
Every generated plot has a PNG download button.

## Example Queries

```
"What's in this dataset?"
"Preprocess the data."
"Show me the UMAP colored by cell type."
"Run differential expression."
"Show a volcano plot for cluster 0."
"Show a heatmap for CD3E, MS4A1, NKG7, LYZ."
"Show a violin plot for MS4A1 across all clusters."
```

## Architecture

```
src/
├── types.py              # DatasetState dataclass
├── agent/
│   ├── core.py           # AgentRunner (gpt-4o-mini)
│   ├── tools.py          # 10 LangChain tools
│   └── prompts.py        # System prompt with dataset context
├── analysis/
│   ├── preprocessing.py  # QC → UMAP → clustering pipeline
│   └── differential.py   # rank_genes_groups wrapper
├── plotting/
│   ├── executor.py       # UMAP, violin, dotplot, heatmap, feature
│   ├── volcano.py        # Volcano plot (matplotlib + adjustText)
│   └── styles.py         # Publication-ready plot style
└── ui/
    ├── app.py            # Streamlit chat loop + upload flow
    └── components.py     # Sidebar widgets
```

## Tech Stack

- Python 3.11, Streamlit, Scanpy, LangChain, OpenAI gpt-4o-mini, Docker
- adjustText for non-overlapping volcano labels
- leidenalg + igraph for Leiden clustering
