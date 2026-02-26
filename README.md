# Nvwa-Lite: Single-Cell RNA-seq Visualization Agent

A web app that lets biologists explore scRNA-seq data through natural language. Powered by OpenAI, LangChain, Scanpy, and Streamlit.

## Try it Online

👉 **[nvwa-lite.streamlit.app](https://nvwa-lite.streamlit.app)** — upload a `.h5ad` file and start chatting. No setup required.

## Run Locally (larger files, more memory)

```bash
# 1. Clone the repo
git clone https://github.com/lisperz/Nvwa-Lite.git && cd Nvwa-Lite

# 2. Add your OpenAI API key
cp .env.example .env
# Edit .env: OPENAI_API_KEY=sk-...

# 3. Run with Docker
docker-compose up

# 4. Open http://localhost:8501
```

Or without Docker:
```bash
uv sync
./scripts/start.sh
```

## Data Format

Upload `.h5ad` files (max 500MB via web, unlimited locally). To convert from Seurat `.rds`:

```r
library(SeuratDisk)
SaveH5Seurat(seurat_obj, filename = "data.h5Seurat")
Convert("data.h5Seurat", dest = "h5ad")
```

## What You Can Do

Ask in plain English — the agent picks the right tool automatically:

| Ask | What happens |
|-----|-------------|
| "What's in this dataset?" | Shows cell/gene counts, obs keys |
| "Preprocess the data" | QC → normalize → PCA → UMAP → Leiden clustering |
| "Show UMAP colored by cell type" | UMAP plot |
| "Violin plot for MS4A1" | Expression distribution across clusters |
| "Run differential expression" | Marker genes per cluster |
| "Volcano plot for cluster 0" | DE volcano plot |
| "Heatmap for CD3E, MS4A1, NKG7" | Gene expression heatmap |

Every plot has a PNG download button.

## Tech Stack

Python 3.11 · Streamlit · Scanpy · LangChain · OpenAI gpt-4o-mini · Docker
