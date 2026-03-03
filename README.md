# Nvwa-Lite

A conversational agent for single-cell RNA-seq data visualization. Upload `.h5ad` files and explore your data through natural language.

## Quick Start

**Docker (recommended):**
```bash
git clone https://github.com/lisperz/Nvwa-Lite.git
cd Nvwa-Lite
cp .env.example .env  # Add your OPENAI_API_KEY
docker-compose up
```

Open http://localhost:8501

**Local:**
```bash
uv sync
./scripts/start.sh
```

## Usage

Upload a `.h5ad` file and ask questions:
- "What's in this dataset?"
- "Preprocess the data"
- "Show UMAP colored by clusters"
- "Find marker genes"
- "Violin plot for CD3E"

## Tech Stack

Python 3.11 · Streamlit · Scanpy · LangChain · OpenAI gpt-4o-mini
