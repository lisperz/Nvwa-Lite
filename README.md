# Nvwa-Lite: Single-Cell RNA-seq Visualization Agent

A Dockerized web app that lets biologists explore scRNA-seq data through natural language. Powered by OpenAI, LangChain, Scanpy, and Streamlit.

## Quick Start

```bash
# 1. Clone the repo
git clone <repo-url> && cd nvwa-lite

# 2. Set your OpenAI API key
cp .env.example .env
# Edit .env and add your key: OPENAI_API_KEY=sk-...

# 3. Run with Docker
docker-compose up

# 4. Open browser
# http://localhost:8501
```

## Local Development

```bash
# Install dependencies
uv sync

# Start the app
./scripts/start.sh

# Stop the app
./scripts/stop.sh
```

## Supported Queries

- "Show me the UMAP plot colored by cell type."
- "Visualize the expression of CD3E on the UMAP."
- "Show a violin plot for MS4A1 across all clusters."
- "Create a dot plot for CD3E, MS4A1, NKG7."
- "What genes are available in this dataset?"

## Architecture

```
src/
├── ui/          # Streamlit interface (chat loop, sidebar widgets)
├── agent/       # LangChain agent (tools, prompts, orchestration)
└── plotting/    # Scanpy plot execution and validation
```

## Tech Stack

- Python 3.11, Streamlit, Scanpy, LangChain, OpenAI gpt-4o, Docker
