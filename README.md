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

## Regression Testing

The test runner exercises the agent in-process (no server required).

**Prerequisites:** `OPENAI_API_KEY` must be set and the project deps installed (`uv sync`).

```bash
# Full suite against a local h5ad file
python scripts/run_tests.py --dataset /absolute/path/to/pbmc3k.h5ad

# Quick sanity check (first 3 cases only)
python scripts/run_tests.py --dataset /absolute/path/to/pbmc3k.h5ad --limit 3

# Custom model and report path
python scripts/run_tests.py \
  --dataset /absolute/path/to/pbmc3k.h5ad \
  --model gpt-4o \
  --report reports/my_run.md
```

The report is written to `reports/regression_report.md` and the console prints a
`PASS/FAIL/WARN/ERROR/SKIP` count summary.  Test cases are defined in
`tests/tests.yaml` (20 cases across happy-path, input-error, resource-scale, and
out-of-scope categories).

## Tech Stack

Python 3.11 · Streamlit · Scanpy · LangChain · OpenAI gpt-4o-mini
