FROM python:3.11-slim

# System deps for matplotlib headless rendering
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgl1 libglib2.0-0 libfontconfig1 curl \
    && rm -rf /var/lib/apt/lists/*

# Install uv
RUN pip install uv

WORKDIR /app

# Install dependencies first (layer caching)
COPY pyproject.toml uv.lock* ./
RUN uv sync --no-dev

# Copy source code
COPY . .

# Pre-cache PBMC3k dataset into the image
RUN uv run python -c "import scanpy as sc; sc.datasets.pbmc3k_processed()"

# Create required directories
RUN mkdir -p logs data/uploads

EXPOSE 8501

HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health || exit 1

CMD ["uv", "run", "streamlit", "run", "src/ui/app.py", \
     "--server.port=8501", "--server.address=0.0.0.0", \
     "--server.headless=true"]
