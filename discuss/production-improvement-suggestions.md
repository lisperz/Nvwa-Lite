# Nvwa-Lite: Production Improvement Suggestions

Based on 5 key production readiness questions, this document maps each concern to concrete improvements for the current MVP.

---

## Q1: Resource Control — Preventing OOM and Cost Overruns

**Current gap:** The MVP runs a single Streamlit process with no memory/CPU limits. A large h5ad upload or a runaway Scanpy computation can OOM the host.

### Recommended improvements

**1. Docker resource limits (immediate)**

Add hard limits to `docker-compose.yml`:

```yaml
services:
  nvwa-lite:
    build: .
    ports:
      - "8501:8501"
    env_file:
      - .env
    volumes:
      - ./logs:/app/logs
    restart: unless-stopped
    deploy:
      resources:
        limits:
          memory: 8g
          cpus: "4.0"
        reservations:
          memory: 2g
```

**2. Pre-execution metadata check in `executor.py`**

Before running any Scanpy plot, check dataset size and reject if it exceeds a configurable threshold:

```python
MAX_CELLS = 500_000  # configurable via env var

def _check_dataset_size(adata: AnnData) -> None:
    if adata.n_obs > MAX_CELLS:
        raise ValueError(
            f"Dataset too large ({adata.n_obs:,} cells). "
            f"Max supported: {MAX_CELLS:,}. Please subsample first."
        )
```

**3. Per-request timeout watchdog**

Wrap agent invocation in `app.py` with a timeout using `concurrent.futures`:

```python
import concurrent.futures

AGENT_TIMEOUT_SECONDS = 120

with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
    future = executor.submit(agent.invoke, prompt, chat_history)
    try:
        response = future.result(timeout=AGENT_TIMEOUT_SECONDS)
    except concurrent.futures.TimeoutError:
        st.error("Request timed out. Please try a simpler query.")
        st.stop()
```

---

## Q2: Secure OpenAI API Usage in Private Deployment

**Current gap:** The API key is passed directly from the UI/env to OpenAI. There is no audit trail, no DLP filtering, and no way to prevent raw biological data from leaking into prompts.

### Recommended improvements

**1. Introduce a FastAPI proxy layer**

Add a thin `src/gateway/` service between the Streamlit UI and OpenAI:

```
Streamlit UI → FastAPI Gateway → OpenAI API
```

The gateway handles:
- API key storage (never exposed to the browser)
- DLP: strip any gene sequences or patient identifiers from prompts before forwarding
- Audit logging: every request/response pair written to PostgreSQL with timestamp, user ID, token count

**2. Enforce "local compute, cloud inference" contract**

The current system already does this correctly — raw h5ad data never leaves the container. Formalize it with a prompt guardrail in `prompts.py`:

```python
SYSTEM_PROMPT_TEMPLATE = """
...
## Data Privacy Rules
- NEVER include raw gene expression values, cell barcodes, or patient identifiers in your responses.
- Only reference gene names, cluster labels, and statistical summaries.
...
"""
```

**3. Azure OpenAI for enterprise clients**

For pharma/biotech clients with compliance requirements (HIPAA, GxP), swap the OpenAI endpoint for Azure OpenAI in `core.py`:

```python
from langchain_openai import AzureChatOpenAI

llm = AzureChatOpenAI(
    azure_deployment=os.environ["AZURE_OPENAI_DEPLOYMENT"],
    azure_endpoint=os.environ["AZURE_OPENAI_ENDPOINT"],
    api_key=os.environ["AZURE_OPENAI_API_KEY"],
    api_version="2024-02-01",
    temperature=0,
)
```

Azure OpenAI provides VPC isolation and guarantees data is not used for model training — critical for regulated industries.

---

## Q3: Redesigning as a Production LLM Agent API Service

**Current gap:** The MVP is a monolithic Streamlit app. For production, the UI and agent logic need to be decoupled so the agent can serve multiple frontends, support async long-running jobs, and scale independently.

### Recommended architecture

```
Browser / CLI / Jupyter
        ↓
  FastAPI REST API  (src/api/)
        ↓
  Task Queue: Redis + Celery  (async jobs)
        ↓
  AgentRunner (LangGraph)  (src/agent/)
        ↓
  Tool Executors  (src/plotting/, src/tools/)
        ↓
  PostgreSQL  (conversation history, job state)
```

**1. Replace Streamlit with FastAPI backend**

Add `src/api/routes.py`:

```python
@router.post("/analyze")
async def analyze(request: AnalyzeRequest) -> AnalyzeResponse:
    task = celery_app.send_task("agent.run", args=[request.query, request.session_id])
    return AnalyzeResponse(task_id=task.id)

@router.get("/analyze/{task_id}")
async def get_result(task_id: str) -> TaskResult:
    result = AsyncResult(task_id, app=celery_app)
    return TaskResult(status=result.status, result=result.result)
```

**2. Upgrade to LangGraph for complex workflows**

Replace the current `AgentRunner` loop with LangGraph for multi-step bioinformatics workflows (e.g., QC → normalization → clustering → visualization):

```python
from langgraph.graph import StateGraph

workflow = StateGraph(AgentState)
workflow.add_node("plan", planning_node)
workflow.add_node("execute_tool", tool_node)
workflow.add_node("summarize", summary_node)
workflow.add_conditional_edges("execute_tool", should_continue)
```

**3. PostgreSQL for persistent conversation history**

Replace the in-memory `st.session_state.chat_history` list with a database-backed store:

```python
# src/storage/conversation.py
class ConversationStore:
    async def append(self, session_id: str, role: str, content: str) -> None: ...
    async def get_history(self, session_id: str) -> list[ChatMessage]: ...
```

**4. Structured output with Pydantic**

All tool inputs/outputs should be validated with Pydantic models to prevent pipeline failures from malformed LLM outputs. The current `PlotResult` dataclass is a good start — extend this pattern to all API request/response schemas.

---

## Q4: Billing and Resource Quota Management

**Current gap:** There is no token tracking, no per-user quota, and no cost visibility. In production, a single user could exhaust the API budget.

### Recommended improvements

**1. LangChain callback for token tracking**

Inject a `UsageTrackingCallback` into the agent in `core.py`:

```python
from langchain_core.callbacks import BaseCallbackHandler

class UsageTrackingCallback(BaseCallbackHandler):
    def __init__(self, session_id: str, usage_store: UsageStore) -> None:
        self.session_id = session_id
        self.usage_store = usage_store

    def on_llm_end(self, response, **kwargs) -> None:
        usage = response.llm_output.get("token_usage", {})
        self.usage_store.record(
            session_id=self.session_id,
            prompt_tokens=usage.get("prompt_tokens", 0),
            completion_tokens=usage.get("completion_tokens", 0),
        )
```

**2. Async write pipeline (Redis → PostgreSQL)**

To avoid adding latency to the agent response, write usage data asynchronously:

```
Agent callback → Redis queue → Background worker → PostgreSQL user_usage table
```

Schema:
```sql
CREATE TABLE user_usage (
    id          SERIAL PRIMARY KEY,
    session_id  TEXT NOT NULL,
    user_id     TEXT,
    model       TEXT,
    prompt_tokens     INT,
    completion_tokens INT,
    created_at  TIMESTAMPTZ DEFAULT NOW()
);
```

**3. FastAPI middleware for request-level tracking**

Capture every API call at the gateway layer regardless of agent internals:

```python
@app.middleware("http")
async def usage_middleware(request: Request, call_next):
    start = time.time()
    response = await call_next(request)
    duration_ms = int((time.time() - start) * 1000)
    await usage_store.record_request(
        user_id=request.headers.get("X-User-ID"),
        endpoint=request.url.path,
        duration_ms=duration_ms,
        status_code=response.status_code,
    )
    return response
```

**4. Per-user quota enforcement**

Before invoking the agent, check the user's remaining quota:

```python
async def check_quota(user_id: str) -> None:
    used = await usage_store.get_monthly_tokens(user_id)
    limit = await quota_store.get_limit(user_id)
    if used >= limit:
        raise QuotaExceededError(f"Monthly token limit ({limit:,}) reached.")
```

---

## Q5: Container Strategy — Docker vs LXC for Production

**Current situation:** The MVP uses a single Docker container, which is correct for this stage.

### Path to production

**Keep Docker for the application layer.** The reasons that made Docker right for the MVP remain valid in production:
- Fast startup/teardown for per-request analysis containers
- Python Docker SDK integration for programmatic container management
- Cloud-native: works directly with AWS ECS, EKS, GCP Cloud Run

**Add Kubernetes (or AWS ECS) for orchestration:**

```yaml
# k8s/deployment.yaml
apiVersion: apps/v1
kind: Deployment
spec:
  replicas: 3
  template:
    spec:
      containers:
        - name: nvwa-lite
          image: nvwa-lite:latest
          resources:
            limits:
              memory: "8Gi"
              cpu: "4"
          env:
            - name: OPENAI_API_KEY
              valueFrom:
                secretKeyRef:
                  name: nvwa-secrets
                  key: openai-api-key
```

**Consider LXC/Singularity only for HPC clients:**

If a pharma client runs analyses on an on-premise HPC cluster that prohibits Docker (common in regulated environments), use Singularity containers instead. Singularity is Docker-compatible (can convert Docker images) but runs without root privileges, which is required on most HPC systems.

**Recommended production container strategy:**

| Scenario | Technology |
|----------|-----------|
| SaaS web app | Docker + Kubernetes/ECS |
| Per-analysis isolation | Docker-in-Docker or sidecar containers |
| On-prem pharma HPC | Singularity (converted from Docker image) |
| Long-running virtual lab | LXC (persistent multi-tool environment) |

---

## Summary: Priority Order for Production Readiness

| Priority | Change | Impact |
|----------|--------|--------|
| P0 | Docker memory/CPU limits + timeout watchdog | Prevents OOM crashes |
| P0 | API key moved to server-side proxy | Security baseline |
| P1 | FastAPI backend + async task queue | Scalability, multi-user support |
| P1 | LangChain usage callbacks + PostgreSQL | Cost visibility |
| P2 | LangGraph for multi-step workflows | Handles complex bioinformatics pipelines |
| P2 | Kubernetes deployment | Horizontal scaling |
| P3 | Azure OpenAI integration | Enterprise/pharma compliance |
| P3 | Singularity support | HPC client compatibility |

