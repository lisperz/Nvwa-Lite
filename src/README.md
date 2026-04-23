# src/

## Folder layout

```
src/
├── agent/          # LLM orchestration: core loop, router, tools, prompts, output_guard, viz_state
├── core/           # shared types + adata schema primitives — zero deps, used everywhere
├── domain/
│   ├── analysis/   # adata → data (DE, QC, composition, clustering, gene lookup, preprocessing)
│   └── plotting/   # adata (+data) → PlotResult images (executor, volcano, comparison, styles)
├── platform/
│   ├── infra/      # external services: auth, storage, db
│   └── observability/  # EventLogger, analytics, dashboard, uns_snapshot
├── session/        # session lifecycle (manager, cleanup_server)
├── tools/          # empty — reserved for T-040 (agent/tools.py split)
└── ui/             # Streamlit surface (app, components, feedback dialogs)
```

**Layer rule:** `domain/analysis` returns data; `domain/plotting` returns PlotResults. Plotting may call analysis; analysis must not import plotting.

**Coming later** (land with their owning tasks, not in T-044): `src/spec_validation/` (T-038), `src/gatekeeper/` (T-004), `src/workflows/` (T-055).

## Where did X go? (T-044 quick reference)

| Before | After |
|---|---|
| `src/types.py` | `src/core/types.py` |
| `src/plotting/validation.py` | `src/core/adata_schema.py` |
| `src/analysis/*` | `src/domain/analysis/*` |
| `src/plotting/*` (non-validation) | `src/domain/plotting/*` |
| `src/auth/service.py` | `src/platform/infra/auth.py` |
| `src/storage/service.py` | `src/platform/infra/storage.py` |
| `src/db/*` | `src/platform/infra/db/*` |
| `src/logging/service.py` | `src/platform/observability/events.py` |
| `src/logging/uns_snapshot.py` | `src/platform/observability/uns_snapshot.py` |
| `src/monitoring/*` | `src/platform/observability/*` |
