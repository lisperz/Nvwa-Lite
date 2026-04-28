"""Domain knowledge for the spec pipeline.

Importing this package triggers @register side-effects from each section
module so src.core.registry.REGISTRY is populated for the agent at import time.

Per local/product/tool_migration_map.md: one module per Yalu Layer 1 section
(or per Yalu subsection if section has >5 tools after collapse). Inspection
is the agent-utility bucket (no Yalu scenario).
"""

from src.domain import inspection  # noqa: F401 — triggers @register side-effect
from src.domain import qc  # noqa: F401 — triggers @register side-effect
