"""Self-registering tool catalog.

Each family module imported below triggers `@register(...)` side-effects that
populate `src.tools.registry.REGISTRY` at import time. T-040 will add
`plotting`, `analysis`, `reasoning`, `subset` to the import list.
"""

from src.tools import inspection  # noqa: F401 — import triggers @register side-effect
