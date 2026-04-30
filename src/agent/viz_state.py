"""
Visualization state tracking for multi-turn plot refinement.

Tracks the last visualization's parameters across chat turns so the agent
can compose constraints instead of dropping them.

Follows the same singleton pattern as _adata / _dataset_state in tools.py.
"""

from dataclasses import dataclass, field, fields
from typing import Optional, List


@dataclass
class VisualizationState:
    """Tracks parameters of the last visualization."""

    plot_type: Optional[str] = None
    color_by: Optional[str] = None
    split_by: Optional[str] = None
    show_labels: bool = False
    show_legend: bool = True
    groupby: Optional[str] = None
    genes: List[str] = field(default_factory=list)
    subset_key: Optional[str] = None
    subset_value: Optional[List[str]] = None
    group1: Optional[str] = None
    group2: Optional[str] = None
    reference: Optional[str] = None
    celltype: Optional[str] = None

    def to_prompt_block(self) -> str:
        """Render state as structured text for prompt injection."""
        if self.plot_type is None:
            return "No previous visualization in this session."

        lines = [f"Plot type: {self.plot_type}"]
        if self.color_by:
            lines.append(f"color_by: {self.color_by}")
        if self.split_by:
            lines.append(f"split_by: {self.split_by}")
        if self.groupby:
            lines.append(f"groupby: {self.groupby}")
        if self.genes:
            lines.append(f"genes: {', '.join(self.genes)}")
        if self.subset_key:
            lines.append(f"subset_key: {self.subset_key}")
        if self.subset_value:
            lines.append(f"subset_value: {', '.join(self.subset_value)}")
        if self.group1:
            lines.append(f"group1: {self.group1}")
        if self.group2:
            lines.append(f"group2: {self.group2}")
        if self.reference:
            lines.append(f"reference: {self.reference}")
        if self.celltype:
            lines.append(f"celltype: {self.celltype}")
        lines.append(f"show_labels: {self.show_labels}")
        lines.append(f"show_legend: {self.show_legend}")
        return "\n".join(lines)

    def clear(self) -> None:
        """Reset all state (called on dataset load)."""
        self.plot_type = None
        self.color_by = None
        self.split_by = None
        self.show_labels = False
        self.show_legend = True
        self.groupby = None
        self.genes = []
        self.subset_key = None
        self.subset_value = None
        self.group1 = None
        self.group2 = None
        self.reference = None
        self.celltype = None


# Module-level singleton (same pattern as tools.py _adata)
_viz_state: Optional[VisualizationState] = None


def bind_viz_state(state: VisualizationState) -> None:
    """Bind a VisualizationState instance for the current session."""
    global _viz_state
    _viz_state = state


def get_viz_state() -> Optional[VisualizationState]:
    """Return the current visualization state."""
    return _viz_state


def update_viz_state(plot_type: str, **params) -> None:
    """Record the parameters of the last successful plot.

    NOTE: signature accepts arbitrary **params for legacy compatibility —
    ``src/agent/tools.py`` callers pass ``row_key``/``col_key``/``cell_type_key``
    that aren't persisted. After M2 (legacy purge), tighten this to explicit
    keyword args so unknown kwargs raise TypeError instead of silent drop.
    """
    if _viz_state is None:
        return
    _viz_state.plot_type = plot_type
    # Reset optional fields before applying new params
    _viz_state.color_by = params.get("color_by")
    _viz_state.split_by = params.get("split_by")
    _viz_state.show_labels = params.get("show_labels", False)
    _viz_state.show_legend = params.get("show_legend", True)
    _viz_state.groupby = params.get("groupby")
    _viz_state.genes = params.get("genes", [])
    _viz_state.subset_key = params.get("subset_key")
    _viz_state.subset_value = params.get("subset_value")
    _viz_state.group1 = params.get("group1")
    _viz_state.group2 = params.get("group2")
    _viz_state.reference = params.get("reference")
    _viz_state.celltype = params.get("celltype")
