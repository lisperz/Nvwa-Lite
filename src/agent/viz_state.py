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
    """Record the parameters of the last successful plot."""
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
