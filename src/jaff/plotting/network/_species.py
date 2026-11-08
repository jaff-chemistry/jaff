# ABOUTME: Draws the species view of net.visualize: the reactions forming one
# ABOUTME: species on its left, those destroying it on its right, with partners
"""Matplotlib renderer of a :class:`SpeciesLayout`."""

from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import FancyBboxPatch

from ._common import FONT_RC
from ._graph import (
    ROLES,
    _edges,
    _reactions,
    draw_legends,
    entry_widths,
    pill,
    role_handles,
    type_handles,
)
from ._layout import (
    drawing_bottom,
    HEADER,
    HEADER_FONT,
    MARGIN,
    SpeciesLayout,
)
from ._markers import label_shapes
from ._rates import RateColours
from ._themes import NetworkTheme, heading_colour, legend_entries, type_colours

if TYPE_CHECKING:
    from ...core.network.network import Network
    from ...core.species.specie import Specie


def involved(net: Network, specie: Specie) -> list[int]:
    """Indices of the reactions that form or destroy *specie* (a net change)."""
    j = specie.index
    change = net.product_matrix[:, j] - net.reactant_matrix[:, j]

    return [i for i in range(len(net.reactions)) if change[i] != 0]


def species_legend_widths(
    net: Network, specie: Specie, theme: NetworkTheme
) -> list[float]:
    """Width in inches of each entry of the legend of the species view of *specie*."""
    colours = type_colours(theme, [r.type for r in net.reactions])
    shown = [net.reactions[i].type for i in involved(net, specie)]
    labels = [t for t, _ in legend_entries(theme, colours, shown)]

    return entry_widths([*labels, *ROLES])


def draw_species(
    net: Network,
    layout: SpeciesLayout,
    theme: NetworkTheme,
    rated: RateColours | None = None,
) -> tuple[Figure, Axes]:
    """Draw *layout* of *net* in *theme*; the axes are in inches, y down.

    *rated* colours the reactions by rate, as in
    :func:`~jaff.plotting.network._graph.draw_graph`.
    """
    types = [r.type for r in net.reactions]
    colours = type_colours(theme, types)
    width, height = layout.size
    with plt.rc_context(FONT_RC):
        fig = plt.figure(figsize=layout.size, facecolor=theme.surface)
        ax = fig.add_axes((0, 0, 1, 1))
        ax.set_xlim(0, width)
        ax.set_ylim(height, 0)
        ax.set_axis_off()
        _lanes(ax, layout, theme)
        _edges(ax, layout, theme, colours, types, rated)
        _reactions(ax, layout, theme, colours, types, rated)
        _pills(ax, net, layout, theme)
        shown = [types[i] for i in sorted(layout.formation + layout.destruction)]
        _legend(ax, layout, theme, colours, shown, rated)

    return fig, ax


def _lanes(ax: Axes, layout: SpeciesLayout, theme: NetworkTheme) -> None:
    """Bands behind the species lanes; headers above every lane that has one."""
    bottom = drawing_bottom(layout)
    for k, lane in enumerate(layout.lanes):
        if k % 2 == 0:
            ax.add_patch(
                FancyBboxPatch(
                    (lane.x0, MARGIN + HEADER / 2),
                    lane.x1 - lane.x0,
                    bottom - MARGIN - HEADER / 2,
                    boxstyle="round,pad=0,rounding_size=0.08",
                    facecolor=theme.band,
                    edgecolor="none",
                    zorder=0,
                )
            )
        if lane.header:
            ax.text(
                lane.centre,
                MARGIN + HEADER / 4,
                lane.header,
                ha="center",
                va="center",
                fontsize=HEADER_FONT,
                color=heading_colour(theme),
            )


def _pills(ax: Axes, net: Network, layout: SpeciesLayout, theme: NetworkTheme) -> None:
    width = layout.pill_width
    for role, positions in (
        ("reactant", layout.reactant_xy),
        ("product", layout.product_xy),
    ):
        for s, (x, y) in positions.items():
            name = net.species[s].name
            pill(ax, x, y, net.species[s], width[s], theme, f"{role}:{name}")
    target = net.species[layout.target]
    x, y = layout.target_xy
    pill(
        ax, x, y, target, width[layout.target], theme, f"target:{target.name}", bold=True
    )


def _legend(
    ax: Axes,
    layout: SpeciesLayout,
    theme: NetworkTheme,
    colours: dict[str, str],
    shown: list[str],
    rated: RateColours | None = None,
) -> None:
    entries = legend_entries(theme, colours, shown)
    handles = type_handles(theme, entries, label_shapes(theme, list(colours)), rated)
    draw_legends(ax, layout, theme, [*handles, *role_handles(theme)], rated)
