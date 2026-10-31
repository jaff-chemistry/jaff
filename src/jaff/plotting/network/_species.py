# ABOUTME: Draws the species view of net.visualize: the reactions forming one
# ABOUTME: species on its left, those destroying it on its right, with partners
"""Matplotlib renderer of a :class:`SpeciesLayout`."""

from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
from matplotlib.artist import Artist
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import FancyBboxPatch, Patch

from ._common import FONT_RC, text_width
from ._graph import LEGEND_ENTRY, ROLES, _edges, _reactions, pill, role_handles
from ._layout import (
    drawing_bottom,
    EQUATION_FONT,
    HEADER,
    HEADER_FONT,
    MARGIN,
    SpeciesLayout,
)
from ._themes import NetworkTheme, legend_entries, type_colours

if TYPE_CHECKING:
    from ...core.network.network import Network
    from ...core.species.specie import Specie


def _shown(net: Network, specie: Specie) -> list[str]:
    """Types of the reactions that form or destroy *specie*, in index order."""
    j = specie.index
    change = net.product_matrix[:, j] - net.reactant_matrix[:, j]

    return [r.type for i, r in enumerate(net.reactions) if change[i] != 0]


def species_legend_widths(
    net: Network, specie: Specie, theme: NetworkTheme
) -> list[float]:
    """Width in inches of each entry of the legend of the species view of *specie*."""
    colours = type_colours(theme, [r.type for r in net.reactions])
    labels = [t for t, _ in legend_entries(theme, colours, _shown(net, specie))]

    return [text_width(t, EQUATION_FONT) + LEGEND_ENTRY for t in [*labels, *ROLES]]


def draw_species(
    net: Network, layout: SpeciesLayout, theme: NetworkTheme
) -> tuple[Figure, Axes]:
    """Draw *layout* of *net* in *theme*; the axes are in inches, y down."""
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
        _edges(ax, layout, theme, colours, types)
        _reactions(ax, layout, theme, colours, types)
        _pills(ax, net, layout, theme)
        shown = [types[i] for i in sorted(layout.formation + layout.destruction)]
        _legend(ax, layout, theme, colours, shown)

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
                color=theme.muted,
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
) -> None:
    handles: list[Artist] = [
        Patch(facecolor=c, edgecolor="none", label=t)
        for t, c in legend_entries(theme, colours, shown)
    ]
    handles += role_handles(theme)
    width, height = layout.size
    ax.legend(
        handles=handles,
        loc="lower center",
        # Anchored in inches, so cropping (tight) leaves it in place.
        bbox_to_anchor=(width / 2, height - MARGIN),
        bbox_transform=ax.transData,
        ncol=layout.legend_columns,
        frameon=False,
        fontsize=EQUATION_FONT,
        labelcolor=theme.text,
        handlelength=1.4,
        columnspacing=1.2,
    )
