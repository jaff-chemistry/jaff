# ABOUTME: Draws the graph view of net.visualize: species pills in atom-count
# ABOUTME: columns, reaction squares in gutters, Bézier edges, headers, legend
"""Matplotlib renderer of a :class:`GraphLayout`."""

from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
from matplotlib import patheffects
from matplotlib.artist import Artist
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch, PathPatch

from ._common import FONT_RC, bezier_path, text_on, text_width
from ._layout import (
    EQUATION_FONT,
    HEADER,
    HEADER_FONT,
    INDEX_FONT,
    LEGEND,
    MARGIN,
    MARKER,
    PILL_HEIGHT,
    SPECIES_FONT,
    GraphLayout,
)
from ._themes import NetworkTheme, legend_types, type_colours

if TYPE_CHECKING:
    from ...core.network.network import Network

ARROW = 0.08
#: Legend space per entry besides its text: handle and column spacing, inches.
LEGEND_ENTRY = 0.4
ROLES = ("reactant", "product")


def draw_graph(
    net: Network, layout: GraphLayout, theme: NetworkTheme
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
        _columns(ax, layout, theme)
        _edges(ax, layout, theme, colours, types)
        _reactions(ax, layout, theme, colours, types)
        _species(ax, net, layout, theme)
        _legend(ax, layout, theme, types)

    return fig, ax


def legend_width(net: Network, theme: NetworkTheme) -> float:
    """Width in inches of the graph legend of *net* in *theme*."""
    labels = [t for t, _ in legend_types(theme, [r.type for r in net.reactions])]

    return sum(text_width(t, EQUATION_FONT) + LEGEND_ENTRY for t in [*labels, *ROLES])


def _columns(ax: Axes, layout: GraphLayout, theme: NetworkTheme) -> None:
    bottom = layout.size[1] - LEGEND - MARGIN
    for lane in (lane for lane in layout.lanes if lane.header is not None):
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
        ax.text(
            lane.centre,
            MARGIN + HEADER / 4,
            lane.header or "",
            ha="center",
            va="center",
            fontsize=HEADER_FONT,
            color=theme.muted,
        )


def _edges(ax, layout, theme, colours, types) -> None:
    for edge in layout.edges:
        path = bezier_path(edge.points)
        typed = colours[types[edge.reaction]]
        if edge.role == "reactant":
            ax.add_patch(
                PathPatch(path, facecolor="none", edgecolor=theme.edge, lw=0.8, zorder=1)
            )
            continue
        line = typed if theme.typed_edges else theme.edge
        head = typed if theme.typed_arrows else line
        ax.add_patch(PathPatch(path, facecolor="none", edgecolor=line, lw=1.1, zorder=1))
        (xa, _), (xb, yb) = edge.points[-2], edge.points[-1]
        step = ARROW if xb > xa else -ARROW
        ax.add_patch(
            FancyArrowPatch(
                (xb - step, yb),
                (xb, yb),
                arrowstyle="-|>",
                mutation_scale=11,
                color=head,
                lw=0,
                shrinkA=0,
                shrinkB=0,
                zorder=2,
            )
        )


def _reactions(ax, layout, theme, colours, types) -> None:
    halo = [patheffects.withStroke(linewidth=3, foreground=theme.halo)]
    for i, (x, y) in layout.reaction_xy.items():
        colour = colours[types[i]]
        ax.add_patch(
            FancyBboxPatch(
                (x - MARKER / 2, y - MARKER / 2),
                MARKER,
                MARKER,
                boxstyle="round,pad=0,rounding_size=0.025",
                facecolor=colour,
                edgecolor="none",
                zorder=3,
                gid=f"reaction:{i}",
            )
        )
        ax.text(
            x,
            y,
            str(i),
            ha="center",
            va="center",
            fontsize=INDEX_FONT,
            fontweight="bold",
            color=text_on(colour, theme),
            zorder=4,
        )
        ax.text(
            x,
            y - MARKER / 2 - 0.04,
            layout.labels[i],
            ha="center",
            va="bottom",
            fontsize=EQUATION_FONT,
            color=theme.text,
            path_effects=halo,
            zorder=4,
        )


def _species(ax, net, layout, theme) -> None:
    for index, (x, y) in layout.species_xy.items():
        specie = net.species[index]
        w = layout.pill_width[index]
        ax.add_patch(
            FancyBboxPatch(
                (x - w / 2, y - PILL_HEIGHT / 2),
                w,
                PILL_HEIGHT,
                boxstyle=f"round,pad=0,rounding_size={PILL_HEIGHT / 2}",
                facecolor=theme.surface,
                edgecolor=theme.ink,
                lw=0.8,
                zorder=3,
                gid=f"specie:{specie.name}",
            )
        )
        ax.text(
            x,
            y,
            f"${specie.latex()}$",
            ha="center",
            va="center",
            fontsize=SPECIES_FONT,
            color=theme.ink,
            zorder=4,
        )


def _legend(ax, layout, theme, types) -> None:
    handles: list[Artist] = [
        Patch(facecolor=c, edgecolor="none", label=t)
        for t, c in legend_types(theme, types)
    ]
    handles += [
        Line2D([], [], color=theme.edge, lw=0.8, label="reactant"),
        Line2D(
            [], [], color=theme.ink, lw=1.1, marker=">", markersize=6, label="product"
        ),
    ]
    width, height = layout.size
    ax.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, MARGIN / height),
        ncol=len(handles),
        frameon=False,
        fontsize=EQUATION_FONT,
        labelcolor=theme.text,
        handlelength=1.4,
        columnspacing=1.2,
    )
