# ABOUTME: Draws the stoichiometry view of net.visualize: equations, reaction-index
# ABOUTME: circles and net-change cells under atom-count bands, with a legend
"""Matplotlib renderer of a :class:`MatrixLayout`."""

from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
from matplotlib.artist import Artist
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyBboxPatch, Patch

from ._common import FONT_RC, text_on, text_width
from ._graph import LEGEND_ENTRY
from ._layout import (
    drawing_bottom,
    CELL,
    EQUATION_FONT,
    HEADER,
    HEADER_FONT,
    INDEX_FONT,
    MARGIN,
    SPECIES_FONT,
    MatrixLayout,
)
from ._themes import NetworkTheme, legend_types, type_colours

if TYPE_CHECKING:
    from ...core.network.network import Network

INSET = 0.02
RADIUS = 0.085


def stoichiometry_legend_widths(net: Network, theme: NetworkTheme) -> list[float]:
    """Width in inches of each entry of the stoichiometry legend of *net*."""
    labels = [t for t, _ in legend_types(theme, [r.type for r in net.reactions])]

    return [
        text_width(t, EQUATION_FONT) + LEGEND_ENTRY
        for t in ["consumed", "produced", *labels]
    ]


def draw_stoichiometry(
    net: Network, layout: MatrixLayout, theme: NetworkTheme
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
        _groups(ax, layout, theme)
        _species(ax, net, layout, theme)
        _rows(ax, layout, theme, colours, types)
        _cells(ax, net, layout, theme)
        _legend(ax, layout, theme, types)

    return fig, ax


def _groups(ax, layout, theme) -> None:
    bottom = drawing_bottom(layout)
    for group in layout.groups:
        ax.add_patch(
            FancyBboxPatch(
                (group.x0, MARGIN + HEADER / 2),
                group.x1 - group.x0,
                bottom - MARGIN - HEADER / 2,
                boxstyle="round,pad=0,rounding_size=0.06",
                facecolor=theme.band,
                edgecolor="none",
                zorder=0,
            )
        )
        ax.text(
            (group.x0 + group.x1) / 2,
            MARGIN + HEADER / 4,
            group.header,
            ha="center",
            va="center",
            fontsize=HEADER_FONT,
            color=theme.muted,
        )


def _species(ax, net, layout, theme) -> None:
    for x, index in zip(layout.column_x, layout.species):
        specie = net.species[index]
        ax.text(
            x,
            layout.label_y,
            f"${specie.latex()}$",
            rotation=45,
            ha="center",
            va="center",
            rotation_mode="anchor",
            fontsize=SPECIES_FONT,
            color=theme.ink,
            gid=f"label:{specie.name}",
        )


def _rows(ax, layout, theme, colours, types) -> None:
    for i, (y, label) in enumerate(zip(layout.row_y, layout.labels)):
        colour = colours[types[i]]
        ax.text(
            layout.equation_x,
            y,
            label,
            ha="right",
            va="center",
            fontsize=EQUATION_FONT,
            color=theme.text,
        )
        ax.add_patch(
            Circle(
                (layout.index_x, y),
                RADIUS,
                facecolor=colour,
                edgecolor="none",
                zorder=3,
                gid=f"reaction:{i}",
            )
        )
        ax.text(
            layout.index_x,
            y,
            str(i),
            ha="center",
            va="center",
            fontsize=INDEX_FONT,
            fontweight="bold",
            color=text_on(colour, theme),
            zorder=4,
        )


def _cells(ax, net, layout, theme) -> None:
    size = CELL - 2 * INSET
    for (i, j), value in _nonzero(layout):
        x, y = layout.column_x[j], layout.row_y[i]
        deep = abs(value) > 1
        if value < 0:
            fill, edge, lw = theme.consumed[deep], "none", 0.0
        elif theme.produced_outline:
            fill, edge, lw = theme.surface, theme.outline, 1.6 if deep else 0.8
        else:
            fill, edge, lw = theme.produced[deep], "none", 0.0
        name = net.species[layout.species[j]].name
        ax.add_patch(
            FancyBboxPatch(
                (x - size / 2, y - size / 2),
                size,
                size,
                boxstyle="round,pad=0,rounding_size=0.03",
                facecolor=fill,
                edgecolor=edge,
                lw=lw,
                zorder=2,
                gid=f"cell:{i}:{name}",
            )
        )
        ax.text(
            x,
            y,
            str(abs(value)),
            ha="center",
            va="center",
            fontsize=EQUATION_FONT,
            color=text_on(fill, theme),
            zorder=3,
            gid="count",
        )


def _nonzero(layout: MatrixLayout):
    rows, cols = layout.values.nonzero()
    for i, j in zip(rows.tolist(), cols.tolist()):
        yield (i, j), int(layout.values[i, j])


def _legend(ax, layout, theme, types) -> None:
    produced = (
        Patch(facecolor=theme.surface, edgecolor=theme.outline, label="produced")
        if theme.produced_outline
        else Patch(facecolor=theme.produced[1], edgecolor="none", label="produced")
    )
    handles: list[Artist] = [
        Patch(facecolor=theme.consumed[1], edgecolor="none", label="consumed"),
        produced,
    ]
    handles += [
        Line2D(
            [],
            [],
            ls="",
            marker="o",
            markersize=6,
            markerfacecolor=c,
            markeredgecolor="none",
            label=t,
        )
        for t, c in legend_types(theme, types)
    ]
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
        handlelength=1.2,
        columnspacing=1.2,
    )
