# ABOUTME: Draws the graph view of net.visualize: species pills in atom-count
# ABOUTME: columns, reaction squares in gutters, Bézier edges, headers, legend
"""Matplotlib renderer of a :class:`GraphLayout`."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Sequence

import matplotlib.pyplot as plt
from matplotlib import patheffects
from matplotlib.artist import Artist
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch, PathPatch

from ._common import FONT_RC, bezier_path, text_on, text_width
from ._layout import (
    LEGEND_ROW,
    drawing_bottom,
    EQUATION_FONT,
    HEADER,
    HEADER_FONT,
    INDEX_FONT,
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
#: Fixed edge widths without rate scaling, points.
REACTANT_WIDTH = 0.8
PRODUCT_WIDTH = 1.1
#: Legend space per entry besides its text: handle and column spacing, inches.
LEGEND_ENTRY = 0.4
ROLES = ("reactant", "product")


def draw_graph(
    net: Network,
    layout: GraphLayout,
    theme: NetworkTheme,
    widths: dict[int, float] | None = None,
    key: Sequence[tuple[float, float]] = (),
) -> tuple[Figure, Axes]:
    """Draw *layout* of *net* in *theme*; the axes are in inches, y down.

    *widths* gives each reaction's edge width (rate scaling, see
    :func:`~jaff.plotting.network._widths.edge_widths`) and *key* the legend's
    ``(rate, width)`` reference entries; without them the edges have fixed widths.
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
        _columns(ax, layout, theme)
        _edges(ax, layout, theme, colours, types, widths)
        _reactions(ax, layout, theme, colours, types, label_halo(widths))
        _species(ax, net, layout, theme)
        _legend(ax, layout, theme, types, key)

    return fig, ax


def legend_widths(net: Network, theme: NetworkTheme) -> list[float]:
    """Width in inches of each entry of the graph legend of *net* in *theme*."""
    labels = [t for t, _ in legend_types(theme, [r.type for r in net.reactions])]

    return entry_widths([*labels, *ROLES])


def entry_widths(labels: Sequence[str]) -> list[float]:
    """Legend entry widths in inches."""
    return [text_width(t, EQUATION_FONT) + LEGEND_ENTRY for t in labels]


def label_halo(widths: dict[int, float] | None) -> float:
    """Halo width behind equation labels: wider when thick arrows may cross them."""
    return 3.0 if widths is None else 3.0 + max(widths.values(), default=0.0)


def key_label(value: float) -> str:
    """Mathtext label of a width-key rate: ``10^-9``, or ``2.3×10^-9`` off a decade."""
    exponent = math.floor(math.log10(value))
    mantissa = value / 10.0**exponent
    if round(mantissa, 1) == 1.0:
        return rf"$k = 10^{{{exponent}}}$"

    return rf"$k = {mantissa:.1f}\times10^{{{exponent}}}$"


def _columns(ax: Axes, layout: GraphLayout, theme: NetworkTheme) -> None:
    bottom = drawing_bottom(layout)
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


def _edges(ax, layout, theme, colours, types, widths=None) -> None:
    """Reactant lines and product arrows; *widths* (by reaction) scales the arrows."""
    for edge in layout.edges:
        path = bezier_path(edge.points)
        typed = colours[types[edge.reaction]]
        gid = f"edge:{edge.role}:{edge.reaction}"
        if edge.role == "reactant":
            ax.add_patch(
                PathPatch(
                    path,
                    facecolor="none",
                    edgecolor=theme.edge,
                    lw=REACTANT_WIDTH,
                    zorder=1,
                    gid=gid,
                )
            )
            continue
        lw = PRODUCT_WIDTH if widths is None else widths[edge.reaction]
        line = typed if theme.typed_edges else theme.edge
        head = typed if theme.typed_arrows else line
        ax.add_patch(
            PathPatch(path, facecolor="none", edgecolor=line, lw=lw, zorder=1, gid=gid)
        )
        (xa, _), (xb, yb) = edge.points[-2], edge.points[-1]
        step = ARROW if xb > xa else -ARROW
        ax.add_patch(
            FancyArrowPatch(
                (xb - step, yb),
                (xb, yb),
                arrowstyle="-|>",
                # Heads grow with thick lines so they stay visible past the line.
                mutation_scale=max(11.0, 11.0 + 3.0 * (lw - PRODUCT_WIDTH)),
                color=head,
                lw=0,
                shrinkA=0,
                shrinkB=0,
                zorder=2,
            )
        )


def _reactions(ax, layout, theme, colours, types, halo_width: float = 3.0) -> None:
    """Reaction squares with their index; equations above, on a *halo_width* halo."""
    halo = [patheffects.withStroke(linewidth=halo_width, foreground=theme.halo)]
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
        pill(ax, x, y, specie, layout.pill_width[index], theme, f"specie:{specie.name}")


def pill(ax, x, y, specie, width, theme, gid, bold=False) -> None:
    """A species pill centred on ``(x, y)``; *bold* draws a heavier outline."""
    ax.add_patch(
        FancyBboxPatch(
            (x - width / 2, y - PILL_HEIGHT / 2),
            width,
            PILL_HEIGHT,
            boxstyle=f"round,pad=0,rounding_size={PILL_HEIGHT / 2}",
            facecolor=theme.surface,
            edgecolor=theme.ink,
            lw=1.8 if bold else 0.8,
            zorder=3,
            gid=gid,
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


def role_handles(theme) -> list[Artist]:
    """Legend handles for reactant lines and product arrows."""
    return [
        Line2D([], [], color=theme.edge, lw=REACTANT_WIDTH, label="reactant"),
        Line2D(
            [],
            [],
            color=theme.ink,
            lw=PRODUCT_WIDTH,
            marker=">",
            markersize=6,
            label="product",
        ),
    ]


def _legend(ax, layout, theme, types, key=()) -> None:
    handles: list[Artist] = [
        Patch(facecolor=c, edgecolor="none", label=t)
        for t, c in legend_types(theme, types)
    ]
    draw_legends(ax, layout, theme, [*handles, *role_handles(theme)], key)


def draw_legends(
    ax: Axes,
    layout,
    theme: NetworkTheme,
    handles: list[Artist],
    key: Sequence[tuple[float, float]],
) -> None:
    """The type legend, with the width key on its own row below it.

    Both are anchored in inches, so cropping (tight) leaves them in place.
    """
    width, height = layout.size
    style = {
        "loc": "lower center",
        "bbox_transform": ax.transData,
        "frameon": False,
        "fontsize": EQUATION_FONT,
        "labelcolor": theme.text,
        "columnspacing": 1.2,
    }
    bottom = height - MARGIN
    if key:
        lines: list[Artist] = [
            Line2D([], [], color=theme.ink, lw=w, label=key_label(value))
            for value, w in key
        ]
        ax.add_artist(
            ax.legend(
                handles=lines,
                bbox_to_anchor=(width / 2, bottom),
                ncol=len(lines),
                handlelength=2.2,
                # Tight padding: the key's row is one LEGEND_ROW tall.
                borderpad=0.1,
                **style,
            )
        )
        bottom -= LEGEND_ROW
    ax.legend(
        handles=handles,
        bbox_to_anchor=(width / 2, bottom),
        ncol=layout.legend_columns,
        handlelength=1.4,
        **style,
    )
