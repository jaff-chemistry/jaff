# ABOUTME: Draws the graph view of net.visualize: species pills in atom-count
# ABOUTME: columns, reaction markers in gutters, Bézier edges and the legend panel
"""Matplotlib renderer of a :class:`GraphLayout`."""

from __future__ import annotations

from typing import TYPE_CHECKING, Sequence

import matplotlib.pyplot as plt
from matplotlib import patheffects
from matplotlib.artist import Artist
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, PathPatch

from ._common import FONT_RC, bezier_path, heading, text_width
from ._layout import (
    PANEL_GAP,
    PANEL_HEAD,
    PANEL_SPACE,
    RATE_BODY,
    TYPE_ROW,
    drawing_bottom,
    EQUATION_FONT,
    HEADER,
    HEADER_FONT,
    MARGIN,
    MARKER,
    PILL_HEIGHT,
    SPECIES_FONT,
    GraphLayout,
)
from ._markers import draw_marker, label_shapes, marker_handle, type_shapes
from ._rates import RateColours, draw_colour_bar
from ._values import draw_values
from ._themes import NetworkTheme, heading_colour, legend_types, type_colours

if TYPE_CHECKING:
    from ...core.network.network import Network

ARROW = 0.08
#: Edge widths, points: reactant lines, product arrows, rate-coloured arrows.
REACTANT_WIDTH = 0.8
PRODUCT_WIDTH = 1.1
RATED_WIDTH = 1.7
#: Halo behind equations, points, so lines passing under stay clear of them.
HALO = 3.0
#: Legend space per entry besides its text: handle and column spacing, inches.
LEGEND_ENTRY = 0.4
ROLES = ("reactant", "product")


def draw_graph(
    net: Network,
    layout: GraphLayout,
    theme: NetworkTheme,
    rated: RateColours | None = None,
) -> tuple[Figure, Axes]:
    """Draw *layout* of *net* in *theme*; the axes are in inches, y down.

    With *rated*, markers and product arrows take each reaction's rate colour
    and a colour bar sits below the legend; otherwise they take the type colours.
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
        _edges(ax, layout, theme, colours, types, rated)
        _reactions(ax, layout, theme, colours, types, rated)
        _species(ax, net, layout, theme)
        _legend(ax, layout, theme, types, rated)

    return fig, ax


def legend_widths(net: Network, theme: NetworkTheme) -> list[float]:
    """Width in inches of each entry of the graph legend of *net* in *theme*."""
    labels = [t for t, _ in legend_types(theme, [r.type for r in net.reactions])]

    return entry_widths([*labels, *ROLES])


def entry_widths(labels: Sequence[str]) -> list[float]:
    """Legend entry widths in inches."""
    return [text_width(t, EQUATION_FONT) + LEGEND_ENTRY for t in labels]


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
            color=heading_colour(theme),
        )


def _edges(ax, layout, theme, colours, types, rated=None) -> None:
    """Reactant lines and product arrows; *rated* colours the arrows by rate."""
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
        if rated is None:
            lw = PRODUCT_WIDTH
            line = typed if theme.typed_edges else theme.edge
            head = typed if theme.typed_arrows else line
        else:
            lw = RATED_WIDTH
            line = head = rated.colours[edge.reaction]
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
                # Heads grow with the line so they stay visible past it.
                mutation_scale=11.0 if rated is None else 12.0,
                color=head,
                lw=0,
                shrinkA=0,
                shrinkB=0,
                zorder=2,
            )
        )


def _reactions(ax, layout, theme, colours, types, rated=None) -> None:
    """Reaction markers with their index, equations above on a halo.

    The shape follows the type; the fill is the type colour, or with *rated*
    the rate colour inside an ink outline.
    """
    halo = [patheffects.withStroke(linewidth=HALO, foreground=theme.halo)]
    shapes = type_shapes(types, len(theme.types))
    for i, (x, y) in layout.reaction_xy.items():
        t = types[i]
        fill = colours[t] if rated is None else rated.colours[i]
        edge = "none" if rated is None else theme.ink
        draw_marker(ax, x, y, shapes[t], fill, theme, i, gid=f"reaction:{i}", edge=edge)
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


def type_handles(theme, entries, shapes, rated) -> list[Artist]:
    """Legend handles of the reaction types: their shape, filled with their colour.

    With rate colouring the colour shows the rate, so the shapes are drawn
    empty (surface fill, ink outline).
    """
    if rated is None:
        return [marker_handle(t, shapes[t], c) for t, c in entries]

    return [marker_handle(t, shapes[t], theme.surface, theme.ink) for t, _ in entries]


def _legend(ax, layout, theme, types, rated=None) -> None:
    entries = legend_types(theme, types)
    handles = type_handles(theme, entries, label_shapes(theme, types), rated)
    draw_legends(ax, layout, theme, [*handles, *role_handles(theme)], rated)


def draw_legends(
    ax: Axes,
    layout,
    theme: NetworkTheme,
    handles: list[Artist],
    rated: RateColours | None,
) -> None:
    """The legend panel under the drawing (see :func:`legend_panel`).

    Top to bottom: the rate colour bar when *rated*, the REACTION TYPE legend of
    *handles*, then the analysis values *rated* lists. Everything is placed in
    inches, so cropping (tight) leaves it in place.
    """
    panel = layout.panel
    top = drawing_bottom(layout) + PANEL_SPACE
    if rated is not None:
        draw_colour_bar(
            ax, layout.size[0], top, theme, rated.cmap, rated.span, rated.quantity
        )
        top += PANEL_HEAD + RATE_BODY + PANEL_GAP
    heading(ax, panel.x0, top, "REACTION TYPE", theme)
    ax.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(panel.x0, top + PANEL_HEAD),
        bbox_transform=ax.transData,
        ncol=layout.legend_columns,
        frameon=False,
        borderpad=0.0,
        borderaxespad=0.0,
        fontsize=EQUATION_FONT,
        labelcolor=theme.text,
        handlelength=1.4,
        columnspacing=1.2,
    )
    top += PANEL_HEAD + layout.legend_rows * TYPE_ROW
    if rated is not None and rated.values:
        draw_values(ax, panel, theme, rated.values, top + PANEL_GAP)
