# ABOUTME: Reaction-type marker shapes of the network figures: assignment by first
# ABOUTME: appearance, drawing a numbered marker, and matching legend handles
"""Marker shapes encoding the reaction type in the graph and species views."""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

from matplotlib.axes import Axes
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyBboxPatch, Patch, RegularPolygon

from ._common import text_on
from ._layout import INDEX_FONT, MARKER
from ._themes import NetworkTheme

#: Shapes given to reaction types in order of first appearance (``TYPE_LIMIT``).
SHAPES = (
    "square",
    "circle",
    "diamond",
    "triangle",
    "hexagon",
    "pentagon",
    "triangle_down",
)
#: Shape of ``unknown`` reactions and of types beyond :data:`SHAPES`.
OTHER_SHAPE = "octagon"
#: Matplotlib marker of each shape, for legend handles.
MPL_MARKERS = {
    "square": "s",
    "circle": "o",
    "diamond": "D",
    "triangle": "^",
    "hexagon": "h",
    "pentagon": "p",
    "triangle_down": "v",
    "octagon": "8",
}

# Regular polygons: (sides, radius, orientation, centre shift, index shift); radius
# and shifts in units of MARKER, shifts down (the axes point y down, so orientation
# pi puts a vertex at the top). Every shape stays within about 0.68 MARKER of the
# centre upwards, clear of the equation above it.
_POLYGONS = {
    "diamond": (4, 0.68, 0.0, 0.0, 0.0),
    "triangle": (3, 0.78, math.pi, 0.12, 0.12),
    "triangle_down": (3, 0.78, 0.0, -0.12, -0.12),
    "hexagon": (6, 0.64, 0.0, 0.0, 0.0),
    "pentagon": (5, 0.68, math.pi, 0.04, 0.1),
    "octagon": (8, 0.62, math.pi / 8, 0.0, 0.0),
}


def type_shapes(types: Sequence[str], limit: int | None = None) -> dict[str, str]:
    """Shape of each reaction type in *types*, by first appearance.

    ``unknown`` and the types after the first *limit* (at most ``len(SHAPES)``)
    get :data:`OTHER_SHAPE`.
    """
    count = len(SHAPES) if limit is None else min(limit, len(SHAPES))
    ordered = [t for t in dict.fromkeys(types) if t != "unknown"]
    shapes = {t: SHAPES[k] if k < count else OTHER_SHAPE for k, t in enumerate(ordered)}
    if "unknown" in types:
        shapes["unknown"] = OTHER_SHAPE

    return shapes


def draw_marker(
    ax: Axes,
    x: float,
    y: float,
    shape: str,
    fill: str,
    theme: NetworkTheme,
    index: int,
    gid: str | None = None,
    edge: str = "none",
) -> Patch:
    """A *shape* marker filled with *fill* at ``(x, y)``, numbered *index*.

    The patch carries *gid* and its shape as its label; the number is drawn in
    the theme's ink or surface, whichever reads better on *fill*.
    """
    style: dict[str, Any] = {
        "facecolor": fill,
        "edgecolor": edge,
        "lw": 0.8,
        "zorder": 3,
        "gid": gid,
    }
    offset = 0.0
    patch: Patch
    if shape == "square":
        patch = FancyBboxPatch(
            (x - MARKER / 2, y - MARKER / 2),
            MARKER,
            MARKER,
            boxstyle="round,pad=0,rounding_size=0.025",
            **style,
        )
    elif shape == "circle":
        patch = Circle((x, y), MARKER * 0.6, **style)
    else:
        sides, size, orientation, centre, lift = _POLYGONS[shape]
        offset = lift * MARKER
        patch = RegularPolygon(
            (x, y + centre * MARKER),
            sides,
            radius=size * MARKER,
            orientation=orientation,
            **style,
        )
    patch.set_label(shape)
    ax.add_patch(patch)
    ax.text(
        x,
        y + offset,
        str(index),
        ha="center",
        va="center",
        fontsize=INDEX_FONT,
        fontweight="bold",
        color=text_on(fill, theme),
        zorder=4,
    )

    return patch


def label_shapes(theme: NetworkTheme, types: Sequence[str]) -> dict[str, str]:
    """Shape per legend label (``other`` included) for a network's *types*."""
    shapes = {
        t.replace("_", " "): s
        for t, s in type_shapes(types, len(theme.types)).items()
        if s != OTHER_SHAPE
    }

    return shapes | {"other": OTHER_SHAPE}


def marker_handle(label: str, shape: str, fill: str, edge: str = "none") -> Line2D:
    """Legend handle: a *shape* marker filled with *fill*, labelled *label*."""
    return Line2D(
        [],
        [],
        ls="none",
        marker=MPL_MARKERS[shape],
        markersize=8,
        markerfacecolor=fill,
        markeredgecolor=edge,
        label=label,
    )
