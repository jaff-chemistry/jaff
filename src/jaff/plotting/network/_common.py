# ABOUTME: Drawing helpers shared by the network figures: text widths, text colour
# ABOUTME: on a fill, Bézier edge paths, and saving / showing / closing a figure
"""Helpers shared by the graph and stoichiometry renderers."""

from __future__ import annotations

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import to_rgb
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties
from matplotlib.path import Path as MplPath
from matplotlib.textpath import TextPath

from ..api import _NON_INTERACTIVE
from ..themes import _BASE_RC
from ._themes import NetworkTheme

#: rcParams the network figures are drawn and saved with (font, maths font).
FONT_RC: dict = _BASE_RC
_FONT = FontProperties(family="DejaVu Sans", math_fontfamily="stix")


def text_width(text: str, size: float) -> float:
    """Width in inches of (mathtext) *text* at *size* points."""
    return TextPath((0, 0), text, size=size, prop=_FONT).get_extents().width / 72


def _luminance(colour: str) -> float:
    r, g, b = to_rgb(colour)

    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def text_on(fill: str, theme: NetworkTheme) -> str:
    """The theme's surface or ink, whichever reads better on *fill*."""
    light, dark = sorted((theme.surface, theme.ink), key=_luminance, reverse=True)

    return dark if _luminance(fill) > 0.5 else light


def bezier_path(points: tuple[tuple[float, float], ...]) -> MplPath:
    """Cubic Béziers through *points*, horizontal at every point."""
    vertices = [points[0]]
    codes = [MplPath.MOVETO]
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        mid = (x0 + x1) / 2
        vertices += [(mid, y0), (mid, y1), (x1, y1)]
        codes += [MplPath.CURVE4] * 3

    return MplPath(vertices, codes)


def finish(
    fig: Figure,
    ax: Axes,
    *,
    theme: NetworkTheme,
    save: str | Path | None,
    show: bool,
    close: bool,
) -> tuple[Figure, Axes]:
    """Save, show and close *fig* as :func:`jaff.plotting.plot` does."""
    if save is not None:
        with plt.rc_context(FONT_RC):
            fig.savefig(
                Path(save), facecolor=theme.surface, bbox_inches="tight", pad_inches=0.05
            )
    shown = show and matplotlib.get_backend().lower() not in _NON_INTERACTIVE
    if shown:
        plt.show()
    # A non-blocking show (interactive mode) leaves the window open: keep it.
    if close and not (shown and plt.isinteractive()):
        plt.close(fig)

    return fig, ax
