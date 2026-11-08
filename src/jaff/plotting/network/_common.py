# ABOUTME: Drawing helpers shared by the network figures: text widths, text colour
# ABOUTME: on a fill, Bézier edge paths, and saving / showing / closing a figure
"""Helpers shared by the graph and stoichiometry renderers."""

from __future__ import annotations

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.artist import Artist
from matplotlib.axes import Axes
from matplotlib.colors import to_rgb
from matplotlib.figure import Figure
from matplotlib.legend import Legend
from matplotlib.font_manager import FontProperties
from matplotlib.path import Path as MplPath
from matplotlib.textpath import TextPath

from ..api import _NON_INTERACTIVE, _VECTOR
from ..themes import _BASE_RC
from ._themes import NetworkTheme

#: rcParams the network figures are drawn and saved with (font, maths font).
FONT_RC: dict = _BASE_RC
#: Resolution of raster files (PNG, JPEG, ...), as in the ``paper`` plot style;
#: vector formats (PDF, SVG, ...) are resolution-free.
SAVE_DPI = 400
#: Space left around the drawn content of a ``tight`` figure, inches.
TIGHT_PAD = 0.03
_FONT = FontProperties(family="DejaVu Sans", math_fontfamily="stix")


def text_width(text: str, size: float) -> float:
    """Width in inches of (mathtext) *text* at *size* points."""
    return TextPath((0, 0), text, size=size, prop=_FONT).get_extents().width / 72


def _luminance(colour: str) -> float:
    """Relative luminance (WCAG) of *colour*."""
    r, g, b = (
        c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
        for c in to_rgb(colour)
    )

    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def _contrast(a: str, b: str) -> float:
    """WCAG contrast ratio of colours *a* and *b* (1 to 21)."""
    high, low = sorted((_luminance(a), _luminance(b)), reverse=True)

    return (high + 0.05) / (low + 0.05)


def text_on(fill: str, theme: NetworkTheme) -> str:
    """The theme's surface or ink, whichever contrasts more with *fill*."""
    return max((theme.ink, theme.surface), key=lambda text: _contrast(fill, text))


def bezier_path(points: tuple[tuple[float, float], ...]) -> MplPath:
    """Cubic Béziers through *points*, horizontal at every point."""
    vertices = [points[0]]
    codes = [MplPath.MOVETO]
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        mid = (x0 + x1) / 2
        vertices += [(mid, y0), (mid, y1), (x1, y1)]
        codes += [MplPath.CURVE4] * 3

    return MplPath(vertices, codes)


def crop(fig: Figure, ax: Axes, pad: float = TIGHT_PAD) -> None:
    """Shrink *fig* to its drawn content plus *pad* inches on each side.

    The axes span the figure in inches (y down), so cropping sets their limits
    and the figure size to the content's bounding box; artists anchored in data
    coordinates, the legends included, keep their place.
    """
    # Every legend, and images such as the rate colour bar.
    legends = [a for a in ax.get_children() if isinstance(a, Legend)]
    artists: list[Artist] = [*ax.patches, *ax.texts, *ax.images, *legends]
    to_data = ax.transData.inverted()
    corners = [to_data.transform(a.get_window_extent().corners()) for a in artists]
    xs = [x for c in corners for x, _ in c]
    ys = [y for c in corners for _, y in c]
    x0, x1 = min(xs) - pad, max(xs) + pad
    y0, y1 = min(ys) - pad, max(ys) + pad
    fig.set_size_inches(x1 - x0, y1 - y0)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)


def finish(
    fig: Figure,
    ax: Axes,
    *,
    theme: NetworkTheme,
    save: str | Path | None,
    show: bool,
    close: bool,
    tight: bool = False,
) -> tuple[Figure, Axes]:
    """Crop (when *tight*), save, show and close *fig* as :func:`jaff.plotting.plot` does."""
    if tight:
        crop(fig, ax)
    if save is not None:
        path = Path(save)
        dpi = None if path.suffix.lower() in _VECTOR else SAVE_DPI
        # A cropped figure is saved as is; otherwise trim to the drawn area.
        bbox = None if tight else "tight"
        with plt.rc_context(FONT_RC):
            fig.savefig(
                path,
                dpi=dpi,
                facecolor=theme.surface,
                bbox_inches=bbox,
                pad_inches=0.05,
            )
    shown = show and matplotlib.get_backend().lower() not in _NON_INTERACTIVE
    if shown:
        plt.show()
    # A non-blocking show (interactive mode) leaves the window open: keep it.
    if close and not (shown and plt.isinteractive()):
        plt.close(fig)

    return fig, ax
