# ABOUTME: Rate colouring for net.visualize: log positions of rate coefficients on
# ABOUTME: a colormap (12-decade floor), decade ticks and the rounded colour bar
"""Colour reactions by rate coefficient in the graph and species views."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from matplotlib import colormaps
from matplotlib.axes import Axes
from matplotlib.colors import Colormap, LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch

from ._layout import COLOUR_BAR_ROWS, EQUATION_FONT, LEGEND_ROW, MARGIN, ValueGroup
from ._themes import NetworkTheme

#: Rates more than this many decades below the largest get the lowest colour.
DECADES = 12
#: Most tick labels on the colour bar.
MAX_TICKS = 6
#: Colour bar title of each quantity.
TITLES = {
    "rate_coefficient": "rate coefficient $k$",
    "rate": r"reaction rate $\mathcal{R}$ (cm$^{-3}$ s$^{-1}$)",
}
#: Colour bar: longest bar, bar thickness, inches.
BAR_LENGTH = 4.0
BAR_HEIGHT = 0.1


@dataclass(frozen=True)
class RateColours:
    """Rate colouring of one figure.

    Attributes
    ----------
    colours : dict of int to str
        Colour of each drawn reaction, by index in ``net.reactions``.
    cmap : Colormap
        The colormap, slow (0) to fast (1).
    span : tuple of float
        ``(log10 low, log10 top)``: the values at the two ends of *cmap*.
    quantity : str
        What is coloured: ``"rate_coefficient"`` or ``"rate"``.
    values : tuple of ValueGroup
        Groups of the values legend under the colour bar (none: no legend).
    """

    colours: dict[int, str]
    cmap: Colormap
    span: tuple[float, float]
    quantity: str = "rate_coefficient"
    values: tuple[ValueGroup, ...] = ()


def rate_positions(values: np.ndarray) -> tuple[np.ndarray, tuple[float, float] | None]:
    """Positions in ``[0, 1]`` of rate *values* on the colormap, and the log range.

    The largest value is at 1 and the smallest at 0, logarithmically between;
    a value more than :data:`DECADES` decades below the largest sits at 0, as
    do zero, negative and non-finite values. Equal values all sit at 1. The
    range is ``(log10 low, log10 top)``, or ``None`` without a positive value.
    """
    values = np.asarray(values, dtype=float)
    good = np.isfinite(values) & (values > 0)
    positions = np.zeros(values.shape)
    if not good.any():
        return positions, None

    logs = np.log10(values[good])
    top = float(logs.max())
    low = max(float(logs.min()), top - DECADES)
    if top == low:
        positions[good] = 1.0
    else:
        positions[good] = np.clip((logs - low) / (top - low), 0.0, 1.0)

    return positions, (low, top)


def decade_ticks(low: float, top: float) -> list[int]:
    """Up to :data:`MAX_TICKS` powers of ten (exponents) in ``[low, top]``.

    Every n-th decade, counted down from the highest so it is always labelled.
    """
    decades = list(range(math.ceil(low), math.floor(top) + 1))
    step = max(1, math.ceil(len(decades) / MAX_TICKS))

    return decades[::-1][::step][::-1]


def rate_colormap(theme: NetworkTheme, cmap: str | Colormap | None = None) -> Colormap:
    """The colormap for rates: *cmap* as given, else the theme's, trimmed."""
    if isinstance(cmap, Colormap):
        return cmap
    if cmap is not None:
        return colormaps[cmap]

    spec = theme.rates
    base = colormaps[spec.name]

    return LinearSegmentedColormap.from_list(
        f"{spec.name}_rates", base(np.linspace(spec.low, spec.high, 256))
    )


def draw_colour_bar(
    ax: Axes,
    size: tuple[float, float],
    theme: NetworkTheme,
    cmap: Colormap,
    span: tuple[float, float],
    quantity: str = "rate_coefficient",
    bottom: float | None = None,
) -> None:
    """A rounded horizontal colour bar of log10 k, centred at the figure's bottom.

    It fills the ``COLOUR_BAR_ROWS`` legend rows above *bottom* (by default the
    bottom margin): the title, the bar, then a power-of-ten label under each tick.
    """
    width, height = size
    length = min(BAR_LENGTH, width - 2 * MARGIN)
    x0 = (width - length) / 2
    if bottom is None:
        bottom = height - MARGIN
    y0 = bottom - COLOUR_BAR_ROWS * LEGEND_ROW + 0.16
    image = ax.imshow(
        np.linspace(0.0, 1.0, 256)[None, :],
        cmap=cmap,
        aspect="auto",
        extent=(x0, x0 + length, y0 + BAR_HEIGHT, y0),
        interpolation="bilinear",
        zorder=2,
    )
    clip = FancyBboxPatch(
        (x0, y0),
        length,
        BAR_HEIGHT,
        boxstyle=f"round,pad=0,rounding_size={BAR_HEIGHT / 2}",
        facecolor="none",
        edgecolor="none",
        transform=ax.transData,
    )
    ax.add_patch(clip)
    image.set_clip_path(clip)
    ax.text(
        x0 + length / 2,
        y0 - 0.04,
        TITLES[quantity],
        ha="center",
        va="bottom",
        fontsize=EQUATION_FONT,
        color=theme.text,
    )
    low, top = span
    for decade in decade_ticks(low, top):
        x = x0 + (decade - low) / (top - low) * length if top > low else x0 + length
        ax.plot(
            [x, x], [y0 + BAR_HEIGHT, y0 + BAR_HEIGHT + 0.035], color=theme.muted, lw=0.6
        )
        ax.text(
            x,
            y0 + BAR_HEIGHT + 0.05,
            rf"$10^{{{decade}}}$",
            ha="center",
            va="top",
            fontsize=EQUATION_FONT - 0.5,
            color=theme.text,
        )
