# ABOUTME: Tick helpers: exponent-labelled log axes (decade / fractional rule)
# ABOUTME: and y tick labels drawn inside the panel with an x strip kept clear

"""Tick helpers.

Exponent-labelled log axes (``10**k`` shown as ``k``): :func:`final_exp_ticks`
puts at most 4 labelled decades on a log axis (every 1, 2, 3, 4, 5, 10, 20,
25, 50 or 100 decades, see :func:`decade_ticks`), unlabelled minor ticks on
decades every :func:`minor_every` (step) decades and no 2..9 sub-ticks.  A
log axis spanning less than :data:`NARROW_DECADES` decades has too few decade
ticks to read well, so :func:`refine_exp_ticks` swaps in ticks at "nice"
fractional exponents (``10**1.05`` shown as ``1.05``) when it finds at least
two distinct labels.

Inside y tick labels (:func:`inside_y_labels`) sit just right of the left
spine; a label that would stick out past the top or bottom of its panel is
blanked (:class:`EdgeTrim`), and :func:`clear_strip` widens the x range on the
left so no curve runs under a label.
"""

from __future__ import annotations

import re

import numpy as np
from matplotlib.axes import Axes
from matplotlib.axis import Axis
from matplotlib.backend_bases import RendererBase
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import (
    FixedLocator,
    Formatter,
    FuncFormatter,
    LogLocator,
    MaxNLocator,
    NullFormatter,
    NullLocator,
)

#: Log axes spanning fewer decades than this get fractional-exponent ticks.
NARROW_DECADES: float = 1.5
_NICE_STEPS = [1, 2, 2.5, 5, 10]
_DECADE_STEPS = (1, 2, 3, 4, 5, 10, 20, 25, 50, 100)


def exp_label(value: float, _pos: int | None = None) -> str:
    """Label a decade tick ``10**k`` as ``k`` (true minus sign); blank otherwise."""
    if value <= 0:
        return ""
    k = np.log10(value)
    if abs(k - round(k)) > 1e-6:
        return ""
    return f"{int(round(k))}".replace("-", "−")


def apply_exp_ticks(axis: Axis) -> None:
    """Put decade major ticks with exponent labels on a log *axis*.

    Not used by the renderer, which sets exponent ticks with
    :func:`final_exp_ticks`.
    """
    axis.set_major_locator(LogLocator(base=10.0, numticks=6))
    axis.set_major_formatter(FuncFormatter(exp_label))
    axis.set_minor_formatter(NullFormatter())


def fractional_labels(exponents: np.ndarray) -> list[str]:
    """Label exponents with just enough decimals to show each exactly.

    Trailing zeros are dropped (``1.10`` -> ``1.1``) and negatives use a true
    minus sign.
    """
    decimals = next(
        (
            d
            for d in range(7)
            if np.allclose(np.round(exponents, d), exponents, rtol=0, atol=1e-9)
        ),
        6,
    )
    labels = []
    for k in exponents:
        text = f"{k:.{decimals}f}"
        if "." in text:
            text = text.rstrip("0").rstrip(".")
        if text in ("-0", ""):
            text = "0"
        labels.append(text.replace("-", "−"))
    return labels


def _nice_exponents(lo: float, hi: float) -> np.ndarray:
    """At least two "nice" exponents inside ``[lo, hi]`` (fewer only if impossible)."""
    tol = 1e-9 * max(1.0, abs(lo), abs(hi))
    ks = np.array([])
    for nbins in (4, 8, 16):
        ks = MaxNLocator(nbins=nbins, steps=_NICE_STEPS).tick_values(lo, hi)
        ks = np.round(ks[(ks >= lo - tol) & (ks <= hi + tol)], 10)
        if len(ks) >= 2:
            break
    return ks


def refine_exp_ticks(axis: Axis) -> None:
    """Fractional-exponent ticks on a log *axis* spanning < 1.5 decades.

    Call once the axis limits are final: the ticks are fixed positions, so
    later limit changes do not re-tick.  Wider (or non-log) axes, and axes too
    narrow for distinct labels, are left alone; minor ticks are switched off on
    a refined axis.
    """
    if axis.get_scale() != "log":
        return
    lo, hi = sorted(axis.get_view_interval())
    if not (np.isfinite(lo) and np.isfinite(hi)) or lo <= 0:
        return
    a, b = np.log10(lo), np.log10(hi)
    if b - a >= NARROW_DECADES:
        return
    ks = _nice_exponents(a, b)
    if len(ks) < 2:
        return
    labels = fractional_labels(ks)
    if len(set(labels)) != len(labels):
        return  # span below label precision: keep the decade ticks

    def label(value: float, _pos: int | None = None) -> str:
        if value <= 0:
            return ""
        hit = np.flatnonzero(np.isclose(ks, np.log10(value), rtol=0, atol=1e-9))
        return labels[hit[0]] if hit.size else ""

    axis.set_major_locator(FixedLocator(10.0**ks))
    axis.set_major_formatter(FuncFormatter(label))
    axis.set_minor_locator(NullLocator())
    axis.set_minor_formatter(NullFormatter())


def log10_label(label: str) -> str:
    """Prefix an axis label with ``log10`` (empty labels stay empty)."""
    return rf"$\log_{{10}}$ {label}" if label else ""


#: Minor-tick spacings (decades) tried for a major step, smallest first.
_MINOR_STEPS = (1, 2, 5, 10, 25, 50)
#: Most minor ticks between two consecutive labelled decades.
MAX_MINORS = 4


def minor_every(step: int) -> int:
    """Decades between minor ticks for a labelled step of *step* decades.

    The smallest spacing in :data:`_MINOR_STEPS` that divides *step* and
    leaves at most :data:`MAX_MINORS` minors between majors (never under a
    decade): 1 for steps up to 5, 2 for 10, 5 for 20 / 25, 10 for 50, 25
    for 100.
    """
    for m in _MINOR_STEPS:
        if step % m == 0 and step // m - 1 <= MAX_MINORS:
            return m
    return step


def decade_ticks(
    lo: float, hi: float, max_labels: int = 4
) -> tuple[list[int], list[int]]:
    """Labelled and minor decade exponents for a log range ``[lo, hi]``.

    Labelled decades are every ``step`` decades, the smallest step in
    :data:`_DECADE_STEPS` giving at most *max_labels* (the steps above 10
    keep very wide axes readable).  Minor ticks fall on decades
    every :func:`minor_every` (step) decades, so a wide axis is not lined
    with one tick per decade.
    """
    k0 = int(np.ceil(np.log10(lo) - 1e-9))
    k1 = int(np.floor(np.log10(hi) + 1e-9))
    every = list(range(k0, k1 + 1))
    major, step = every, 1
    for step in _DECADE_STEPS:
        major = [k for k in every if k % step == 0]
        if len(major) <= max_labels:
            break
    m = minor_every(step)
    return major, [k for k in every if k % m == 0]


def final_exp_ticks(axis: Axis) -> None:
    """Exponent ticks on a log *axis* from its final view limits.

    Narrow axes get fractional exponents (:func:`refine_exp_ticks`); wider
    ones the decade rule of :func:`decade_ticks`.  Non-log axes are left alone.
    """
    if axis.get_scale() != "log":
        return
    lo, hi = sorted(axis.get_view_interval())
    if not (np.isfinite(lo) and np.isfinite(hi)) or lo <= 0:
        return
    if np.log10(hi) - np.log10(lo) < NARROW_DECADES:
        refine_exp_ticks(axis)
        if isinstance(axis.get_major_locator(), FixedLocator):
            return
    major, minor = decade_ticks(lo, hi)
    axis.set_major_locator(FixedLocator([10.0**k for k in major]))
    axis.set_major_formatter(FuncFormatter(exp_label))
    axis.set_minor_locator(FixedLocator([10.0**k for k in minor]))
    axis.set_minor_formatter(NullFormatter())


class EdgeTrim(Formatter):
    """Wrap a y formatter: blank labels that would stick out of the panel.

    ``bottom`` (axes fraction) extends the blanked zone at the bottom, over
    a strip holding inside x tick labels.  A label level with an inside
    panel title is blanked too when it would come within
    :data:`TITLE_CLEAR_PT` of the title (measured on the drawn title).
    """

    def __init__(
        self,
        ax: Axes,
        inner: Formatter,
        fontsize: float,
        tick_pt: float = 0.0,
        scale: float = 1.0,
    ) -> None:
        self.ax, self.inner, self.fontsize, self.tick_pt = ax, inner, fontsize, tick_pt
        self.scale = scale
        self.bottom = 0.0

    def set_locs(self, locs: list[float]) -> None:
        """Forward the tick locations to the wrapped formatter."""
        self.inner.set_locs(locs)

    def __call__(self, value: float, pos: int | None = None) -> str:
        """The wrapped label, or ``""`` near an edge or an inside title."""
        to_axes = self.ax.transData + self.ax.transAxes.inverted()
        frac = to_axes.transform((1.0, value))[1]
        h_pt = self.ax.get_window_extent().height * 72 / self.ax.figure.dpi
        half = (0.5 * self.fontsize + 1.0) / max(h_pt, 1e-9)
        if frac < self.bottom + half or frac > 1 - half:
            return ""
        text = self.inner(value, pos)
        return "" if text and self._hits_title(text, value) else text

    def _hits_title(self, text: str, value: float) -> bool:
        """Whether *text* at *value* would crowd an inside title."""
        title = self.ax.title
        if not title.get_text():
            return False
        r = self.ax.figure._get_renderer()  # any backend, also while saving
        box, panel = title.get_window_extent(r), self.ax.get_window_extent(r)
        if box.y1 > panel.y1 + 1:  # title above the panel
            return False
        px = self.ax.figure.dpi / 72
        y = self.ax.transData.transform((1.0, value))[1]
        h = 0.5 * self.fontsize * px
        if y + h < box.y0 or y - h > box.y1:
            return False
        width, _, _ = r.get_text_width_height_descent(
            text, FontProperties(size=self.fontsize), ismath="$" in text
        )
        right = panel.x0 + (self.tick_pt + INSIDE_GAP_PT * self.scale) * px + width
        return right + TITLE_CLEAR_PT * px > box.x0  # titles don't scale


class XEdgeTrim(Formatter):
    """Wrap an x formatter: blank labels that would cross a vertical edge.

    *sides* holds ``"left"`` and/or ``"right"``: the right spine where the
    next panel touches this one, both for inside x labels.  ``left`` (axes
    fraction) extends the blanked zone on the left, over a strip holding
    inside y tick labels.  The label width is estimated from its character
    count (mathtext markup stripped).
    """

    def __init__(
        self, ax: Axes, inner: Formatter, fontsize: float, sides: tuple[str, ...]
    ) -> None:
        self.ax, self.inner, self.fontsize, self.sides = ax, inner, fontsize, sides
        self.left = 0.0

    def set_locs(self, locs: list[float]) -> None:
        """Forward the tick locations to the wrapped formatter."""
        self.inner.set_locs(locs)

    def __call__(self, value: float, pos: int | None = None) -> str:
        """The wrapped label, or ``""`` when it would reach a trimmed edge."""
        text = self.inner(value, pos)
        chars = len(re.sub(r"\\mathdefault|[{}$^\\]", "", text))
        to_axes = self.ax.transData + self.ax.transAxes.inverted()
        frac = to_axes.transform((value, 1.0))[0]
        w_pt = self.ax.get_window_extent().width * 72 / self.ax.figure.dpi
        half = (0.3 * self.fontsize * chars + 1.0) / max(w_pt, 1e-9)
        if "right" in self.sides and frac > 1 - half:
            return ""
        if "left" in self.sides and frac < self.left + half:
            return ""
        return text


#: Points between the left spine and an inside label, and after the label.
INSIDE_GAP_PT: float = 3.0
#: Least clearance (points) between an inside y label and an inside title.
TITLE_CLEAR_PT: float = 3.0


def inside_x_labels(ax: Axes, tick_pt: float, scale: float = 1.0) -> None:
    """Draw *ax*'s bottom x tick labels inside, just above the inward ticks."""
    ax.tick_params(axis="x", which="major", pad=-(tick_pt + INSIDE_GAP_PT * scale))
    for label in ax.get_xticklabels():
        label.set_verticalalignment("bottom")


def x_label_strip(
    ax: Axes, renderer: RendererBase, tick_pt: float, scale: float = 1.0
) -> float:
    """Tallest inside x label plus tick and gaps, as a fraction of the axes height."""
    height = ax.get_window_extent(renderer).height
    if height <= 0:
        return 0.0
    gaps = (tick_pt + 2 * INSIDE_GAP_PT * scale) * ax.figure.dpi / 72
    heights = [
        t.get_window_extent(renderer).height for t in ax.get_xticklabels() if t.get_text()
    ]
    return (max(heights) + gaps) / height if heights else 0.0


def inside_y_labels(
    ax: Axes, fontsize: float, tick_pt: float, scale: float = 1.0
) -> None:
    """Draw *ax*'s y tick labels inside, just right of the inward ticks.

    *tick_pt* is the major tick length: the label starts that far plus
    :data:`INSIDE_GAP_PT` (times *scale*, the style's column-count scale
    factor, see :func:`~jaff.plotting.themes.scaled`) from the spine, so the
    tick never touches it (a tick running into a "2" reads as "−2").
    """
    ax.tick_params(axis="y", which="major", pad=-(tick_pt + INSIDE_GAP_PT * scale))
    for label in ax.get_yticklabels():
        label.set_horizontalalignment("left")
    inner = ax.yaxis.get_major_formatter()
    if not isinstance(inner, EdgeTrim):
        ax.yaxis.set_major_formatter(EdgeTrim(ax, inner, fontsize, tick_pt, scale))


def label_strip(
    ax: Axes, renderer: RendererBase, tick_pt: float, scale: float = 1.0
) -> float:
    """Widest inside y label plus tick and gaps, as a fraction of the axes width."""
    width = ax.get_window_extent(renderer).width
    if width <= 0:
        return 0.0
    gaps = (tick_pt + 2 * INSIDE_GAP_PT * scale) * ax.figure.dpi / 72
    widths = [
        t.get_window_extent(renderer).width
        for t in ax.get_yticklabels()
        if t.get_text() and t.get_visible()
    ]
    return (max(widths) + gaps) / width if widths else 0.0


def clear_strip(ax: Axes, frac: float, axis: str = "x") -> None:
    """Widen *ax*'s range at the low end so its first *frac* is empty.

    ``axis="x"`` clears a strip on the left (inside y labels), ``"y"`` one at
    the bottom (inside x labels); log axes are widened in log space.
    """
    if frac <= 0 or frac >= 0.9:
        return
    get, put, scale = (
        (ax.get_xlim, ax.set_xlim, ax.get_xscale())
        if axis == "x"
        else (ax.get_ylim, ax.set_ylim, ax.get_yscale())
    )
    lo, hi = get()
    log = scale == "log"
    a, b = (np.log10(lo), np.log10(hi)) if log else (lo, hi)
    c = (a - frac * b) / (1 - frac)
    put(10.0**c if log else c, hi)
