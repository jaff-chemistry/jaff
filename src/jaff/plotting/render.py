# ABOUTME: Pure-matplotlib renderer: FigureSpec -> (Figure, 2-D ndarray of Axes)
# ABOUTME: Applies the theme, cleans data, draws lines/markers/bars, legends

"""Render a validated :class:`~jaff.plotting.spec.FigureSpec`.

:func:`render` has no side effects beyond creating the figure: it never saves
or shows, and it applies the style's rcParams only inside
:func:`matplotlib.pyplot.rc_context`, so global matplotlib state is
untouched.

Order of work: draw every panel's data (band bars and hollow markers are
only queued); place outer / sup labels and each panel's x tick-label role;
put titles inside (or above when too wide) and leave headroom above the
data for inside text; settle ticks (final exponent ticks, inside x labels on
middle panels, inside y labels); legends; the manual tight layout
(:mod:`jaff.plotting._layout`); widen y ranges downward to clear inside x
labels (re-ticking y); widen x ranges to clear the inside y label strip
(final x ticks); a second layout pass when the first widening moved outside
y labels; trim x labels that would cross an edge; draw the queued band bars;
finally hollow markers, spaced from the final layout.  Markers and the fixed
exponent ticks are therefore chosen from the limits at render time: changing
limits afterwards does not re-space or re-tick them.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.collections import LineCollection
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from ..io import JaffLogger
from . import _layout, _legend, _marks, _ticks
from .spec import AxesSpec, FigureSpec, SeriesSpec
from .themes import (
    REFERENCE_LW_FACTOR,
    REFERENCE_ZORDER,
    Palette,
    Style,
    assign,
    label_looks,
    resolve_palette,
    resolve_style,
    scaled,
    split_style,
)

#: Headroom above the data for inside panel text, as a fraction of the data
#: span (log axes: of at least two decades).
HEADROOM: float = 0.25
#: Autoscale margin of a log y axis (fraction of its span).
LOG_MARGIN: float = 0.08
#: Most decades a log y axis spans without an explicit ``ylim``: data further
#: below its top (e.g. an ``exp(-E/T)`` rate underflowing at low T) runs off
#: the bottom instead of stretching the axis.
LOG_MAX_DECADES: float = 30.0
#: Points between an inside title and the tag under it.
TAG_GAP_PT: float = 3.0
#: A title wider than this fraction of the panel goes above it instead.
TITLE_FIT: float = 0.85


def render(spec: FigureSpec) -> tuple[Figure, np.ndarray]:
    """Draw *spec* and return ``(fig, axes)``.

    Parameters
    ----------
    spec : FigureSpec
        Validated figure specification.

    Returns
    -------
    tuple[matplotlib.figure.Figure, numpy.ndarray]
        The figure and a 2-D array of axes shaped like ``spec.grid``; cells
        beyond ``len(spec.axes)`` are hidden.

    Raises
    ------
    ValueError
        If the style or palette name is unknown.

    Notes
    -----
    Unlike :func:`~jaff.plotting.plot`, the figure stays open in pyplot;
    close it with :func:`matplotlib.pyplot.close` when done.

    Examples
    --------
    >>> from jaff.plotting import FigureSpec, render
    >>> spec = FigureSpec.from_dict(
    ...     {"grid": (1, 2), "axes": [{"series": [{"x": [1, 2], "y": [3, 4]}]}]}
    ... )
    >>> fig, axes = render(spec)
    >>> axes.shape, axes[0, 1].axison
    ((1, 2), False)
    """
    style = effective_style(spec)
    palette = resolve_palette(spec.palette)
    m, n = spec.grid
    pw, ph = style.panel_size
    with plt.rc_context(style.rc):
        fig, axes = plt.subplots(
            m,
            n,
            squeeze=False,
            figsize=spec.figsize or (n * pw, m * ph),
            sharex=spec.sharex,
            sharey=spec.sharey,
        )
        fig.patch.set_facecolor(style.surface)
        grid = _trim_grid(axes, len(spec.axes))
        flat = list(grid.ravel())
        used, free = flat[: len(spec.axes)], flat[len(spec.axes) :]
        for ax in free:
            ax.set_axis_off()

        modes = _legend.fallback(spec.axes, _legend.modes(spec.axes, spec.legend), style)
        shared = _shared_looks(spec, modes, style, palette)
        solo = _solo_color(style, palette, shared)
        pending: list[_marks.PendingMarks] = []
        bars: list[_PendingBars] = []
        entries: list[_legend.Entries] = []
        for j, (ax, aspec) in enumerate(zip(used, spec.axes)):
            looks = shared if modes[j] == _legend.SHARED else None
            marks, ent = _draw_axes(ax, aspec, style, palette, j, bars, looks, solo)
            pending += marks
            entries.append(ent)
        if spec.title:
            fig.suptitle(spec.title)
        for group in _y_groups(len(used), spec.sharey):
            _cap_log_range([used[j] for j in group], [spec.axes[j] for j in group])

        groups = _x_groups(used, spec.sharex)
        roles = _x_roles(len(used), n, groups)
        _outer_labels(fig, grid, spec, roles)
        width = pw if spec.figsize is None else 0.8 * spec.figsize[0] / n
        texts = [
            _panel_texts(ax, aspec, j, modes[j], entries[j], style, width)
            for j, (ax, aspec) in enumerate(zip(used, spec.axes))
        ]
        for group in _y_groups(len(used), spec.sharey):
            if any(texts[j] for j in group):
                _headroom([used[j] for j in group], [spec.axes[j] for j in group])
        inside_x = [r == INSIDE and a.ylim is None for r, a in zip(roles, spec.axes)]
        _settle_ticks(used, spec, inside_x, style)
        placed = _legend.apply(
            fig, used, free, spec.axes, entries, modes, style, spec.legend_title
        )
        _layout.apply(fig, grid, len(used), placed, style, spec.figsize)
        bottoms = _bottom_strips(fig, used, spec, inside_x, style)
        lefts = _clear_label_strips(fig, used, spec, groups, style)
        if any(bottoms) and not style.inside_ticks:
            # outside y labels may have changed with the widened y ranges
            _layout.apply(fig, grid, len(used), placed, style, spec.figsize)
        _trim_x_labels(used, spec, inside_x, lefts, style)
        _draw_bars(bars)
        _marks.place(fig, pending, style)

    return fig, axes


def _block(grid: tuple[int, int], n_used: int) -> tuple[int, int]:
    """Rows and columns of the grid that hold used panels.

    Panels fill row-major, so only trailing rows (and, with a single row,
    trailing columns) can be empty.
    """
    m, n = grid
    rows = max(1, -(-n_used // n))
    return rows, n if rows > 1 else max(1, min(n, n_used))


def effective_style(spec: FigureSpec) -> Style:
    """The spec's style, scaled for its column count.

    With ``k = max(1, used columns / style.ref_cols)`` every visual size is
    scaled by ``k`` (:func:`jaff.plotting.themes.scaled`); an explicit
    ``figsize`` keeps ``k = 1``.  Used by rendering and saving alike.
    """
    style = resolve_style(spec.style)
    if spec.figsize is not None:
        return style
    _, cols = _block(spec.grid, len(spec.axes))
    return scaled(style, max(1.0, cols / style.ref_cols))


def _trim_grid(axes: np.ndarray, n_used: int) -> np.ndarray:
    """The block of the grid holding used panels (:func:`_block`); the rest is hidden.

    Hidden cells inside the block stay (they can host the shared legend).
    """
    n = axes.shape[1]
    rows, cols = _block(axes.shape, n_used)
    for i, ax in enumerate(axes.ravel()):
        if i // n >= rows or i % n >= cols:
            ax.set_axis_off()
            ax.set_visible(False)
    return axes[:rows, :cols]


def _solo_color(style: Style, palette: Palette, shared: dict[str, tuple]) -> str:
    """Colour for a lone unlabelled series, never one the shared legend uses."""
    taken = {mcolors.to_hex(color) for color, _ in shared.values()}
    candidates = [style.solo_color, *palette.cycle]
    return next(
        (c for c in candidates if mcolors.to_hex(c) not in taken), style.solo_color
    )


#: ``(axes, series, colour, linewidth)`` for a band-bar series awaiting final limits.
_PendingBars = tuple[Axes, SeriesSpec, str, float]


def _style_axes(ax: Axes, style: Style) -> None:
    """Surface, ink spines on all sides, inward major+minor ticks on all sides."""
    ax.set_facecolor(style.surface)
    for spine in ax.spines.values():
        spine.set_color(style.ink)
    ax.minorticks_on()
    t = style.ticks
    ax.tick_params(
        which="both",
        direction="in",
        colors=style.ink,
        top=True,
        bottom=True,
        left=True,
        right=True,
    )
    ax.tick_params(which="major", length=t["major_length"], width=t["major_width"])
    ax.tick_params(which="minor", length=t["minor_length"], width=t["minor_width"])
    ax.tick_params(axis="x", labelsize=style.rc.get("xtick.labelsize"))
    ax.tick_params(axis="y", labelsize=style.rc.get("ytick.labelsize"))


def plottable_mask(
    kind: str,
    x: np.ndarray,
    y: np.ndarray,
    xscale: str,
    yscale: str,
    x_lo: np.ndarray | None = None,
    x_hi: np.ndarray | None = None,
) -> np.ndarray:
    """Boolean mask of the points that can be drawn on the given scales.

    Points must be finite (bars: non-NaN edges, an open ``inf`` edge being
    allowed) and positive on a log axis.
    """
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(y)
    if yscale == "log":
        keep &= y > 0
    if kind == "bars":
        x_lo, x_hi = np.asarray(x_lo, dtype=float), np.asarray(x_hi, dtype=float)
        keep &= ~np.isnan(x_lo) & ~np.isnan(x_hi)
        if xscale == "log":
            keep &= ((x_lo > 0) | np.isneginf(x_lo)) & ((x_hi > 0) | np.isposinf(x_hi))
        return keep
    x = np.asarray(x, dtype=float)
    keep &= np.isfinite(x)
    if xscale == "log":
        keep &= x > 0
    return keep


def min_points(kind: str) -> int:
    """Fewest plottable points a series needs: one bar, two line points."""
    return 1 if kind == "bars" else 2


def _clean(s: SeriesSpec, aspec: AxesSpec, path: str) -> SeriesSpec | None:
    """Drop non-finite points, and non-positive ones on log axes.

    Returns ``None`` (with a warning) when too few points remain: fewer than
    two for a line, none for bars.
    """
    keep = plottable_mask(s.kind, s.x, s.y, aspec.xscale, aspec.yscale, s.x_lo, s.x_hi)
    if keep.sum() < min_points(s.kind):
        JaffLogger().get_logger().warning(
            f"{path}: {s.label or 'series'} has {keep.sum()} plottable point(s); skipped."
        )
        return None
    changes: dict[str, Any] = {}
    if s.kind == "bars":
        changes = {"x_lo": s.x_lo[keep], "x_hi": s.x_hi[keep]}
    return replace(s, x=s.x[keep], y=s.y[keep], **changes)


def _handle(kind: str, color: str, look: Any, lw: float, style: Style) -> Line2D:
    """Legend proxy matching how a series is drawn."""
    if kind == "bars":
        return Line2D([], [], color=color, lw=2 * lw)
    ls, marker = split_style(look)
    return Line2D(
        [],
        [],
        color=color,
        linestyle=ls if ls is not None else "",
        lw=lw,
        marker=marker,
        markersize=style.marker["size"] * 1.4,
        markerfacecolor=style.surface,
        markeredgecolor=color,
        markeredgewidth=lw,
    )


def _shared_looks(
    spec: FigureSpec, modes: list[str], style: Style, palette: Palette
) -> dict[str, tuple]:
    """Figure-wide ``(colour, style)`` per label of the shared-legend panels."""
    labels = [
        s.label
        for aspec, mode in zip(spec.axes, modes)
        if mode == _legend.SHARED
        for s in aspec.series
        if s.label and s.role != "reference"
    ]
    return label_looks(labels, style, palette)


def _width_order(s: SeriesSpec, style: Style) -> tuple[float, float | None]:
    """Line width and z-order: explicit, else the role's default."""
    ref = s.role == "reference"
    lw = s.lw if s.lw is not None else style.lw * (REFERENCE_LW_FACTOR if ref else 1.0)
    zorder = s.zorder if s.zorder is not None else (REFERENCE_ZORDER if ref else None)
    return lw, zorder


def _draw_series(
    ax: Axes,
    s: SeriesSpec,
    color: str,
    look: Any,
    lw: float,
    zorder: float | None,
    line_gap_pt: float,
) -> _marks.PendingMarks | None:
    """Draw a line (and/or queue its hollow markers); return queued markers."""
    ls, marker = split_style(look)
    if ls is None:
        # Invisible line: claims the data limits; markers come after layout.
        ax.plot(s.x, s.y, visible=False)
        return _marks.PendingMarks(ax, s.x, s.y, color, marker, lw, zorder)
    ax.plot(s.x, s.y, color=color, linestyle=ls, lw=lw, zorder=zorder)
    if marker is None:
        return None
    return _marks.PendingMarks(ax, s.x, s.y, color, marker, lw, zorder, line_gap_pt)


def _draw_axes(
    ax: Axes,
    aspec: AxesSpec,
    style: Style,
    palette: Palette,
    index: int,
    bars: list[_PendingBars],
    looks: dict[str, tuple] | None = None,
    solo: str | None = None,
) -> tuple[list[_marks.PendingMarks], _legend.Entries]:
    """Draw one panel's data and labels; return queued markers and legend entries.

    *looks* fixes ``(colour, style)`` by label (panels keyed by the shared
    legend); *solo* colours an unlabelled series alone on the panel.
    """
    ax.set_xscale(aspec.xscale)
    ax.set_yscale(aspec.yscale)
    if aspec.yscale == "log":
        ax.margins(y=LOG_MARGIN)  # keeps edge decade labels clear of the frame
    _style_axes(ax, style)

    cleaned = [
        c
        for i, s in enumerate(aspec.series)
        if (c := _clean(s, aspec, f"axes[{index}].series[{i}]")) is not None
    ]
    if solo and len(cleaned) == 1:
        lone = cleaned[0]
        if not (
            lone.label or lone.color or lone.role == "reference" or lone.kind == "bars"
        ):
            cleaned[0] = replace(lone, color=solo)
    pending: list[_marks.PendingMarks] = []
    entries: _legend.Entries = []
    for s, (color, look) in zip(cleaned, assign(cleaned, style, palette, looks)):
        lw, zorder = _width_order(s, style)
        if s.kind == "bars":
            _queue_bars(ax, s, color, lw, bars, f"axes[{index}]")
        elif (
            marks := _draw_series(
                ax, s, color, look, lw, zorder, style.marker["line_gap_pt"]
            )
        ) is not None:
            pending.append(marks)
        if s.label:
            entries.append((s.label, _handle(s.kind, color, look, lw, style)))
    if aspec.xlim is not None:
        ax.set_xlim(aspec.xlim)
    if aspec.ylim is not None:
        ax.set_autoscaley_on(False)  # skip the pending autoscale (may underflow)
        ax.set_ylim(aspec.ylim)
    ax.set_xlabel(_axis_label(aspec.xscale, aspec.xlabel, style))
    ax.set_ylabel(_axis_label(aspec.yscale, aspec.ylabel, style))
    return pending, entries


def _queue_bars(
    ax: Axes, s: SeriesSpec, color: str, lw: float, bars: list[_PendingBars], path: str
) -> None:
    """Feed a bar series' finite edges to the data limits and queue it.

    Open (infinite) edges never touch the limits; they are resolved against the
    final limits by :func:`_draw_bars`.  A series with no finite edge is skipped.
    """
    edges = np.concatenate([s.x_lo, s.x_hi])
    ys = np.concatenate([s.y, s.y])
    ok = np.isfinite(edges)
    if not ok.any():
        JaffLogger().get_logger().warning(
            f"{path}: bars have no finite edges to size the axis; skipped."
        )
        return
    ax.update_datalim(np.column_stack([edges[ok], ys[ok]]))
    ax.autoscale_view()
    bars.append((ax, s, color, lw))


def _draw_bars(bars: list[_PendingBars]) -> None:
    """Draw queued band bars once every panel's limits are final.

    Infinite edges clip to the axes' x-limits; the collections are added
    without autoscaling so they never feed back into the limits.
    """
    for ax, s, color, lw in bars:
        lo, hi = ax.get_xlim()
        x_lo = np.where(np.isneginf(s.x_lo), lo, s.x_lo)
        x_hi = np.where(np.isposinf(s.x_hi), hi, s.x_hi)
        segs = [[(a, y), (b, y)] for a, b, y in zip(x_lo, x_hi, s.y)]
        lc = LineCollection(segs, colors=color, linewidths=2 * lw, zorder=s.zorder)
        ax.add_collection(lc, autolim=False)


def _x_groups(used: list[Axes], sharex: bool) -> list[list[int]]:
    """Panels sharing one x range: all (``sharex``), else equal scale + limits."""
    if sharex:
        return [list(range(len(used)))] if used else []
    groups: dict[tuple, list[int]] = {}
    for i, ax in enumerate(used):
        lo, hi = ax.get_xlim()
        key = (
            ax.get_xscale(),
            float(np.format_float_positional(lo, 9, fractional=False)),
            float(np.format_float_positional(hi, 9, fractional=False)),
        )
        groups.setdefault(key, []).append(i)
    return list(groups.values())


def _common(labels: list[str]) -> str:
    """The label every panel shows, or ``""`` when they differ or are empty."""
    return labels[0] if labels and all(lb == labels[0] for lb in labels) else ""


#: X-axis role of a panel within its column (see :func:`_x_roles`).
TOP, BOTTOM, SAME, INSIDE = "top", "bottom", "same", "inside"


def _x_roles(n_used: int, ncols: int, groups: list[list[int]]) -> list[str]:
    """Where each panel's x tick labels go, by its place in its column.

    Topmost used panel: on the top spine.  Bottommost (also a lone panel):
    at the bottom.  Middle: none when it shares the bottommost panel's x
    range, else inside above the bottom spine.
    """
    group_of = {i: k for k, g in enumerate(groups) for i in g}
    roles = []
    for i in range(n_used):
        column = list(range(i % ncols, n_used, ncols))
        if i == column[-1]:
            roles.append(BOTTOM)
        elif i == column[0]:
            roles.append(TOP)
        else:
            roles.append(SAME if group_of[i] == group_of[column[-1]] else INSIDE)
    return roles


def _outer_labels(
    fig: Figure, axes: np.ndarray, spec: FigureSpec, roles: list[str]
) -> None:
    """Place x tick labels by column role; draw shared axis labels once.

    Compares the labels as displayed.  When every used panel shows the same
    y label and the used panels span at least two rows, it becomes one
    ``fig.supylabel``; likewise a common x label over at least two columns
    becomes ``fig.supxlabel``.  Otherwise the x label follows the tick labels
    (top spine for the topmost panel, bottom for the bottommost) and a middle
    panel repeating its column's bottom label drops it; a y label is dropped
    when the panel to the left shows the same one.
    """
    n_used, ncols = len(spec.axes), spec.grid[1]
    flat = axes.ravel()
    xl = [a.get_xlabel() for a in flat[:n_used]]
    yl = [a.get_ylabel() for a in flat[:n_used]]
    size = plt.rcParams["axes.labelsize"]
    sup_x = _common(xl) if min(n_used, ncols) >= 2 else ""
    sup_y = _common(yl) if -(-n_used // ncols) >= 2 else ""
    if sup_x:
        fig.supxlabel(sup_x, fontsize=size)
    if sup_y:
        fig.supylabel(sup_y, fontsize=size)
    for i, role in enumerate(roles):
        ax, bottom = flat[i], list(range(i % ncols, n_used, ncols))[-1]
        ax.tick_params(
            axis="x",
            which="both",
            labeltop=role == TOP,
            labelbottom=role in (BOTTOM, INSIDE),
        )
        if role == TOP:
            ax.xaxis.set_label_position("top")
        if sup_x or (role in (SAME, INSIDE) and xl[i] == xl[bottom]):
            ax.set_xlabel("")
        if sup_y or (i % ncols and yl[i - 1] == yl[i]):
            ax.set_ylabel("")


def _panel_texts(
    ax: Axes,
    aspec: AxesSpec,
    index: int,
    mode: str,
    entries: _legend.Entries,
    style: Style,
    panel_w: float,
) -> bool:
    """Title inside (or above when too wide) plus tag / letter inside.

    Returns whether any text sits inside the panel (it then needs headroom).
    """
    inside_title, title_h = False, 0.0
    if aspec.title:
        w, h = _text_size(ax, aspec.title, style.tag["fontsize"])
        if w <= TITLE_FIT * panel_w:
            inside_title, title_h = True, h
            ax.set_title(
                aspec.title,
                x=0.95,
                y=0.94,
                pad=0,
                ha="right",
                va="top",
                fontsize=style.tag["fontsize"],
                bbox=_legend.box(style),
            )
        else:
            # Auto-positioned by matplotlib: above top-spine tick labels / x label.
            ax.set_title(aspec.title)
    tag = aspec.tag
    if not tag and mode == _legend.GROUP and entries and not aspec.title:
        tag = _legend.letter(index)
    if tag:
        drop_pt = title_h * 72 + TAG_GAP_PT if inside_title else 0.0  # unscaled
        _legend.add_tag(ax, tag, style, drop_pt=drop_pt)
    return inside_title or bool(tag)


def _text_size(ax: Axes, text: str, fontsize: float) -> tuple[float, float]:
    """``(width, height)`` in inches of *text* at *fontsize*, without keeping it.

    Measured on a throw-away artist so the panel title is set only once:
    ``set_title`` with an explicit ``y`` permanently turns off matplotlib's
    automatic title placement above top-spine labels.
    """
    probe = ax.text(0, 0, text, fontsize=fontsize, transform=ax.transAxes)
    box = probe.get_window_extent(ax.figure.canvas.get_renderer())
    probe.remove()
    return box.width / ax.figure.dpi, box.height / ax.figure.dpi


def _y_groups(n_used: int, sharey: bool) -> list[list[int]]:
    """Panels sharing one y range: all with ``sharey``, else each alone."""
    if sharey:
        return [list(range(n_used))] if n_used else []
    return [[j] for j in range(n_used)]


def _span(ax: Axes, log: bool) -> float | None:
    """Data span of *ax* in y; ``None`` if unusable.

    On a log axis the span is in decades, at most :data:`LOG_MAX_DECADES`
    (the range :func:`_cap_log_range` shows); taking the difference of the
    logs keeps it finite when ``y1 / y0`` would overflow (subnormal ``y0``).
    """
    y0, y1 = ax.dataLim.y0, ax.dataLim.y1
    if not (np.isfinite(y0) and np.isfinite(y1)):
        return None
    if log:
        if not (y0 > 0 and y1 > 0):
            return None
        return min(float(np.log10(y1) - np.log10(y0)), LOG_MAX_DECADES)
    return (y1 - y0) or abs(y1) or 1.0


def _cap_log_range(group: list[Axes], specs: list[AxesSpec]) -> None:
    """Limit a log y range to :data:`LOG_MAX_DECADES` below the data's top.

    Applied once per group of panels sharing a y range, before anything reads
    the limits.  Without it, data reaching subnormal values (``1e-320``) makes
    the autoscale margin underflow to 0, which matplotlib rejects, leaving
    limits that hide every curve.  Panels with an explicit ``ylim`` are left
    alone, as are ranges within the cap.
    """
    if any(s.ylim is not None for s in specs) or group[0].get_yscale() != "log":
        return
    lows = [a.dataLim.minposy for a in group]
    highs = [a.dataLim.y1 for a in group]
    lo = min((v for v in lows if np.isfinite(v) and v > 0), default=None)
    hi = max((v for v in highs if np.isfinite(v) and v > 0), default=None)
    if lo is None or hi is None:
        return
    top = float(np.log10(hi))
    if top - float(np.log10(lo)) <= LOG_MAX_DECADES:
        return
    ax = group[0]
    # Off first: set_ylim would otherwise run the pending autoscale, which is
    # what underflows.
    ax.set_autoscaley_on(False)
    ax.set_ylim(
        10.0 ** (top - LOG_MAX_DECADES), 10.0 ** (top + LOG_MARGIN * LOG_MAX_DECADES)
    )


def _headroom(group: list[Axes], specs: list[AxesSpec]) -> None:
    """Raise a y range above the data so inside text stays clear of curves.

    Applied once per group of panels sharing a y range, with the group's
    largest data span (re-applying per panel would compound under sharey).
    """
    if any(s.ylim is not None for s in specs):
        return
    ax = group[0]
    log = ax.get_yscale() == "log"
    spans = [s for a in group if (s := _span(a, log)) is not None]
    if not spans:
        return
    lo, hi = ax.get_ylim()
    if log:
        ax.set_ylim(lo, hi * 10 ** (HEADROOM * max(max(spans), 2.0)))
    else:
        ax.set_ylim(lo, hi + HEADROOM * max(spans))


def _settle_ticks(
    used: list[Axes], spec: FigureSpec, inside_x: list[bool], style: Style
) -> None:
    """Final exponent ticks, then inside y (and middle-panel x) tick labels.

    A panel with an explicit ``xlim`` keeps outside y labels: its x range
    cannot be widened to clear a strip for inside ones (``inside_x`` already
    excludes panels with an explicit ``ylim``).  X ticks are redone later.
    """
    if style.exp_ticks:
        for ax in used:
            _ticks.final_exp_ticks(ax.yaxis)
            _ticks.final_exp_ticks(ax.xaxis)
    for ax, inside in zip(used, inside_x):
        if inside:
            _ticks.inside_x_labels(ax, style.ticks["major_length"], style.scale)
    if style.inside_ticks:
        size = plt.rcParams["ytick.labelsize"]
        size = size if isinstance(size, (int, float)) else plt.rcParams["font.size"]
        tick = style.ticks["major_length"]
        for ax, aspec in zip(used, spec.axes):
            if aspec.xlim is None:
                _ticks.inside_y_labels(ax, float(size), tick, style.scale)


def _clear_label_strips(
    fig: Figure,
    used: list[Axes],
    spec: FigureSpec,
    groups: list[list[int]],
    style: Style,
) -> dict[int, float]:
    """Widen x ranges so no curve runs under an inside y label.

    Runs on the final layout (the strip is a fraction of the final panel
    width); panels sharing an x range get the widest strip of their group.
    One redraw first, so the labels read are those of the final size (edge
    trimming depends on the panel height).  Returns the cleared fraction per
    panel index (for x label trimming).
    """
    if not style.inside_ticks:
        return {}
    fig.draw_without_rendering()
    r = fig.canvas.get_renderer()
    tick = style.ticks["major_length"]
    cleared: dict[int, float] = {}
    for group in groups:
        inside = [i for i in group if spec.axes[i].xlim is None]
        frac = max(
            (_ticks.label_strip(used[i], r, tick, style.scale) for i in inside),
            default=0.0,
        )
        for i in inside:
            cleared[i] = frac
        targets = [used[i] for i in inside]
        for ax in targets[:1] if spec.sharex else targets:
            _ticks.clear_strip(ax, frac)
    if style.exp_ticks:
        for ax in used:
            _ticks.final_exp_ticks(ax.xaxis)
    return cleared


def _bottom_strips(
    fig: Figure,
    used: list[Axes],
    spec: FigureSpec,
    inside_x: list[bool],
    style: Style,
) -> list[float]:
    """Widen y ranges downward so no curve runs under inside x labels.

    Mirrors the inside-y strip: one group per shared y range, the tallest
    strip of the group; y ticks are re-chosen afterwards and inside y labels
    keep clear of the strip.  Returns the cleared fraction per panel.
    """
    out = [0.0] * len(used)
    if not any(inside_x):
        return out
    fig.draw_without_rendering()
    r = fig.canvas.get_renderer()
    tick = style.ticks["major_length"]
    for group in _y_groups(len(used), spec.sharey):
        members = [j for j in group if inside_x[j]]
        frac = max(
            (_ticks.x_label_strip(used[j], r, tick, style.scale) for j in members),
            default=0,
        )
        if not members or frac <= 0:
            continue
        _ticks.clear_strip(used[group[0]], frac, axis="y")
        for j in group:
            out[j] = frac
            _retick_y(used[j], spec.axes[j], frac, style)
    return out


def _retick_y(ax: Axes, aspec: AxesSpec, bottom: float, style: Style) -> None:
    """Re-choose y ticks after a limit change; inside y labels skip *bottom*."""
    if style.exp_ticks:
        _ticks.final_exp_ticks(ax.yaxis)
    if style.inside_ticks and aspec.xlim is None:
        size = plt.rcParams["ytick.labelsize"]
        size = size if isinstance(size, (int, float)) else plt.rcParams["font.size"]
        _ticks.inside_y_labels(ax, float(size), style.ticks["major_length"], style.scale)
        trim = ax.yaxis.get_major_formatter()
        if isinstance(trim, _ticks.EdgeTrim):
            trim.bottom = bottom


#: Panels closer than this (inches) count as touching for x label trimming.
TOUCH_GAP: float = 0.02


def _trim_x_labels(
    used: list[Axes],
    spec: FigureSpec,
    inside_x: list[bool],
    lefts: dict[int, float],
    style: Style,
) -> None:
    """Blank x labels that would cross an edge.

    Inside x labels are trimmed at both side spines (and kept off the inside
    y label strip); outside ones at a right spine another panel touches
    (top and bottom labels share the formatter).
    """
    ncols = spec.grid[1]
    size = plt.rcParams["xtick.labelsize"]
    size = float(size if isinstance(size, (int, float)) else plt.rcParams["font.size"])
    for i, ax in enumerate(used):
        touching = style.gap < TOUCH_GAP and (i + 1) % ncols and i + 1 < len(used)
        sides = ("left", "right") if inside_x[i] else ("right",) if touching else ()
        if not sides:
            continue
        trim = _ticks.XEdgeTrim(ax, ax.xaxis.get_major_formatter(), size, sides)
        trim.left = lefts.get(i, 0.0) if inside_x[i] else 0.0
        ax.xaxis.set_major_formatter(trim)


def _axis_label(scale: str, label: str, style: Style) -> str:
    """The axis label, prefixed with ``log10`` on exponent-labelled log axes."""
    if scale == "log" and style.exp_ticks:
        return _ticks.log10_label(label)
    return label
