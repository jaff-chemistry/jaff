# ABOUTME: Manual tight grid layout: equal absolute gaps, margins sized to labels
# ABOUTME: Positions axes, sup labels and the legends kept in the right margin

"""Tight uniform grid layout, done by hand in inches.

After one draw, each used panel's decorations outside its box are measured:
y tick labels + y label on the left/right, x tick labels + x label below, an
outside title above.  Columns are then spaced ``left decor + panel + right
decor + gap`` and rows ``top decor + panel + bottom decor + gap``, so the
empty space between neighbouring panels' decorations is the same ``gap``
across rows and columns.  Bottom decorations of a panel with no used panel
below it hang into the hidden cell beneath and do not widen the row gap.

Without ``figsize`` every panel box is ``style.panel_size`` and the figure is
sized to fit; with ``figsize`` the panel size is what is left over.  Margin
legends are stacked top-down right of the grid (wrapping into further
columns when taller than the grid); a legend hosted in a hidden cell is
lowered below the decorations hanging into that cell.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from matplotlib.axes import Axes
from matplotlib.backend_bases import RendererBase
from matplotlib.figure import Figure, SubFigure
from matplotlib.legend import Legend
from matplotlib.text import Text

from ._legend import Placement
from .themes import Style


def renderer(fig: Figure | SubFigure) -> RendererBase:
    """The renderer that measures *fig*'s text, on any backend.

    ``fig.canvas.get_renderer()`` exists only on Agg canvases (``pdf``,
    ``svg``, ``ps`` have none); the figure's own renderer works everywhere.
    ``Figure._get_renderer`` is private (missing from the stubs) but it is what
    matplotlib itself uses, and :mod:`jaff.plotting._ticks` already relies on it.
    """
    return fig._get_renderer()  # ty: ignore[unresolved-attribute]


#: Outer figure padding, inches.
PAD: float = 0.04
#: Vertical space between stacked margin legends, inches.
STACK_GAP: float = 0.12
#: Smallest panel side when a fixed figsize leaves too little room, inches.
MIN_PANEL: float = 0.2


@dataclass
class Decor:
    """Decoration extents outside one panel box, in inches."""

    left: float = 0.0
    right: float = 0.0
    bottom: float = 0.0
    top: float = 0.0


def _decor(ax: Axes, r: RendererBase) -> Decor:
    """Measure how far *ax*'s labels reach outside its box."""
    dpi = ax.figure.dpi
    bb = ax.get_window_extent(r)
    d = Decor()
    ytb = ax.yaxis.get_tightbbox(r)
    if ytb is not None:
        d.left = max(0.0, bb.x0 - ytb.x0) / dpi
        d.right = max(0.0, ytb.x1 - bb.x1) / dpi
    xtb = ax.xaxis.get_tightbbox(r)
    if xtb is not None:
        d.bottom = max(0.0, bb.y0 - xtb.y0) / dpi
        d.top = max(0.0, xtb.y1 - bb.y1) / dpi
    if ax.get_title():
        tb = ax.title.get_window_extent(r)
        d.top = max(d.top, (tb.y1 - bb.y1) / dpi)
        d.left = max(d.left, (bb.x0 - tb.x0) / dpi)  # a wide outside title
        d.right = max(d.right, (tb.x1 - bb.x1) / dpi)  # overhangs its panel
    return d


def _size(artist: Legend | Text, r: RendererBase, dpi: float) -> tuple[float, float]:
    """``(width, height)`` of *artist* in inches."""
    box = artist.get_window_extent(r)
    return box.width / dpi, box.height / dpi


def _stack(sizes: list[tuple[float, float]], height: float) -> list[list[int]]:
    """Split legends, in order, into columns no taller than *height*."""
    cols: list[list[int]] = []
    used = height + 1.0
    for i, (_, h) in enumerate(sizes):
        if cols and used + STACK_GAP + h <= height:
            cols[-1].append(i)
            used += STACK_GAP + h
        else:
            cols.append([i])
            used = h
    return cols


@dataclass
class _Grid:
    """Per-column / per-row decoration extents of the used panels."""

    left: list[float]
    right: list[float]
    top: list[float]
    bottom: list[float]
    hang: dict[int, float]  # flat index of a hidden cell -> decor hanging into it


def _grid(axes: np.ndarray, n_used: int, decors: list[Decor]) -> _Grid:
    """Collect decoration extents per column and row."""
    m, n = axes.shape
    g = _Grid([0.0] * n, [0.0] * n, [0.0] * m, [0.0] * m, {})
    for i, d in enumerate(decors):
        r, c = divmod(i, n)
        g.left[c] = max(g.left[c], d.left)
        g.right[c] = max(g.right[c], d.right)
        g.top[r] = max(g.top[r], d.top)
        below = i + n
        if below < n_used or r == m - 1:
            g.bottom[r] = max(g.bottom[r], d.bottom)
        else:
            g.hang[below] = d.bottom
    return g


@dataclass
class _Frame:
    """Solved layout: figure size, panel size and cell origins, in inches."""

    width: float
    height: float
    pw: float
    ph: float
    col_x: list[float]
    row_top: list[float]
    grid_h: float
    host_drop: float


def _sup_texts(fig: Figure) -> tuple[Text | None, Text | None, Text | None]:
    """The figure's sup x label, sup y label and suptitle (``None`` when unset).

    matplotlib has no public getter for these artists (``get_suptitle()`` gives
    the string); the private attributes are read here only.
    """
    return (
        fig._supxlabel,  # ty: ignore[unresolved-attribute]
        fig._supylabel,  # ty: ignore[unresolved-attribute]
        fig._suptitle,  # ty: ignore[unresolved-attribute]
    )


def _extras(fig: Figure, r: RendererBase) -> tuple[float, float, float]:
    """Space taken by the sup y label (left), sup x label (bottom), suptitle (top)."""
    dpi = fig.dpi
    sup_x, sup_y, title = _sup_texts(fig)
    left = _size(sup_y, r, dpi)[0] + PAD if sup_y is not None else 0.0
    bottom = _size(sup_x, r, dpi)[1] + PAD if sup_x is not None else 0.0
    top = _size(title, r, dpi)[1] + PAD if title is not None else 0.0
    return left, bottom, top


def _origins(
    g: _Grid, pw: float, ph: float, left: float, top: float, gap: float
) -> tuple[list[float], list[float]]:
    """Left edge of each column and top edge of each row of panel boxes (inches)."""
    col_x, x = [], PAD + left
    for c in range(len(g.left)):
        x += g.left[c]
        col_x.append(x)
        x += pw + g.right[c] + gap
    row_top, y = [], PAD + top
    for row in range(len(g.top)):
        y += g.top[row]
        row_top.append(y)
        y += ph + g.bottom[row] + gap
    return col_x, row_top


def _solve(
    g: _Grid,
    shape: tuple[int, int],
    placed: Placement,
    sizes: list[tuple[float, float]],
    host_size: tuple[float, float] | None,
    extras: tuple[float, float, float],
    style: Style,
    figsize: tuple[float, float] | None,
) -> tuple[_Frame, list[list[int]], list[float]]:
    """Panel size, figure size and cell origins; plus the margin legend columns."""
    m, n = shape
    gap, (sup_l, sup_b, sup_t) = style.gap, extras
    y_fixed = 2 * PAD + sup_t + sup_b + sum(g.top) + sum(g.bottom) + gap * (m - 1)
    ph = (
        style.panel_size[1]
        if figsize is None
        else max((figsize[1] - y_fixed) / m, MIN_PANEL)
    )
    grid_h = m * ph + sum(g.top[1:]) + sum(g.bottom[:-1]) + gap * (m - 1)

    columns = _stack(sizes, grid_h) if sizes else []
    col_w = [max(sizes[i][0] for i in col) for col in columns]
    leg_gap = max(gap, 0.1)
    right = leg_gap * len(col_w) + sum(col_w)
    x_fixed = 2 * PAD + sup_l + sum(g.left) + sum(g.right) + gap * (n - 1) + right
    pw = (
        style.panel_size[0]
        if figsize is None
        else max((figsize[0] - x_fixed) / n, MIN_PANEL)
    )

    host_drop, bottom = 0.0, 0.0
    if placed.host is not None and host_size is not None:
        host_drop = g.hang.get(placed.host, 0.0) + (gap if placed.host in g.hang else 0.0)
        right += max(0.0, host_size[0] - pw - g.right[n - 1])
        bottom = max(0.0, host_drop + host_size[1] - ph - g.bottom[m - 1])
    if figsize is None:
        width = (
            2 * PAD + sup_l + sum(g.left) + sum(g.right) + gap * (n - 1) + right + n * pw
        )
        height = y_fixed + m * ph + bottom
        if columns:  # a margin column taller than the grid must still fit
            tallest = max(
                sum(sizes[i][1] for i in col) + STACK_GAP * (len(col) - 1)
                for col in columns
            )
            height = max(height, 2 * PAD + sup_t + g.top[0] + tallest)
    else:
        width, height = figsize

    col_x, row_top = _origins(g, pw, ph, sup_l, sup_t, gap)
    frame = _Frame(width, height, pw, ph, col_x, row_top, grid_h, host_drop)
    return frame, columns, col_w


def _place_texts(fig: Figure, f: _Frame, g: _Grid, n_used: int, n: int) -> None:
    """Sup labels centred on the grid, suptitle on top."""
    W, H = f.width, f.height
    grid_left, grid_right = f.col_x[0], f.col_x[-1] + f.pw + g.right[-1]
    used_right = max((f.col_x[i % n] + f.pw for i in range(n_used)), default=grid_right)
    sup_x, sup_y, title = _sup_texts(fig)
    if sup_y is not None:
        y = 1 - (f.row_top[0] + f.grid_h / 2) / H
        sup_y.set(x=PAD / W, y=y, ha="left", va="center")
    if sup_x is not None:
        x = (grid_left + used_right) / 2 / W
        grid_bottom = f.row_top[-1] + f.ph + g.bottom[-1]
        sup_x.set(x=x, y=1 - grid_bottom / H, ha="center", va="top")
    if title is not None:
        x = (grid_left + grid_right) / 2 / W
        title.set(x=x, y=1 - PAD / H, ha="center", va="top")


def _place_legends(
    fig: Figure,
    flat: np.ndarray,
    placed: Placement,
    f: _Frame,
    g: _Grid,
    columns: list[list[int]],
    col_w: list[float],
    style: Style,
) -> None:
    """Hosted legend below the labels hanging into its cell; margin stack."""
    W, H = f.width, f.height
    if placed.host is not None and placed.shared is not None:
        anchor = (0.0, 1.0 - f.host_drop / f.ph)
        placed.shared.set_bbox_to_anchor(anchor, transform=flat[placed.host].transAxes)
    margin = placed.margin()
    leg_gap = max(style.gap, 0.1)
    x = f.col_x[-1] + f.pw + g.right[-1] + leg_gap
    for col, width in zip(columns, col_w):
        y = f.row_top[0]
        for i in col:
            item = margin[i]
            if isinstance(item, Legend):
                item.set_bbox_to_anchor((x / W, 1 - y / H), transform=fig.transFigure)
            else:  # the legend_title heading over a stack of groups
                item.set(x=x / W, y=1 - y / H, transform=fig.transFigure)
            y += _height(item, fig) + STACK_GAP
        x += width + leg_gap


def _height(artist: Legend | Text, fig: Figure) -> float:
    """Legend (or heading) height in inches."""
    return _size(artist, renderer(fig), fig.dpi)[1]


#: Measure -> solve -> place passes at most; decorations usually settle in two.
MAX_PASSES: int = 3
#: Decorations changing less than this (inches) between passes count as settled.
SETTLED: float = 0.005


def _measure(
    fig: Figure, axes: np.ndarray, n_used: int, placed: Placement
) -> tuple[_Grid, list[tuple[float, float]], tuple[float, float] | None, tuple]:
    """One draw; decoration extents, margin sizes, hosted legend size, extras."""
    fig.draw_without_rendering()
    r = renderer(fig)
    flat = axes.ravel()
    g = _grid(axes, n_used, [_decor(ax, r) for ax in flat[:n_used]])
    sizes = [_size(item, r, fig.dpi) for item in placed.margin()]
    host = placed.shared if placed.host is not None else None
    host_size = _size(host, r, fig.dpi) if host is not None else None
    return g, sizes, host_size, _extras(fig, r)


def _settled(a: _Grid, b: _Grid) -> bool:
    """Whether two measurements agree within :data:`SETTLED` inches."""
    pairs = zip(a.left + a.right + a.top + a.bottom, b.left + b.right + b.top + b.bottom)
    return all(abs(x - y) < SETTLED for x, y in pairs)


def apply(
    fig: Figure,
    axes: np.ndarray,
    n_used: int,
    placed: Placement,
    style: Style,
    figsize: tuple[float, float] | None,
) -> None:
    """Size the figure and position every axes, sup label and margin legend.

    Tick labels depend on the panel size (a taller panel may get more
    ticks), so the layout is repeated on the placed geometry until the
    measured decorations stop changing (at most :data:`MAX_PASSES` passes).
    """
    m, n = axes.shape
    flat = axes.ravel()
    previous: _Grid | None = None
    for _ in range(MAX_PASSES):
        g, sizes, host_size, extras = _measure(fig, axes, n_used, placed)
        if previous is not None and _settled(previous, g):
            break
        f, columns, col_w = _solve(
            g, (m, n), placed, sizes, host_size, extras, style, figsize
        )
        fig.set_size_inches(f.width, f.height)
        for i, ax in enumerate(flat):
            row, c = divmod(i, n)
            top = f.row_top[row] + f.ph
            ax.set_position(
                [
                    f.col_x[c] / f.width,
                    1 - top / f.height,
                    f.pw / f.width,
                    f.ph / f.height,
                ]
            )
        _place_texts(fig, f, g, n_used, n)
        _place_legends(fig, flat, placed, f, g, columns, col_w, style)
        previous = g
