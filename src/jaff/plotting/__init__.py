# ABOUTME: jaff.plotting public API: spec-driven matplotlib figures
# ABOUTME: plot / plot_rates / plot_xsecs builders plus theme registries

"""Spec-driven, publication-style plotting for JAFF.

Pure matplotlib.  A figure is described as **data**: a ``grid`` of axes,
each holding a list of series (see :func:`plot` for the figure-spec format,
and :class:`FigureSpec`, :class:`AxesSpec` and :class:`SeriesSpec` for its
keys).  :func:`plot_rates` and :func:`plot_xsecs` build that spec from
reactions.  Themes supply the look (a *style*, see :data:`STYLES`) and the
colours (a *palette*, see :data:`PALETTES`); every series can override them.

The ``Reaction.plot_rate_coefficient``, ``Reaction.plot_xsecs``,
``Reactions.plot_rates`` and ``Reactions.plot_xsecs`` methods forward their
keyword arguments to these functions.

Notes
-----
**Panels.**  In :func:`plot_rates` / :func:`plot_xsecs` each top-level item
is one panel and a nested list is overlaid on one axes.  Without ``grid``
everything shares a single axes, unless there are more curves than the
style's ``legend_max`` (default 12): then the curves are split in order over
an automatic grid.  Panels fill the grid row by row; leftover cells are
hidden (the bottom-right one hosts the shared legend).  An axis label common
to every panel is drawn once for the whole figure (``supxlabel`` /
``supylabel``) when the panels span two or more columns / rows.

**Closing.**  Figures are closed in pyplot after saving / showing
(``close=True``), so loops don't accumulate open figures; the returned
``(fig, axes)`` still work.  Pass ``close=False`` to keep the figure open.
After an interactive, non-blocking show (``plt.isinteractive()``) the figure
is left open.

**Look.**  White surface, black spines on all four sides, inward major and
minor ticks on all sides, DejaVu Sans with STIX maths.  Panels sit on a tight
grid: the same gap (``style.gap``, inches) between columns and rows, margins
sized to the outer labels.  In each column the top panel carries its x tick
labels on the top spine and the bottom panel at the bottom; a middle panel
shows none when it shares the bottom panel's x range, else small labels
inside above its bottom spine.  Rows and columns without a used panel take no
space.  Each panel box is ``style.panel_size``, unless ``figsize`` fixes the
figure.  In the default ``"paper"`` style:

- y tick labels are drawn **inside** the panel, right of the left spine
  (labels at the top/bottom edge are blanked, and the x range gets an empty
  strip on the left so no curve runs under a label);
- panel **titles** sit inside, upper right, with headroom above the data (a
  title too wide for the panel goes above it);
- log axes are labelled by exponent (``−3`` for 10⁻³) with at most ~4
  labelled decades and minor ticks on decades (at most 4 between labels);
  axes under 1.5 decades get fractional exponents (``1.05``); the axis label
  gains a ``log₁₀`` prefix.

Marker series are hollow markers spaced evenly along the curve on screen; a
line+marker style (``"-o"``, ``"--D"``) adds sparse markers to a line.  A
series with ``"role": "reference"`` is a thick dark-grey line drawn beneath
the rest.

**Legends.**  The default ``legend="auto"`` (same as ``"figure"``) draws
**one shared legend** for the whole figure, de-duplicated by label, in the
bottom-right hidden grid cell when the grid is unfinished, else in the right
margin.  A label keeps one colour and style across every panel, so the shared
legend keys them all.  ``legend="axes"`` (opt-in, for panels with different
labels) draws one legend **group per panel**, stacked in the right margin and
headed by the panel's title, tag or letter ``(a)``, ``(b)``, ... (an untitled
panel shows its letter inside).  ``"none"`` draws no legend; an axes-level
``legend`` (``"axes"`` / ``"none"``) overrides the figure policy for that
panel.  A legend with more than the style's ``legend_max`` entries is not
drawn (with a warning); a shared legend over the cap that keys several panels
falls back to per-panel groups.  ``legend_title`` adds a bold heading to the
shared legend (or above the stack of groups).  Panel names are explicit: use
the axes ``title`` or ``tag``; a lone series is not turned into a tag.

**Themes.**  A theme has two independent halves, each chosen by name or a
plain dict: the *style* (everything except colour: rcParams, panel size and
gap, axes look, tick geometry, exponent and inside tick labels, the
line-style / marker ``cycle``, line width, marker size and spacing, legend and
tag formatting, the legend size cap, save options) and the *palette*
(colours).  Built-in styles:

| Name | Look |
| --- | --- |
| ``"paper"`` | Default.  Print size (1.9 × 1.6 in panels, 0.08 in gaps); inside y tick labels; exponent log axes. |
| ``"talk"`` | Large fonts and ticks, thick lines (5 × 4 in panels, 0.25 in gaps); ``10^k`` log labels. |

Built-in palettes:

| Name | Colours |
| --- | --- |
| ``"jaff"`` | Default: Okabe-Ito without yellow plus distinct GOW hues (14) |
| ``"okabe-ito"`` | Colour-blind safe (8) |
| ``"gow"`` | Distinct species hues (15) |
| ``"logo"`` | Logo colours (4) |
| ``"sunset"`` | Navy to peach (4) |
| ``"vivid"`` | Red, blue, green, purple, grey (5) |
| ``"muted"``, ``"deep"`` | seaborn-style cycles (10 each) |

A style dict is a partial override of a base style (``"base"``, default
``"paper"``); see :func:`register_style` for its fields.  A palette list is a
custom colour cycle; a palette dict maps **labels** to colours (other labels
use the default ``"jaff"`` cycle).  Add named themes with
:func:`register_style` / :func:`register_palette`.  Rendering applies the
style inside :func:`matplotlib.pyplot.rc_context`; global ``rcParams`` are
never changed.

**Cycling.**  Series without an explicit ``color`` / ``style`` are assigned
**colour first, then style**: series *i* gets ``palette[i % n]`` and
``cycle[i // n]``, so dashes and markers appear only once the colours run out
(built-in cycle: four linestyles, then markers ``o X P D s ^``).  A style
entry that is a marker draws hollow markers spaced evenly along the curve; a
fmt string such as ``"-o"`` or ``"--D"`` draws the line plus sparse hollow
markers.  The index counts only series without an explicit colour (resp.
style, or a label-mapped colour); bar series don't use the style cycle, and
the style cycle stays on its last entry once exhausted.  When a shared figure
legend keys several panels, the look is assigned per **label** across the
whole figure, so one label has one colour and style in every panel;
unlabelled series continue the cycle after the legend's looks, so they never
resemble a legend entry.  A series with ``"role": "reference"`` is outside
both cycles: dark grey ``#4d4d4d``, 1.9 × the style line width, drawn beneath
the other series.

Examples
--------
>>> from jaff import Network
>>> from jaff.plotting import plot_rates, plot_xsecs
>>> net = Network("GOW")
>>> r = list(net.reactions.with_type("two_body"))
>>> fig, axes = plot_rates(r[:3])  # one axes, three curves
>>> fig, axes = plot_rates([r[0], [r[1], r[2]]], grid=(1, 2))  # nested = overlay
>>> fig, axes = plot_rates(r[:4], grid=(2, 2), style="talk", palette="sunset")

Cross sections need a network with photo reactions:

>>> net = Network("h_photoionization", radiation={"bands": [1, 13.6, 100, "inf"]})
>>> photo = [p for p in net.reactions.photo_reactions() if p.photo.xsec]
>>> fig, axes = plot_xsecs(photo, bands=True, save="xsecs.pdf")
"""

from .api import plot, plot_rates, plot_xsecs
from .render import render
from .spec import AxesSpec, FigureSpec, SeriesSpec
from .themes import PALETTES, STYLES, register_palette, register_style

__all__ = [
    "plot",
    "plot_rates",
    "plot_xsecs",
    "render",
    "FigureSpec",
    "AxesSpec",
    "SeriesSpec",
    "STYLES",
    "PALETTES",
    "register_style",
    "register_palette",
]
