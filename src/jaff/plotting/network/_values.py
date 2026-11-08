# ABOUTME: Analysis values in the legend panel of rate-coloured figures: number
# ABOUTME: formatting, label templates, groups from net.analysis and their rows
"""The analysis values a rate-coloured figure used, listed in its legend panel."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Any

from matplotlib.axes import Axes
from sympy import Expr, Float, Symbol, latex

from ...core.analysis._evaluate import _function
from ._common import heading
from ._layout import (
    EQUATION_FONT,
    PANEL_GAP,
    PANEL_HEAD,
    VALUE_BLOCK_GAP,
    VALUE_ROW,
    LegendPanel,
    ValueGroup,
)
from ._themes import NetworkTheme

if TYPE_CHECKING:
    from ...core.analysis import Analysis

#: Title of each table's group, in drawing order.
TITLES = {
    "parameters": "PARAMETERS",
    "densities": r"NUMBER DENSITIES (cm$^{-3}$)",
    "functions": "FUNCTIONS",
}
#: Significant figures of a value.
FIGURES = 3


def value_tex(value: float) -> str:
    """*value* as mathtext without dollars: 3 significant figures.

    Plain in ``[1e-2, 1e3)``, else ``m\\times10^{n}`` (``10^{n}`` when ``m`` is 1).
    """
    if value == 0:
        return "0"
    sign = "-" if value < 0 else ""
    size = abs(value)
    exponent = math.floor(math.log10(size))
    mantissa = float(f"{size / 10**exponent:.{FIGURES - 1}f}")
    if mantissa >= 10:
        mantissa, exponent = mantissa / 10, exponent + 1
    if -2 <= exponent < 3:
        return f"{sign}{float(f'{size:.{FIGURES}g}'):g}"
    if mantissa == 1:
        return f"{sign}10^{{{exponent}}}"

    return rf"{sign}{mantissa:g}\times10^{{{exponent}}}"


def entry_text(label: str, tex: str) -> str:
    """A legend entry: *label* with the value *tex* at its ``{}``.

    Inside maths the value goes in bare, outside it in dollars; a label without
    ``{}`` reads ``label = value``, in one maths span when the label is one.
    """
    if "{}" not in label:
        if label.startswith("$") and label.endswith("$") and label.count("$") == 2:
            return f"{label[:-1]} = {tex}$"
        return f"{label} = ${tex}$"
    at = label.index("{}")
    inside = label[:at].replace(r"\$", "").count("$") % 2 == 1

    return label[:at] + (tex if inside else f"${tex}$") + label[at + 2 :]


def _function_tex(name: str, value: Any, parameters: dict[str, float]) -> str:
    """A function value: its number for the parameters, else its LaTeX."""
    expr: Expr = _function(name, value)
    numeric = expr.xreplace(
        {
            s: Float(parameters[s.name])
            for s in expr.free_symbols
            if isinstance(s, Symbol) and s.name in parameters
        }
    )
    if numeric.is_number:
        return value_tex(float(numeric.evalf()))

    return latex(expr)


def value_groups(analysis: Analysis) -> tuple[ValueGroup, ...]:
    """The labelled entries of *analysis*, one group per table that has any."""
    if not analysis.enabled:
        return ()
    params, dens, funcs = analysis.tables()
    numbers = params.values()
    groups = []
    for key, table in (("parameters", params), ("densities", dens), ("functions", funcs)):
        entries = []
        for name, entry in table.get().items():
            if entry.label is None or entry.value is None:
                continue
            if key == "functions":
                tex = _function_tex(name, entry.value, numbers)
            else:
                tex = value_tex(entry.value)
            entries.append(entry_text(entry.label, tex))
        if entries:
            groups.append(ValueGroup(TITLES[key], tuple(entries)))

    return tuple(groups)


def draw_values(
    ax: Axes,
    panel: LegendPanel,
    theme: NetworkTheme,
    groups: tuple[ValueGroup, ...],
    top: float,
) -> None:
    """The value rows of the legend *panel*, from *top* down (y down, inches).

    Each group is a heading over its entries, in its block's columns; blocks of
    one line sit side by side from the panel's left edge.
    """
    for line in panel.lines:
        x = panel.x0
        for block in line:
            group = groups[block.group]
            heading(ax, x, top, group.title, theme)
            for k, entry in enumerate(group.entries):
                row, column = divmod(k, block.columns)
                ax.text(
                    x + column * block.column_width,
                    top + PANEL_HEAD + (row + 0.5) * VALUE_ROW,
                    entry,
                    ha="left",
                    va="center",
                    fontsize=EQUATION_FONT,
                    color=theme.text,
                    zorder=3,
                )
            x += block.width + VALUE_BLOCK_GAP
        top += PANEL_HEAD + max(b.rows for b in line) * VALUE_ROW + PANEL_GAP
