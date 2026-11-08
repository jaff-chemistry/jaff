# ABOUTME: Geometry for net.visualize: layered graph and stoichiometry-matrix layouts
# ABOUTME: Positions in inches; text widths come from a measure callable (no matplotlib)
"""Layout of the network figures, in inches, with the y axis pointing down."""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, NamedTuple

import numpy as np

from ...chemistry._periodic_table import periodic_table

if TYPE_CHECKING:
    from ...core.network.network import Network
    from ...core.reaction.reaction import Reaction
    from ...core.species.specie import Specie

#: Width in inches of a (mathtext) string at a font size in points.
Measure = Callable[[str, float], float]

#: Equation text for pseudo-species; other pseudo-species show their name.
PSEUDO_LABELS: dict[str, str] = {
    "_PHOTON": r"\gamma",
    "_CR": r"{\rm CR}",
    "_CRP": r"{\rm CR}",
    "_CRPHOT": r"{\rm CR}\gamma",
}

SPECIES_FONT = 9.0
EQUATION_FONT = 9.0
HEADER_FONT = 8.5
INDEX_FONT = 7.0
ROW_PITCH = 0.6
LANE_BELOW = 0.2
PILL_HEIGHT = 0.22
PILL_PAD = 0.1
MARKER = 0.17
GUTTER_PAD = 0.25
MARGIN = 0.25
HEADER = 0.35
LEGEND = 0.45
SWEEPS = 4
CELL = 0.26
INDEX_COLUMN = 0.3
LABEL_GAP = 0.1
GROUP_GAP = 0.06
EMPTY_WIDTH = 1.0
#: Height of each legend row after the first, inches.
LEGEND_ROW = 0.22
#: Legend panel under the graph and species views, inches: space above it, gap
#: between its rows, heading height, colour bar body (bar and tick labels) and
#: its longest length, value row pitch, gap between value columns and blocks.
PANEL_SPACE = 0.25
PANEL_GAP = 0.14
PANEL_HEAD = 0.22
RATE_BODY = 0.36
#: Row pitch of the reaction-type legend in the panel (9 pt, 0.5 label spacing).
TYPE_ROW = 0.19
BAR_LENGTH = 4.0
VALUE_ROW = 0.18
VALUE_COLUMNS = 2
VALUE_COLUMN_GAP = 0.25
VALUE_BLOCK_GAP = 0.6


class ValueGroup(NamedTuple):
    """One headed group of analysis values: its heading and entry texts."""

    title: str
    entries: tuple[str, ...]


@dataclass(frozen=True)
class ValueBlock:
    """Place of one value group in the legend panel.

    Attributes
    ----------
    group : int
        Index of the group.
    columns : int
        Entry columns.
    column_width : float
        Width of one column (widest entry and the gap), inches.
    width : float
        Width of the block (columns or heading, the wider), inches.
    rows : int
        Entry rows.
    """

    group: int
    columns: int
    column_width: float
    width: float
    rows: int


@dataclass(frozen=True)
class LegendPanel:
    """Geometry of the legend panel under a graph or species view.

    Attributes
    ----------
    x0 : float
        Left edge of the type legend and value rows, inches.
    rate : bool
        Whether the rate colour bar row is drawn.
    lines : tuple of tuple of ValueBlock
        Value blocks, line by line.
    height : float
        Height of the panel, the space above it included, inches.
    """

    x0: float = MARGIN
    rate: bool = False
    lines: tuple[tuple[ValueBlock, ...], ...] = ()
    height: float = PANEL_SPACE + PANEL_HEAD + TYPE_ROW


def _value_block(
    index: int, group: ValueGroup, measure: Measure, columns: int
) -> ValueBlock:
    column = max(measure(e, EQUATION_FONT) for e in group.entries) + VALUE_COLUMN_GAP
    columns = max(1, min(columns, len(group.entries)))
    width = max(columns * column - VALUE_COLUMN_GAP, measure(group.title, HEADER_FONT))

    return ValueBlock(
        index, columns, column, width, math.ceil(len(group.entries) / columns)
    )


def _value_lines(
    groups: Sequence[ValueGroup], available: float, measure: Measure
) -> tuple[tuple[ValueBlock, ...], ...]:
    """Value blocks packed into lines no wider than *available* inches."""
    lines: list[list[ValueBlock]] = []
    used = 0.0
    for i, group in enumerate(groups):
        block = _value_block(i, group, measure, VALUE_COLUMNS)
        if block.width > available:
            block = _value_block(i, group, measure, 1)
        if lines and used + VALUE_BLOCK_GAP + block.width <= available:
            lines[-1].append(block)
            used += VALUE_BLOCK_GAP + block.width
        else:
            lines.append([block])
            used = block.width

    return tuple(tuple(line) for line in lines)


def line_width(line: Sequence[ValueBlock]) -> float:
    """Width of one line of value blocks, inches."""
    return sum(b.width for b in line) + VALUE_BLOCK_GAP * (len(line) - 1)


def _table_width(widths: Sequence[float], columns: int) -> float:
    """Width of a legend of *widths* in *columns* (filled column by column)."""
    if not widths:
        return 0.0
    per, large = divmod(len(widths), columns)
    sizes = [per + 1] * large + [per] * (columns - large)
    starts = [sum(sizes[:c]) for c in range(columns)]

    return sum(max(widths[a : a + n]) for a, n in zip(starts, sizes) if n)


def legend_panel(
    widths: Sequence[float],
    width: float,
    measure: Measure,
    colour_bar: bool = False,
    values: Sequence[ValueGroup] = (),
) -> tuple[int, int, LegendPanel]:
    """Type-legend columns and rows, and the panel under a drawing *width* wide.

    Rows top to bottom: the rate colour bar (with *colour_bar*), the reaction
    types (entries *widths* inches wide, as many columns as fit) and the value
    groups, side by side and wrapping. The type and value rows share a left
    edge, the widest of them centred; the panel never needs more than *width*.
    """
    available = width - 2 * MARGIN
    columns = legend_columns(widths, available)
    rows = max(math.ceil(len(widths) / columns), 1)
    lines = _value_lines(values, available, measure)
    content = max([_table_width(widths, columns), *(line_width(ln) for ln in lines)])
    height = PANEL_SPACE + PANEL_HEAD + rows * TYPE_ROW
    if colour_bar:
        height += PANEL_HEAD + RATE_BODY + PANEL_GAP
    for line in lines:
        height += PANEL_GAP + PANEL_HEAD + max(b.rows for b in line) * VALUE_ROW

    return (
        columns,
        rows,
        LegendPanel(max(MARGIN, (width - content) / 2), colour_bar, lines, height),
    )


def legend_columns(widths: Sequence[float], available: float) -> int:
    """Legend columns for the fewest rows that fit in *available* inches.

    For ``r`` rows the entries are spread evenly over ``ceil(n / r)`` columns.
    Matplotlib fills a legend column by column, giving the first ``n % k`` of
    ``k`` columns one entry more than the rest; each column is as wide as its
    widest entry. One column when nothing narrower fits (and for no entries).
    """
    n = len(widths)
    for rows in range(1, n + 1):
        k = math.ceil(n / rows)
        per, large = divmod(n, k)
        sizes = [per + 1] * large + [per] * (k - large)
        starts = [sum(sizes[:c]) for c in range(k)]
        if sum(max(widths[a : a + size]) for a, size in zip(starts, sizes)) <= available:
            return k

    return 1


def _legend(widths: Sequence[float], width: float) -> tuple[int, int]:
    """Columns and rows of a legend under a drawing *width* inches wide."""
    columns = legend_columns(widths, width - 2 * MARGIN)

    return columns, max(math.ceil(len(widths) / columns), 1)


def drawing_bottom(layout: GraphLayout | MatrixLayout | SpeciesLayout) -> float:
    """Lowest height (inches, y down) the drawing may reach: the legend's top."""
    if isinstance(layout, MatrixLayout):
        return layout.size[1] - MARGIN - LEGEND - (layout.legend_rows - 1) * LEGEND_ROW

    return layout.size[1] - MARGIN - layout.panel.height


def atom_count(specie: Specie) -> int:
    """Number of atoms in *specie*: its element tokens (charges, ``e-``, tags not)."""
    table = periodic_table()

    # Charges, electrons, grains and name tags are table entries without protons.
    return sum(
        token in table and table[token].protons is not None for token in specie.exploded
    )


def _column_key(specie: Specie) -> int:
    """Atom-count group: atom-free species share the 1-atom group."""
    return max(atom_count(specie), 1)


def _header(key: int, counts: set[int]) -> str:
    if key == 1:
        return "0–1 ATOMS" if 0 in counts else "1 ATOM"

    return f"{key} ATOMS"


def specie_label(specie: Specie) -> str:
    """Mathtext (without ``$``) for a species or pseudo-species."""
    if specie.is_special:
        return PSEUDO_LABELS.get(specie.name, rf"{{\rm {specie.name.lstrip('_')}}}")

    return specie.latex()


def equation_label(reaction: Reaction) -> str:
    """Mathtext equation ``$A + B \\rightarrow C$`` of *reaction*."""
    left = " + ".join(specie_label(s) for s in reaction.reactants)
    right = " + ".join(specie_label(s) for s in reaction.products)

    return f"${left} \\rightarrow {right}$"


def _groups(net: Network) -> tuple[list[Specie], dict[int, int], list[str]]:
    """Core species, their atom-count column, and the column headers."""
    core = [s for s in net.species if s.is_core]
    keys = sorted({_column_key(s) for s in core})
    column = {s.index: keys.index(_column_key(s)) for s in core}
    headers = [
        _header(k, {atom_count(s) for s in core if _column_key(s) == k}) for k in keys
    ]

    return core, column, headers


def _links(reaction: Reaction) -> list[tuple[int, str]]:
    """``(specie index, role)`` once per core species and side of *reaction*."""
    links: list[tuple[int, str]] = []
    for role, side in (
        ("reactant", reaction.reactants.core),
        ("product", reaction.products.core),
    ):
        for specie in side:
            if (specie.index, role) not in links:
                links.append((specie.index, role))

    return links


@dataclass(frozen=True)
class Lane:
    """A vertical strip: a species column (with a header) or a reaction gutter."""

    x0: float
    x1: float
    header: str | None

    @property
    def centre(self) -> float:
        """Horizontal centre of the lane."""
        return (self.x0 + self.x1) / 2


@dataclass(frozen=True)
class Edge:
    """A species-reaction connection; *points* are the Bézier waypoints."""

    specie: int
    reaction: int
    role: str
    points: tuple[tuple[float, float], ...]


@dataclass(frozen=True)
class GraphLayout:
    """Geometry of the graph view (inches, y down).

    ``lanes`` run left to right: species columns (with a header) and reaction
    gutters (header ``None``). ``species_lane`` / ``reaction_lane`` give each
    node's lane.
    """

    size: tuple[float, float]
    lanes: tuple[Lane, ...]
    species_lane: dict[int, int]
    reaction_lane: dict[int, int]
    species_xy: dict[int, tuple[float, float]]
    reaction_xy: dict[int, tuple[float, float]]
    pill_width: dict[int, float]
    labels: dict[int, str]
    edges: tuple[Edge, ...]
    legend_columns: int = 1
    legend_rows: int = 1
    panel: LegendPanel = LegendPanel()


def _row_y(row: float) -> float:
    return MARGIN + HEADER + (row + 0.5) * ROW_PITCH


def _row_of(y: float) -> float:
    return (y - MARGIN - HEADER) / ROW_PITCH - 0.5


def _channel(y0: float, y1: float) -> float:
    """The routing lane below the row nearest the mean height of *y0* and *y1*.

    Equations sit above their marker, so an edge running under a row never
    crosses one.
    """
    return _row_y(round(_row_of((y0 + y1) / 2))) + LANE_BELOW


def _sweep(
    members: list[list[tuple[str, int]]], links: list[list[tuple[int, str]]]
) -> None:
    """Order each lane by the mean row of its neighbours (barycentre sweeps)."""
    neighbours: dict[tuple[str, int], list[tuple[str, int]]] = defaultdict(list)
    for i, reaction_links in enumerate(links):
        for s, _ in reaction_links:
            neighbours[("r", i)].append(("s", s))
            neighbours[("s", s)].append(("r", i))

    row = {node: k for lane in members for k, node in enumerate(lane)}
    order = list(range(len(members)))
    for _ in range(SWEEPS):
        for sweep in (order, order[::-1]):
            for lane in (members[i] for i in sweep):

                def key(node: tuple[str, int]) -> tuple[float, int]:
                    near = neighbours[node]
                    bary = sum(row[n] for n in near) / len(near) if near else row[node]
                    return bary, node[1]

                lane.sort(key=key)
                row.update({node: k for k, node in enumerate(lane)})


def _sub_columns(
    entries: list[tuple[int, int]], target: int
) -> list[list[tuple[str, int]]]:
    """Split ``(reach, reaction)`` entries into sub-columns of at most *target*.

    Sorted by reach, so reactions reaching further right land further right.
    """
    if not entries:
        return []
    entries = sorted(entries)
    size = math.ceil(len(entries) / math.ceil(len(entries) / target))

    return [
        [("r", i) for _, i in entries[k : k + size]] for k in range(0, len(entries), size)
    ]


def _lanes(
    core: list[Specie],
    column: dict[int, int],
    headers: list[str],
    links: list[list[tuple[int, str]]],
) -> tuple[list[list[tuple[str, int]]], list[str | None]]:
    """Lane members (nodes) and headers, left to right.

    A reaction whose species share one column goes in that column's side gutter
    (left of the first column, right of any other); one spanning columns goes in
    the gutter between the middle of its span. Crowded gutters are split into
    sub-columns no taller than the tallest species column.
    """
    ncols = len(headers)
    own: dict[int, list[tuple[int, int]]] = defaultdict(list)
    between: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for i, reaction_links in enumerate(links):
        cols = [column[s] for s, _ in reaction_links] or [0]
        lo, hi = min(cols), max(cols)
        if lo == hi:
            own[lo].append((0, i))
        else:
            g = lo + (hi - 1 - lo) // 2
            between[g].append(((hi - g - 1) - (g - lo), i))

    target = max(sum(column[s.index] == c for s in core) for c in range(ncols))
    members: list[list[tuple[str, int]]] = []
    headers_: list[str | None] = []

    def gutters(subs: list[list[tuple[str, int]]]) -> None:
        members.extend(subs)
        headers_.extend([None] * len(subs))

    for c in range(ncols):
        if c == 0:
            gutters(_sub_columns(own[0], target))
        members.append([("s", s.index) for s in core if column[s.index] == c])
        headers_.append(headers[c])
        if c > 0:
            gutters(_sub_columns(own[c], target))
        if c < ncols - 1:
            gutters(_sub_columns(between[c], target) or [[]])

    return members, headers_


def graph_layout(
    net: Network,
    measure: Measure,
    legend: Sequence[float] = (),
    colour_bar: bool = False,
    values: Sequence[ValueGroup] = (),
) -> GraphLayout:
    """Lay out *net* as a layered species-reaction graph.

    Species stand in columns by atom count. Reactions sit in gutters: a reaction
    within one column beside that column, one spanning columns between them,
    and a crowded gutter splits into sub-columns. Each species sits at the mean
    row of its reactions; edges crossing lanes run in the lane below a row.

    Parameters
    ----------
    net : Network
        The network to draw.
    measure : callable
        ``measure(text, size)``: width in inches of mathtext *text* at *size* points.
    legend : sequence of float, optional
        Width in inches of each legend entry; a legend wider than the drawing
        wraps onto more rows (the figure grows taller, not wider).
    colour_bar : bool, optional
        Leave room in the legend panel for the rate colour bar.
    values : sequence of ValueGroup, optional
        Groups of analysis values listed at the bottom of the legend panel.

    Returns
    -------
    GraphLayout
        Node positions, edges and the figure size, in inches.
    """
    core, column, headers = _groups(net)
    reactions = list(net.reactions)
    labels = {i: equation_label(r) for i, r in enumerate(reactions)}
    if not core:
        width = 2 * MARGIN + EMPTY_WIDTH
        columns, legend_rows, panel = legend_panel(
            legend, width, measure, colour_bar, values
        )
        size = (width, 2 * MARGIN + HEADER + panel.height)
        return GraphLayout(
            size, (), {}, {}, {}, {}, {}, labels, (), columns, legend_rows, panel
        )

    links = [_links(r) for r in reactions]
    members, lane_headers = _lanes(core, column, headers, links)
    _sweep(members, links)
    lane_of = {node: k for k, lane in enumerate(members) for node in lane}

    pill = {s.index: measure(f"${s.latex()}$", SPECIES_FONT) + 2 * PILL_PAD for s in core}
    lanes: list[Lane] = []
    x = MARGIN
    for lane, header in zip(members, lane_headers):
        if header is not None:
            width = max([pill[n[1]] for n in lane], default=PILL_HEIGHT)
        else:
            texts = [measure(labels[n[1]], EQUATION_FONT) for n in lane]
            width = max([*texts, MARKER]) + 2 * GUTTER_PAD
        lanes.append(Lane(x, x + width, header))
        x += width

    reaction_xy: dict[int, tuple[float, float]] = {}
    for k, lane in enumerate(members):
        if lane_headers[k] is None:
            for row, (_, i) in enumerate(lane):
                reaction_xy[i] = (lanes[k].centre, _row_y(row))

    # Each species at the mean row of its reactions, at least one row apart.
    species_xy: dict[int, tuple[float, float]] = {}
    for k, lane in enumerate(members):
        if lane_headers[k] is None:
            continue
        previous = -1
        for _, s in lane:
            ys = [
                reaction_xy[i][1]
                for i, ls in enumerate(links)
                if any(t == s for t, _ in ls)
            ]
            wanted = round(_row_of(sum(ys) / len(ys))) if ys else 0
            previous = max(wanted, previous + 1)
            species_xy[s] = (lanes[k].centre, _row_y(previous))

    edges: list[Edge] = []
    for i, reaction_links in enumerate(links):
        rx, ry = reaction_xy[i]
        rlane = lane_of[("r", i)]
        for s, role in reaction_links:
            sx, sy = species_xy[s]
            slane = lane_of[("s", s)]
            side = 1 if sx > rx else -1
            at_species = (sx - side * pill[s] / 2, sy)
            at_reaction = (rx + side * MARKER / 2, ry)
            start, end = (
                (at_species, at_reaction)
                if role == "reactant"
                else (at_reaction, at_species)
            )
            points = [start]
            if abs(slane - rlane) > 1:
                lo, hi = sorted((slane, rlane))
                channel = _channel(sy, ry)
                via = [(lanes[lo + 1].x0, channel), (lanes[hi - 1].x1, channel)]
                points += via if start[0] < end[0] else via[::-1]
            points.append(end)
            edges.append(Edge(s, i, role, tuple(points)))

    bottom = max(y for _, y in [*species_xy.values(), *reaction_xy.values()])
    columns, legend_rows, panel = legend_panel(
        legend, x + MARGIN, measure, colour_bar, values
    )
    size = (x + MARGIN, bottom + ROW_PITCH / 2 + panel.height + MARGIN)

    return GraphLayout(
        size,
        tuple(lanes),
        {s.index: lane_of[("s", s.index)] for s in core},
        {i: lane_of[("r", i)] for i in range(len(reactions))},
        species_xy,
        reaction_xy,
        pill,
        labels,
        tuple(edges),
        columns,
        legend_rows,
        panel,
    )


@dataclass(frozen=True)
class Group:
    """An atom-count group of matrix columns."""

    x0: float
    x1: float
    header: str


@dataclass(frozen=True)
class MatrixLayout:
    """Geometry of the stoichiometry view (inches, y down).

    ``values[i, j]`` is the net change of species ``species[j]`` (a catalogue
    index) in reaction ``i``.
    """

    size: tuple[float, float]
    species: tuple[int, ...]
    values: np.ndarray = field(compare=False)
    labels: tuple[str, ...]
    groups: tuple[Group, ...]
    column_x: tuple[float, ...]
    row_y: tuple[float, ...]
    equation_x: float
    index_x: float
    label_y: float
    legend_columns: int = 1
    legend_rows: int = 1


def matrix_layout(
    net: Network, measure: Measure, legend: Sequence[float] = ()
) -> MatrixLayout:
    """Lay out the net stoichiometry matrix of *net*.

    Parameters
    ----------
    net : Network
        The network to draw.
    measure : callable
        ``measure(text, size)``: width in inches of mathtext *text* at *size* points.
    legend : sequence of float, optional
        Width in inches of each legend entry; a legend wider than the drawing
        wraps onto more rows (the figure grows taller, not wider).

    Returns
    -------
    MatrixLayout
        Cell, label and group positions and the figure size, in inches.
    """
    core, column, headers = _groups(net)
    order = sorted(core, key=lambda s: (column[s.index], s.index))
    species = tuple(s.index for s in order)
    values = (net.product_matrix - net.reactant_matrix)[:, list(species)]
    labels = tuple(equation_label(r) for r in net.reactions)

    equation_w = max([measure(t, EQUATION_FONT) for t in labels], default=0.0)
    species_w = max([measure(f"${s.latex()}$", SPECIES_FONT) for s in order], default=0.0)
    # A label rotated 45 degrees, centred on its column: half its extent above and below.
    label_h = (species_w + SPECIES_FONT / 72) * math.sqrt(0.5)
    label_y = MARGIN + HEADER + label_h / 2
    top = MARGIN + HEADER + label_h + LABEL_GAP

    equation_x = MARGIN + equation_w
    index_x = equation_x + LABEL_GAP + INDEX_COLUMN / 2
    x = equation_x + LABEL_GAP + INDEX_COLUMN + LABEL_GAP
    column_x: list[float] = []
    groups: list[Group] = []
    for c, header in enumerate(headers):
        x0 = x
        for s in order:
            if column[s.index] == c:
                column_x.append(x + CELL / 2)
                x += CELL
        groups.append(Group(x0, x, header))
        x += GROUP_GAP

    row_y = tuple(top + (i + 0.5) * CELL for i in range(len(labels)))
    columns, legend_rows = _legend(legend, x + MARGIN)
    size = (
        x + MARGIN,
        top + len(labels) * CELL + LEGEND + (legend_rows - 1) * LEGEND_ROW + MARGIN,
    )

    return MatrixLayout(
        size,
        species,
        values,
        labels,
        tuple(groups),
        tuple(column_x),
        row_y,
        equation_x,
        index_x,
        label_y,
        columns,
        legend_rows,
    )


@dataclass(frozen=True)
class SpeciesLayout:
    """Geometry of the species view (inches, y down).

    Five lanes, left to right: other reactants of the reactions forming the
    species, those reactions, the species itself, the reactions destroying it,
    and their other products.
    """

    size: tuple[float, float]
    lanes: tuple[Lane, ...]
    target: int
    target_xy: tuple[float, float]
    formation: tuple[int, ...]
    destruction: tuple[int, ...]
    reaction_xy: dict[int, tuple[float, float]]
    reactant_xy: dict[int, tuple[float, float]]
    product_xy: dict[int, tuple[float, float]]
    pill_width: dict[int, float]
    labels: dict[int, str]
    edges: tuple[Edge, ...]
    legend_columns: int = 1
    legend_rows: int = 1
    panel: LegendPanel = LegendPanel()


#: Space each side of a pill in the species view's outer lanes, inches.
SIDE_PAD = 0.15


def _stack(
    partners: dict[int, list[int]], reaction_y: dict[int, float], top: float
) -> dict[int, float]:
    """Heights of partner species: each at the mean of its reactions, a row apart."""
    mean = {s: sum(reaction_y[i] for i in rs) / len(rs) for s, rs in partners.items()}
    heights: dict[int, float] = {}
    previous = top - ROW_PITCH / 2
    for s in sorted(partners, key=lambda s: (mean[s], s)):
        previous = max(mean[s], previous + ROW_PITCH)
        heights[s] = previous

    return heights


def species_layout(
    net: Network,
    specie: Specie,
    measure: Measure,
    legend: Sequence[float] = (),
    colour_bar: bool = False,
    values: Sequence[ValueGroup] = (),
) -> SpeciesLayout:
    """Lay out the reactions forming and destroying *specie*.

    A reaction forms the species when its net change in it is positive and
    destroys it when negative; a reaction with the species on both sides and no
    net change is left out. Reactions stack in their lanes, centred on each
    other; each partner species sits at the mean height of its reactions.

    Parameters
    ----------
    net : Network
        The network the species belongs to.
    specie : Specie
        The species to draw.
    measure : callable
        ``measure(text, size)``: width in inches of mathtext *text* at *size* points.
    legend : sequence of float, optional
        Width in inches of each legend entry; a legend wider than the drawing
        wraps onto more rows (the figure grows taller, not wider).
    colour_bar : bool, optional
        Leave room in the legend panel for the rate colour bar.
    values : sequence of ValueGroup, optional
        Groups of analysis values listed at the bottom of the legend panel.

    Returns
    -------
    SpeciesLayout
        Node positions, edges and the figure size, in inches.
    """
    reactions = list(net.reactions)
    j = specie.index
    change = net.product_matrix[:, j] - net.reactant_matrix[:, j]
    formation = tuple(i for i in range(len(reactions)) if change[i] > 0)
    destruction = tuple(i for i in range(len(reactions)) if change[i] < 0)
    labels = {i: equation_label(reactions[i]) for i in formation + destruction}

    reactants: dict[int, list[int]] = defaultdict(list)
    for i in formation:
        for s in dict.fromkeys(reactions[i].reactants.core):
            if s.index != j:
                reactants[s.index].append(i)
    products: dict[int, list[int]] = defaultdict(list)
    for i in destruction:
        for s in dict.fromkeys(reactions[i].products.core):
            if s.index != j:
                products[s.index].append(i)

    rows = max(len(formation), len(destruction), len(reactants), len(products), 1)
    top = MARGIN + HEADER
    reaction_y: dict[int, float] = {}
    for column in (formation, destruction):
        offset = (rows - len(column)) / 2
        for k, i in enumerate(column):
            reaction_y[i] = _row_y(k + offset)

    pill = {
        s: measure(f"${net.species[s].latex()}$", SPECIES_FONT) + 2 * PILL_PAD
        for s in {j, *reactants, *products}
    }

    def equations(column: tuple[int, ...]) -> float:
        texts = [measure(labels[i], EQUATION_FONT) for i in column]
        return max([*texts, MARKER]) + 2 * GUTTER_PAD

    widths = [
        max([pill[s] for s in reactants], default=PILL_HEIGHT) + 2 * SIDE_PAD,
        equations(formation),
        pill[j] + 2 * GUTTER_PAD,
        equations(destruction),
        max([pill[s] for s in products], default=PILL_HEIGHT) + 2 * SIDE_PAD,
    ]
    headers = [
        "REACTANTS",
        "FORMATION CHANNEL",
        None,
        "DESTRUCTION CHANNEL",
        "PRODUCTS",
    ]
    lanes: list[Lane] = []
    x = MARGIN
    for width, header in zip(widths, headers):
        # A lane is never narrower than its header, so headers never overlap.
        if header is not None:
            width = max(width, measure(header, HEADER_FONT) + 2 * LABEL_GAP)
        lanes.append(Lane(x, x + width, header))
        x += width

    reaction_xy = {
        i: (lanes[1 if i in formation else 3].centre, y) for i, y in reaction_y.items()
    }
    target_xy = (lanes[2].centre, _row_y((rows - 1) / 2))
    reactant_xy = {
        s: (lanes[0].centre, y) for s, y in _stack(reactants, reaction_y, top).items()
    }
    product_xy = {
        s: (lanes[4].centre, y) for s, y in _stack(products, reaction_y, top).items()
    }

    tx, ty = target_xy
    edges: list[Edge] = []
    for i in formation:
        rx, ry = reaction_xy[i]
        for s in dict.fromkeys(t for t, rs in reactants.items() if i in rs):
            sx, sy = reactant_xy[s]
            edges.append(
                Edge(s, i, "reactant", ((sx + pill[s] / 2, sy), (rx - MARKER / 2, ry)))
            )
        edges.append(
            Edge(j, i, "product", ((rx + MARKER / 2, ry), (tx - pill[j] / 2, ty)))
        )
    for i in destruction:
        rx, ry = reaction_xy[i]
        edges.append(
            Edge(j, i, "reactant", ((tx + pill[j] / 2, ty), (rx - MARKER / 2, ry)))
        )
        for s in dict.fromkeys(t for t, rs in products.items() if i in rs):
            sx, sy = product_xy[s]
            edges.append(
                Edge(s, i, "product", ((rx + MARKER / 2, ry), (sx - pill[s] / 2, sy)))
            )

    bottom = max(
        y
        for _, y in [
            target_xy,
            *reaction_xy.values(),
            *reactant_xy.values(),
            *product_xy.values(),
        ]
    )
    columns, legend_rows, panel = legend_panel(
        legend, x + MARGIN, measure, colour_bar, values
    )
    size = (x + MARGIN, bottom + ROW_PITCH / 2 + panel.height + MARGIN)

    return SpeciesLayout(
        size,
        tuple(lanes),
        j,
        target_xy,
        formation,
        destruction,
        reaction_xy,
        reactant_xy,
        product_xy,
        pill,
        labels,
        tuple(edges),
        columns,
        legend_rows,
        panel,
    )
