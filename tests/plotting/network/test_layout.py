# ABOUTME: Tests for the network-visualisation layout: atom counts, labels, lanes,
# ABOUTME: rows, edges and sizes of the graph and stoichiometry-matrix layouts

import numpy as np
import pytest

from jaff.core.species import Specie
from jaff.plotting.network._layout import (
    HEADER,
    LANE_BELOW,
    MARGIN,
    MARKER,
    ROW_PITCH,
    atom_count,
    equation_label,
    graph_layout,
    matrix_layout,
)


def fake_measure(text: str, size: float) -> float:
    """Width in inches: 0.6 em per character, independent of the font."""
    return 0.6 * size / 72 * len(text)


# Columns by atom count: 1 = H, H+, O, e-; 2 = H2, H2+, CO, CO+; 3 = H2O, HCN, HNC
LINES = [
    "H + H -> H2 [10,1000] 1e-10",  # 0: columns 0..1: between-gutter 0|1
    "H2 + H+ -> H2+ + H [10,1000] 1e-9",  # 1: columns 0..1: between-gutter 0|1
    "H2 + O -> H2O [10,1000] 1e-10",  # 2: columns 0..2: left middle gutter 0|1
    "H2+ + CO -> H2 + CO+ [10,1000] 1e-9",  # 3: column 1 only: its side gutter
    "HCN + H2O -> HNC + H2O [10,1000] 1e-9",  # 4: last column only: its side gutter
    "H+ + e- -> H [10,1000] 1e-12",  # 5: column 0 only: its side gutter (left)
]

# Columns: 1 = H, O; 2 = H2, O2, OH; 3 = H2O.  Eight reactions span the first
# gutter, the tallest column has 3 species: the gutter splits into 3 sub-columns.
CROWDED = [
    "H + H -> H2 [10,1000] 1e-10",
    "O + O -> O2 [10,1000] 1e-10",
    "H + O -> OH [10,1000] 1e-10",
    "H2 + O -> OH + H [10,1000] 1e-10",
    "O2 + H -> OH + O [10,1000] 1e-10",
    "OH + H -> H2 + O [10,1000] 1e-10",
    "OH + O -> O2 + H [10,1000] 1e-10",
    "H + OH -> H2O [10,1000] 1e-10",  # 7: reaches column 2
]


@pytest.fixture
def net(make_network):
    return make_network(LINES)


@pytest.mark.parametrize(
    "name, n",
    [
        ("H", 1),
        ("H2+", 2),
        ("e-", 0),
        ("H2O", 3),
        ("c-C3H2+", 5),
        ("CO_DUST", 2),
        ("GRAIN0", 0),
    ],
)
def test_atom_count_counts_element_tokens(name, n):
    assert atom_count(Specie(name)) == n


def test_equation_label_is_mathtext(make_network):
    net = make_network(["H + H -> H2 [10,1000] 1e-10"])

    assert (
        equation_label(net.reactions[0]) == r"${\rm H} + {\rm H} \rightarrow {\rm H_{2}}$"
    )


def test_equation_label_names_photons(make_network):
    net = make_network(["H2 -> H + H  []  PHOTO, 11.2"])

    assert r"\gamma" in equation_label(net.reactions[0])


def _column_of(layout, lane):
    """Header of the species column a lane belongs to (None for gutters)."""
    return layout.lanes[lane].header


def _column_lanes(layout):
    return [k for k, lane in enumerate(layout.lanes) if lane.header is not None]


def _on_a_row(y):
    """Whether *y* is a row centre."""
    row = (y - MARGIN - HEADER) / ROW_PITCH - 0.5
    return row == pytest.approx(round(row))


class TestGraphLayout:
    def test_species_columns_by_atom_count(self, net):
        layout = graph_layout(net, fake_measure)
        column = {
            net.species[i].name: _column_of(layout, lane)
            for i, lane in layout.species_lane.items()
        }

        assert column["H"] == column["H+"] == column["O"] == column["e-"] == "0–1 ATOMS"
        assert column["H2"] == column["CO"] == column["CO+"] == "2 ATOMS"
        assert column["H2O"] == column["HCN"] == "3 ATOMS"

    def test_headers(self, net):
        layout = graph_layout(net, fake_measure)

        assert [layout.lanes[k].header for k in _column_lanes(layout)] == [
            "0–1 ATOMS",
            "2 ATOMS",
            "3 ATOMS",
        ]

    def test_single_column_reactions_get_a_side_gutter(self, net):
        layout = graph_layout(net, fake_measure)
        first, middle, last = _column_lanes(layout)

        assert layout.reaction_lane[5] == first - 1  # left of the first column
        assert layout.reaction_lane[3] == middle + 1  # right of its column
        assert layout.reaction_lane[4] == last + 1  # right of the last column

    def test_spanning_reactions_sit_between_their_columns(self, net):
        layout = graph_layout(net, fake_measure)
        first, middle, _ = _column_lanes(layout)

        assert all(first < layout.reaction_lane[i] < middle for i in (0, 1, 2))

    def test_single_column_network(self, make_network):
        layout = graph_layout(
            make_network(["H+ + e- -> H [10,1000] 1e-12"]), fake_measure
        )

        assert [lane.header for lane in layout.lanes] == [None, "0–1 ATOMS"]
        assert layout.reaction_lane == {0: 0}

    def test_crowded_gutter_splits_into_sub_columns(self, make_network):
        layout = graph_layout(make_network(CROWDED), fake_measure)
        first, middle, _ = _column_lanes(layout)
        gutters = range(first + 1, middle)
        rows = [sum(lane == g for lane in layout.reaction_lane.values()) for g in gutters]

        assert len(gutters) == 3
        assert max(rows) <= 3

    def test_sub_columns_ordered_by_reach(self, make_network):
        layout = graph_layout(make_network(CROWDED), fake_measure)
        _, middle, _ = _column_lanes(layout)

        assert layout.reaction_lane[7] == middle - 1  # reaches furthest right

    def test_no_two_nodes_share_a_slot(self, net):
        layout = graph_layout(net, fake_measure)
        slots = list(layout.species_xy.values()) + list(layout.reaction_xy.values())

        assert len(set(slots)) == len(slots)

    def test_pseudo_species_are_not_nodes(self, make_network):
        net = make_network(["H2 -> H + H  []  PHOTO, 11.2"])
        layout = graph_layout(net, fake_measure)

        assert {net.species[i].name for i in layout.species_xy} == {"H", "H2"}

    def test_one_edge_per_species_role(self, net):
        layout = graph_layout(net, fake_measure)
        h = net.species["H"].index
        edges = [(e.specie, e.role) for e in layout.edges if e.reaction == 0]

        assert edges.count((h, "reactant")) == 1

    def test_reactant_edge_ends_at_the_marker(self, net):
        layout = graph_layout(net, fake_measure)
        edge = next(e for e in layout.edges if e.reaction == 0 and e.role == "reactant")
        rx, ry = layout.reaction_xy[0]

        assert edge.points[-1] == pytest.approx((rx - MARKER / 2, ry))

    def test_long_edges_use_the_lane_below_a_row(self, net):
        layout = graph_layout(net, fake_measure)
        h2o = net.species["H2O"].index
        h2 = net.species["H2"].index
        far = next(e for e in layout.edges if e.reaction == 2 and e.specie == h2o)
        near = next(e for e in layout.edges if e.reaction == 2 and e.specie == h2)

        # H2O (column 2) from the 0|1 gutter crosses column 1: a detour.
        assert len(far.points) == 4
        assert far.points[1][1] == far.points[2][1]
        assert _on_a_row(far.points[1][1] - LANE_BELOW)
        assert len(near.points) == 2

    def test_species_sit_beside_their_reactions(self, make_network):
        net = make_network(
            [
                "H + H -> H2 [10,1000] 1e-10",
                "O + O -> O2 [10,1000] 1e-10",
                "C + C -> C2 [10,1000] 1e-10",
            ]
        )
        layout = graph_layout(net, fake_measure)

        for i, reaction in enumerate(net.reactions):
            for specie in [*reaction.reactants.core, *reaction.products.core]:
                assert layout.species_xy[specie.index][1] == layout.reaction_xy[i][1]

    def test_species_on_rows_and_a_row_apart(self, net):
        layout = graph_layout(net, fake_measure)
        for lane in _column_lanes(layout):
            ys = sorted(
                layout.species_xy[s][1]
                for s, k in layout.species_lane.items()
                if k == lane
            )
            assert all(_on_a_row(y) for y in ys)
            assert all(b - a >= ROW_PITCH - 1e-9 for a, b in zip(ys, ys[1:]))

    def test_wide_enough_for_the_legend(self, net):
        layout = graph_layout(net, fake_measure, legend_width=30.0)

        assert layout.size[0] >= 30.0 + 2 * MARGIN

    def test_deterministic(self, net):
        assert graph_layout(net, fake_measure) == graph_layout(net, fake_measure)

    def test_size_grows_with_reactions(self, make_network):
        small = graph_layout(make_network(LINES[:2]), fake_measure)
        large = graph_layout(make_network(LINES), fake_measure)

        assert large.size[1] >= small.size[1]
        assert large.size[0] > small.size[0]

    def test_empty_network(self, make_network):
        layout = graph_layout(make_network(["# no reactions"]), fake_measure)

        assert layout.lanes == () and layout.edges == ()
        assert layout.size[0] > 0 and layout.size[1] > 0


class TestMatrixLayout:
    def test_values_are_net_change_in_column_order(self, net):
        layout = matrix_layout(net, fake_measure)
        expected = (net.product_matrix - net.reactant_matrix)[:, list(layout.species)]

        np.testing.assert_array_equal(layout.values, expected)

    def test_columns_grouped_by_atom_count(self, net):
        layout = matrix_layout(net, fake_measure)
        counts = [max(atom_count(net.species[i]), 1) for i in layout.species]

        assert counts == sorted(counts)
        assert [g.header for g in layout.groups] == ["0–1 ATOMS", "2 ATOMS", "3 ATOMS"]

    def test_catalogue_order_within_a_group(self, net):
        layout = matrix_layout(net, fake_measure)
        group = [i for i in layout.species if max(atom_count(net.species[i]), 1) == 1]

        assert group == sorted(group)

    def test_one_row_per_reaction(self, net):
        layout = matrix_layout(net, fake_measure)

        assert layout.values.shape[0] == len(layout.row_y) == net.reactions.count
        assert np.diff(layout.row_y) == pytest.approx(
            [layout.row_y[1] - layout.row_y[0]] * 5
        )

    def test_labels_clear_the_cells(self, net):
        layout = matrix_layout(net, fake_measure)

        assert layout.label_y < layout.row_y[0]

    def test_equations_left_of_index_left_of_cells(self, net):
        layout = matrix_layout(net, fake_measure)

        assert layout.equation_x < layout.index_x < layout.column_x[0]
