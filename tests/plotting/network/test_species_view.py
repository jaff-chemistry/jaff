# ABOUTME: Tests for the species view (net.visualize.species): formation and
# ABOUTME: destruction of one species, layout geometry and rendering

import matplotlib.pyplot as plt
import pytest
from matplotlib.colors import to_hex

from jaff.plotting.network import NETWORK_THEMES, NetworkVisualizer
from jaff.plotting.network._layout import (
    MARGIN,
    MARKER,
    ROW_PITCH,
    species_layout,
)
from tests.plotting.network.test_layout import fake_measure

# H3+ is formed by reaction 0 and destroyed by 1-3; He in 4 is a catalyst only.
LINES = [
    "H2 + H2+ -> H + H3+ [10,1000] 2e-9",
    "H3+ + e- -> H + H2 [10,1000] 2e-8",
    "H3+ + O -> H2 + OH+ [10,1000] 8e-10",
    "H3+ + CO -> H2 + HCO+ [10,1000] 2e-9",
    "H2 + He -> H + H + He [10,1000] 1e-15",
]


@pytest.fixture
def net(make_network):
    return make_network(LINES)


def _names(net, indices):
    return {net.species[i].name for i in indices}


class TestSpeciesLayout:
    def test_formation_and_destruction(self, net):
        layout = species_layout(net, net.species["H3+"], fake_measure)

        assert layout.formation == (0,)
        assert layout.destruction == (1, 2, 3)

    def test_catalyst_only_reactions_are_left_out(self, net):
        layout = species_layout(net, net.species["He"], fake_measure)

        assert layout.formation == () and layout.destruction == ()

    def test_partner_columns(self, net):
        layout = species_layout(net, net.species["H3+"], fake_measure)

        assert _names(net, layout.reactant_xy) == {"H2", "H2+"}
        assert _names(net, layout.product_xy) == {"H", "H2", "OH+", "HCO+"}

    def test_lanes(self, net):
        layout = species_layout(net, net.species["H3+"], fake_measure)

        assert [lane.header for lane in layout.lanes] == [
            "REACTANTS",
            "FORMS",
            None,
            "DESTROYS",
            "PRODUCTS",
        ]
        assert layout.reaction_xy[0][0] == pytest.approx(layout.lanes[1].centre)
        assert all(
            layout.reaction_xy[i][0] == pytest.approx(layout.lanes[3].centre)
            for i in (1, 2, 3)
        )

    def test_target_is_centred_on_the_rows(self, net):
        layout = species_layout(net, net.species["H3+"], fake_measure)
        ys = [layout.reaction_xy[i][1] for i in layout.destruction]

        assert layout.target_xy[0] == pytest.approx(layout.lanes[2].centre)
        assert layout.target_xy[1] == pytest.approx((min(ys) + max(ys)) / 2)

    def test_partners_a_row_apart(self, net):
        layout = species_layout(net, net.species["H3+"], fake_measure)
        ys = sorted(y for _, y in layout.product_xy.values())

        assert all(b - a >= ROW_PITCH - 1e-9 for a, b in zip(ys, ys[1:]))

    def test_partner_beside_its_only_reaction(self, net):
        layout = species_layout(net, net.species["H3+"], fake_measure)
        h = net.species["H"].index  # a product of reaction 1 only

        assert layout.product_xy[h][1] == pytest.approx(layout.reaction_xy[1][1])

    def test_edges(self, net):
        layout = species_layout(net, net.species["H3+"], fake_measure)
        target = net.species["H3+"].index
        roles = {(e.specie, e.reaction): e.role for e in layout.edges}

        assert roles[(target, 0)] == "product"  # reaction 0 forms H3+
        assert roles[(target, 1)] == "reactant"  # H3+ feeds reaction 1
        assert roles[(net.species["H2+"].index, 0)] == "reactant"
        assert roles[(net.species["OH+"].index, 2)] == "product"

    def test_edges_meet_the_markers(self, net):
        layout = species_layout(net, net.species["H3+"], fake_measure)
        target = net.species["H3+"].index
        into = next(e for e in layout.edges if e.specie == target and e.reaction == 1)
        rx, ry = layout.reaction_xy[1]

        assert into.points[-1] == pytest.approx((rx - MARKER / 2, ry))

    def test_species_without_reactions(self, make_network):
        net = make_network(["H + H -> H2 [10,1000] 1e-10"])
        layout = species_layout(net, net.species["H"], fake_measure)

        assert layout.formation == ()
        assert layout.destruction == (0,)

    def test_wide_enough_for_the_legend(self, net):
        layout = species_layout(net, net.species["H3+"], fake_measure, legend_width=20.0)

        assert layout.size[0] >= 20.0 + 2 * MARGIN


class TestSpeciesView:
    @pytest.fixture
    def drawn(self, net):
        fig, ax = NetworkVisualizer(net).species("H3+", show=False, close=False)
        return fig, ax

    def _gids(self, ax, prefix):
        return {p.get_gid() for p in ax.patches if (p.get_gid() or "").startswith(prefix)}

    def test_pills(self, drawn):
        _, ax = drawn

        assert self._gids(ax, "target:") == {"target:H3+"}
        assert self._gids(ax, "reactant:") == {"reactant:H2", "reactant:H2+"}
        assert self._gids(ax, "product:") == {
            "product:H",
            "product:H2",
            "product:OH+",
            "product:HCO+",
        }

    def test_markers(self, drawn):
        _, ax = drawn

        assert self._gids(ax, "reaction:") == {f"reaction:{i}" for i in range(4)}

    def test_legend_lists_only_the_types_shown(self, make_network):
        net = make_network(
            ["H2 -> H + H  []  PHOTO, 11.2", "H + H -> H2 [10,1000] 1e-10"]
        )
        _, ax = NetworkVisualizer(net).species("H", show=False, close=False)
        labels = [t.get_text() for t in ax.get_legend().get_texts()]
        handles = ax.get_legend().legend_handles
        colour = {t: to_hex(h.get_facecolor()) for t, h in zip(labels[:2], handles)}

        assert labels[:2] == ["photo", "other"]
        assert colour["photo"] == NETWORK_THEMES["vivid"].types[0]

    def test_unknown_species_raises(self, net):
        with pytest.raises(KeyError):
            NetworkVisualizer(net).species("XYZ", show=False)

    def test_save_and_close(self, net, tmp_path):
        NetworkVisualizer(net).species("H3+", save=tmp_path / "h3p.svg", show=False)

        assert (tmp_path / "h3p.svg").stat().st_size > 0
        assert plt.get_fignums() == []

    def test_through_the_network(self, net):
        fig, _ = net.visualize.species("H3+", show=False, close=False)

        assert fig is not None
