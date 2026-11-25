# ABOUTME: Tests for the stoichiometry view (net.visualize.stoichiometry): cells,
# ABOUTME: unsigned counts, index circles, centre-anchored rotated species labels

import pytest
from matplotlib.colors import to_hex

from jaff import Network
from jaff.plotting.network import NETWORK_THEMES, NetworkVisualizer
from tests.plotting.network.test_layout import LINES


@pytest.fixture
def drawn(make_network):
    net = make_network(LINES)
    fig, ax = NetworkVisualizer(net).stoichiometry(show=False, close=False)
    return net, fig, ax


def _cells(ax):
    return {p.get_gid(): p for p in ax.patches if (p.get_gid() or "").startswith("cell:")}


def test_one_cell_per_nonzero_entry(drawn):
    net, _, ax = drawn
    net_change = net.product_matrix - net.reactant_matrix
    nonzero = {f"cell:{i}:{net.species[j].name}" for i, j in zip(*net_change.nonzero())}

    assert set(_cells(ax)) == nonzero


def test_counts_are_unsigned(drawn):
    _, _, ax = drawn
    counts = {t.get_text() for t in ax.texts if t.get_gid() == "count"}

    assert counts <= {"1", "2", "3", "4"}
    assert "2" in counts  # H + H -> H2 consumes two H


def test_consumed_cells_use_theme_colours(drawn):
    net, _, ax = drawn
    h = net.species["H"].name
    cell = _cells(ax)[f"cell:0:{h}"]  # reaction 0 consumes 2 H

    assert to_hex(cell.get_facecolor()) == NETWORK_THEMES["vivid"].consumed[1]


def test_one_index_circle_per_reaction(drawn):
    net, _, ax = drawn
    gids = {
        p.get_gid() for p in ax.patches if (p.get_gid() or "").startswith("reaction:")
    }

    assert gids == {f"reaction:{i}" for i in range(net.reactions.count)}


def test_species_labels_centred_and_rotated(drawn):
    _, _, ax = drawn
    labels = [t for t in ax.texts if (t.get_gid() or "").startswith("label:")]

    assert labels
    assert all(t.get_rotation() == pytest.approx(45) for t in labels)
    assert all(t.get_ha() == "center" and t.get_va() == "center" for t in labels)
    assert all(t.get_rotation_mode() == "anchor" for t in labels)
    assert len({t.get_position()[1] for t in labels}) == 1


def test_save_and_bundled_network(tmp_path):
    NetworkVisualizer(Network("GOW")).stoichiometry(save=tmp_path / "s.pdf", show=False)

    assert (tmp_path / "s.pdf").stat().st_size > 0


def test_legend_is_centred(drawn):
    _, _, ax = drawn
    legend = ax.get_legend()

    assert legend._loc == 8  # "lower center"
    assert legend.get_bbox_to_anchor().bounds[0] == pytest.approx(ax.bbox.width / 2)
