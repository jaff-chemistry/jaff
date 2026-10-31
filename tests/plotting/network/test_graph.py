# ABOUTME: Tests for the graph view (net.visualize.graph): one pill per species,
# ABOUTME: one marker per reaction, theme colours, saving and closing

import matplotlib.pyplot as plt
import pytest
from matplotlib.colors import to_hex

from jaff import Network
from jaff.plotting.network import NETWORK_THEMES, NetworkVisualizer
from tests.plotting.network.test_layout import LINES


def _gid(ax, prefix):
    return {p.get_gid(): p for p in ax.patches if (p.get_gid() or "").startswith(prefix)}


@pytest.fixture
def drawn(make_network):
    net = make_network(LINES)
    fig, ax = NetworkVisualizer(net).graph(show=False, close=False)
    return net, fig, ax


def test_one_pill_per_core_species(drawn):
    net, _, ax = drawn

    assert set(_gid(ax, "specie:")) == {
        f"specie:{s.name}" for s in net.species if s.is_core
    }


def test_one_marker_per_reaction(drawn):
    net, _, ax = drawn

    assert set(_gid(ax, "reaction:")) == {
        f"reaction:{i}" for i in range(net.reactions.count)
    }


def test_marker_colour_follows_type(drawn):
    net, _, ax = drawn
    theme = NETWORK_THEMES["vivid"]
    marker = _gid(ax, "reaction:")["reaction:0"]
    expected = theme.other if net.reactions[0].type == "unknown" else theme.types[0]

    assert to_hex(marker.get_facecolor()) == expected


def test_figure_size_is_the_layout_size(drawn):
    _, fig, _ = drawn

    assert fig.get_size_inches()[1] > 1


def test_equations_are_drawn(drawn):
    _, _, ax = drawn
    texts = [t.get_text() for t in ax.texts]

    assert r"${\rm H} + {\rm H} \rightarrow {\rm H_{2}}$" in texts


def test_theme_by_name(make_network):
    net = make_network(LINES)
    fig, _ = NetworkVisualizer(net).graph(theme="dark", show=False, close=False)

    assert to_hex(fig.get_facecolor()) == NETWORK_THEMES["dark"].surface


def test_save_writes_a_file(make_network, tmp_path):
    net = make_network(LINES)
    NetworkVisualizer(net).graph(save=tmp_path / "net.svg", show=False)

    assert (tmp_path / "net.svg").stat().st_size > 0


def test_close_leaves_no_open_figure(make_network):
    net = make_network(LINES)
    NetworkVisualizer(net).graph(show=False)

    assert plt.get_fignums() == []


def test_empty_network_draws(make_network):
    fig, _ = NetworkVisualizer(make_network(["# none"])).graph(show=False, close=False)

    assert fig is not None


@pytest.mark.parametrize("name", ["GOW", "networks/COthin/react_COthin.jet"])
def test_bundled_networks_draw(name, tmp_path):
    NetworkVisualizer(Network(name)).graph(save=tmp_path / "g.png", show=False)

    assert (tmp_path / "g.png").exists()


def test_arrowheads_are_readable(drawn):
    from matplotlib.patches import FancyArrowPatch

    _, _, ax = drawn
    heads = [p for p in ax.patches if isinstance(p, FancyArrowPatch)]

    assert heads
    assert all(p.get_mutation_scale() >= 11 for p in heads)


def test_legend_is_centred(drawn):
    _, _, ax = drawn
    legend = ax.get_legend()

    assert legend._loc == 8  # "lower center"
    assert legend.get_bbox_to_anchor().bounds[0] == pytest.approx(ax.bbox.width / 2)


def test_raster_saves_at_print_resolution(make_network, tmp_path):
    from matplotlib.image import imread

    from jaff.plotting.network._common import SAVE_DPI

    net = make_network(LINES)
    fig, _ = NetworkVisualizer(net).graph(save=tmp_path / "g.png", show=False)
    height = imread(tmp_path / "g.png").shape[0]

    assert SAVE_DPI == 400
    assert height == pytest.approx(fig.get_size_inches()[1] * SAVE_DPI, rel=0.05)
