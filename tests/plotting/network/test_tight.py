# ABOUTME: Tests for the tight option of the net.visualize views: the figure is
# ABOUTME: cropped to its drawn content, with the legend kept in place

import pytest

from jaff.plotting.network import NetworkVisualizer
from jaff.plotting.network._common import TIGHT_PAD
from tests.plotting.network.test_layout import LINES

VIEWS = {
    "graph": lambda v, **kw: v.graph(**kw),
    "stoichiometry": lambda v, **kw: v.stoichiometry(**kw),
    "species": lambda v, **kw: v.species("H2", **kw),
}


def _content_inches(fig, ax):
    """Bounding box (inches) of everything drawn, the axes background excluded."""
    renderer = fig.canvas.get_renderer()
    artists = [*ax.patches, *ax.texts, ax.get_legend()]
    boxes = [a.get_window_extent(renderer) for a in artists]
    x0 = min(b.x0 for b in boxes) / fig.dpi
    x1 = max(b.x1 for b in boxes) / fig.dpi
    y0 = min(b.y0 for b in boxes) / fig.dpi
    y1 = max(b.y1 for b in boxes) / fig.dpi
    return x0, y0, x1, y1


@pytest.mark.parametrize("view", VIEWS)
def test_tight_figure_is_smaller(make_network, view):
    net = make_network(LINES)
    loose, _ = VIEWS[view](NetworkVisualizer(net), show=False, close=False)
    tight, _ = VIEWS[view](NetworkVisualizer(net), tight=True, show=False, close=False)

    assert tight.get_size_inches()[0] < loose.get_size_inches()[0]
    assert tight.get_size_inches()[1] < loose.get_size_inches()[1]


@pytest.mark.parametrize("view", VIEWS)
def test_tight_figure_hugs_its_content(make_network, view):
    net = make_network(LINES)
    fig, ax = VIEWS[view](NetworkVisualizer(net), tight=True, show=False, close=False)
    width, height = fig.get_size_inches()
    x0, y0, x1, y1 = _content_inches(fig, ax)

    assert x0 == pytest.approx(TIGHT_PAD, abs=0.02)
    assert y0 == pytest.approx(TIGHT_PAD, abs=0.02)
    assert width - x1 == pytest.approx(TIGHT_PAD, abs=0.02)
    assert height - y1 == pytest.approx(TIGHT_PAD, abs=0.02)


def test_tight_save(make_network, tmp_path):
    from matplotlib.image import imread

    from jaff.plotting.network._common import SAVE_DPI

    net = make_network(LINES)
    fig, _ = NetworkVisualizer(net).graph(tight=True, save=tmp_path / "g.png", show=False)

    assert imread(tmp_path / "g.png").shape[0] == pytest.approx(
        fig.get_size_inches()[1] * SAVE_DPI, abs=2
    )
