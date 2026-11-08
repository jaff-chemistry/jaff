# ABOUTME: Tests for the legend panel under the graph and species views: rate bar,
# ABOUTME: reaction-type legend and value groups as headed rows, order and geometry

import matplotlib.pyplot as plt
import pytest
from matplotlib.colors import to_hex

from jaff.plotting.network import NETWORK_THEMES, NetworkVisualizer
from jaff.plotting.network._layout import (
    BAR_LENGTH,
    PANEL_HEAD,
    ValueGroup,
    legend_panel,
)
from jaff.plotting.network._themes import heading_colour

LINES = [
    "H2 + H2+ -> H + H3+ [10,1000] 2e-9",
    "H3+ + e- -> H + H2 [10,1000] 2e-8*(tgas/300)**(-0.5)",
    "H3+ + O -> H2 + OH+ [10,1000] 8e-10",
    "H2 -> H2+ + e- [] 1.2*crate",
]
K = "rate_coefficient"
HEADINGS = (
    "RATE COEFFICIENT $k$",
    "REACTION TYPE",
    "PARAMETERS",
    r"NUMBER DENSITIES (cm$^{-3}$)",
)


@pytest.fixture
def net(make_network):
    net = make_network(LINES)
    net.analysis.enable(
        parameters={
            "tgas": {"value": 300.0, "label": r"$T_\mathrm{gas} = {}$ K"},
            "crate": {"value": 1e-17, "label": True},
        },
        densities={
            "H2": {"value": 100.0, "label": True},
            "H": {"value": 1.0, "label": True},
        },
    )

    return net


def _drawn(net, view="graph", **kwargs):
    args = ("H3+",) if view == "species" else ()
    fig, ax = getattr(NetworkVisualizer(net), view)(
        *args, show=False, close=False, **kwargs
    )
    fig.canvas.draw()

    return fig, ax


def _text(ax, text):
    (found,) = [t for t in ax.texts if t.get_text() == text]

    return found


def _extent(artist):
    return artist.get_window_extent()


# --------------------------------------------------------------------------- #
# Rows                                                                         #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("view", ["graph", "species"])
def test_without_colour_by_only_the_type_row(net, view):
    fig, ax = _drawn(net, view)
    texts = [t.get_text() for t in ax.texts]
    plt.close(fig)

    assert "REACTION TYPE" in texts
    assert not {"RATE COEFFICIENT $k$", "PARAMETERS"} & set(texts)
    assert len(ax.images) == 0


@pytest.mark.parametrize("view", ["graph", "species"])
def test_rows_run_rate_types_values_top_to_bottom(net, view):
    fig, ax = _drawn(net, view, colour_by=K)
    rate = _extent(_text(ax, "RATE COEFFICIENT $k$"))
    bar = _extent(ax.images[0])
    types = _extent(_text(ax, "REACTION TYPE"))
    legend = _extent(ax.get_legend())
    values = _extent(_text(ax, "PARAMETERS"))
    plt.close(fig)

    # display y points up
    assert rate.y0 > bar.y1 > types.y1
    assert types.y0 > legend.y1 - 1
    assert legend.y0 > values.y1
    assert values.y0 > 0


def test_rate_heading_follows_the_quantity(net):
    net.analysis.enable(densities=[1.0] * len(net.species))
    fig, ax = _drawn(net, colour_by="rate")
    texts = [t.get_text() for t in ax.texts]
    plt.close(fig)

    assert r"REACTION RATE $\mathcal{R}$ (cm$^{-3}$ s$^{-1}$)" in texts


def test_bar_is_centred_and_full_length(net):
    fig, ax = _drawn(net, "species", colour_by=K)
    x0, x1, _, _ = ax.images[0].get_extent()
    width = fig.get_size_inches()[0]
    plt.close(fig)

    assert (x0 + x1) / 2 == pytest.approx(width / 2)
    assert x1 - x0 == pytest.approx(min(BAR_LENGTH, width - 0.5))


def test_types_and_values_share_a_left_edge(net):
    fig, ax = _drawn(net, "species", colour_by=K)
    types = _text(ax, "REACTION TYPE").get_position()[0]
    values = _text(ax, "PARAMETERS").get_position()[0]
    legend_x = ax.transData.inverted().transform(
        (ax.get_legend().get_bbox_to_anchor().bounds[0], 0)
    )[0]
    plt.close(fig)

    assert types == pytest.approx(values)
    assert legend_x == pytest.approx(types)


def test_value_groups_sit_side_by_side(net):
    fig, ax = _drawn(net, "species", colour_by=K)
    params = _text(ax, "PARAMETERS").get_position()
    dens = _text(ax, r"NUMBER DENSITIES (cm$^{-3}$)").get_position()
    plt.close(fig)

    assert dens[1] == pytest.approx(params[1])
    assert dens[0] > params[0]


def test_headings_are_bold_in_the_text_colour(net):
    theme = NETWORK_THEMES["vivid"]
    fig, ax = _drawn(net, colour_by=K)
    headings = [t for t in ax.texts if t.get_text() in HEADINGS]
    plt.close(fig)

    assert len(headings) == 4
    for t in headings:
        assert to_hex(t.get_color()) == to_hex(theme.text)
        assert t.get_fontweight() == "bold"


def test_panel_never_widens_the_figure(net):
    plain, _ = _drawn(net, "species")
    coloured, _ = _drawn(net, "species", colour_by=K)
    plt.close("all")

    assert coloured.get_size_inches()[0] == pytest.approx(plain.get_size_inches()[0])


def test_no_frames_or_rules(net):
    fig, ax = _drawn(net, colour_by=K)
    plt.close(fig)

    assert not [p for p in ax.patches if p.get_gid() == "values-box"]


# --------------------------------------------------------------------------- #
# Geometry                                                                     #
# --------------------------------------------------------------------------- #


def _measure(text, size):
    return 0.1 * len(text)


GROUPS = (
    ValueGroup("PARAMETERS", ("a = 1", "b = 2", "c = 3")),
    ValueGroup("DENSITIES", ("n(H) = 1", "n(H2) = 2")),
)


def test_groups_share_a_line_when_they_fit():
    _, _, panel = legend_panel([0.5] * 3, 20.0, _measure, values=GROUPS)

    assert [len(line) for line in panel.lines] == [2]


def test_groups_wrap_when_they_do_not_fit():
    _, _, panel = legend_panel([0.5] * 3, 2.6, _measure, values=GROUPS)

    assert [len(line) for line in panel.lines] == [1, 1]


def test_a_block_too_wide_alone_gets_one_column():
    _, _, panel = legend_panel([0.5], 1.6, _measure, values=GROUPS[:1])

    assert panel.lines[0][0].columns == 1


def test_panel_height_counts_its_rows():
    _, _, types = legend_panel([0.5] * 3, 20.0, _measure)
    _, _, rated = legend_panel([0.5] * 3, 20.0, _measure, colour_bar=True)
    _, _, full = legend_panel([0.5] * 3, 20.0, _measure, colour_bar=True, values=GROUPS)

    assert types.height < rated.height < full.height
    assert full.height - rated.height > PANEL_HEAD


# --------------------------------------------------------------------------- #
# Lane headers                                                                 #
# --------------------------------------------------------------------------- #


def test_heading_colour_is_between_muted_and_text():
    theme = NETWORK_THEMES["vivid"]
    colour = heading_colour(theme)

    assert to_hex(colour) not in (to_hex(theme.muted), to_hex(theme.text))


def test_lane_headers_use_the_heading_colour(net):
    fig, ax = _drawn(net, "species")
    header = _text(ax, "REACTANTS")
    plt.close(fig)

    assert to_hex(header.get_color()) == to_hex(heading_colour(NETWORK_THEMES["vivid"]))
