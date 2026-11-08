# ABOUTME: Tests for the analysis values of rate-coloured figures: value formatting,
# ABOUTME: label templates, groups from net.analysis and their rows in the legend

import matplotlib.pyplot as plt
import pytest
import sympy as sp

from jaff.plotting.network import NetworkVisualizer
from jaff.plotting.network._values import entry_text, value_groups, value_tex

LINES = [
    "H2 + H2+ -> H + H3+ [10,1000] 2e-9",
    "H3+ + e- -> H + H2 [10,1000] 2e-8*(tgas/300)**(-0.5)",
    "H2 -> H2+ + e- [] 1.2*crate*fshield(av)",
]
K = "rate_coefficient"


@pytest.fixture
def net(make_network):
    net = make_network(LINES, funcfile=False)
    net.analysis.enable(
        parameters={
            "tgas": {"value": 300.0, "label": r"$T_\mathrm{gas} = {}$ K"},
            "crate": 1e-17,
            "av": 2.0,
        },
        densities={"H2": {"value": 100.0, "label": True}},
        functions={"fshield": {"value": sp.exp(-sp.Symbol("av")), "label": True}},
    )

    return net


# --------------------------------------------------------------------------- #
# Formatting                                                                   #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    ("value", "tex"),
    [
        (20.0, "20"),
        (0.02814, "0.0281"),
        (1.0, "1"),
        (999.0, "999"),
        (0.0, "0"),
        (2e-16, r"2\times10^{-16}"),
        (1e18, "10^{18}"),
        (1234.5, r"1.23\times10^{3}"),
        (-3e-5, r"-3\times10^{-5}"),
        (9.999e-6, "10^{-5}"),
    ],
)
def test_value_tex(value, tex):
    assert value_tex(value) == tex


def test_placeholder_inside_maths_takes_the_bare_value():
    assert entry_text(r"$T = {}$ K", "20") == r"$T = 20$ K"


def test_placeholder_outside_maths_takes_the_value_in_maths():
    assert entry_text("{} K", r"2\times10^{3}") == r"$2\times10^{3}$ K"


def test_maths_label_without_placeholder_joins_the_value_in_maths():
    assert entry_text("$A_V$", "10") == "$A_V = 10$"


def test_text_label_without_placeholder_reads_label_equals_value():
    assert entry_text("visual extinction", "10") == "visual extinction = $10$"


# --------------------------------------------------------------------------- #
# Groups                                                                       #
# --------------------------------------------------------------------------- #


def test_groups_hold_the_labelled_entries_only(net):
    groups = value_groups(net.analysis)

    assert [g.title for g in groups] == [
        "PARAMETERS",
        r"NUMBER DENSITIES (cm$^{-3}$)",
        "FUNCTIONS",
    ]
    assert groups[0].entries == (r"$T_\mathrm{gas} = 300$ K",)
    assert groups[1].entries == (r"$n({\rm H_{2}}) = 100$",)


def test_a_function_shows_its_value_for_the_parameters(net):
    (_, _, functions) = value_groups(net.analysis)

    assert functions.entries == ("fshield = $0.135$",)


def test_groups_without_labels_are_left_out(net):
    net.analysis.densities.set("H2", 100.0, label=False)
    net.analysis.functions.set("fshield", 0.5, label=False)

    assert [g.title for g in value_groups(net.analysis)] == ["PARAMETERS"]


def test_disabled_analysis_has_no_groups(net):
    net.analysis.disable()

    assert value_groups(net.analysis) == ()


# --------------------------------------------------------------------------- #
# Figures                                                                      #
# --------------------------------------------------------------------------- #


def _texts(ax):
    return [t.get_text() for t in ax.texts]


@pytest.mark.parametrize("view", ["graph", "species"])
def test_coloured_figure_lists_the_values(net, view):
    args = ("H3+",) if view == "species" else ()
    _, ax = getattr(NetworkVisualizer(net), view)(*args, show=False, colour_by=K)

    assert "PARAMETERS" in _texts(ax)
    assert r"$T_\mathrm{gas} = 300$ K" in _texts(ax)
    assert "FUNCTIONS" in _texts(ax)


def test_no_values_without_colour_by(net):
    _, ax = NetworkVisualizer(net).graph(show=False)

    assert "PARAMETERS" not in _texts(ax)


def test_no_values_without_labels(net):
    for table in (
        net.analysis.parameters,
        net.analysis.densities,
        net.analysis.functions,
    ):
        for name in table.names:
            entry = table.get(name)
            if entry.label is not None:
                table.set(name, entry.value, label=False)
    _, ax = NetworkVisualizer(net).graph(show=False, colour_by=K)

    assert not {"PARAMETERS", "FUNCTIONS"} & set(_texts(ax))


def test_values_make_the_figure_taller(net):
    with_values, _ = NetworkVisualizer(net).graph(show=False, colour_by=K)
    net.analysis.parameters.set("tgas", 300.0, label=False)
    net.analysis.densities.set("H2", 100.0, label=False)
    net.analysis.functions.set("fshield", 0.5, label=False)
    without, _ = NetworkVisualizer(net).graph(show=False, colour_by=K)

    assert with_values.get_size_inches()[1] > without.get_size_inches()[1]


def test_values_sit_inside_the_figure_below_the_type_legend(net):
    fig, ax = NetworkVisualizer(net).graph(show=False, close=False, colour_by=K)
    fig.canvas.draw()
    legend = ax.get_legend().get_window_extent()
    entries = [
        t.get_window_extent()
        for t in ax.texts
        if t.get_text() in ("PARAMETERS", r"$T_\mathrm{gas} = 300$ K")
    ]
    plt.close(fig)

    for e in entries:
        assert e.y1 < legend.y0 and e.y0 > 0
        assert 0 < e.x0 and e.x1 < fig.bbox.width


def test_tight_crop_keeps_the_values(net):
    fig, ax = NetworkVisualizer(net).graph(
        show=False, close=False, tight=True, colour_by=K
    )
    fig.canvas.draw()
    (entry,) = [t for t in ax.texts if t.get_text() == r"$T_\mathrm{gas} = 300$ K"]
    box = entry.get_window_extent()
    plt.close(fig)

    assert box.y0 >= -1 and box.x0 >= -1 and box.x1 <= fig.bbox.width + 1
