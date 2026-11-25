# ABOUTME: Tests for jaff.plotting.render — artists produced from a FigureSpec
# ABOUTME: Grid, hidden cells, scales, cycling, markers, cleaning, bars, rc scoping

import matplotlib as mpl
import matplotlib.colors as mcolors
import numpy as np
import pytest
from matplotlib.collections import LineCollection

from jaff.plotting._marks import mark_indices
from jaff.plotting._ticks import exp_label
from jaff.plotting import render
from jaff.plotting.spec import FigureSpec

X = np.logspace(0, 2, 50)


def _render(axes, **kw):
    return render(FigureSpec.from_dict({"axes": axes, **kw}))


def _visible_lines(ax):
    return [ln for ln in ax.lines if ln.get_visible()]


def test_grid_shape_and_hidden_cells():
    fig, axes = _render([{"series": []}] * 3, grid=(2, 2))
    assert axes.shape == (2, 2)
    assert [a.axison for a in axes.ravel()] == [True, True, True, False]


def test_axes_box_from_panel_size():
    # panel_size is the axes box (spec §14.1); margins and gaps come on top
    fig, axes = _render([{"series": []}] * 2, grid=(1, 2), style={"panel_size": (3, 2)})
    fig.canvas.draw()
    for ax in axes.ravel():
        box = ax.get_window_extent()
        assert (box.width / fig.dpi, box.height / fig.dpi) == pytest.approx(
            (3, 2), rel=0.01
        )
    assert fig.get_size_inches()[0] >= 6


def test_scales_applied():
    _, axes = _render([{"series": [{"x": X, "y": X}], "xscale": "log"}])
    ax = axes[0, 0]
    assert ax.get_xscale() == "log" and ax.get_yscale() == "linear"


def test_colour_then_style_cycle():
    series = [{"x": X, "y": X * (i + 1)} for i in range(3)]
    _, axes = _render(
        [{"series": series}], palette=["#ff0000", "#00ff00"], style={"cycle": ["-", "--"]}
    )
    lines = _visible_lines(axes[0, 0])
    assert [mcolors.to_hex(ln.get_color()) for ln in lines] == [
        "#ff0000",
        "#00ff00",
        "#ff0000",
    ]
    assert [ln.get_linestyle() for ln in lines] == ["-", "-", "--"]


def test_marker_series_is_hollow_and_sparse():
    _, axes = _render([{"series": [{"x": X, "y": X, "style": "o", "color": "k"}]}])
    (marks,) = _visible_lines(axes[0, 0])
    assert marks.get_marker() == "o" and marks.get_linestyle() == "None"
    assert mcolors.to_hex(marks.get_markerfacecolor()) == "#ffffff"
    assert 2 <= len(marks.get_xdata()) < len(X)


def test_nonfinite_and_nonpositive_dropped_on_log():
    y = np.array([1.0, np.nan, -1.0, 0.0, 2.0, 3.0])
    x = np.arange(1.0, 7.0)
    _, axes = _render([{"series": [{"x": x, "y": y}], "yscale": "log"}])
    (line,) = _visible_lines(axes[0, 0])
    assert list(line.get_ydata()) == [1.0, 2.0, 3.0]


def test_too_few_points_skipped():
    _, axes = _render([{"series": [{"x": [1.0, 2.0], "y": [np.nan, 1.0]}]}])
    assert _visible_lines(axes[0, 0]) == []


def test_bars_with_open_edge_clip_to_axis():
    series = [
        {"x": X, "y": X},
        {
            "x": [5.0, 50.0],
            "y": [2.0, 3.0],
            "kind": "bars",
            "x_lo": [2.0, 20.0],
            "x_hi": [10.0, np.inf],
        },
    ]
    _, axes = _render([{"series": series, "xscale": "log"}])
    ax = axes[0, 0]
    (bars,) = [c for c in ax.collections if isinstance(c, LineCollection)]
    segs = bars.get_segments()
    assert segs[0][:, 0].tolist() == [2.0, 10.0]
    assert np.isfinite(segs[1][:, 0]).all()
    assert 100.0 <= segs[1][1, 0] <= ax.get_xlim()[1]


def test_exp_ticks_label_prefix_and_formatter():
    _, axes = _render([{"series": [{"x": X, "y": X}], "xscale": "log", "xlabel": "T"}])
    ax = axes[0, 0]
    assert ax.get_xlabel() == r"$\log_{10}$ T"
    assert ax.xaxis.get_major_formatter()(100.0) == "2"


def test_exp_ticks_off_in_talk():
    _, axes = _render(
        [{"series": [{"x": X, "y": X}], "xscale": "log", "xlabel": "T"}], style="talk"
    )
    assert axes[0, 0].get_xlabel() == "T"


def test_outer_labels_only():
    # Mixed labels (no figure-wide sup label) fall back to the outer-panel rule.
    def ax(xl, yl):
        return {"series": [{"x": X, "y": X}], "xlabel": xl, "ylabel": yl}

    _, axes = _render(
        [ax("x", "y"), ax("u", "y"), ax("x", "v"), ax("u", "v")], grid=(2, 2)
    )
    # §15: the top row keeps its x label, on the top spine
    assert [a.get_xlabel() for a in axes.ravel()] == ["x", "u", "x", "u"]
    assert [a.xaxis.get_label_position() for a in axes.ravel()] == [
        "top",
        "top",
        "bottom",
        "bottom",
    ]
    assert [a.get_ylabel() for a in axes.ravel()] == ["y", "", "v", ""]


def test_axes_title():
    _, axes = _render([{"series": [], "title": "CO"}])
    assert axes[0, 0].get_title() == "CO"


def test_rcparams_untouched():
    before = mpl.rcParams["font.size"]
    _render([{"series": [{"x": X, "y": X}]}], style={"rc": {"font.size": 31}})
    assert mpl.rcParams["font.size"] == before


def test_unknown_style_name():
    with pytest.raises(ValueError, match="unknown style"):
        _render([{"series": []}], style="nope")


def test_mark_indices_even_spacing():
    disp = np.column_stack([np.arange(0.0, 101.0), np.zeros(101)])
    assert mark_indices(disp, 25.0).tolist() == [0, 25, 50, 75, 100]


@pytest.mark.parametrize(
    ("value", "label"), [(1e-3, "−3"), (1.0, "0"), (1e5, "5"), (3.0, ""), (0.0, "")]
)
def test_exp_label(value, label):
    assert exp_label(value) == label


def _bar_segments(ax):
    (bars,) = [c for c in ax.collections if isinstance(c, LineCollection)]
    return bars.get_segments()


def test_open_bar_edge_clips_to_explicit_xlim():
    line = {"x": np.linspace(0.0, 10.0, 20), "y": np.linspace(1.0, 2.0, 20)}
    bar = {"x": [7.0], "y": [1.5], "kind": "bars", "x_lo": [5.0], "x_hi": [np.inf]}
    _, axes = _render([{"series": [line, bar], "xlim": [0, 100]}])
    ax = axes[0, 0]
    assert _bar_segments(ax)[0][1, 0] == 100.0
    assert ax.get_xlim() == (0.0, 100.0)


def test_bars_do_not_widen_limits():
    line = {"x": np.linspace(0.0, 10.0, 20), "y": np.linspace(1.0, 2.0, 20)}
    bar = {"x": [7.0], "y": [1.5], "kind": "bars", "x_lo": [5.0], "x_hi": [np.inf]}
    _, with_bar = _render([{"series": [line, bar]}])
    _, without = _render([{"series": [line]}])
    assert with_bar[0, 0].get_xlim() == without[0, 0].get_xlim()


def test_bars_only_without_finite_edges_skipped():
    bar = {"x": [1.0], "y": [1.0], "kind": "bars", "x_lo": [-np.inf], "x_hi": [np.inf]}
    _, axes = _render([{"series": [bar]}])
    assert not [c for c in axes[0, 0].collections if isinstance(c, LineCollection)]


def test_open_left_edge_kept_on_log_axis():
    bar = {"x": [3.0], "y": [2.0], "kind": "bars", "x_lo": [-np.inf], "x_hi": [5.0]}
    _, axes = _render([{"series": [{"x": X, "y": X}, bar], "xscale": "log"}])
    ax = axes[0, 0]
    seg = _bar_segments(ax)[0]
    assert seg[0, 0] == ax.get_xlim()[0] and seg[1, 0] == 5.0


def test_outer_labels_three_columns():
    ax = {"series": [{"x": X, "y": X}], "ylabel": "y"}
    _, axes = _render([ax] * 3, grid=(1, 3))
    assert [a.get_ylabel() for a in axes.ravel()] == ["y", "", ""]


def _bar_line_panels(order):
    short = {
        "series": [
            {"x": np.linspace(0.0, 10.0, 20), "y": np.linspace(1.0, 2.0, 20)},
            {"x": [7.0], "y": [1.5], "kind": "bars", "x_lo": [5.0], "x_hi": [np.inf]},
        ]
    }
    wide = {
        "series": [{"x": np.linspace(0.0, 100.0, 20), "y": np.linspace(1.0, 2.0, 20)}]
    }
    panels = [short, wide] if order == "bars_first" else [wide, short]
    return _render(panels, grid=(2, 1), sharex=True)


@pytest.mark.parametrize("order", ["bars_first", "bars_last"])
def test_sharex_bars_do_not_clip_other_panels(order):
    _, axes = _bar_line_panels(order)
    assert axes[0, 0].get_xlim() == axes[1, 0].get_xlim()
    assert axes[0, 0].get_xlim()[1] >= 100.0
    bars_ax = axes[0, 0] if order == "bars_first" else axes[1, 0]
    assert _bar_segments(bars_ax)[0][1, 0] == bars_ax.get_xlim()[1]


def test_sharex_limits_independent_of_panel_order():
    _, a = _bar_line_panels("bars_first")
    _, b = _bar_line_panels("bars_last")
    assert a[0, 0].get_xlim() == b[0, 0].get_xlim()


def _tick_texts(axis):
    return [t.get_text() for t in axis.get_ticklabels() if t.get_text()]


def test_narrow_log_axis_gets_fractional_exponent_ticks():
    x = np.linspace(11.0, 14.0, 30)
    fig, axes = _render([{"series": [{"x": x, "y": x}], "xscale": "log"}])
    fig.canvas.draw()
    ax = axes[0, 0]
    labels = _tick_texts(ax.xaxis)
    assert len(labels) >= 2
    assert all(float(t.replace("−", "-")) for t in labels)
    assert len(set(labels)) == len(labels)
    assert not len(ax.xaxis.get_minorticklocs())


def test_narrow_negative_log_axis_uses_true_minus():
    x = np.linspace(0.011, 0.014, 30)
    fig, axes = _render([{"series": [{"x": x, "y": x}], "xscale": "log"}])
    fig.canvas.draw()
    labels = _tick_texts(axes[0, 0].xaxis)
    assert len(labels) >= 2 and all(t.startswith("−") for t in labels)


def test_wide_log_axis_keeps_decade_labels():
    x = np.logspace(0, 4, 30)
    fig, axes = _render([{"series": [{"x": x, "y": x}], "xscale": "log"}])
    fig.canvas.draw()
    labels = _tick_texts(axes[0, 0].xaxis)
    assert set(labels) <= {"−1", "0", "1", "2", "3", "4", "5"} and "2" in labels


def test_narrow_axis_untouched_without_exp_ticks():
    x = np.linspace(11.0, 14.0, 30)
    _, axes = _render(
        [{"series": [{"x": x, "y": x}], "xscale": "log"}], style={"exp_ticks": False}
    )
    assert type(axes[0, 0].xaxis.get_major_locator()).__name__ == "LogLocator"


@pytest.mark.parametrize(
    ("ticks", "labels"),
    [
        ([1.05, 1.1, 1.15], ["1.05", "1.1", "1.15"]),
        ([-2.5, -2.0, -1.5], ["−2.5", "−2", "−1.5"]),
        ([0.0, 0.2, 0.4], ["0", "0.2", "0.4"]),
    ],
)
def test_fractional_exponent_labels(ticks, labels):
    from jaff.plotting._ticks import fractional_labels

    assert fractional_labels(np.array(ticks)) == labels


def _sup(fig, which):
    artist = getattr(fig, f"_sup{which}label")
    return artist.get_text() if artist is not None else None


def test_shared_labels_2x2_use_sup_labels():
    ax = {"series": [{"x": X, "y": X}], "xlabel": "x", "ylabel": "y"}
    fig, axes = _render([ax] * 4, grid=(2, 2))
    assert _sup(fig, "y") == "y" and _sup(fig, "x") == "x"
    assert [a.get_ylabel() for a in axes.ravel()] == [""] * 4
    assert [a.get_xlabel() for a in axes.ravel()] == [""] * 4


def test_shared_labels_one_row_only_supxlabel():
    ax = {"series": [{"x": X, "y": X}], "xlabel": "x", "ylabel": "y"}
    fig, axes = _render([ax] * 2, grid=(1, 2))
    assert _sup(fig, "y") is None and _sup(fig, "x") == "x"
    assert [a.get_ylabel() for a in axes.ravel()] == ["y", ""]
    assert [a.get_xlabel() for a in axes.ravel()] == ["", ""]


def test_shared_labels_fontsize_matches_axes_labels():
    ax = {"series": [{"x": X, "y": X}], "xlabel": "x", "ylabel": "y"}
    fig, _ = _render([ax] * 4, grid=(2, 2))
    assert fig._supylabel.get_fontsize() == pytest.approx(11)


def test_mixed_labels_no_sup_labels():
    a = {"series": [{"x": X, "y": X}], "xlabel": "x", "ylabel": "y"}
    b = {"series": [{"x": X, "y": X}], "xlabel": "u", "ylabel": "v"}
    fig, axes = _render([a, b, a, b], grid=(2, 2))
    assert _sup(fig, "y") is None and _sup(fig, "x") is None
    assert [a.get_ylabel() for a in axes.ravel()] == ["y", "v", "y", "v"]
    # §15: the top row keeps its x label, on the top spine
    assert [a.get_xlabel() for a in axes.ravel()] == ["x", "u", "x", "u"]
    assert [a.xaxis.get_label_position() for a in axes.ravel()] == [
        "top",
        "top",
        "bottom",
        "bottom",
    ]


def test_sup_xlabel_keeps_sharex_tick_labels():
    ax = {"series": [{"x": X, "y": X}], "xlabel": "x"}
    fig, axes = _render([ax] * 3, grid=(2, 2), sharex=True)
    top_right = axes[0, 1]
    assert _sup(fig, "x") == "x"
    assert top_right.xaxis.get_tick_params()["labelbottom"]


def test_nearly_constant_log_axis_labels_distinct():
    y = 1e5 * (1 + 1e-8 * np.linspace(0.0, 1.0, 50))
    fig, axes = _render([{"series": [{"x": X, "y": y}], "yscale": "log"}])
    fig.canvas.draw()
    labels = _tick_texts(axes[0, 0].yaxis)
    assert len(set(labels)) == len(labels)


def test_render_function_does_not_shadow_its_module(monkeypatch):
    """jaff.plotting.render is the function; the module is jaff.plotting._render."""
    import importlib
    import types

    import jaff.plotting

    module = importlib.import_module("jaff.plotting._render")

    assert isinstance(module, types.ModuleType)
    assert jaff.plotting.render is module.render
    # Dotted-path targets resolve through attributes, which shadowing broke.
    monkeypatch.setattr("jaff.plotting._render.LOG_MAX_DECADES", 10.0)
    assert module.LOG_MAX_DECADES == 10.0
