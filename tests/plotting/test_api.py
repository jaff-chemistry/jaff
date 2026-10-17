# ABOUTME: Tests for jaff.plotting.api builders and Reaction/Reactions wrappers
# ABOUTME: Nesting -> panels, per-panel scales, save by suffix, None on nothing

from types import SimpleNamespace

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest
import sympy

from jaff.plotting import plot, plot_rates, plot_xsecs

TGAS = sympy.Symbol("tgas")


def _r(k):
    return SimpleNamespace(
        rate=k * (TGAS / 300) ** 0.5, tmin=10, tmax=1000, get_latex=lambda: f"$R_{k}$"
    )


def _photo(name):
    energy = np.linspace(10.0, 20.0, 11)
    return SimpleNamespace(
        xsecs_dict={
            "photon_energy": energy,
            "photodecay": np.full(11, 1e-18),
            "photo_absorption": np.full(11, 2e-18),
        },
        get_latex=lambda: name,
    )


def _labels(ax):
    """Entries keying *ax*: its own legend group, else the figure's shared legend."""
    legend = ax.get_legend()
    if legend is None:
        fig = ax.figure
        hosted = [a.get_legend() for a in fig.axes if not a.axison and a.get_legend()]
        shared = list(fig.legends) + hosted
        legend = shared[0] if shared else None
    return [t.get_text() for t in legend.get_texts()] if legend else []


class TestPlot:
    def test_returns_fig_and_axes(self):
        fig, axes = plot({"axes": [{"series": [{"x": [1, 2], "y": [1, 2]}]}]}, show=False)
        assert axes.shape == (1, 1)

    def test_no_series_returns_none(self):
        assert plot({"axes": [{"series": []}]}, show=False) is None

    @pytest.mark.parametrize("suffix", [".png", ".pdf", ".svg"])
    def test_save_by_suffix(self, tmp_path, suffix):
        path = tmp_path / f"out{suffix}"
        plot({"axes": [{"series": [{"x": [1, 2], "y": [1, 2]}]}]}, save=path, show=False)
        assert path.stat().st_size > 0


class TestPlotRates:
    def test_no_grid_overlays_everything(self):
        fig, axes = plot_rates([_r(1), [_r(2), _r(3)]], show=False)
        assert axes.shape == (1, 1)
        assert _labels(axes[0, 0]) == ["$R_1$", "$R_2$", "$R_3$"]

    def test_nesting_with_grid(self):
        _, axes = plot_rates([_r(1), [_r(2), _r(3)]], grid=(1, 2), show=False)
        left, right = axes.ravel()
        # one shared legend keys the figure; the lone curve is named by its title
        assert left.get_title() == "$R_1$"
        assert left.get_legend() is None and right.get_legend() is None
        assert _labels(right) == ["$R_2$", "$R_3$"]

    def test_per_panel_scales(self):
        _, axes = plot_rates(
            [_r(1), _r(2)], grid=(1, 2), yscale=["log", "linear"], show=False
        )
        assert [a.get_yscale() for a in axes.ravel()] == ["log", "linear"]

    def test_scale_list_length_checked(self):
        with pytest.raises(ValueError, match="yscale: got 1 values for 2 panels"):
            plot_rates([_r(1), _r(2)], grid=(1, 2), yscale=["log"], show=False)

    def test_labels_override(self):
        _, axes = plot_rates([_r(1), _r(2)], labels=["a", "b"], show=False)
        assert _labels(axes[0, 0]) == ["a", "b"]

    def test_unevaluable_panel_dropped(self):
        bad = SimpleNamespace(rate=sympy.Symbol("ntot"), tmin=None, tmax=None)
        _, axes = plot_rates([bad, _r(1)], grid=(1, 2), show=False)
        assert [a.axison for a in axes.ravel()] == [True, False]

    def test_nothing_evaluable_returns_none(self):
        bad = SimpleNamespace(rate=sympy.Symbol("ntot"), tmin=None, tmax=None)
        assert plot_rates(bad, show=False) is None

    def test_figure_title_not_leaked_from_panels(self):
        fig, _ = plot_rates([_r(1)], grid=(1, 1), show=False)
        assert fig._suptitle is None


class TestPlotXsecs:
    def test_single_reaction_title_and_process_legend(self):
        _, axes = plot_xsecs(_photo("$C$"), show=False)
        ax = axes[0, 0]
        assert ax.get_title() == "$C$"
        assert _labels(ax) == ["Photoabsorption", "Photodecay"]

    def test_several_reactions_prefixed(self):
        _, axes = plot_xsecs(
            [_photo("$C$"), _photo("$O$")], processes="photodecay", show=False
        )
        assert _labels(axes[0, 0]) == ["$C$ — Photodecay", "$O$ — Photodecay"]

    def test_grid(self):
        _, axes = plot_xsecs([_photo("$C$"), _photo("$O$")], grid=(2, 1), show=False)
        assert [a.get_title() for a in axes.ravel()] == ["$C$", "$O$"]

    def test_nothing_returns_none(self):
        assert plot_xsecs(SimpleNamespace(xsecs_dict=None), show=False) is None


class TestReactionWrappers:
    def test_reaction_plot_rate_coefficient(self, make_network):
        net = make_network(["H + H -> H2 [10,1000] 1e-10*(tgas/300)**0.5"])
        fig, axes = net.reactions[0].plot_rate_coefficient(show=False)
        assert axes[0, 0].get_title() == net.reactions[0].get_latex()

    def test_reactions_plot_rates(self, make_network):
        net = make_network(["H + H -> H2 [10,1000] 1e-10", "H + O -> OH [10,1000] 2e-10"])
        _, axes = net.reactions.plot_rates(grid=(1, 2), show=False)
        assert axes.shape == (1, 2)


class TestSvgText:
    def test_svg_keeps_text_as_text(self, tmp_path):
        path = tmp_path / "out.svg"
        plot(
            {"axes": [{"title": "hello", "series": [{"x": [1, 2], "y": [1, 2]}]}]},
            save=path,
            show=False,
        )
        assert "<text" in path.read_text()


NETWORKS = Path(__file__).resolve().parents[2] / "networks"


@pytest.fixture(scope="module")
def uclchem():
    from jaff import Network

    return Network(str(NETWORKS / "uclchem_small_gas/uclchem_small_gas_network.jet"))


class TestRealNetworks:
    def test_uclchem_unbounded_tmin_plots(self, uclchem):
        # The first catalogue entries have no tmin (UCLCHEM's placeholder 0) and
        # k = 0 (desorption): sampled fine, but nothing survives on a log y axis.
        zeros = list(uclchem.reactions)[:5]
        assert all(r.tmin is None for r in zeros)
        assert plot_rates(zeros, yscale="linear", show=False) is not None
        assert zeros[0].plot_rate_coefficient(yscale="linear", show=False) is not None
        assert plot_rates(zeros, show=False) is None

    def test_uclchem_catalogue_has_no_blank_panels(self, uclchem):
        fig, axes = uclchem.reactions.plot_rates(show=False)
        used = [ax for ax in axes.ravel() if ax.axison]
        assert len(used) > 1
        assert all(any(ln.get_visible() for ln in ax.lines) for ax in used)

    def test_catalogue_is_a_group(self, make_network):
        net = make_network(["H + H -> H2 [10,1000] 1e-10", "H + O -> OH [10,1000] 2e-10"])
        _, axes = plot_rates(net.reactions, grid=(1, 2), show=False)
        assert axes.shape == (1, 2)


class TestGrouping:
    def test_tuple_of_reactions_overlays(self):
        _, axes = plot_rates((_r(1), _r(2)), show=False)
        assert _labels(axes[0, 0]) == ["$R_1$", "$R_2$"]

    def test_xy_pair_is_a_leaf(self):
        _, axes = plot_rates(([1, 2], [3, 4]), show=False)
        assert len(axes[0, 0].lines) == 1

    def test_ndarray_and_generator_groups(self):
        _, axes = plot_rates(np.array([_r(1), _r(2)], dtype=object), show=False)
        assert len(axes[0, 0].lines) == 2
        _, axes = plot_rates((_r(k) for k in (1, 2, 3)), show=False)
        assert len(axes[0, 0].lines) == 3

    def test_xsec_tuple_is_group(self):
        _, axes = plot_xsecs(
            (_photo("$C$"), _photo("$O$")), processes="photodecay", show=False
        )
        assert _labels(axes[0, 0]) == ["$C$ — Photodecay", "$O$ — Photodecay"]


class TestMisc:
    def test_show_under_agg_no_warning(self):
        with warnings.catch_warnings(record=True) as rec:
            warnings.simplefilter("always")
            plot({"axes": [{"series": [{"x": [1, 2], "y": [1, 2]}]}]}, show=True)
        assert not [w for w in rec if issubclass(w.category, UserWarning)]

    def test_too_many_panels_for_grid(self):
        with pytest.raises(ValueError, match="grid .* has 1 cells but 2 panels"):
            plot_rates([_r(1), _r(2)], grid=(1, 1), show=False)

    def test_prefix_counts_contributing_reactions(self):
        none = SimpleNamespace(xsecs_dict=None, get_latex=lambda: "N")
        _, axes = plot_xsecs([none, _photo("$C$")], show=False)
        ax = axes[0, 0]
        assert ax.get_title() == "$C$"
        assert _labels(ax) == ["Photoabsorption", "Photodecay"]

    def test_labels_align_across_dropped_panel(self):
        bad = SimpleNamespace(rate=sympy.Symbol("ntot"), tmin=None, tmax=None)
        _, axes = plot_rates(
            [[_r(1)], [bad], [_r(2)]], labels=["a", "b", "c"], grid=(1, 3), show=False
        )
        flat = axes.ravel()
        assert [a.get_title() for a in flat[:2]] == ["a", "c"]
        assert not flat[2].axison

    def test_save_png_uses_style_dpi(self, tmp_path):
        from PIL import Image

        path = tmp_path / "o.png"
        plot({"axes": [{"series": [{"x": [1, 2], "y": [1, 2]}]}]}, save=path, show=False)
        with Image.open(path) as img:
            assert img.info.get("dpi") is not None


def test_tuple_of_xy_pairs_overlays():
    x = np.linspace(1.0, 2.0, 5)
    _, axes = plot_rates(((x, x), (x, 2 * x)), labels=["a", "b"], show=False)
    assert _labels(axes[0, 0]) == ["a", "b"]


class TestAutoGrid:
    @pytest.mark.parametrize(
        ("n", "grid"),
        [(1, (1, 1)), (2, (1, 2)), (3, (1, 3)), (4, (2, 2)), (5, (2, 3)), (48, (7, 7))],
    )
    def test_auto_grid(self, n, grid):
        from jaff.plotting.api import _auto_grid

        assert _auto_grid(n) == grid


class TestLegendCapChunking:
    def test_rates_chunked_into_panels(self):
        _, axes = plot_rates([_r(k) for k in range(1, 31)], show=False)
        assert axes.shape == (1, 3)
        counts = [len(_labels(ax)) for ax in axes.ravel()]
        assert counts == [12, 12, 6]
        assert _labels(axes[0, 0])[0] == "$R_1$" and _labels(axes[0, 1])[0] == "$R_13$"

    def test_style_legend_max_raises_cap(self):
        _, axes = plot_rates(
            [_r(k) for k in range(1, 14)], style={"legend_max": 20}, show=False
        )
        assert axes.shape == (1, 1) and len(_labels(axes[0, 0])) == 13

    def test_explicit_grid_not_chunked(self, jaff_warnings):
        _, axes = plot_rates([[_r(k) for k in range(1, 14)]], grid=(1, 1), show=False)
        assert axes.shape == (1, 1) and axes[0, 0].get_legend() is None
        assert any("legend" in m for m in jaff_warnings)

    def test_unevaluable_not_counted(self):
        bad = SimpleNamespace(rate=sympy.Symbol("ntot"), tmin=None, tmax=None)
        _, axes = plot_rates([bad] + [_r(k) for k in range(1, 13)], show=False)
        assert axes.shape == (1, 1) and len(_labels(axes[0, 0])) == 12

    def test_lone_curve_chunk_titled(self):
        _, axes = plot_rates([_r(k) for k in range(1, 14)], show=False)
        assert axes.shape == (1, 2) and axes[0, 1].get_title() == "$R_13$"

    def test_scale_list_without_grid_rejected_when_chunked(self):
        with pytest.raises(ValueError, match="yscale.*grid"):
            plot_rates([_r(k) for k in range(1, 14)], yscale=["log"], show=False)

    def test_scale_scalar_applies_to_every_chunk(self):
        _, axes = plot_rates([_r(k) for k in range(1, 14)], yscale="linear", show=False)
        assert [a.get_yscale() for a in axes.ravel()] == ["linear", "linear"]

    def test_xsecs_chunked_keep_reactions_together(self):
        photos = [_photo(f"$P_{k}$") for k in range(7)]
        _, axes = plot_xsecs(photos, show=False)
        # 7 reactions x 2 processes = 14 curves -> 6 reactions (12) + 1 reaction (2).
        assert axes.shape == (1, 2)
        assert len(_labels(axes[0, 0])) == 12
        assert axes[0, 1].get_title() == "$P_6$"
        assert _labels(axes[0, 1]) == ["Photoabsorption", "Photodecay"]


class TestClose:
    def test_plot_closes_figure_by_default(self):
        fig, _ = plot({"axes": [{"series": [{"x": [1, 2], "y": [1, 2]}]}]}, show=False)
        assert fig.number not in plt.get_fignums()

    def test_close_false_keeps_figure(self):
        fig, _ = plot(
            {"axes": [{"series": [{"x": [1, 2], "y": [1, 2]}]}]}, show=False, close=False
        )
        assert fig.number in plt.get_fignums()

    def test_closed_figure_still_saves(self, tmp_path):
        fig, axes = plot({"axes": [{"series": [{"x": [1, 2], "y": [1, 2]}]}]}, show=False)
        axes[0, 0].set_title("edited")
        fig.savefig(tmp_path / "after.png")
        assert (tmp_path / "after.png").stat().st_size > 0

    @pytest.mark.parametrize("close", [True, False])
    def test_builders_forward_close(self, close):
        fig, _ = plot_rates([_r(1)], show=False, close=close)
        assert (fig.number in plt.get_fignums()) is not close
        fig, _ = plot_xsecs(_photo("$C$"), show=False, close=close)
        assert (fig.number in plt.get_fignums()) is not close


def _zero(k):
    return SimpleNamespace(
        rate=0 * TGAS, tmin=10, tmax=1000, get_latex=lambda: f"$Z_{k}$"
    )


class TestUnplottableSeries:
    def test_zero_rates_dropped_before_chunking(self, jaff_warnings):
        items = [_zero(k) for k in range(5)] + [_r(k) for k in range(1, 31)]
        _, axes = plot_rates(items, show=False)
        assert axes.shape == (1, 3)
        assert [len(_labels(ax)) for ax in axes.ravel()] == [12, 12, 6]
        assert all(any(ln.get_visible() for ln in ax.lines) for ax in axes.ravel())
        summary = [m for m in jaff_warnings if "no plottable points" in m]
        assert len(summary) == 1 and "5 series" in summary[0]
        assert not [m for m in jaff_warnings if "plottable point(s)" in m]

    def test_zero_rates_kept_on_linear_axis(self):
        _, axes = plot_rates([_zero(1), _r(1)], yscale="linear", show=False)
        assert len(_labels(axes[0, 0])) == 2

    def test_all_unplottable_returns_none(self):
        assert plot_rates([_zero(1), _zero(2)], show=False) is None

    def test_zero_xsec_reaction_dropped_before_chunking(self):
        zero = _photo("$Z$")
        zero.xsecs_dict = {
            k: v * 0 if k != "photon_energy" else v for k, v in zero.xsecs_dict.items()
        }
        photos = [zero] + [_photo(f"$P_{k}$") for k in range(6)]
        _, axes = plot_xsecs(photos, show=False)
        assert axes.shape == (1, 1) and len(_labels(axes[0, 0])) == 12


class TestInteractiveClose:
    @pytest.fixture
    def interactive(self, monkeypatch):
        import jaff.plotting.api as api

        monkeypatch.setattr(api.matplotlib, "get_backend", lambda: "qtagg")
        monkeypatch.setattr(api.plt, "show", lambda *a, **k: None)

    @pytest.mark.parametrize(("ion", "kept"), [(True, True), (False, False)])
    def test_close_skipped_after_nonblocking_show(
        self, interactive, monkeypatch, ion, kept
    ):
        import jaff.plotting.api as api

        monkeypatch.setattr(api.plt, "isinteractive", lambda: ion)
        fig, _ = plot({"axes": [{"series": [{"x": [1, 2], "y": [1, 2]}]}]}, show=True)
        assert (fig.number in plt.get_fignums()) is kept

    def test_close_when_not_shown_even_if_interactive(self, interactive, monkeypatch):
        import jaff.plotting.api as api

        monkeypatch.setattr(api.plt, "isinteractive", lambda: True)
        fig, _ = plot({"axes": [{"series": [{"x": [1, 2], "y": [1, 2]}]}]}, show=False)
        assert fig.number not in plt.get_fignums()


class TestReviewMinor:
    def test_pack_heavy_items_own_runs(self):
        from jaff.plotting.api import _pack

        assert _pack([3, 20, 2], 12) == [[0], [1], [2]]

    def test_figure_legend_falls_back_when_chunked(self, jaff_warnings):
        fig, axes = plot_rates([_r(k) for k in range(1, 14)], legend="figure", show=False)
        assert not fig.legends
        assert len(_labels(axes[0, 0])) == 12
        assert sum("legend" in m and "figure" in m for m in jaff_warnings) == 1

    @pytest.mark.parametrize("builder", [plot_rates, plot_xsecs])
    def test_empty_input_no_panel_warning(self, builder, jaff_warnings):
        assert builder([], show=False) is None
        assert not jaff_warnings


class TestWarningSummary:
    def test_unevaluable_rates_summarised(self, jaff_warnings):
        bad = [SimpleNamespace(rate=sympy.Symbol("ntot"), tmin=None, tmax=None)] * 3
        _, axes = plot_rates(bad + [_r(1), _r(2)], show=False)
        assert len(_labels(axes[0, 0])) == 2
        hits = [m for m in jaff_warnings if "cannot be evaluated" in m]
        assert hits == [
            "plot_rates: 3 rates cannot be evaluated vs temperature; skipped."
        ]
        assert not [m for m in jaff_warnings if "Cannot evaluate" in m]


def test_grid_with_unevaluable_panels_warns_once_per_cause(jaff_warnings):
    bad = [SimpleNamespace(rate=sympy.Symbol("ntot"), tmin=None, tmax=None)] * 5
    plot_rates([*bad, _r(1), _r(2)], grid=(3, 3), show=False)
    assert not any("has nothing to plot" in m for m in jaff_warnings)
    assert sum("cannot be evaluated" in m for m in jaff_warnings) == 1


def test_reactions_plot_xsecs_grid_no_per_panel_spam(jaff_warnings):
    plot_xsecs(
        [SimpleNamespace(xsecs_dict=None)] * 4 + [_photo("$C$")], grid=(3, 3), show=False
    )
    assert not any("has nothing to plot" in m for m in jaff_warnings)
