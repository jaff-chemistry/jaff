# ABOUTME: Tests for jaff.plotting.adapters — reactions to series dicts
# ABOUTME: Rate evaluation over T, skips, xsec processes/units, band bars

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import sympy
from matplotlib.collections import LineCollection

from jaff.plotting.adapters import normalize_processes, rate_series, xsec_series
from tests.xsec_props import make_xsec_props, stub_photo

TGAS = sympy.Symbol("tgas")


def _reaction(rate, tmin=None, tmax=None):
    return SimpleNamespace(
        rate_coefficient=rate, tmin=tmin, tmax=tmax, get_latex=lambda: "$A \\to B$"
    )


def _photo(processes=("photodecay",), band_xsecs=None):
    energy = np.linspace(10.0, 20.0, 11)
    xs = make_xsec_props(
        photon_energy=energy, **{p: np.full(11, 1e-18) for p in processes}
    )
    return SimpleNamespace(
        photo=SimpleNamespace(
            xsec=xs,
            band_xsecs=band_xsecs
            if band_xsecs is not None
            else pd.DataFrame(columns=["lower", "upper", "eavg", "xsec", "xsec_frac"]),
        ),
        get_latex=lambda: "$C \\to C^+$",
    )


class TestRateSeries:
    def test_reaction_uses_own_bounds_and_latex(self):
        s = rate_series(_reaction(1e-10 * (TGAS / 300) ** 0.5, 10, 1000), npoints=5)
        assert s["x"][0] == pytest.approx(10) and s["x"][-1] == pytest.approx(1000)
        assert s["y"][-1] == pytest.approx(1e-10 * (1000 / 300) ** 0.5)
        assert s["label"] == "$A \\to B$"

    def test_defaults_when_unbounded(self):
        s = rate_series(_reaction(sympy.Float(1e-9)), npoints=3)
        assert s["x"][0] == pytest.approx(2.73) and s["x"][-1] == pytest.approx(1e6)

    def test_caller_bounds_override(self):
        s = rate_series(_reaction(TGAS, 10, 1000), tmin=1, tmax=2, npoints=2)
        assert s["x"].tolist() == pytest.approx([1, 2])

    def test_unevaluable_returns_none(self):
        assert rate_series(_reaction(TGAS * sympy.Symbol("ntot"))) is None

    def test_bare_expression_needs_bounds(self):
        with pytest.raises(ValueError, match="tmin and tmax"):
            rate_series(TGAS**2)

    def test_bare_expression(self):
        s = rate_series(TGAS**2, tmin=1, tmax=10, npoints=2)
        assert s["y"].tolist() == pytest.approx([1, 100]) and s["label"] == "tgas**2"

    def test_xy_tuple(self):
        s = rate_series(([1, 2], [3, 4]), label="data")
        assert s["y"].tolist() == [3, 4] and s["label"] == "data"

    def test_bad_type(self):
        with pytest.raises(TypeError, match="Cannot plot rate input"):
            rate_series(42)


class TestRateRobustness:
    @pytest.mark.parametrize("tmin,tmax", [(0, 10), (100, 10), (-1, 5)])
    def test_bad_bounds(self, tmin, tmax):
        with pytest.raises(ValueError, match="0 < tmin < tmax"):
            rate_series(_reaction(TGAS), tmin=tmin, tmax=tmax)

    def test_none_rate(self):
        assert rate_series(_reaction(None, 1, 10)) is None

    def test_zoo_rate(self):
        assert rate_series(_reaction(sympy.zoo * TGAS, 1, 10)) is None


class TestXsecSeries:
    def test_one_series_per_process_with_data(self):
        out = xsec_series(_photo(("photodecay",)))
        assert [s["label"] for s in out] == ["Photodecay"]

    def test_units(self):
        (s,) = xsec_series(_photo(), energy_unit="eV", xsec_unit="Mb")
        assert s["y"][0] == pytest.approx(1.0)  # 1e-18 cm^2 = 1 Mb

    def test_prefix(self):
        (s,) = xsec_series(_photo(), prefix=True)
        assert s["label"] == "$C \\to C^+$ — Photodecay"

    def test_no_xsecs(self):
        assert xsec_series(SimpleNamespace(photo=None)) == []

    def test_empty_band_xsecs_adds_nothing(self):
        out = xsec_series(_photo(), bands=True)
        assert [s.get("kind", "line") for s in out] == ["line"]

    def test_band_bars(self):
        bands = pd.DataFrame(
            {
                "lower": [11.2, 13.6],
                "upper": [13.6, np.inf],
                "eavg": [12.0, 20.0],
                "xsec": [2e-18, np.nan],
                "xsec_frac": [1.0, 0.0],
            }
        )
        out = xsec_series(_photo(band_xsecs=bands), bands=True)
        bar = out[-1]
        assert bar["kind"] == "bars" and bar["label"] == "Band average"
        assert bar["y"].tolist() == pytest.approx([2.0])
        assert bar["x_lo"].tolist() == [11.2] and bar["x_hi"].tolist() == [13.6]

    def test_invalid_process(self):
        with pytest.raises(KeyError, match="Invalid cross-section"):
            normalize_processes("photoionisation")

    def test_process_list(self):
        assert normalize_processes(None) == ["photo_absorption", "photodecay"]
        assert normalize_processes("photodecay") == ["photodecay"]


class TestXsecExtras:
    @staticmethod
    def _bands():
        return pd.DataFrame(
            {
                "lower": [11.2, 13.6],
                "upper": [13.6, np.inf],
                "eavg": [12.0, 20.0],
                "xsec": [2e-18, 3e-18],
                "xsec_frac": [1.0, 0.0],
            }
        )

    def test_band_bars_nm_open_edge(self):
        out = xsec_series(_photo(band_xsecs=self._bands()), bands=True, energy_unit="nm")
        bar = out[-1]
        assert np.all(bar["x_lo"] <= bar["x_hi"])
        assert bar["x_lo"][1] == -np.inf
        assert bar["x_hi"][1] == pytest.approx(91.16, abs=0.05)
        assert bar["x_lo"][0] < bar["x_hi"][0]

    def test_band_bars_ev_open_edge_stays_high(self):
        bar = xsec_series(_photo(band_xsecs=self._bands()), bands=True)[-1]
        assert bar["x_hi"][1] == np.inf and bar["x_lo"][1] == pytest.approx(13.6)

    def test_band_bars_render_nm(self):
        from jaff.plotting.render import render
        from jaff.plotting.spec import FigureSpec

        bar = xsec_series(_photo(band_xsecs=self._bands()), bands=True, energy_unit="nm")[
            -1
        ]
        line = {"x": [80.0, 130.0], "y": [1.0, 2.0]}
        spec = FigureSpec.from_dict({"axes": [{"series": [line, bar], "xscale": "log"}]})
        _, axes = render(spec)
        ax = axes[0, 0]
        (bars,) = [c for c in ax.collections if isinstance(c, LineCollection)]
        segs = bars.get_segments()
        xlim = ax.get_xlim()
        assert segs[0][:, 0].min() == pytest.approx(91.16, abs=0.05)
        assert segs[0][:, 0].max() == pytest.approx(110.7, abs=0.1)
        assert segs[1][:, 0].min() == pytest.approx(xlim[0])
        assert segs[1][:, 0].max() == pytest.approx(91.16, abs=0.05)

    def test_band_bars_ev_edges(self):
        bar = xsec_series(_photo(band_xsecs=self._bands()), bands=True)[-1]
        assert bar["x_lo"].tolist() == pytest.approx([11.2, 13.6])
        assert bar["x_hi"][0] == pytest.approx(13.6)

    def test_band_prefix(self):
        out = xsec_series(_photo(band_xsecs=self._bands()), bands=True, prefix=True)
        assert out[-1]["label"] == "$C \\to C^+$ — Band average"

    def test_verner_sampling_nm(self):
        E = sympy.Symbol("ph_E")
        expr = sympy.Piecewise((0, E < 10), (1e-18 * (10 / E) ** 3, E < 100), (0, True))
        r = SimpleNamespace(
            photo=stub_photo(
                xsec=make_xsec_props(photodecay_expr=expr, database="verner")
            ),
            get_latex=lambda: "x",
        )
        (s,) = xsec_series(r, energy_unit="nm")
        assert len(s["x"]) > 10 and np.all(np.isfinite(s["y"]))
        assert s["x"].max() > s["x"].min()
        assert s["x"].max() == pytest.approx(123.98, abs=0.05)

    def test_missing_process_returns_empty(self):
        assert xsec_series(_photo(("photodecay",)), processes=["photo_absorption"]) == []

    def test_no_energy_returns_empty(self):
        r = SimpleNamespace(photo=stub_photo(xsec=make_xsec_props()))
        assert xsec_series(r) == []


class TestReactionBoundsSanitised:
    @pytest.mark.parametrize("bad", [0.0, -5.0, None, float("nan"), float("inf")])
    def test_bad_reaction_tmin_falls_back(self, bad):
        s = rate_series(_reaction(TGAS, bad, 1000), npoints=3)
        assert s["x"][0] == pytest.approx(2.73)

    @pytest.mark.parametrize("bad", [0.0, -1.0, float("nan"), float("inf")])
    def test_bad_reaction_tmax_falls_back(self, bad):
        s = rate_series(_reaction(TGAS, 10, bad), npoints=3)
        assert s["x"][-1] == pytest.approx(1e6)

    def test_caller_tmin_zero_still_raises(self):
        with pytest.raises(ValueError, match="0 < tmin < tmax"):
            rate_series(_reaction(TGAS, 0.0, 1000), tmin=0, tmax=10)
