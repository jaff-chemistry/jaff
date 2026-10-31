# ABOUTME: photoCrossSection<N> aux function: a custom σ(ph_E) replaces the database
# ABOUTME: cross section unless overridden (global + per-reaction override_custom_cross_section)

from types import SimpleNamespace

import pytest
import sympy as sp

from tests.jaffgen_run import run_jaffgen
from jaff.errors import ParserError
from jaff.physics.photo_reactions.config import RadiationConfig

BANDS = [6.0, 13.6, 100.0]
PHOTO = "C -> C+ + E  []  PHOTO, 11.3\n"
KEY = "C._PHOTON__C+.e-"
XSEC = "@function photoCrossSection0()\n@return 1e-18\n"


def _photo(net):
    return next(r for r in net.reactions if r.type == "photo")


def _toml(tmp_path, body):
    cfg = tmp_path / "jaff.toml"
    cfg.write_text(f'[network.reactions.reaction."{KEY}".photo]\n{body}\n')

    return str(cfg)


@pytest.fixture
def xsec_network(make_network, tmp_path):
    """Network from *lines* with *jfunc* as its aux file; kwargs go to ``Network``."""

    def _make(lines=PHOTO, jfunc=XSEC, **kwargs):
        (tmp_path / "p.jfunc").write_text(jfunc)
        return make_network(lines, name="p.jet", **kwargs)

    return _make


class TestFunctionConflicts:
    def test_jfunc_with_both_raises(self, xsec_network):
        jfunc = XSEC + "@function chemRateCoeff0()\n@return 1e-10\n"
        with pytest.raises(ParserError, match="photoCrossSection"):
            xsec_network(jfunc=jfunc)


class TestCustomCrossSection:
    def test_replaces_database_xsec(self, xsec_network):
        rxn = _photo(xsec_network())

        assert rxn.function("photoCrossSection") is not None
        assert rxn.photo.xsec.database == "custom"
        assert float(rxn.photo.xsec.photodecay_expr) == pytest.approx(
            1e-18, rel=1e-6, abs=0
        )
        assert rxn.photo.xsec.photon_energy is None
        assert rxn.photo.xsec.process.absorption is False
        assert rxn.photo.xsec.process.decay_type == "ionization"

    def test_dissociation_decay_type(self, xsec_network):
        rxn = _photo(xsec_network("CO -> C + O  []  PHOTO, 11\n"))
        assert rxn.photo.xsec.process.decay_type == "dissociation"

    def test_band_and_total_averages_use_it(self, xsec_network):
        props = RadiationConfig(bands=BANDS)
        rxn = _photo(xsec_network(radiation=props))

        assert float(rxn.photo.xsec_total) == pytest.approx(1e-18, rel=1e-6, abs=0)
        bands = [float(band.xsec) for band in rxn.photo.bands]
        assert bands == pytest.approx([1e-18, 1e-18], rel=1e-6, abs=0)
        assert rxn.rate_coefficient.has(rxn.net.radiation.den[0])

    def test_energy_dependent_xsec_uses_photon_energy_symbol(self, xsec_network):
        """σ = 1e-18 above 13.6 eV only: zero in the lower band, full in the upper."""
        jfunc = (
            "@function photoCrossSection0()\n"
            "@return Piecewise((1e-18, ph_E >= 13.6), (0, True))\n"
        )
        net = xsec_network(jfunc=jfunc, radiation=RadiationConfig(bands=BANDS))
        rxn = _photo(net)

        assert rxn.photo.xsec.photodecay_expr.free_symbols == {net.radiation.ph_E}
        bands = [float(band.xsec) for band in rxn.photo.bands]
        # abs: far below the 1e-18 cross section, above noise in the zero band.
        assert bands == pytest.approx([0.0, 1e-18], rel=1e-6, abs=1e-30)

    def test_extra_symbols_stay_symbolic(self, xsec_network):
        jfunc = "@function photoCrossSection0(tgas)\n@return 1e-18*tgas\n"
        net = xsec_network(jfunc=jfunc, radiation=RadiationConfig(bands=BANDS))
        rxn = _photo(net)
        tgas = net.symbols.tgas

        assert float(rxn.photo.xsec_total / tgas) == pytest.approx(1e-18, rel=1e-6, abs=0)
        assert all(band.xsec.has(tgas) for band in rxn.photo.bands)
        assert rxn.rate_coefficient.has(tgas)

    def test_density_symbols_are_standardized(self, xsec_network):
        """A density in a custom σ reaches the band data and equations as nden."""
        jfunc = "@function photoCrossSection0(n_C)\n@return 1e-18*n_C\n"
        net = xsec_network(jfunc=jfunc, radiation=RadiationConfig(bands=BANDS))
        rxn = _photo(net)
        raw = sp.Symbol("n_C")
        nden = net.symbols.ndens[net.species["C"].index]

        exprs = [rxn.photo.xsec_total, *net.sradodes(order=2)]
        for band in rxn.photo.bands:
            exprs += [band.k, band.xsec]
        assert not any(sp.sympify(e).has(raw) for e in exprs)
        assert all(band.k.has(nden) for band in rxn.photo.bands)

    def test_per_reaction_pi_database_raises(self, xsec_network, tmp_path):
        with pytest.raises(ParserError, match="pi_database"):
            xsec_network(config=_toml(tmp_path, 'pi_database = "verner"'))


def _global_photo(**photo):
    """``[network]`` keys setting the global photo settings ([network.reactions.photo])."""
    return {"reactions": {"photo": photo}}


class TestOverride:
    def test_global_flag_uses_database(self, xsec_network):
        rxn = _photo(
            xsec_network(
                radiation=RadiationConfig(bands=BANDS),
                config=_global_photo(override_custom_cross_section=True),
            )
        )

        assert rxn.photo.xsec.database != "custom"

    def test_per_reaction_flag_needs_radiation(self, xsec_network, tmp_path):
        # Like the global flag, a reaction's own flag has no effect without
        # radiation: there is no cross-section rate to switch.
        cfg = _toml(tmp_path, "override_custom_cross_section = true")
        rxn = _photo(xsec_network(config=cfg))

        assert rxn.config.photo.override_custom_cross_section is True
        assert rxn.photo.xsec.database == "custom"

    def test_per_reaction_flag_applies_with_radiation(self, xsec_network, tmp_path):
        cfg = _toml(tmp_path, "override_custom_cross_section = true")
        rxn = _photo(xsec_network(config=cfg, radiation=RadiationConfig(bands=BANDS)))

        assert rxn.photo.xsec.database != "custom"

    def test_per_reaction_flag_beats_global(self, xsec_network, tmp_path):
        cfg = _toml(tmp_path, "override_custom_cross_section = false")
        rxn = _photo(
            xsec_network(
                config=cfg,
                radiation=RadiationConfig(bands=BANDS),
                **_global_photo(override_custom_cross_section=True),
            )
        )

        assert rxn.photo.xsec.database == "custom"

    def test_override_allows_pi_database(self, xsec_network, tmp_path):
        body = 'override_custom_cross_section = true\npi_database = "verner"'
        rxn = _photo(
            xsec_network(
                config=_toml(tmp_path, body), radiation=RadiationConfig(bands=BANDS)
            )
        )

        assert rxn.photo.xsec.database == "verner"

    def test_non_photo_reaction_raises(self, make_network, tmp_path):
        cfg = tmp_path / "jaff.toml"
        cfg.write_text(
            '[network.reactions.reaction."H.H__H2".photo]\noverride_custom_cross_section = true\n'
        )
        with pytest.raises(ParserError, match="photo"):
            make_network("H + H -> H2  []  1e-10\n", name="t.jet", config=str(cfg))

    def test_non_bool_raises(self, xsec_network, tmp_path):
        with pytest.raises(ParserError, match="override_custom_cross_section"):
            xsec_network(config=_toml(tmp_path, 'override_custom_cross_section = "yes"'))


def test_jaffgen_global_flag_reaches_radiation(tmp_path):
    cfg = tmp_path / "jaffgen.toml"
    cfg.write_text(
        "[network.radiation]\nbands = [6, 13.6, 100]\n"
        "[network.reactions.photo]\noverride_custom_cross_section = true\n"
    )
    gen = run_jaffgen(
        SimpleNamespace(
            network="h_photoionization",
            config=str(cfg),
            label=None,
            funcfile=None,
            duplicate_policy=None,
            expand_nuclei=None,
            strict=None,
            network_config=None,
            outdir=str(tmp_path / "out"),
            indir=None,
            files=None,
            template="microphysics",
        )
    )

    assert gen.net.radiation is not None
    assert gen.net.config.reactions.photo.override_custom_cross_section is True
