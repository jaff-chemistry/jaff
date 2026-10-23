# ABOUTME: spectralYield<N> aux function: radiation energy a reaction emits per event per
# ABOUTME: unit photon energy, integrated per group and added to the density equations

import math

import astropy.units as u
import pytest
import sympy as sp

from jaff.core.network import NetworkSymbols
from jaff.physics.photo_reactions.config import RadiationConfig

EV2ERG = u.eV.to(u.erg)
BANDS = [6.0, 13.6, 100.0]
# Radiative recombination: a non-photo reaction that emits.
RECOMB = "H+ + E -> H  []  2.5e-13\n"
FLAT = "@function spectralYield0()\n@return 1e-12\n"


@pytest.fixture
def yield_network(make_network, tmp_path):
    """Network from *lines* with *jfunc* as its aux file; kwargs go to ``Network``."""

    def _make(lines=RECOMB, jfunc=FLAT, **kwargs):
        (tmp_path / "y.jfunc").write_text(jfunc)
        return make_network(lines, name="y.jet", **kwargs)

    return _make


def _floats(values):
    return [float(v) for v in values]


def test_yield_is_carried_on_rad(yield_network):
    rxn = yield_network().reactions[0]

    assert float(rxn.rad.spectral_yield) == pytest.approx(1e-12, rel=1e-6, abs=0)


def test_reaction_without_yield(yield_network):
    rxn = yield_network(jfunc="").reactions[0]

    assert rxn.rad.spectral_yield is None
    assert rxn.rad.group_yields is None


def test_radiation_emitters_lists_reactions_with_a_yield(yield_network):
    net = yield_network(RECOMB + "H + H -> H2  []  1e-17\n")

    emitters = net.reactions.radiation_emitters()
    assert [r.verbatim for r in emitters] == ["H+ + e- -> H"]
    assert net.reactions.radiation_emitters() is emitters


def test_radiation_emitters_empty_without_yields(yield_network):
    assert len(yield_network(jfunc="").reactions.radiation_emitters()) == 0


def test_no_radiation_no_group_yields(yield_network):
    assert yield_network().reactions[0].rad.group_yields is None


def test_energy_mode_integrates_yield(yield_network):
    props = RadiationConfig(bands=BANDS, mode="u", profile_index=0.5)
    rxn = yield_network(radiation=props).reactions[0]

    expected = [1e-12 * (13.6 - 6.0), 1e-12 * (100.0 - 13.6)]
    assert _floats(rxn.rad.group_yields) == pytest.approx(expected, rel=1e-6, abs=0)


def test_photon_mode_counts_photons(yield_network):
    rxn = yield_network(radiation=RadiationConfig(bands=BANDS)).reactions[0]

    expected = [
        1e-12 * math.log(13.6 / 6.0) / EV2ERG,
        1e-12 * math.log(100.0 / 13.6) / EV2ERG,
    ]
    assert _floats(rxn.rad.group_yields) == pytest.approx(expected, rel=1e-6, abs=0)


def test_line_lands_in_its_group_only(yield_network):
    jfunc = (
        "@function spectralYield0()\n"
        "@return Piecewise((1e-12, (ph_E >= 20) & (ph_E <= 21)), (0, True))\n"
    )
    props = RadiationConfig(bands=BANDS, mode="u", profile_index=0.5)
    rxn = yield_network(jfunc=jfunc, radiation=props).reactions[0]

    # abs: far below the 1e-12 yield, above numerical noise in the zero group.
    assert _floats(rxn.rad.group_yields) == pytest.approx(
        [0.0, 1e-12], rel=1e-6, abs=1e-24
    )


def test_extra_symbols_stay_symbolic(yield_network):
    jfunc = "@function spectralYield0(tgas)\n@return 1e-16*tgas\n"
    props = RadiationConfig(bands=BANDS, mode="u", profile_index=0.5)
    rxn = yield_network(jfunc=jfunc, radiation=props).reactions[0]

    assert all(g.has(NetworkSymbols.tgas) for g in rxn.rad.group_yields)


def test_yield_symbols_reported_in_network_variables(yield_network):
    jfunc = "@function spectralYield0(tgas)\n@return 1e-16*tgas\n"
    net = yield_network(jfunc=jfunc, radiation=RadiationConfig(bands=BANDS))

    assert NetworkSymbols.tgas in net.symbols.variables


def test_symbolic_yield_without_closed_form_raises(yield_network):
    jfunc = "@function spectralYield0(tgas)\n@return exp(-tgas*ph_E)/(ph_E + log(ph_E))\n"
    props = RadiationConfig(bands=BANDS, mode="u", profile_index=0.5)
    with pytest.raises(ValueError, match="closed-form"):
        yield_network(jfunc=jfunc, radiation=props)


def test_yield_enters_density_equations_only(yield_network):
    net = yield_network(radiation=RadiationConfig(bands=BANDS))
    rxn = net.reactions[0]
    odes = net.sradodes(order=2)
    nb = net.radiation.nbands

    for i in range(nb):
        assert sp.simplify(odes[i] - rxn.rad.group_yields[i] * rxn.rate) == 0
        assert odes[nb + i].is_zero
