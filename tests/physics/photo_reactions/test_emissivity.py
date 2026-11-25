# ABOUTME: globalEmissivity aux function: network-wide radiation emitted per volume per time
# ABOUTME: per photon energy, integrated per group and added to the density equations

import math

import astropy.units as u
import pytest
import sympy as sp

from jaff.physics.photo_reactions.config import RadiationConfig

EV2ERG = u.eV.to(u.erg)
BANDS = [6.0, 13.6, 100.0]
NET = "H + H -> H2  []  1e-17\n"
FLAT = "@function globalEmissivity()\n@return 1e-20\n"
U_MODE = {"mode": "u", "profile_index": 0.5}


@pytest.fixture
def emissive_network(make_network, tmp_path):
    """Network with *jfunc* as its aux file; kwargs go to ``Network``."""

    def _make(jfunc=FLAT, **kwargs):
        (tmp_path / "e.jfunc").write_text(jfunc)
        return make_network(NET, name="e.jet", **kwargs)

    return _make


def _floats(values):
    return [float(v) for v in values]


def test_without_radiation_ignored(emissive_network):
    assert emissive_network().radiation is None


def test_without_function_is_none(emissive_network):
    rad = emissive_network(jfunc="", radiation=RadiationConfig(bands=BANDS)).radiation

    assert rad.emissivity is None
    assert all(g.emissivity is None for g in rad.groups)
    assert rad.total_emissivity is None


def test_plain_emissivity_name_is_not_reserved(emissive_network):
    jfunc = "@function emissivity()\n@return 1e-20\n"
    rad = emissive_network(jfunc=jfunc, radiation=RadiationConfig(bands=BANDS)).radiation

    assert rad.emissivity is None


def test_spectral_expression_stored(emissive_network):
    rad = emissive_network(radiation=RadiationConfig(bands=BANDS)).radiation

    assert float(rad.emissivity) == pytest.approx(1e-20, rel=1e-6, abs=0)


def test_energy_mode_group_integrals(emissive_network):
    rad = emissive_network(radiation=RadiationConfig(bands=BANDS, **U_MODE)).radiation

    expected = [1e-20 * (13.6 - 6.0), 1e-20 * (100.0 - 13.6)]
    assert _floats(g.emissivity for g in rad.groups) == pytest.approx(
        expected, rel=1e-6, abs=0
    )
    assert float(rad.total_emissivity) == pytest.approx(
        1e-20 * (100.0 - 6.0), rel=1e-6, abs=0
    )


def test_photon_mode_group_integrals(emissive_network):
    rad = emissive_network(radiation=RadiationConfig(bands=BANDS)).radiation

    expected = [
        1e-20 * math.log(13.6 / 6.0) / EV2ERG,
        1e-20 * math.log(100.0 / 13.6) / EV2ERG,
    ]
    assert _floats(g.emissivity for g in rad.groups) == pytest.approx(
        expected, rel=1e-6, abs=0
    )


def test_density_symbols_are_standardized(emissive_network):
    jfunc = "@function globalEmissivity(n_H, tgas)\n@return 1e-30*n_H*tgas\n"
    net = emissive_network(jfunc=jfunc, radiation=RadiationConfig(bands=BANDS))

    nden = net.symbols.ndens[net.species["H"].index]
    for group in net.radiation.groups:
        assert group.emissivity.has(nden)
        assert group.emissivity.has(net.symbols.tgas)
        assert not group.emissivity.has(sp.Symbol("n_H"))


def test_symbols_reported_in_network_variables(emissive_network):
    jfunc = "@function globalEmissivity(tgas)\n@return 1e-30*tgas\n"
    net = emissive_network(jfunc=jfunc, radiation=RadiationConfig(bands=BANDS))

    assert net.symbols.tgas in net.symbols.variables


def test_enters_density_equations_only(emissive_network):
    net = emissive_network(radiation=RadiationConfig(bands=BANDS))
    odes = net.radiation.odes(order=2)
    nb = net.radiation.nbands

    for i, group in enumerate(net.radiation.groups):
        assert sp.simplify(odes[i] - group.emissivity) == 0
        assert odes[nb + i].is_zero


def test_radiation_computes_its_own_emissivity(emissive_network):
    from jaff import Network

    net = emissive_network(radiation=RadiationConfig(bands=BANDS))
    groups = net.radiation.groups
    before = [g.emissivity for g in groups]
    for g in groups:
        g.emissivity = None

    net.radiation.compute()

    assert [g.emissivity for g in groups] == before
    assert not hasattr(Network, "_Network__set_emissivity")
    assert not hasattr(Network, "_Network__set_absorption_coefficient")
