# ABOUTME: Regression tests for Radiation.odes: how per-band photo-reaction data is
# ABOUTME: assembled into the radiation density/flux moment equations and laid out

import pytest
import sympy as sp

from jaff.core.network import NetworkSymbols
from jaff.physics.photo_reactions.config import RadiationConfig

BANDS = [6.0, 13.6, 100.0]
PHOTO = "C -> C+ + E  []  PHOTO, 11.3\nCO -> C + O  []  PHOTO, 11\n"
NON_PHOTO = "H + H -> H2  []  1e-17\n"


@pytest.fixture
def radiation_network(make_network):
    def _make(lines=PHOTO, name="rad.jet"):
        return make_network(lines, name=name, radiation=RadiationConfig(bands=BANDS))

    return _make


def _flux(net, rxn, k):
    """Reaction flux for band coefficient *k*: k times every reactant density."""
    nden = net.symbols.ndens
    for reactant in rxn.reactants.core:
        k *= nden[net.species[str(reactant)].index]

    return k


def _expected(net):
    """Density and flux equations per band, built from each reaction's band data."""
    groups = net.radiation.groups
    rflux = sp.IndexedBase("rflux", shape=(net.radiation.nbands,))
    flux_map = {g.sym: rflux[i] for i, g in enumerate(groups)}
    loss = [sp.Float(0.0)] * len(groups)
    for rxn in net.reactions:
        if rxn.photo is None or not rxn.photo.bands:
            continue
        for i, band in enumerate(rxn.photo.bands):
            loss[i] -= _flux(net, rxn, band.k)

    density = loss
    fluxes = [loss[i].xreplace(flux_map) for i in range(len(groups))]

    return density, fluxes


def _same(a, b):
    return sp.simplify(sp.sympify(a) - sp.sympify(b)) == 0


def test_density_equation_removes_photo_reaction_flux(radiation_network):
    net = radiation_network()
    density, _ = _expected(net)
    odes = net.radiation.odes(order=2)
    nb = net.radiation.nbands

    assert len(odes) == 2 * nb
    for i in range(nb):
        assert _same(odes[i], density[i])
        assert odes[i] != 0


def test_flux_equation_is_density_loss_with_rflux(radiation_network):
    net = radiation_network()
    _, fluxes = _expected(net)
    odes = net.radiation.odes(order=2)
    nb = net.radiation.nbands

    for i in range(nb):
        assert _same(odes[nb + i], fluxes[i])
        assert not odes[nb + i].has(net.radiation.groups[i].sym)


def test_no_reaction_radiation_source_term(radiation_network):
    """Without a spectralYield<N>, reactions emit nothing (the old dRad is gone)."""
    net = radiation_network()

    assert not hasattr(net, "dRad_dt_extra")
    assert not hasattr(net.reactions, "dRad")
    for rxn in net.reactions:
        assert rxn.rad.spectral_yield is None
        assert rxn.rad.group_yields is None
        assert rxn.function("deltaRad") is None
        assert all(not hasattr(band, "delta_rad") for band in rxn.photo.bands)


def test_spectral_profiles_use_photon_energy_symbol(radiation_network):
    rad = radiation_network().radiation

    assert rad.ph_E is NetworkSymbols.ph_E
    assert rad.nph_profile.free_symbols == {NetworkSymbols.ph_E}
    assert all(g.ph_E is NetworkSymbols.ph_E for g in rad.groups)


def test_non_photo_reactions_do_not_contribute(radiation_network):
    net = radiation_network(PHOTO + NON_PHOTO)
    nden = net.symbols.ndens
    non_photo = [nden[net.species[s].index] for s in ("H", "H2")]

    for ode in net.radiation.odes(order=2):
        assert not any(ode.has(n) for n in non_photo)


@pytest.mark.parametrize("order", [0, 1, 2, 3])
def test_layout_places_each_band_by_ordered_index(radiation_network, order):
    net = radiation_network()
    density, fluxes = _expected(net)
    odes = net.radiation.odes(order=order)

    for i in range(net.radiation.nbands):
        ei, fi = net.radiation.ordered_index(i, order)
        assert _same(odes[ei], density[i])
        assert _same(odes[fi], fluxes[i])


def test_photo_reaction_without_band_data_raises(radiation_network):
    net = radiation_network()
    rxn = next(r for r in net.reactions if r.photo is not None)
    rxn.photo.bands = None

    with pytest.raises(RuntimeError, match="no per-band radiation data"):
        net.radiation.odes()


def test_rejects_unknown_order(radiation_network):
    with pytest.raises(ValueError, match="Invalid order"):
        radiation_network().radiation.odes(order=4)


def test_flux_equations_use_the_radiation_flux_symbol(radiation_network):
    net = radiation_network()
    odes = net.radiation.odes()
    fluxes = {
        s for e in odes for s in e.atoms(sp.Indexed) if s.base == net.radiation.flux
    }

    assert fluxes
