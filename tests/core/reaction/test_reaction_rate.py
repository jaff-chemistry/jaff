# ABOUTME: Tests for Reaction.rate, the shared
# ABOUTME: "rate coefficient times reactant densities" used by every equation builder

import pytest
import sympy as sp

NETWORK = (
    "H + H -> H2  []  1e-17\nC -> C+ + E  []  PHOTO, 11.3\nH + CR -> H+ + E  []  1e-17\n"
)


def _rxn(net, serialized, rtype="unknown"):
    return net.reactions[serialized, rtype]


def test_rate_is_coefficient_times_reactant_densities(make_network):
    net = make_network(NETWORK, name="r.jet")
    nden = net.symbols.ndens
    h, c = nden[net.species["H"].index], nden[net.species["C"].index]

    hh = _rxn(net, "H.H__H2")
    assert sp.simplify(hh.rate.expanded - hh.rate_coefficient * h**2) == 0
    # Pseudo-species (_PHOTON) carry no density.
    photo = _rxn(net, "C._PHOTON__C+.e-", "photo")
    assert sp.simplify(photo.rate.expanded - photo.rate_coefficient * c) == 0


def test_rate_coefficient_is_read_only(make_network):
    rxn = _rxn(make_network(NETWORK, name="r.jet"), "H.H__H2")

    with pytest.raises(AttributeError):
        rxn.rate_coefficient = sp.Float(2.0e-17)


def test_rate_follows_the_segments(make_network):
    net = make_network(NETWORK, name="r.jet")
    rxn = _rxn(net, "H.H__H2")
    before = rxn.rate.expanded

    rxn.rate_segments.map_rates(lambda rate: 2 * rate)

    assert rxn.rate.expanded != before
    h = net.symbols.ndens[net.species["H"].index]
    assert rxn.rate_coefficient == pytest.approx(2.0e-17)
    assert rxn.rate.expanded == rxn.rate_coefficient * h**2


def test_rate_is_cached(make_network):
    rxn = _rxn(make_network(NETWORK, name="r.jet"), "H.H__H2")

    assert rxn.rate.expanded is rxn.rate.expanded


def test_photo_band_rate_is_the_rate_coefficient(make_network):
    net = make_network(
        NETWORK, name="r.jet", radiation={"bands": [6.0, 13.6, 100.0]}, funcfile=False
    )
    photo = _rxn(net, "C._PHOTON__C+.e-", "photo")

    assert photo.photo.band_rate is not None
    assert photo.rate_coefficient is photo.photo.band_rate
    assert photo.rate_segments[0].rate is not photo.photo.band_rate


def test_photo_without_radiation_uses_its_segments(make_network):
    photo = _rxn(make_network(NETWORK, name="r.jet"), "C._PHOTON__C+.e-", "photo")

    assert photo.photo.band_rate is None
    assert photo.rate_coefficient == photo.rate_segments[0].rate


def test_fluxes_use_rate(make_network):
    net = make_network(NETWORK, name="r.jet")

    assert net.sfluxes() == [rxn.rate.expanded for rxn in net.reactions]


def test_catalogue_rates_are_reaction_rates(make_network):
    net = make_network(NETWORK, name="r.jet")

    assert list(net.reactions.rates()) == [rxn.rate.expanded for rxn in net.reactions]


def test_catalogue_rate_coefficients_are_coefficients(make_network):
    net = make_network(NETWORK, name="r.jet")

    assert list(net.reactions.rate_coefficients()) == [
        rxn.rate_coefficient for rxn in net.reactions
    ]
