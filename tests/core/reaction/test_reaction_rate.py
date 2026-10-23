# ABOUTME: Tests for Reaction.rate, the shared
# ABOUTME: "rate coefficient times reactant densities" used by every equation builder

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
    assert sp.simplify(hh.rate - hh.rate_coefficient * h**2) == 0
    # Pseudo-species (_PHOTON) carry no density.
    photo = _rxn(net, "C._PHOTON__C+.e-", "photo")
    assert sp.simplify(photo.rate - photo.rate_coefficient * c) == 0


def test_reassigning_rate_coefficient_refreshes_rate(make_network):
    net = make_network(NETWORK, name="r.jet")
    rxn = _rxn(net, "H.H__H2")
    before = rxn.rate

    rxn.rate_coefficient = sp.Float(2.0e-17)

    assert rxn.rate != before
    h = net.symbols.ndens[net.species["H"].index]
    assert rxn.rate == sp.Float(2.0e-17) * h**2


def test_fluxes_use_rate(make_network):
    net = make_network(NETWORK, name="r.jet")

    assert net.sfluxes() == [rxn.rate for rxn in net.reactions]


def test_catalogue_rates_are_reaction_rates(make_network):
    net = make_network(NETWORK, name="r.jet")

    assert list(net.reactions.rates()) == [rxn.rate for rxn in net.reactions]


def test_catalogue_rate_coefficients_are_coefficients(make_network):
    net = make_network(NETWORK, name="r.jet")

    assert list(net.reactions.rate_coefficients()) == [
        rxn.rate_coefficient for rxn in net.reactions
    ]
