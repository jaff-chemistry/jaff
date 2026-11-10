# ABOUTME: Tests for Rate (reaction.rate): the expanded and reduced forms of a reaction
# ABOUTME: rate, the krate symbols, and the form argument of the rate and heating builders

import random

import pytest
import sympy as sp
from sympy.core.function import AppliedUndef

from jaff.core.network import NetworkSymbols

from jaff.core.reaction import Rate

NETWORK = (
    "H + H -> H2  []  1e-17\n"
    "C -> C+ + E  []  PHOTO, 11.3\n"
    "H + CR -> H+ + E  []  1e-17\n"
    "H+ + E -> H  []  2e-10*(tgas/300)**(-0.7)\n"
)


@pytest.fixture
def net(make_network):
    return make_network(NETWORK, name="r.jet")


def _rxn(net, serialized, rtype="unknown"):
    return net.reactions[serialized, rtype]


def _definitions(net):
    return {r.rate.symbol: r.rate.coefficient for r in net.reactions}


def _same(a, b):
    """*a* and *b* agree numerically at random values of their symbols and functions."""
    rng = random.Random(7)
    leaves = set()
    for expr in (a, b):
        leaves |= expr.atoms(sp.Indexed, sp.Symbol, AppliedUndef)
    values = {
        leaf: rng.uniform(1.0, 3.0) * (100.0 if leaf == NetworkSymbols.tgas else 1.0)
        for leaf in leaves
        if not isinstance(leaf, sp.Symbol) or leaf == NetworkSymbols.tgas
    }

    def number(expr):
        return float(sp.sympify(expr).xreplace(values))

    return number(a) == pytest.approx(number(b), rel=1e-12)


# --------------------------------------------------------------------------- #
# Rate                                                                         #
# --------------------------------------------------------------------------- #


def test_rate_is_a_rate(net):
    rxn = _rxn(net, "H.H__H2")

    assert isinstance(rxn.rate, Rate)
    assert rxn.rate.reaction is rxn


def test_symbol_is_krate_at_the_catalogue_position(net):
    krate = net.symbols.krates

    assert [r.rate.symbol for r in net.reactions] == [
        krate[i] for i in range(len(net.reactions))
    ]
    assert krate.shape == (len(net.reactions),)


def test_densities_are_the_core_reactant_densities(net):
    nden = net.symbols.ndens
    h, c = nden[net.species["H"].index], nden[net.species["C"].index]

    assert _rxn(net, "H.H__H2").rate.densities == h**2
    # Pseudo-species (_PHOTON, _CR) carry no density.
    assert _rxn(net, "C._PHOTON__C+.e-", "photo").rate.densities == c


def test_forms(net):
    rate = _rxn(net, "H+.e-__H").rate

    assert _same(rate.expanded, rate.coefficient * rate.densities)
    assert rate.reduced == rate.symbol * rate.densities
    assert rate.form("expanded") == rate.expanded
    assert rate.form("reduced") == rate.reduced


def test_coefficient_is_the_merged_segments(net):
    rxn = _rxn(net, "H+.e-__H")

    assert rxn.rate.coefficient == rxn.rate_segments.equivalent_rate


def test_coefficient_is_the_only_way_to_the_coefficient(net):
    assert not hasattr(_rxn(net, "H.H__H2"), "rate_coefficient")


def test_unknown_form_raises(net):
    with pytest.raises(ValueError, match="expanded, reduced"):
        _rxn(net, "H.H__H2").rate.form("inline")


def test_expanded_follows_the_coefficient(net):
    rxn = _rxn(net, "H.H__H2")
    rate = rxn.rate
    before = rate.expanded

    rxn.rate_segments.map_rates(lambda r: 2 * r)

    assert rxn.rate is rate
    assert rate.expanded != before
    assert _same(rate.expanded, rxn.rate.coefficient * rate.densities)


def test_expanded_is_cached(net):
    rate = _rxn(net, "H.H__H2").rate

    assert rate.expanded is rate.expanded


# --------------------------------------------------------------------------- #
# Builders                                                                     #
# --------------------------------------------------------------------------- #


def test_catalogue_rates_take_the_form(net):
    assert list(net.reactions.rates()) == [r.rate.expanded for r in net.reactions]
    assert list(net.reactions.rates("reduced")) == [r.rate.reduced for r in net.reactions]


@pytest.mark.parametrize("name", ["dEdt_chemical", "dEdt_tot"])
def test_reduced_heating_gives_the_expanded_one(net, name):
    # Set before the reduced form is first built (the expanded one is built on load).
    for r in net.reactions:
        r.dE = 1e-12
    reduced = getattr(net.thermodynamics, name)("reduced").volumetric
    expected = sum(1e-12 * r.rate.expanded for r in net.reactions)

    assert reduced.has(net.symbols.krates)
    assert _same(reduced.xreplace(_definitions(net)), expected)


@pytest.mark.parametrize("name", ["dTdt_chemical", "dTdt_tot"])
def test_reduced_temperature_rate_gives_the_expanded_one(net, name):
    thermo = net.thermodynamics
    reduced = getattr(thermo, name)("reduced")

    assert _same(reduced.xreplace(_definitions(net)), getattr(thermo, name)())


def test_thermodynamics_caches_each_form(net):
    thermo = net.thermodynamics

    assert thermo.dTdt_tot() is thermo.dTdt_tot("expanded")
    assert thermo.dTdt_tot("reduced") is thermo.dTdt_tot("reduced")
    assert thermo.dEdt_tot("reduced") is not thermo.dEdt_tot()


def test_thermodynamics_rejects_an_unknown_form(net):
    with pytest.raises(ValueError, match="expanded, reduced"):
        net.thermodynamics.dTdt_tot("inline")
