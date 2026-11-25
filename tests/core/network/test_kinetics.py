# ABOUTME: Tests for Kinetics (net.kinetics): the species ODEs and the chemical rates of
# ABOUTME: change of weighted density sums (weighted_rate, dntot_dt, drho_dt)

import random

import pytest
import sympy as sp
from sympy.core.function import AppliedUndef

from jaff.core.network import Kinetics, NetworkSymbols

TWO_SPECIES = "@format:idx,R,R,P,rate\n1,H,H,H2,1\n2,H2,H,H,1\n"
NUMBER_CONSERVING = "@format:idx,R,R,P,P,rate\n1,H2,O,OH,H,1\n"
IONISATION = "H -> H+ + E  []  1e-9\nH+ + E -> H  []  2e-10*(tgas/300)**(-0.7)\n"


@pytest.fixture
def net(make_network):
    return make_network(TWO_SPECIES, funcfile=False)


def _fluxes(net):
    return list(net.reactions.rates())


def _definitions(net):
    return {r.rate.symbol: r.rate.coefficient for r in net.reactions}


def _same(a, b):
    """*a* and *b* agree numerically at random values of their symbols and functions."""
    rng = random.Random(7)
    leaves = set()
    for expr in (a, b):
        leaves |= sp.sympify(expr).atoms(sp.Indexed, sp.Symbol, AppliedUndef)
    values = {
        leaf: rng.uniform(1.0, 3.0) * (100.0 if leaf == NetworkSymbols.tgas else 1.0)
        for leaf in leaves
        if not isinstance(leaf, sp.Symbol) or leaf == NetworkSymbols.tgas
    }

    def number(expr):
        return float(sp.sympify(expr).xreplace(values))

    return number(a) == pytest.approx(number(b), rel=1e-12)


# --------------------------------------------------------------------------- #
# The component                                                                #
# --------------------------------------------------------------------------- #


def test_network_has_one_kinetics(net):
    assert isinstance(net.kinetics, Kinetics)
    assert net.kinetics is net.kinetics
    assert net.kinetics.net is net


@pytest.mark.parametrize(
    ("owner", "name"),
    [
        ("net", "sfluxes"),
        ("net", "sodes"),
        ("net", "sradodes"),
        ("reactions", "odes"),
        ("symbols", "weighted_rate"),
        ("symbols", "dntot_dt"),
        ("symbols", "drho_dt"),
    ],
)
def test_old_equation_builders_are_gone(net, owner, name):
    obj = net if owner == "net" else getattr(net, owner)

    assert not hasattr(obj, name)


# --------------------------------------------------------------------------- #
# Species ODEs                                                                 #
# --------------------------------------------------------------------------- #


def test_odes_balance(make_network):
    net = make_network(IONISATION, name="e.jet", funcfile=False)
    odes = net.kinetics.odes()
    h, hp, e = (net.species[name].index for name in ("H", "H+", "e-"))
    rates = net.reactions.rates()

    assert len(odes) == net.species.core.count
    assert odes[h] == rates[1] - rates[0]
    assert odes[hp] == rates[0] - rates[1]
    assert odes[e] == odes[hp]


def test_reduced_odes_give_the_expanded_ones(make_network):
    net = make_network(IONISATION, name="e.jet", funcfile=False)
    reduced = net.kinetics.odes("reduced")

    assert all(row.has(net.symbols.krates) for row in reduced if row != 0)
    for row, expected in zip(reduced, net.kinetics.odes(), strict=True):
        assert _same(row.xreplace(_definitions(net)), expected)


def test_odes_reject_an_unknown_form(net):
    with pytest.raises(ValueError, match="expanded, reduced"):
        net.kinetics.odes("inline")


# --------------------------------------------------------------------------- #
# Weighted rates                                                               #
# --------------------------------------------------------------------------- #


def test_weighted_rate_matches_hand_sum(net):
    w_h, w_h2 = sp.symbols("w_h w_h2")
    weights = [0, 0]
    weights[net.species["H"].index], weights[net.species["H2"].index] = w_h, w_h2
    flux = _fluxes(net)
    # r1: H+H -> H2 : ν_H=-2, ν_H2=+1 ; r2: H2+H -> H : ν_H=0, ν_H2=-1
    expected = (-2 * w_h + w_h2) * flux[0] + (-w_h2) * flux[1]

    assert sp.simplify(net.kinetics.weighted_rate(weights) - expected) == 0


def test_weighted_rate_rejects_wrong_length(net):
    with pytest.raises(ValueError, match="one weight per species"):
        net.kinetics.weighted_rate([1])


def test_weighted_rate_rtol_only_drops_rounding(net):
    # rtol drops a coefficient only when it is negligible against the
    # magnitudes that cancel to produce it; real changes survive.
    masses = [s.mass for s in net.species]

    assert net.kinetics.weighted_rate(masses, rtol=1e-12) == net.kinetics.drho_dt()
    assert net.kinetics.weighted_rate(masses, rtol=1e-12) != 0


def test_reduced_weighted_rate_gives_the_expanded_one(make_network):
    net = make_network(IONISATION, name="e.jet", funcfile=False)
    weights = [1] * net.species.count
    reduced = net.kinetics.weighted_rate(weights, form="reduced")

    assert reduced.has(net.symbols.krates)
    assert _same(reduced.xreplace(_definitions(net)), net.kinetics.weighted_rate(weights))


# --------------------------------------------------------------------------- #
# dntot_dt / drho_dt                                                           #
# --------------------------------------------------------------------------- #


def test_dntot_dt_counts_particle_change(net):
    # TWO_SPECIES: H + H -> H2 (Δn = -1) ; H2 + H -> H (Δn = -1)
    flux = _fluxes(net)

    assert sp.simplify(net.kinetics.dntot_dt() - (-flux[0] - flux[1])) == 0


def test_dntot_dt_zero_when_number_conserving(make_network):
    # H2 + O -> OH + H : Δn = 0, so the reaction drops out exactly
    net = make_network(NUMBER_CONSERVING, funcfile=False)

    assert net.kinetics.dntot_dt() == 0


def test_drho_dt_is_mass_weighted_stoichiometric_sum(net):
    # Never forced to zero: Σ_r (Σ_i m_i ν_ri) F_r for every reaction.
    # TWO_SPECIES: r1 H + H -> H2 (Δm = m_H2 - 2 m_H), r2 H2 + H -> H (Δm = -m_H2)
    flux = _fluxes(net)
    m_h, m_h2 = net.species["H"].mass, net.species["H2"].mass
    expected = (m_h2 - 2 * m_h) * flux[0] + (-m_h2) * flux[1]

    assert sp.simplify(net.kinetics.drho_dt() - expected) == 0


def test_drho_dt_equals_weighted_rate_of_masses(net):
    masses = [s.mass for s in net.species]

    assert net.kinetics.drho_dt() == net.kinetics.weighted_rate(masses)


def test_drho_dt_drops_float_rounding_residues(make_network):
    # H + e- -> H+ + e- + e- conserves mass; with float masses its Δm is a
    # rounding residue (~1e-41 g), which must not leak into ρ̇.
    net = make_network("@format:idx,R,R,P,P,P,rate\n1,H,e-,H+,e-,e-,1\n", funcfile=False)
    h, e, hp = (net.species[s] for s in ("H", "e-", "H+"))

    assert hp.mass + 2 * e.mass - h.mass - e.mass != 0  # residue really exists
    assert net.kinetics.drho_dt() == 0


def test_drho_dt_does_not_trust_check_mass(make_network, monkeypatch):
    # H -> H+ loses one electron mass; even if the reaction were deemed
    # "conserved", the per-cell mass rate must keep its -m_e term.
    from jaff.core.reaction import Reaction

    monkeypatch.setattr(Reaction, "check_mass", lambda self: True)
    net = make_network("@format:idx,R,P,rate\n1,H,H+,1\n", funcfile=False)
    m_e = 9.10938371e-28
    flux = _fluxes(net)
    expected = (net.species["H+"].mass - net.species["H"].mass) * flux[0]

    assert sp.simplify(net.kinetics.drho_dt() - expected) == 0
    assert abs(net.species["H+"].mass - net.species["H"].mass + m_e) < 1e-40


@pytest.mark.parametrize("name", ["dntot_dt", "drho_dt"])
def test_reduced_density_rates_give_the_expanded_ones(make_network, name):
    net = make_network(IONISATION, name="e.jet", funcfile=False)
    reduced = getattr(net.kinetics, name)("reduced")

    assert _same(reduced.xreplace(_definitions(net)), getattr(net.kinetics, name)())
