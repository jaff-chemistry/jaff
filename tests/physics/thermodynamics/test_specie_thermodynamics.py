# ABOUTME: Tests for specie.thermodynamics (SpecieThermodynamics): each species' adiabatic
# ABOUTME: index under the network's EOS, attached by Thermodynamics; None without one

import pytest
import sympy as sp

from jaff.physics.constants import k_B
from jaff.physics.thermodynamics import SpecieThermodynamics, Thermodynamics
from jaff.physics.thermodynamics.config import ThermodynamicsConfig
from jaff.physics.thermodynamics.eos.config import IdealEosConfig

GAMMA = 5.0 / 3.0
TWO_SPECIES = "@format:idx,R,R,P,rate\n1,H,H,H2,1\n2,H2,H,H,1\n"
MULTI = {"type": "multi_gamma", "default_gamma": GAMMA, "gamma_map": {"H2": 1.4}}


def _gammas(net) -> dict[str, sp.Expr]:
    return {s.name: s.thermodynamics.gamma for s in net.species}


def test_ideal_gives_its_gamma_to_every_species(make_network):
    net = make_network(TWO_SPECIES, funcfile=False)

    assert _gammas(net) == {"H": pytest.approx(GAMMA), "H2": pytest.approx(GAMMA)}


def test_multi_gamma_uses_the_map_then_the_default(make_network):
    net = make_network(TWO_SPECIES, funcfile=False, thermodynamics={"eos": MULTI})

    assert _gammas(net) == {"H": pytest.approx(GAMMA), "H2": pytest.approx(1.4)}


def test_gamma_is_a_constant_expression(make_network):
    net = make_network(TWO_SPECIES, funcfile=False, thermodynamics={"eos": MULTI})
    thermo = net.species["H2"].thermodynamics

    assert isinstance(thermo, SpecieThermodynamics)
    assert isinstance(thermo.gamma, sp.Expr)
    assert thermo.is_constant


def test_eos_without_adiabatic_indices_attaches_none(make_network):
    net = make_network(
        TWO_SPECIES, funcfile=False, thermodynamics={"eos": {"type": "fermi_degenerate"}}
    )

    assert [s.thermodynamics for s in net.species] == [None, None]


def test_multi_gamma_energy_is_built_from_the_species_gammas(make_network):
    net = make_network(TWO_SPECIES, funcfile=False, thermodynamics={"eos": MULTI})
    sym = net.symbols
    expected = sum(
        sym.ndens[s.index] * k_B.cgs.value * sym.tgas / (s.thermodynamics.gamma - 1)
        for s in net.species
    )

    assert sp.simplify(net.thermodynamics.eos.volumetric - expected) == 0


def test_a_thermodynamics_attaches_from_its_own_section(make_network):
    net = make_network(TWO_SPECIES, funcfile=False)
    Thermodynamics(net, ThermodynamicsConfig(eos=IdealEosConfig(gamma=1.4)))

    assert _gammas(net) == {"H": pytest.approx(1.4), "H2": pytest.approx(1.4)}


def test_eos_has_no_gammas(make_network):
    net = make_network(TWO_SPECIES, funcfile=False)

    assert not hasattr(net.thermodynamics.eos, "gammas")
