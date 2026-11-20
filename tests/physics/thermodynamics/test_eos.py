# ABOUTME: Unit tests for the symbolic ideal-gas EOS used by the Jacobian energy column
# ABOUTME: Volumetric, per-mass and per-particle internal energies on a real network


import pytest
import sympy as sp

from jaff import Network
from jaff.physics.thermodynamics import InternalEnergy
from jaff.physics.thermodynamics.eos.config import IdealEosConfig
from jaff.physics.constants import k_B
from tests._paths import FIXTURES_DIR

GAMMA = 5.0 / 3.0
TGAS = sp.symbols("tgas")


@pytest.fixture(scope="module")
def net() -> Network:
    return Network(
        str(FIXTURES_DIR / "react_cie_hepp.jet"),
        thermodynamics={"eos": IdealEosConfig(gamma=GAMMA)},
    )


def _volumetric(net: Network) -> sp.Expr:
    return net.symbols.ntot * k_B.cgs.value * TGAS / (GAMMA - 1.0)


def test_volumetric(net: Network) -> None:
    assert sp.simplify(net.thermodynamics.eos.volumetric - _volumetric(net)) == 0


def test_specific_per_mass(net: Network) -> None:
    expected = _volumetric(net) / net.symbols.rho
    assert sp.simplify(net.thermodynamics.eos.specific - expected) == 0


def test_specific_per_particle(net: Network) -> None:
    expected = k_B.cgs.value * TGAS / (GAMMA - 1.0)
    assert sp.simplify(net.thermodynamics.eos.per_particle - expected) == 0


class TestGammas:
    """``eos.gammas``: the adiabatic index of every species (ideal and multi_gamma)."""

    MULTI = {"type": "multi_gamma", "default_gamma": GAMMA, "gamma_map": {"H": 1.4}}

    def test_ideal_gives_its_gamma_to_every_species(self, net: Network) -> None:
        gammas = net.thermodynamics.eos.gammas

        assert gammas == [pytest.approx(GAMMA)] * net.species.count

    def test_multi_gamma_uses_the_map_then_the_default(self) -> None:
        net = Network(
            str(FIXTURES_DIR / "react_cie_hepp.jet"), thermodynamics={"eos": self.MULTI}
        )
        gammas = dict(zip(net.species.names(), net.thermodynamics.eos.gammas))

        assert gammas == {
            name: pytest.approx(1.4 if name == "H" else GAMMA)
            for name in net.species.names()
        }

    def test_match_the_volumetric_energy(self) -> None:
        net = Network(
            str(FIXTURES_DIR / "react_cie_hepp.jet"), thermodynamics={"eos": self.MULTI}
        )
        eos = net.thermodynamics.eos
        expected = sum(
            net.symbols.ndens[i] * k_B.cgs.value * TGAS / (g - 1.0)
            for i, g in enumerate(eos.gammas)
        )

        assert sp.simplify(eos.volumetric - expected) == 0

    def test_an_energy_without_them_has_none(self, net: Network) -> None:
        with pytest.raises(AttributeError, match="adiabatic indices"):
            _ = InternalEnergy(sp.Symbol("E"), net).gammas

    def test_a_sum_of_energies_has_none(self, net: Network) -> None:
        eos = net.thermodynamics.eos
        with pytest.raises(AttributeError, match="adiabatic indices"):
            _ = (eos + eos).gammas
