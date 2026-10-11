# ABOUTME: Unit tests for the symbolic ideal-gas EOS used by the Jacobian energy column
# ABOUTME: Volumetric, per-mass and per-particle internal energies on a real network

from pathlib import Path

import pytest
import sympy as sp

from jaff import Network
from jaff.physics import EosProps
from jaff.physics.constants import k_B

GAMMA = 5.0 / 3.0
TGAS = sp.symbols("tgas")


@pytest.fixture(scope="module")
def net() -> Network:
    return Network(
        str(Path(__file__).parent / "fixtures" / "react_cie_hepp.jet"),
        eos_props=EosProps("ideal", gamma=GAMMA),
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
