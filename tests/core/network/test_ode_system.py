# ABOUTME: Tests for OdeSystem: the state layout per thermal mode and radiation order,
# ABOUTME: the rhs assembled from the components, and validation of the choices

import pytest
import sympy as sp

from jaff import Network
from jaff.core.network import OdeSystem, StateVariable
from tests._paths import FIXTURES_DIR

BANDS = {"bands": [1, 13.6, 100, "inf"]}


@pytest.fixture(scope="module")
def dedt_network():
    """A single-reaction network that also carries an internal-energy rate."""
    return Network(str(FIXTURES_DIR / "test_jac_dedt.dat"))


@pytest.fixture(scope="module")
def rad_network():
    """Three species, three radiation bands."""
    return Network("h_photoionization", radiation=BANDS)


def _kinds(system: OdeSystem) -> list[str]:
    return [var.kind for var in system.state]


class TestState:
    def test_default_is_species_only(self, dedt_network) -> None:
        system = dedt_network.ode_system()
        ndens = dedt_network.symbols.ndens
        assert system.state == tuple(
            StateVariable("species", ndens[i], i) for i in range(3)
        )

    @pytest.mark.parametrize(
        ("thermal", "kind"), [("dedt", "energy"), ("dtdt", "temperature")]
    )
    def test_thermal_slot_follows_species(self, dedt_network, thermal, kind) -> None:
        system = dedt_network.ode_system(thermal=thermal)
        assert _kinds(system) == ["species"] * 3 + [kind]
        assert system.state[3] == StateVariable(kind, dedt_network.symbols.tgas, 0)

    @pytest.mark.parametrize("order", [0, 1, 2, 3])
    def test_radiation_slots_follow_order(self, rad_network, order) -> None:
        system = rad_network.ode_system(thermal="dtdt", radiation=True, rad_order=order)
        rad = rad_network.radiation
        n = rad_network.species.count + 1
        assert len(system.state) == n + 2 * rad.nbands
        for i in range(rad.nbands):
            ei, fi = rad.ordered_index(i, order)
            assert system.state[n + ei] == StateVariable(
                "radiation_density", rad.den[i], i
            )
            assert system.state[n + fi] == StateVariable("radiation_flux", rad.flux[i], i)

    def test_options_are_kept(self, rad_network) -> None:
        system = rad_network.ode_system(
            thermal="dedt", energy="specific", radiation=True, rad_order=2
        )
        assert (system.thermal, system.energy, system.radiation, system.rad_order) == (
            "dedt",
            "specific",
            True,
            2,
        )


class TestRhs:
    def test_is_column_matrix_in_state_order(self, rad_network) -> None:
        system = rad_network.ode_system(thermal="dtdt", radiation=True, rad_order=3)
        assert isinstance(system.rhs, sp.ImmutableMatrix)
        assert system.rhs.shape == (len(system.state), 1)

    def test_rows_come_from_the_components(self, rad_network) -> None:
        system = rad_network.ode_system(thermal="dtdt", radiation=True, rad_order=1)
        expected = [
            *rad_network.reactions.odes(rad_network.species),
            rad_network.thermodynamics.dTdt_tot,
            *rad_network.radiation.odes(1),
        ]
        assert list(system.rhs) == expected

    @pytest.mark.parametrize(
        "energy", ["volumetric", "specific", "per_particle", "molar"]
    )
    def test_dedt_row_is_energy_rate_in_form(self, dedt_network, energy) -> None:
        system = dedt_network.ode_system(thermal="dedt", energy=energy)
        assert system.rhs[3] == dedt_network.thermodynamics.dEdt_tot.normaliser(energy)

    def test_none_has_species_rows_only(self, dedt_network) -> None:
        system = dedt_network.ode_system()
        assert list(system.rhs) == dedt_network.reactions.odes(dedt_network.species)

    def test_is_cached(self, dedt_network) -> None:
        system = dedt_network.ode_system(thermal="dedt")
        assert system.rhs is system.rhs


class TestValidation:
    def test_unknown_thermal(self, dedt_network) -> None:
        with pytest.raises(ValueError, match="thermal"):
            dedt_network.ode_system(thermal="bogus")

    def test_unknown_energy_with_dedt(self, dedt_network) -> None:
        with pytest.raises(ValueError, match="bogus"):
            dedt_network.ode_system(thermal="dedt", energy="bogus")

    def test_radiation_without_field(self, dedt_network) -> None:
        with pytest.raises(ValueError, match="radiation"):
            dedt_network.ode_system(radiation=True)

    def test_bad_rad_order(self, rad_network) -> None:
        with pytest.raises(ValueError, match="order"):
            rad_network.ode_system(radiation=True, rad_order=4)

    def test_rad_order_ignored_without_radiation(self, rad_network) -> None:
        assert len(rad_network.ode_system(rad_order=7).state) == 3
