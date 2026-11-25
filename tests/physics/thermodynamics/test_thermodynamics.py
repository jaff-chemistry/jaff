# ABOUTME: Tests for Thermodynamics: cached EOS, the chemical and extra rates,
# ABOUTME: quotient-rule dE/dt forms and the dT/dt expression

import pytest
import sympy as sp

from jaff.physics.thermodynamics.config import ThermodynamicsConfig
from jaff.physics.thermodynamics.eos.config import IdealEosConfig
from jaff.physics.constants import N_A, k_B
from jaff.physics.thermodynamics import InternalEnergy, Thermodynamics

GAMMA = 5.0 / 3.0
TWO_SPECIES = "@format:idx,R,R,P,rate\n1,H,H,H2,1\n2,H2,H,H,1\n"


@pytest.fixture
def net(make_network):
    return make_network(
        TWO_SPECIES, funcfile=False, thermodynamics={"eos": IdealEosConfig(gamma=GAMMA)}
    )


def test_eos_is_cached_internal_energy(net):
    eos = net.thermodynamics.eos
    assert isinstance(eos, InternalEnergy)
    assert net.thermodynamics.eos is eos


def test_eos_uses_network_eos_config(net):
    sym = net.symbols
    expected = sym.ntot * k_B.cgs.value * sym.tgas / (GAMMA - 1.0)
    assert sp.simplify(net.thermodynamics.eos.volumetric - expected) == 0


def test_network_has_no_eos_method(net):
    assert not hasattr(net, "eos")


def test_dedt_specific_uses_quotient_rule(net):
    th, sym = net.thermodynamics, net.symbols
    e, edot, rho = th.eos.volumetric, th.dEdt_tot().volumetric, sym.rho
    expected = edot / rho - e * net.kinetics.drho_dt() / rho**2
    assert sp.simplify(th.dEdt_tot().specific - expected) == 0


def test_dedt_per_particle_uses_quotient_rule(net):
    th, sym = net.thermodynamics, net.symbols
    e, edot, n = th.eos.volumetric, th.dEdt_tot().volumetric, sym.ntot
    expected = edot / n - e * net.kinetics.dntot_dt() / n**2
    assert sp.simplify(th.dEdt_tot().per_particle - expected) == 0


def test_dedt_molar_is_avogadro_times_per_particle(net):
    th, sym = net.thermodynamics, net.symbols
    e, edot, n = th.eos.volumetric, th.dEdt_tot().volumetric, sym.ntot
    expected = N_A.cgs.value * (edot / n - e * net.kinetics.dntot_dt() / n**2)
    assert sp.simplify(th.dEdt_tot().molar - expected) == 0


def test_dtdt_tot_ideal_gas(net):
    th, sym = net.thermodynamics, net.symbols
    cv = k_B.cgs.value / (GAMMA - 1.0)
    expected = (th.dEdt_tot().volumetric - cv * sym.tgas * net.kinetics.dntot_dt()) / (
        cv * sym.ntot
    )
    assert sp.simplify(th.dTdt_tot() - expected) == 0


def _numerically_equal(a: sp.Expr, b: sp.Expr) -> bool:
    """Compare two expressions at fixed positive values of all their atoms.

    Avoids false negatives from SymPy keeping equal floats at different precisions.
    """
    atoms = sorted((a - b).atoms(sp.Indexed), key=str)
    values: dict = {x: 1.0 + 0.37 * i for i, x in enumerate(atoms)}
    a_num, b_num = a.xreplace(values), b.xreplace(values)
    symbols = sorted(a_num.free_symbols | b_num.free_symbols, key=str)
    values = {x: 2.0 + 0.53 * i for i, x in enumerate(symbols)}
    a_val, b_val = float(a_num.xreplace(values)), float(b_num.xreplace(values))
    return abs(a_val - b_val) <= 1e-12 * max(abs(a_val), abs(b_val))


LAMBDA = sp.Symbol("lam")


@pytest.fixture
def thermo_with_extra(make_network, tmp_path):
    """Thermodynamics with a non-zero heating/cooling rate ``lam`` (a .jfunc)."""
    (tmp_path / "heated.jfunc").write_text(
        "@function heatingCoolingRate()\n@return lam\n"
    )
    net = make_network(
        TWO_SPECIES,
        name="heated.dat",
        thermodynamics={"eos": IdealEosConfig(gamma=GAMMA)},
    )

    return net.thermodynamics


def test_dtdt_extra_has_no_composition_term(thermo_with_extra):
    # Heating/cooling does not change particle number: Ṫ_extra = Λ / (∂E/∂T)
    cv = k_B.cgs.value / (GAMMA - 1.0)
    expected = LAMBDA / (cv * thermo_with_extra.net.symbols.ntot)
    assert _numerically_equal(thermo_with_extra.dTdt_extra, expected)


def test_dtdt_chemical_carries_composition_term(thermo_with_extra):
    th, sym = thermo_with_extra, thermo_with_extra.net.symbols
    cv = k_B.cgs.value / (GAMMA - 1.0)
    edot_chem = th.dEdt_chemical().volumetric
    expected = (edot_chem - cv * sym.tgas * th.net.kinetics.dntot_dt()) / (cv * sym.ntot)
    assert sp.simplify(th.dTdt_chemical() - expected) == 0


def test_dtdt_parts_sum_to_total(thermo_with_extra):
    th = thermo_with_extra
    assert sp.simplify(th.dTdt_chemical() + th.dTdt_extra - th.dTdt_tot()) == 0


@pytest.mark.parametrize("name", ["dEdt_extra", "dTdt_extra"])
def test_extra_rates_are_lazy_and_cached(net, name):
    th = Thermodynamics(net, net.config.thermodynamics)
    assert name not in vars(th)
    assert getattr(th, name) is getattr(th, name)


@pytest.mark.parametrize(
    "name", ["dEdt_chemical", "dEdt_tot", "dTdt_chemical", "dTdt_tot"]
)
def test_reaction_rates_are_lazy_and_cached(net, name):
    th = Thermodynamics(net, net.config.thermodynamics)
    assert not th._by_form
    assert getattr(th, name)() is getattr(th, name)()


class TestConfig:
    """Thermodynamics takes its own section, like every configured component."""

    def test_network_passes_its_section(self, net):
        assert net.thermodynamics.config is net.config.thermodynamics

    def test_eos_comes_from_the_section_given(self, net):
        th = Thermodynamics(net, ThermodynamicsConfig(eos=IdealEosConfig(gamma=1.4)))
        sym = net.symbols
        expected = sym.ntot * k_B.cgs.value * sym.tgas / (1.4 - 1.0)

        assert sp.simplify(th.eos.volumetric - expected) == 0


ENERGY_FORMS = ["volumetric", "specific", "per_particle", "molar"]


class TestTemperatureDerivatives:
    """``temperature_derivatives(energy)``: how T responds to the evolved energy and n."""

    @staticmethod
    def _solved_temperature(net, energy: str) -> tuple[sp.Symbol, sp.Expr]:
        """T(e, n) solved from the EOS, with ``e`` a plain symbol."""
        e = sp.Symbol("e_state")
        expr = net.thermodynamics.eos.normaliser(energy)
        (tgas,) = sp.solve(sp.Eq(e, expr), net.symbols.tgas)
        return e, tgas

    @pytest.mark.parametrize("energy", ENERGY_FORMS)
    def test_match_derivatives_of_solved_temperature(self, net, energy) -> None:
        sym = net.symbols
        e, t_of_e = self._solved_temperature(net, energy)
        dT_de, dT_dn = net.thermodynamics.temperature_derivatives(energy)
        n = net.species.count
        assert len(dT_dn) == n

        point = {sym.ndens[i]: 1.0e3 * (i + 1) for i in range(n)}
        t0 = 250.0
        e0 = float(
            net.thermodynamics.eos.normaliser(energy).xreplace({**point, sym.tgas: t0})
        )
        at_t = {**point, sym.tgas: t0}
        at_e = {**point, e: e0}
        assert float(dT_de.xreplace(at_t)) == pytest.approx(
            float(sp.diff(t_of_e, e).xreplace(at_e)), rel=1e-12
        )
        for j in range(n):
            assert float(dT_dn[j].xreplace(at_t)) == pytest.approx(
                float(sp.diff(t_of_e, sym.ndens[j]).xreplace(at_e)), rel=1e-12
            )

    def test_unknown_energy_form(self, net) -> None:
        with pytest.raises(ValueError, match="bogus"):
            net.thermodynamics.temperature_derivatives("bogus")
