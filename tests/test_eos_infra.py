# ABOUTME: Robustness tests for the EOS infrastructure (EosProps, EosFactory)
# ABOUTME: Validation, InternalEnergy forms, builders, Thermodynamics.eos wiring

import gc
import weakref
from pathlib import Path
from types import SimpleNamespace

import pytest
import sympy as sp

from jaff import Network
from jaff.physics.constants import N_A, k_B
from jaff.physics import EosFactory, EosProps
from jaff.physics.thermodynamics import InternalEnergy

GAMMA = 5.0 / 3.0
TGAS = sp.Symbol("tgas")
E = sp.Symbol("E")


def _stub_net(tag: str) -> SimpleNamespace:
    symbols = SimpleNamespace(
        rho=sp.Symbol(f"rho_{tag}"), ntot=sp.Symbol(f"ntot_{tag}"), tgas=sp.Symbol("tgas")
    )
    return SimpleNamespace(symbols=symbols)


def _stub_species_net() -> SimpleNamespace:
    species = [SimpleNamespace(index=0, name="H"), SimpleNamespace(index=1, name="H2")]
    ndens = sp.IndexedBase("nden", shape=(2,))
    symbols = SimpleNamespace(
        ndens=ndens,
        ntot=ndens[0] + ndens[1],
        rho=sp.Symbol("rho"),
        tgas=sp.Symbol("tgas"),
    )
    return SimpleNamespace(species=species, symbols=symbols)


class TestEosIdentity:
    def test_same_expr_different_networks_bind_own_network(self) -> None:
        net_a, net_b = _stub_net("a"), _stub_net("b")
        expr = sp.Symbol("e")
        eos_a = InternalEnergy(expr, net_a)
        eos_b = InternalEnergy(expr, net_b)
        assert eos_b.specific == expr / net_b.symbols.rho
        assert eos_a.specific == expr / net_a.symbols.rho


class TestEosFactoryLifetime:
    def test_factory_not_kept_alive_after_generation(self) -> None:
        factory = EosFactory(_stub_net("a"), EosProps("ideal", gamma=5.0 / 3.0))
        factory.ideal()
        ref = weakref.ref(factory)
        del factory
        gc.collect()
        assert ref() is None


class TestEosPropsValidation:
    def test_missing_required_key_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="gamma_map"):
            EosProps("multi_gamma", default_gamma=GAMMA)

    def test_ideal_gamma_defaults_to_monoatomic(self) -> None:
        assert EosProps("ideal").gamma == 1.6666666666667

    def test_unknown_key_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="gama"):
            EosProps("ideal", gamma=5.0 / 3.0, gama=2.0)

    def test_unknown_type_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="bad"):
            EosProps("bad")

    def test_empty_gamma_map_accepted(self) -> None:
        props = EosProps("multi_gamma", default_gamma=5.0 / 3.0, gamma_map={})
        assert props.gamma_map == {}

    def test_non_numeric_gamma_raises_type_error(self) -> None:
        with pytest.raises(TypeError, match="gamma"):
            EosProps("ideal", gamma="5/3")

    def test_non_dict_gamma_map_raises_type_error(self) -> None:
        with pytest.raises(TypeError, match="gamma_map"):
            EosProps("multi_gamma", default_gamma=5.0 / 3.0, gamma_map=[1.4])

    @pytest.mark.parametrize("gamma", [1.0, 0.0, 0.5])
    def test_gamma_not_above_one_raises(self, gamma: float) -> None:
        with pytest.raises(ValueError, match="gamma"):
            EosProps("ideal", gamma=gamma)

    def test_gamma_map_value_of_one_raises(self) -> None:
        with pytest.raises(ValueError, match="H2"):
            EosProps("multi_gamma", default_gamma=5.0 / 3.0, gamma_map={"H2": 1.0})

    def test_default_gamma_of_one_raises(self) -> None:
        with pytest.raises(ValueError, match="default_gamma"):
            EosProps("multi_gamma", default_gamma=1.0, gamma_map={})


class TestEosFactoryDispatch:
    def test_relativistic_fermi_degenerate_not_implemented(self) -> None:
        factory = EosFactory(_stub_net("a"), EosProps("relativistic_fermi_degenerate"))
        with pytest.raises(NotImplementedError):
            factory.relativistic_fermi_degenerate()

    def test_invalid_type_message_separates_sentences(self) -> None:
        props = EosProps("ideal", gamma=5.0 / 3.0)
        props.type = "bad"
        with pytest.raises(ValueError, match=r"Invalid eos: 'bad'\. Valid eos types"):
            EosFactory(_stub_net("a"), props).generate()

    def test_generate_builds_eos(self) -> None:
        net = _stub_net("a")
        eos = EosFactory(net, EosProps("ideal", gamma=GAMMA)).generate()
        assert isinstance(eos, InternalEnergy)

    def test_builder_registry_matches_props_types(self) -> None:
        assert set(EosFactory._BUILDERS) == set(EosProps._REQUIRED)


class TestEosForms:
    @pytest.fixture
    def eos(self) -> InternalEnergy:
        return InternalEnergy(sp.Symbol("E"), _stub_net("a"))

    def test_volumetric_is_wrapped_expr(self, eos: InternalEnergy) -> None:
        assert eos.volumetric == sp.Symbol("E")

    def test_specific_divides_by_rho(self, eos: InternalEnergy) -> None:
        assert eos.specific == sp.Symbol("E") / sp.Symbol("rho_a")

    def test_per_particle_divides_by_ntot(self, eos: InternalEnergy) -> None:
        assert eos.per_particle == sp.Symbol("E") / sp.Symbol("ntot_a")

    def test_molar_scales_per_particle_by_avogadro(self, eos: InternalEnergy) -> None:
        ratio = sp.simplify(eos.molar / (eos.per_particle * N_A.cgs.value))
        assert float(ratio) == pytest.approx(1.0, rel=1e-14)


class TestEosNormaliser:
    @pytest.mark.parametrize(
        "form, expected",
        [
            ("volumetric", E),
            ("specific", E / sp.Symbol("rho_a")),
            ("per_particle", E / sp.Symbol("ntot_a")),
            ("molar", E * N_A.cgs.value / sp.Symbol("ntot_a")),
        ],
    )
    def test_normaliser_per_form(self, form: str, expected: sp.Expr) -> None:
        eos = InternalEnergy(sp.Symbol("E"), _stub_net("a"))
        assert sp.simplify(eos.normaliser(form) - expected) == 0

    def test_unknown_form_raises(self) -> None:
        eos = InternalEnergy(sp.Symbol("E"), _stub_net("a"))
        with pytest.raises(ValueError, match="bogus"):
            eos.normaliser("bogus")


class TestEosBuilders:
    def test_ideal_volumetric_energy(self) -> None:
        net = _stub_net("a")
        eos = EosFactory(net, EosProps("ideal", gamma=GAMMA)).ideal()
        expected = net.symbols.ntot * k_B.cgs.value * TGAS / (GAMMA - 1.0)
        assert sp.simplify(eos.volumetric - expected) == 0

    def test_multi_gamma_returns_eos(self) -> None:
        props = EosProps("multi_gamma", default_gamma=GAMMA, gamma_map={})
        assert isinstance(
            EosFactory(_stub_species_net(), props).multi_gamma(), InternalEnergy
        )

    def test_multi_gamma_per_species_energy(self) -> None:
        net = _stub_species_net()
        props = EosProps("multi_gamma", default_gamma=GAMMA, gamma_map={"H2": 1.4})
        eos = EosFactory(net, props).multi_gamma()
        kt = k_B.cgs.value * TGAS
        nden = net.symbols.ndens
        expected = nden[0] * kt / (GAMMA - 1.0) + nden[1] * kt / (1.4 - 1.0)
        assert sp.simplify(eos.volumetric - expected) == 0

    def test_fermi_degenerate_not_implemented(self) -> None:
        factory = EosFactory(_stub_net("a"), EosProps("fermi_degenerate"))
        with pytest.raises(NotImplementedError):
            factory.fermi_degenerate()


class TestNetworkEos:
    def test_without_props_defaults_to_ideal(self) -> None:
        fresh = Network(str(Path(__file__).parent / "fixtures" / "react_cie_hepp.jet"))
        sym = fresh.symbols
        expected = sym.ntot * k_B.cgs.value * TGAS / (1.6666666666667 - 1.0)
        assert sp.simplify(fresh.thermodynamics.eos.volumetric - expected) == 0

    def test_constructor_props_used(self) -> None:
        fresh = Network(
            str(Path(__file__).parent / "fixtures" / "react_cie_hepp.jet"),
            eos_props=EosProps("ideal", gamma=1.4),
        )
        expected = fresh.symbols.ntot * k_B.cgs.value * TGAS / (1.4 - 1.0)
        assert sp.simplify(fresh.thermodynamics.eos.volumetric - expected) == 0
