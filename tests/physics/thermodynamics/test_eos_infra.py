# ABOUTME: Robustness tests for the EOS infrastructure (EOS config models, EosFactory)
# ABOUTME: InternalEnergy forms, builders, Thermodynamics.eos wiring from the config

import gc
import weakref
from types import SimpleNamespace

import pytest
import sympy as sp

from jaff import Network
from jaff.physics.constants import N_A, k_B
from jaff.physics import EosFactory
from jaff.physics.thermodynamics.eos.config import (
    EOS_TYPES,
    FermiDegenerateEosConfig,
    IdealEosConfig,
    MultiGammaEosConfig,
    RelativisticFermiDegenerateEosConfig,
)
from jaff.physics.thermodynamics import InternalEnergy
from tests._paths import FIXTURES_DIR

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
        factory = EosFactory(_stub_net("a"), IdealEosConfig(gamma=5.0 / 3.0))
        factory.ideal()
        ref = weakref.ref(factory)
        del factory
        gc.collect()
        assert ref() is None


class TestEosFactoryDispatch:
    def test_relativistic_fermi_degenerate_not_implemented(self) -> None:
        factory = EosFactory(_stub_net("a"), RelativisticFermiDegenerateEosConfig())
        with pytest.raises(NotImplementedError):
            factory.relativistic_fermi_degenerate()

    def test_generate_builds_eos(self) -> None:
        net = _stub_net("a")
        eos = EosFactory(net, IdealEosConfig(gamma=GAMMA)).generate()
        assert isinstance(eos, InternalEnergy)

    def test_builder_registry_matches_props_types(self) -> None:
        assert set(EosFactory._BUILDERS) == set(EOS_TYPES)


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
        eos = EosFactory(net, IdealEosConfig(gamma=GAMMA)).ideal()
        expected = net.symbols.ntot * k_B.cgs.value * TGAS / (GAMMA - 1.0)
        assert sp.simplify(eos.volumetric - expected) == 0

    def test_multi_gamma_returns_eos(self) -> None:
        props = MultiGammaEosConfig(default_gamma=GAMMA, gamma_map={})
        assert isinstance(
            EosFactory(_stub_species_net(), props).multi_gamma(), InternalEnergy
        )

    def test_multi_gamma_per_species_energy(self) -> None:
        net = _stub_species_net()
        props = MultiGammaEosConfig(default_gamma=GAMMA, gamma_map={"H2": 1.4})
        eos = EosFactory(net, props).multi_gamma()
        kt = k_B.cgs.value * TGAS
        nden = net.symbols.ndens
        expected = nden[0] * kt / (GAMMA - 1.0) + nden[1] * kt / (1.4 - 1.0)
        assert sp.simplify(eos.volumetric - expected) == 0

    def test_fermi_degenerate_not_implemented(self) -> None:
        factory = EosFactory(_stub_net("a"), FermiDegenerateEosConfig())
        with pytest.raises(NotImplementedError):
            factory.fermi_degenerate()


class TestNetworkEos:
    def test_without_props_defaults_to_ideal(self) -> None:
        fresh = Network(str(FIXTURES_DIR / "react_cie_hepp.jet"))
        sym = fresh.symbols
        expected = sym.ntot * k_B.cgs.value * TGAS / (1.6666666666667 - 1.0)
        assert sp.simplify(fresh.thermodynamics.eos.volumetric - expected) == 0

    def test_constructor_props_used(self) -> None:
        fresh = Network(
            str(FIXTURES_DIR / "react_cie_hepp.jet"),
            thermodynamics={"eos": IdealEosConfig(gamma=1.4)},
        )
        expected = fresh.symbols.ntot * k_B.cgs.value * TGAS / (1.4 - 1.0)
        assert sp.simplify(fresh.thermodynamics.eos.volumetric - expected) == 0

    def test_jaff_toml_eos_used(self, tmp_path) -> None:
        net_file = tmp_path / "react_cie_hepp.jet"
        net_file.write_text((FIXTURES_DIR / "react_cie_hepp.jet").read_text())
        (tmp_path / "jaff.toml").write_text("[network.thermodynamics.eos]\ngamma = 1.3\n")

        fresh = Network(str(net_file))

        expected = fresh.symbols.ntot * k_B.cgs.value * TGAS / (1.3 - 1.0)
        assert sp.simplify(fresh.thermodynamics.eos.volumetric - expected) == 0

    def test_constructor_eos_replaces_jaff_toml_eos_whole(self, tmp_path) -> None:
        net_file = tmp_path / "react_cie_hepp.jet"
        net_file.write_text((FIXTURES_DIR / "react_cie_hepp.jet").read_text())
        (tmp_path / "jaff.toml").write_text(
            '[network.thermodynamics.eos]\ntype = "multi_gamma"\n'
            "default_gamma = 1.4\ngamma_map = { H2 = 1.4 }\n"
        )

        fresh = Network(str(net_file), thermodynamics={"eos": IdealEosConfig(gamma=1.2)})

        assert fresh.config.thermodynamics.eos == IdealEosConfig(gamma=1.2)
