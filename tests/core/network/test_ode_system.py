# ABOUTME: Tests for OdeSystem: the state layout per thermal mode and radiation order,
# ABOUTME: the rhs and Jacobian in both rate forms, and validation of the choices

import random

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
            thermal="dedt", energy="specific", radiation=True, rad_order=2, form="reduced"
        )
        assert (
            system.thermal,
            system.energy,
            system.radiation,
            system.rad_order,
            system.form,
        ) == ("dedt", "specific", True, 2, "reduced")


class TestRhs:
    def test_is_column_matrix_in_state_order(self, rad_network) -> None:
        system = rad_network.ode_system(thermal="dtdt", radiation=True, rad_order=3)
        assert isinstance(system.rhs, sp.ImmutableMatrix)
        assert system.rhs.shape == (len(system.state), 1)

    def test_rows_come_from_the_components(self, rad_network) -> None:
        system = rad_network.ode_system(thermal="dtdt", radiation=True, rad_order=1)
        expected = [
            *rad_network.kinetics.odes(),
            rad_network.thermodynamics.dTdt_tot(),
            *rad_network.radiation.odes(1),
        ]
        assert list(system.rhs) == expected

    @pytest.mark.parametrize(
        "energy", ["volumetric", "specific", "per_particle", "molar"]
    )
    def test_dedt_row_is_energy_rate_in_form(self, dedt_network, energy) -> None:
        system = dedt_network.ode_system(thermal="dedt", energy=energy)
        assert system.rhs[3] == dedt_network.thermodynamics.dEdt_tot().normaliser(energy)

    def test_none_has_species_rows_only(self, dedt_network) -> None:
        system = dedt_network.ode_system()
        assert list(system.rhs) == dedt_network.kinetics.odes()

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


TGAS_NETWORK = [
    "H+ + e- -> H [10,1e4] 2e-10*tgas**(-0.7)",
    "H + H -> H2 [10,1e4] 1e-17*(tgas/300)**0.5",
    "H2 + e- -> H + H + e- [10,1e4] 3e-9*exp(-5e4/tgas)",
]


@pytest.fixture(scope="module")
def tgas_network(tmp_path_factory):
    """Temperature-dependent rates on an ideal-gas EOS (T solvable from e)."""
    path = tmp_path_factory.mktemp("ode") / "tgas.dat"
    path.write_text("\n".join(TGAS_NETWORK) + "\n", encoding="utf-8")
    return Network(str(path), funcfile=False)


def _symbols(system: OdeSystem) -> list:
    return [var.symbol for var in system.state]


class TestJacobian:
    @pytest.mark.parametrize("thermal", ["none", "dtdt"])
    def test_is_plain_jacobian_without_energy(self, tgas_network, thermal) -> None:
        system = tgas_network.ode_system(thermal=thermal)
        expected = system.rhs.jacobian(_symbols(system))
        assert system.jacobian == sp.ImmutableSparseMatrix(expected)

    def test_radiation_columns(self, rad_network) -> None:
        system = rad_network.ode_system(thermal="dtdt", radiation=True, rad_order=2)
        expected = system.rhs.jacobian(_symbols(system))
        assert system.jacobian == sp.ImmutableSparseMatrix(expected)

    def test_is_sparse_immutable_and_cached(self, rad_network) -> None:
        system = rad_network.ode_system(radiation=True)
        jac = system.jacobian
        assert isinstance(jac, sp.ImmutableSparseMatrix)
        n = len(system.state)
        assert jac.shape == (n, n)
        assert 0 not in jac.todok().values()
        assert system.jacobian is jac

    @pytest.mark.parametrize(
        "energy", ["volumetric", "specific", "per_particle", "molar"]
    )
    def test_dedt_matches_rhs_written_in_energy(self, tgas_network, energy) -> None:
        """Independent check: substitute T(e, n) and differentiate by (n, e)."""
        net = tgas_network
        sym = net.symbols
        system = net.ode_system(thermal="dedt", energy=energy)
        e = sp.Symbol("e_state")
        e_expr = net.thermodynamics.eos.normaliser(energy)
        (t_of_e,) = sp.solve(sp.Eq(e, e_expr), sym.tgas)
        rhs_in_e = system.rhs.xreplace({sym.tgas: t_of_e})
        state_in_e = [e if var.kind == "energy" else var.symbol for var in system.state]
        expected = rhs_in_e.jacobian(state_in_e)

        n = net.species.count
        point = {sym.ndens[i]: 1.0e2 * (i + 2) for i in range(n)}
        t0 = 800.0
        e0 = float(e_expr.xreplace({**point, sym.tgas: t0}))
        got = system.jacobian.xreplace({**point, sym.tgas: t0})
        want = expected.xreplace({**point, e: e0})
        size = len(system.state)
        for i in range(size):
            # Exact zeros may come out as round-off of the row's scale
            scale = max(abs(float(want[i, j])) for j in range(size))
            for j in range(size):
                assert float(got[i, j]) == pytest.approx(
                    float(want[i, j]), rel=1e-9, abs=1e-12 * scale
                ), (i, j)


# --------------------------------------------------------------------------- #
# Reduced form: krate[r] symbols in the rows, the chain rule in the Jacobian   #
# --------------------------------------------------------------------------- #

FORM_CASES = [
    pytest.param("tgas", {}, id="none"),
    pytest.param("tgas", {"thermal": "dtdt"}, id="dtdt"),
    *(
        pytest.param("tgas", {"thermal": "dedt", "energy": e}, id=f"dedt-{e}")
        for e in ("volumetric", "specific", "per_particle", "molar")
    ),
    pytest.param("rad", {"radiation": True}, id="radiation"),
    pytest.param(
        "rad", {"thermal": "dtdt", "radiation": True, "rad_order": 2}, id="rad-dtdt"
    ),
]


@pytest.fixture
def networks(tgas_network, rad_network):
    return {"tgas": tgas_network, "rad": rad_network}


def _state_values(system: OdeSystem, seed: int) -> dict:
    """Random positive values of every state symbol (tgas in 100-300 K)."""
    rng = random.Random(seed)
    tgas = system.network.symbols.tgas
    values = {var.symbol: rng.uniform(1.0, 3.0) for var in system.state}
    values[tgas] = rng.uniform(100.0, 300.0)

    return values


def _numbers(matrix: sp.MatrixBase, values: dict) -> list[list[float]]:
    rows, cols = matrix.shape
    evaluated = matrix.xreplace(values)

    return [[float(evaluated[i, j]) for j in range(cols)] for i in range(rows)]


def _assert_close(got: list[list[float]], want: list[list[float]]) -> None:
    for i, (got_row, want_row) in enumerate(zip(got, want, strict=True)):
        # Exact zeros may come out as round-off of the row's scale
        scale = max(abs(x) for x in want_row)
        for j, (g, w) in enumerate(zip(got_row, want_row, strict=True)):
            assert g == pytest.approx(w, rel=1e-9, abs=1e-12 * scale), (i, j)


class TestReducedForm:
    def test_expanded_is_the_default(self, rad_network) -> None:
        system = rad_network.ode_system(radiation=True)

        assert system.form == "expanded"
        assert system.rate_coefficients == {}
        assert not system.rhs.has(rad_network.symbols.krates)

    def test_rate_coefficients_define_every_krate(self, rad_network) -> None:
        system = rad_network.ode_system(form="reduced")

        assert system.rate_coefficients == {
            r.rate.symbol: r.rate.coefficient for r in rad_network.reactions
        }
        assert list(system.rate_coefficients) == [
            rad_network.symbols.krates[i] for i in range(len(rad_network.reactions))
        ]

    def test_rows_hold_the_krate_symbols(self, rad_network) -> None:
        system = rad_network.ode_system(thermal="dtdt", radiation=True, form="reduced")
        species = rad_network.species.count

        assert all(row.has(rad_network.symbols.krates) for row in system.rhs[:species])
        assert system.rhs[species].has(rad_network.symbols.krates)

    def test_unknown_form_raises(self, rad_network) -> None:
        with pytest.raises(ValueError, match="expanded, reduced"):
            rad_network.ode_system(form="inline")

    def test_temperature_column_reuses_the_rate_symbols(self, tgas_network) -> None:
        # dk/dT = krate[r] * dln(k)/dT: the Arrhenius rates are not written out again.
        system = tgas_network.ode_system(thermal="dtdt", form="reduced")
        t = tgas_network.species.count
        column = [system.jacobian[i, t] for i in range(t)]

        assert any(entry.has(tgas_network.symbols.krates) for entry in column)
        assert not any(entry.has(sp.exp) for entry in column)

    def test_photo_rates_depend_on_the_radiation_state(self, rad_network) -> None:
        # The radiation case exercises the chain rule through a non-thermal column.
        system = rad_network.ode_system(radiation=True, form="reduced")
        state = {var.symbol for var in system.state}

        assert any(
            coefficient.free_symbols & state - {rad_network.symbols.tgas}
            for coefficient in system.rate_coefficients.values()
        )

    @pytest.mark.parametrize(("which", "options"), FORM_CASES)
    def test_rhs_gives_the_expanded_one(self, networks, which, options) -> None:
        net = networks[which]
        reduced = net.ode_system(**options, form="reduced")
        expanded = net.ode_system(**options)
        rhs = reduced.rhs.xreplace(reduced.rate_coefficients)

        for seed in range(3):
            values = _state_values(reduced, seed)
            _assert_close(_numbers(rhs, values), _numbers(expanded.rhs, values))

    @pytest.mark.parametrize(("which", "options"), FORM_CASES)
    def test_jacobian_gives_the_expanded_one(self, networks, which, options) -> None:
        net = networks[which]
        reduced = net.ode_system(**options, form="reduced")
        expanded = net.ode_system(**options)
        jacobian = reduced.jacobian.xreplace(reduced.rate_coefficients)

        for seed in range(3):
            values = _state_values(reduced, seed)
            _assert_close(_numbers(jacobian, values), _numbers(expanded.jacobian, values))

    @pytest.mark.parametrize(("which", "options"), FORM_CASES)
    def test_jacobian_has_the_same_sparsity(self, networks, which, options) -> None:
        net = networks[which]
        reduced = net.ode_system(**options, form="reduced")
        expanded = net.ode_system(**options)

        assert set(reduced.jacobian.todok()) == set(expanded.jacobian.todok())


PARAM_NETWORK = [
    "H -> H+ + e- [10,1e4] 1e-17*crate*exp(-2.5*av)",
    "H+ + e- -> H [10,1e4] 2e-10*tgas**(-0.7)",
    "H + H -> H2 [10,1e4] 3e-17*(tdust/100)**0.5",
]


@pytest.fixture(scope="module")
def param_network(tmp_path_factory):
    """Rates carrying free parameters: ``crate``, ``av``, ``tdust`` and ``tgas``."""
    path = tmp_path_factory.mktemp("ode") / "params.dat"
    path.write_text("\n".join(PARAM_NETWORK) + "\n", encoding="utf-8")
    return Network(str(path), funcfile=False)


def _names(system: OdeSystem) -> list[str]:
    return [str(p) for p in system.parameters]


class TestParameters:
    def test_free_symbols_sorted_by_name(self, param_network) -> None:
        system = param_network.ode_system()
        assert _names(system) == ["av", "crate", "tdust", "tgas"]

    def test_temperature_is_state_not_parameter_when_evolved(self, param_network) -> None:
        for thermal in ("dtdt", "dedt"):
            system = param_network.ode_system(thermal=thermal)
            assert _names(system) == ["av", "crate", "tdust"]

    def test_reduced_form_has_the_same_parameters(self, param_network) -> None:
        for thermal in ("none", "dtdt"):
            expanded = param_network.ode_system(thermal=thermal)
            reduced = param_network.ode_system(thermal=thermal, form="reduced")
            assert _names(reduced) == _names(expanded)

    def test_are_plain_symbols(self, param_network) -> None:
        parameters = param_network.ode_system(thermal="dtdt").parameters
        assert all(isinstance(p, sp.Symbol) for p in parameters)

    def test_are_the_network_variables_minus_the_state(self) -> None:
        net = Network("GOW")
        thermal = _names(net.ode_system(thermal="dtdt"))

        assert thermal == [
            "av", "chi", "crate", "d2g", "gradv",
            "ncol_C", "ncol_CO", "ncol_H", "ncol_H2", "vdisp",
        ]  # fmt: skip
        # network-wide: gradv, used only by the energy equation, is listed without it
        assert _names(net.ode_system()) == sorted([*thermal, "tgas"])


def _assembled(system: OdeSystem) -> sp.Matrix:
    """A + a c^T from the system's split Jacobian."""
    a_matrix, a, c = system.jacobian_split
    return sp.Matrix(a_matrix) + sp.Matrix(a) * sp.Matrix(c).T


SPLIT_CASES = [
    pytest.param("tgas", "specific", "expanded", id="tgas-specific"),
    pytest.param("tgas", "volumetric", "expanded", id="tgas-volumetric"),
    pytest.param("tgas", "specific", "reduced", id="tgas-specific-reduced"),
    pytest.param("rad", "specific", "expanded", id="rad-specific"),
]


class TestJacobianSplit:
    """``jacobian_split``: J = A + a c^T, A at fixed T, a = drhs/dT, c = dT/dstate."""

    @pytest.mark.parametrize(("which", "energy", "form"), SPLIT_CASES)
    def test_reassembles_the_jacobian(self, networks, which, energy, form) -> None:
        radiation = which == "rad"  # its photo rates read the radiation state
        system = networks[which].ode_system(
            thermal="dedt", energy=energy, form=form, radiation=radiation
        )
        if form == "reduced":
            definitions = system.rate_coefficients
            assembled = _assembled(system).xreplace(definitions)
            jacobian = system.jacobian.xreplace(definitions)
        else:
            assembled, jacobian = _assembled(system), system.jacobian

        for seed in range(3):
            values = _state_values(system, seed)
            _assert_close(_numbers(assembled, values), _numbers(jacobian, values))

    def test_fixed_temperature_part_has_no_energy_column(self, tgas_network) -> None:
        system = tgas_network.ode_system(thermal="dedt", energy="specific")
        a_matrix, _, _ = system.jacobian_split
        energy = tgas_network.species.count

        assert all(j != energy for _, j in a_matrix.todok())

    def test_vectors_span_the_state(self, rad_network) -> None:
        system = rad_network.ode_system(thermal="dedt", energy="specific", radiation=True)
        _, a, c = system.jacobian_split
        kinds = [var.kind for var in system.state]

        assert len(a) == len(c) == len(system.state)
        # T depends on the densities and the energy only
        assert all(
            c[j] == 0 for j, kind in enumerate(kinds) if kind.startswith("radiation")
        )

    def test_fixed_temperature_part_is_sparser(self) -> None:
        net = Network("COthin")
        system = net.ode_system(thermal="dedt", energy="specific")
        a_matrix, _, _ = system.jacobian_split

        assert len(a_matrix.todok()) < len(system.jacobian.todok())

    @pytest.mark.parametrize("thermal", ["none", "dtdt"])
    def test_only_for_an_evolved_energy(self, tgas_network, thermal) -> None:
        with pytest.raises(ValueError, match="dedt"):
            tgas_network.ode_system(thermal=thermal).jacobian_split


class TestPowerOfRate:
    """A rate raised to a state-dependent power: split before differentiating, the
    Jacobian holds no exp(+E/T) (which overflows at low T) and keeps its value."""

    RATE = "(1.18e-10*exp(-69500/tgas)*sqrt(tgas))**(1/(1 + tgas/1000))"

    @pytest.fixture
    def network(self, make_network):
        return make_network([f"H2 + H2 -> H + H + H2 [10,1e5] {self.RATE}"])

    @staticmethod
    def _expanded(system, expr):
        return expr.xreplace(system.rate_coefficients)

    @staticmethod
    def _growing(expr, tgas=90.0):
        t = sp.Symbol("tgas")
        return [
            a
            for a in expr.atoms(sp.exp)
            if a.args[0].free_symbols <= {t} and float(a.args[0].subs(t, tgas)) > 700
        ]

    @pytest.mark.parametrize("form", ["expanded", "reduced"])
    def test_jacobian_has_no_growing_exponential(self, network, form) -> None:
        system = network.ode_system(thermal="dtdt", form=form)
        entries = [self._expanded(system, e) for e in system.jacobian.values()]

        assert not [e for e in entries if self._growing(e)]

    @pytest.mark.parametrize("form", ["expanded", "reduced"])
    def test_jacobian_keeps_its_value(self, network, form) -> None:
        system = network.ode_system(thermal="dtdt", form=form)
        tgas, nden = sp.Symbol("tgas"), network.symbols.ndens
        rows = [r.xreplace(system.rate_coefficients) for r in system.rhs]
        plain = [
            r.diff(tgas)
            for r in network.kinetics.odes() + [network.thermodynamics.dTdt_tot()]
        ]
        point = {tgas: 400.0, nden[0]: 1e2, nden[1]: 3e1}
        point.update({x: 1.0 for e in plain for x in e.free_symbols if x not in point})
        column = len(system.state) - 1

        for i, want in enumerate(plain):
            got = self._expanded(system, system.jacobian[i, column])
            assert float(got.xreplace(point)) == pytest.approx(
                float(want.xreplace(point).evalf(30)), rel=1e-12
            )
        assert len(rows) == len(plain)

    def test_rhs_keeps_its_value(self, network) -> None:
        system = network.ode_system(thermal="dtdt")
        tgas, nden = sp.Symbol("tgas"), network.symbols.ndens
        point = {tgas: 400.0, nden[0]: 1e2, nden[1]: 3e1}
        for got, want in zip(system.rhs, network.kinetics.odes(), strict=False):
            point.update({x: 1.0 for x in want.free_symbols if x not in point})
            assert float(got.xreplace(point)) == pytest.approx(
                float(want.xreplace(point)), rel=1e-12
            )


def test_rate_coefficient_derivatives_are_reported_first(make_network, monkeypatch):
    import jaff.core.network.ode_system as ode_system_module

    descriptions = []

    def track(sequence, *args, description="Working...", **kwargs):
        descriptions.append(description)
        yield from sequence

    monkeypatch.setattr(ode_system_module.jaff_progress, "track", track)
    network = make_network(["H + H -> H2 [10,1e4] 1e-10*tgas**0.7*exp(-5e4/tgas)"])
    _ = network.ode_system(thermal="dtdt", form="reduced").jacobian
    rates = descriptions.index("Differentiating the rate coefficients")

    assert rates < descriptions.index("Differentiating the rhs")
