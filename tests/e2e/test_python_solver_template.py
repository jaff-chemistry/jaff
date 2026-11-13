# ABOUTME: The python_solve_ivp template as a standalone solver: generated for real
# ABOUTME: networks, its Jacobian matches the rhs and solve() integrates with Radau

import importlib
import re
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace

import numpy as np
import pytest

from jaff import Network
from tests._paths import REPO_ROOT
from tests.jaffgen_run import run_jaffgen

NETWORKS = {
    "demo2": str(REPO_ROOT / "networks" / "demos" / "demo2.jet"),
    "COthin": "COthin",
    "GOW": "GOW",
}

# Parameter values by name (1 for any other); unknown functions return 1.
VALUES = {
    "av": 1.0,
    "chi": 1.0,
    "crate": 1e-16,
    "d2g": 0.01,
    "gradv": 1e-14,
    "tdust": 20.0,
    "vdisp": 1e5,
}
NCOL = 1e14
# COthin's kl11 = 10**(-27.029 + 3.801 logT - 29487/T) underflows below ~200 K
TEMPERATURE = 1000.0


def _generate(network: str, outdir: Path) -> Path:
    args = SimpleNamespace(
        network=network,
        config=None,
        label=None,
        funcfile=None,
        duplicate_policy=None,
        expand_nuclei=None,
        strict=None,
        network_config=None,
        outdir=str(outdir),
        indir=None,
        files=None,
        template="python_solve_ivp",
    )

    return run_jaffgen(args).output_dir


def _load(folder: Path) -> SimpleNamespace:
    """Import the generated ``chemistry`` and ``solver`` modules of *folder* afresh."""
    for name in ("chemistry", "solver"):
        sys.modules.pop(name, None)
    sys.path.insert(0, str(folder))
    try:
        chemistry = importlib.import_module("chemistry")
        solver = importlib.import_module("solver")
    finally:
        sys.path.remove(str(folder))
        for name in ("chemistry", "solver"):
            sys.modules.pop(name, None)

    return SimpleNamespace(chemistry=chemistry, solver=solver)


def _stub_unknown_functions(chemistry: ModuleType, net: Network, folder: Path) -> None:
    """Unknown functions return 1, so their ``<name>_partial_<n>`` derivatives 0."""
    source = (folder / "chemistry.py").read_text(encoding="utf-8")
    for name in net.symbols.interp_functions | net.symbols.undefined_functions:
        setattr(chemistry, name, lambda *args: 1.0)
        for partial in set(re.findall(rf"\b{name}_partial_\w+(?=\()", source)):
            setattr(chemistry, partial, lambda *args: 0.0)


def _params(chemistry: ModuleType) -> dict[str, float]:
    return {
        p: NCOL if p.startswith("ncol_") else VALUES.get(p, 1.0)
        for p in chemistry.PARAMETERS
    }


@pytest.fixture(scope="module", params=list(NETWORKS))
def generated(request, tmp_path_factory) -> SimpleNamespace:
    network = NETWORKS[request.param]
    net = Network(network)
    folder = _generate(network, tmp_path_factory.mktemp(request.param))
    modules = _load(folder)
    _stub_unknown_functions(modules.chemistry, net, folder)

    return SimpleNamespace(net=net, folder=folder, **vars(modules))


def _state(nspec: int) -> tuple[np.ndarray, float]:
    return 1.0e2 * (1.0 + 0.37 * np.arange(nspec)), TEMPERATURE


def test_writes_the_solver_files(generated):
    assert (generated.folder / "chemistry.py").is_file()
    assert (generated.folder / "solver.py").is_file()


def test_species_map_has_the_network_names(generated):
    assert generated.chemistry.SPECIES == {
        name: i for i, name in enumerate(generated.net.species.names())
    }
    assert generated.chemistry.NSPEC == generated.net.species.count


def test_parameters_are_the_dtdt_system_ones(generated):
    expected = generated.net.ode_system(thermal="dtdt").parameters
    assert generated.chemistry.PARAMETERS == tuple(map(str, expected))


def test_rhs_has_a_row_per_species_and_the_temperature(generated):
    chem = generated.chemistry
    nden, tgas = _state(chem.NSPEC)
    dydt = chem.rhs(nden, tgas, _params(chem))

    assert dydt.shape == (chem.NSPEC + 1,)
    assert np.all(np.isfinite(dydt))


def test_internal_energy_is_positive(generated):
    chem = generated.chemistry
    nden, tgas = _state(chem.NSPEC)

    assert chem.internal_energy(nden, tgas) > 0.0


def test_jacobian_matches_finite_differences(generated):
    chem = generated.chemistry
    params = _params(chem)
    nden, tgas = _state(chem.NSPEC)
    y = np.append(nden, tgas)

    def f(state: np.ndarray) -> np.ndarray:
        return chem.rhs(state[:-1], state[-1], params)

    jac = chem.jacobian(nden, tgas, params)
    assert jac.shape == (chem.NSPEC + 1, chem.NSPEC + 1)

    fd = np.empty_like(jac)
    for j in range(y.size):
        h = 1e-6 * abs(y[j])
        up, down = y.copy(), y.copy()
        up[j] += h
        down[j] -= h
        fd[:, j] = (f(up) - f(down)) / (2 * h)

    for row in range(y.size):
        scale = max(np.abs(fd[row]).max(), np.abs(jac[row]).max(), 1e-300)
        assert jac[row] == pytest.approx(fd[row], rel=1e-5, abs=1e-7 * scale)


def _neutral_start(net: Network) -> dict[str, float]:
    """Every species at 1 cm^-3, electrons set to the total ion charge."""
    n0 = {s.name: 1.0 for s in net.species}
    electron = net.species.e_idx()
    if electron is not None:
        charge = sum(s.charge for s in net.species if s.charge > 0)
        n0[net.species.names()[electron]] = float(charge)

    return n0


def _has_unknown_functions(net: Network) -> bool:
    return bool(net.symbols.interp_functions | net.symbols.undefined_functions)


@pytest.mark.parametrize("evolve", [False, True])
def test_solve_integrates_with_radau(generated, evolve):
    if evolve and _has_unknown_functions(generated.net):
        # Stubbed tables (= 1) make GOW cool at ~1e14 K/s: T -> 0 within 1e-11 s.
        pytest.skip("temperature evolution needs the real tables (benchmarks)")
    chem, solver = generated.chemistry, generated.solver
    n0 = _neutral_start(generated.net)

    result = solver.solve(
        n0, TEMPERATURE, 1.0e8, _params(chem), evolve_temperature=evolve
    )

    assert result.solution.success
    assert set(result.n) == set(chem.SPECIES)
    assert all(np.all(np.isfinite(v)) for v in result.n.values())
    assert np.all(np.isfinite(result.temperature))
    if not evolve:
        assert np.all(result.temperature == TEMPERATURE)


def test_solve_uses_tight_radau_defaults(generated):
    chem, solver = generated.chemistry, generated.solver
    n0 = _neutral_start(generated.net)
    result = solver.solve(n0, TEMPERATURE, 1.0e3, _params(chem))

    assert result.solution.success
    assert solver.DEFAULT_RTOL == pytest.approx(1e-12)
    assert solver.DEFAULT_ATOL == pytest.approx(1e-20)


def test_missing_and_unknown_parameters_are_named(generated):
    chem, solver = generated.chemistry, generated.solver
    params = _params(chem)
    missing = chem.PARAMETERS[0]
    del params[missing]
    params["bogus"] = 1.0

    with pytest.raises(ValueError, match=f"{missing}.*bogus|bogus.*{missing}"):
        solver.solve({}, 100.0, 1.0, params)


def test_unknown_species_is_named(generated):
    chem, solver = generated.chemistry, generated.solver

    with pytest.raises(ValueError, match="Xq"):
        solver.solve({"Xq": 1.0}, 100.0, 1.0, _params(chem))
