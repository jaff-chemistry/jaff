# ABOUTME: The rosenbrock_cpp template as a header-only integrator for host codes:
# ABOUTME: generated for real networks it builds, its Jacobian matches and solve() works
"""End-to-end tests of the generated C++ (``rosenbrock_cpp``).

Each network is generated; small C++ drivers include ``solver.hpp`` and call the
``chem`` API. The state is the species number densities [cm^-3] followed by the
specific internal energy e [erg g^-1]. Unknown functions (GOW's tables) are made
to return 1, and their partial derivatives 0, by rewriting the generated stubs.
"""

import re
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from jaff import Network
from tests._paths import FIXTURES_DIR, REPO_ROOT
from tests.python_solver import generate

pytestmark = pytest.mark.skipif(
    not sys.platform.startswith("linux") or shutil.which("g++") is None,
    reason="needs Linux with g++",
)

NETWORKS = {
    "demo2": str(REPO_ROOT / "networks" / "demos" / "demo2.jet"),
    "COthin": "COthin",
    "GOW": str(FIXTURES_DIR / "GOW_benchmark" / "GOW.jet"),
}

# Parameter values by name (1 for any other); as in the other template tests
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
TEMPERATURE = 1000.0
CXXFLAGS = ["-std=c++17", "-O2", "-Wall"]


def _stub_unknown_functions(folder: Path) -> None:
    """Make the generated stubs return 1 (functions) and 0 (partial derivatives)."""
    path = folder / "user_functions.hpp"
    lines = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.lstrip().startswith("inline double "):
            value = "0.0" if "_partial_" in line.split("(")[0] else "1.0"
            line = re.sub(r"\{.*\}\s*$", f"{{ return {value}; }}", line)
        lines.append(line)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _value(name: str) -> float:
    return NCOL if name.startswith("ncol_") else VALUES.get(name, 1.0)


def _parameters(net: Network) -> str:
    """C++ statements setting every parameter of ``chem::Parameters p``."""
    names = net.ode_system(thermal="dedt", energy="specific").parameters
    return "\n".join(f"    p.{name} = {_value(str(name))!r};" for name in names)


def _neutral_start(net: Network) -> str:
    """Every species at 1 cm^-3, electrons equal to the total ion charge."""
    electron = net.species.e_idx()
    charge = sum(s.charge for s in net.species if s.charge > 0)
    values = [float(charge) if i == electron else 1.0 for i in range(net.species.count)]

    return "\n".join(f"    n[{i}] = {v!r};" for i, v in enumerate(values))


def _run(folder: Path, name: str, body: str) -> subprocess.CompletedProcess:
    """Compile a driver whose main() is *body* against the generated headers; run it."""
    source = folder / f"{name}.cpp"
    source.write_text(
        '#include <cstdio>\n#include <string>\n#include "solver.hpp"\n\n'
        "int main() {\n" + body + "\n}\n",
        encoding="utf-8",
    )
    exe = folder / name
    subprocess.run(
        ["g++", *CXXFLAGS, f"-I{folder}", str(source), "-o", str(exe)],
        check=True,
        capture_output=True,
        text=True,
    )

    return subprocess.run([str(exe)], cwd=folder, capture_output=True, text=True)


def _numbers(text: str) -> np.ndarray:
    return np.array([float(v) for v in text.split()])


@pytest.fixture(scope="module", params=list(NETWORKS))
def generated(request, tmp_path_factory) -> SimpleNamespace:
    network = NETWORKS[request.param]
    net = Network(network)
    folder = generate(network, tmp_path_factory.mktemp(request.param), "rosenbrock_cpp")
    _stub_unknown_functions(folder)

    return SimpleNamespace(net=net, folder=folder)


def test_writes_the_solver_headers(generated):
    for name in ("chemistry.hpp", "user_functions.hpp", "solver.hpp", "tableau.hpp"):
        assert (generated.folder / name).is_file()
    for name in ("parameters.hpp", "rhs.hpp", "jac.hpp", "state.hpp", "eos.hpp"):
        assert not (generated.folder / name).exists()


def test_species_and_parameters_by_name(generated):
    names = generated.net.species.names()
    body = """
    std::printf("%d %d\\n", chem::nspecs, chem::neq);
    for (int i = 0; i < chem::nspecs; ++i) std::printf("%s %d\\n",
        std::string(chem::species[i]).c_str(), chem::species_index(chem::species[i]));
    chem::Parameters p;
    for (auto name : chem::unset_parameters(p)) std::printf("unset %s\\n",
        std::string(name).c_str());
"""
    proc = _run(generated.folder, "names", body)
    lines = proc.stdout.splitlines()

    assert lines[0] == f"{len(names)} {len(names) + 1}"
    assert lines[1 : 1 + len(names)] == [f"{n} {i}" for i, n in enumerate(names)]
    params = generated.net.ode_system(thermal="dedt", energy="specific").parameters
    assert lines[1 + len(names) :] == [f"unset {p}" for p in params]


def test_energy_and_temperature_round_trip(generated):
    body = f"""
    double n[chem::nspecs];
{_neutral_start(generated.net)}
    const double e = chem::energy(n, {TEMPERATURE!r});
    std::printf("%.17e %.17e\\n", e, chem::temperature(n, e));
"""
    e, temperature = _numbers(_run(generated.folder, "round_trip", body).stdout)

    assert e > 0.0
    assert temperature == pytest.approx(TEMPERATURE, rel=1e-13)


def test_jacobian_matches_finite_differences(generated):
    nspec = generated.net.species.count
    body = f"""
    chem::Parameters p;
{_parameters(generated.net)}
    double y[chem::neq], up[chem::neq], down[chem::neq];
    double fu[chem::neq], fd[chem::neq], f[chem::neq];
    for (int i = 0; i < chem::nspecs; ++i) y[i] = 1e2 * (1.0 + 0.37 * i);
    y[chem::nspecs] = chem::energy(y, {TEMPERATURE!r});
    chem::Matrix jac{{}};
    chem::rhs_jac(y, p, f, jac);
    for (int i = 0; i < chem::neq; ++i) {{
        for (int j = 0; j < chem::neq; ++j) std::printf("%.17e ", jac[i][j]);
        std::printf("\\n");
    }}
    for (int i = 0; i < chem::neq; ++i) {{
        for (int j = 0; j < chem::neq; ++j) {{
            const double h = 1e-6 * std::abs(y[j]);
            for (int k = 0; k < chem::neq; ++k) up[k] = down[k] = y[k];
            up[j] += h;
            down[j] -= h;
            chem::rhs(up, p, fu);
            chem::rhs(down, p, fd);
            std::printf("%.17e ", (fu[i] - fd[i]) / (2 * h));
        }}
        std::printf("\\n");
    }}
"""
    proc = _run(generated.folder, "check_jacobian", body)
    jac, fd = _numbers(proc.stdout).reshape(2, nspec + 1, nspec + 1)

    for row in range(nspec + 1):
        scale = max(np.abs(fd[row]).max(), np.abs(jac[row]).max(), 1e-300)
        assert jac[row] == pytest.approx(fd[row], rel=1e-5, abs=1e-7 * scale)


def _has_unknown_functions(net: Network) -> bool:
    return bool(net.symbols.applied_functions)


@pytest.mark.parametrize("evolve", [False, True])
def test_solve_advances_a_cell(generated, evolve):
    if evolve and _has_unknown_functions(generated.net):
        # Stubbed tables (= 1) make GOW cool at ~1e14 K/s: T -> 0 within 1e-11 s.
        pytest.skip("energy evolution needs the real tables (benchmarks)")
    flag = "true" if evolve else "false"
    body = f"""
    chem::Parameters p;
{_parameters(generated.net)}
    double n[chem::nspecs];
{_neutral_start(generated.net)}
    double e = chem::energy(n, {TEMPERATURE!r});
    const double e0 = e;
    chem::Options options;
    options.evolve_energy = {flag};
    const chem::Status status = chem::solve(n, e, 1e8, p, options);
    std::printf("%d %.17e %.17e", static_cast<int>(status), e0, e);
    for (int i = 0; i < chem::nspecs; ++i) std::printf(" %.17e", n[i]);
"""
    proc = _run(generated.folder, "run_solve", body)
    status, e0, e, *densities = _numbers(proc.stdout)

    assert status == 0, proc.stdout
    assert np.all(np.isfinite(densities))
    if not evolve:
        # the specific energy is held; T may still move as the composition changes
        assert e == e0


def test_unset_parameters_are_reported(generated):
    body = """
    chem::Parameters p;
    double n[chem::nspecs];
    for (int i = 0; i < chem::nspecs; ++i) n[i] = 1.0;
    double e = chem::energy(n, 100.0);
    std::printf("%d\\n", static_cast<int>(chem::solve(n, e, 1.0, p)));
"""
    status = _run(generated.folder, "unset", body).stdout.split()[0]

    assert int(status) == int(_status(generated.folder, "unset_parameters"))


def test_failed_solve_leaves_the_input_unchanged(generated):
    body = f"""
    chem::Parameters p;
{_parameters(generated.net)}
    double n[chem::nspecs];
{_neutral_start(generated.net)}
    double before[chem::nspecs];
    for (int i = 0; i < chem::nspecs; ++i) before[i] = n[i];
    double e = chem::energy(n, {TEMPERATURE!r});
    const double e_before = e;
    chem::Options options;
    options.max_steps = 1;
    options.h_init = 1e-6;
    const chem::Status status = chem::solve(n, e, 1e8, p, options);
    bool same = e == e_before;
    for (int i = 0; i < chem::nspecs; ++i) same = same && n[i] == before[i];
    std::printf("%d %d\\n", static_cast<int>(status), same ? 1 : 0);
"""
    status, same = _numbers(_run(generated.folder, "fail", body).stdout)

    assert status == _status(generated.folder, "too_many_steps")
    assert same == 1


def _status(folder: Path, name: str) -> int:
    """Value of ``chem::Status::<name>``."""
    body = f'    std::printf("%d\\n", static_cast<int>(chem::Status::{name}));'
    return int(_run(folder, f"status_{name}", body).stdout)


def test_sparse_lu_is_the_default(generated):
    body = """
    std::printf("%d\\n", chem::Options{}.lu == chem::LinearSolver::sparse ? 1 : 0);
"""
    assert _run(generated.folder, "default_lu", body).stdout.split() == ["1"]


def test_pattern_covers_every_jacobian_entry(generated):
    body = f"""
    chem::Parameters p;
{_parameters(generated.net)}
    double y[chem::neq], f[chem::neq];
    for (int i = 0; i < chem::nspecs; ++i) y[i] = 1e2 * (1.0 + 0.37 * i);
    y[chem::nspecs] = chem::energy(y, {TEMPERATURE!r});
    chem::Matrix jac{{}};
    chem::rhs_jac(y, p, f, jac);
    bool in_pattern[chem::neq][chem::neq] = {{}};
    for (const auto &entry : chem::jacobian_pattern) in_pattern[entry[0]][entry[1]] = true;
    int missing = 0;
    for (int i = 0; i < chem::neq; ++i)
        for (int j = 0; j < chem::neq; ++j)
            if (jac[i][j] != 0.0 && !in_pattern[i][j]) ++missing;
    std::printf("%d\\n", missing);
"""
    assert _run(generated.folder, "pattern", body).stdout.split() == ["0"]


@pytest.mark.parametrize("evolve", [False, True])
def test_sparse_and_dense_lu_agree(generated, evolve):
    if evolve and _has_unknown_functions(generated.net):
        pytest.skip("energy evolution needs the real tables (benchmarks)")
    flag = "true" if evolve else "false"
    body = f"""
    chem::Parameters p;
{_parameters(generated.net)}
    for (const chem::LinearSolver lu : {{chem::LinearSolver::sparse,
                                         chem::LinearSolver::dense}}) {{
        double n[chem::nspecs];
{_neutral_start(generated.net)}
        double e = chem::energy(n, {TEMPERATURE!r});
        chem::Options options;
        options.evolve_energy = {flag};
        options.lu = lu;
        const chem::Status status = chem::solve(n, e, 1e8, p, options);
        std::printf("%d %.17e", static_cast<int>(status), e);
        for (int i = 0; i < chem::nspecs; ++i) std::printf(" %.17e", n[i]);
        std::printf("\\n");
    }}
"""
    out = _run(generated.folder, "sparse_dense", body).stdout.splitlines()
    sparse, dense = (_numbers(line) for line in out)

    assert sparse[0] == 0 and dense[0] == 0, out
    assert sparse[1:] == pytest.approx(dense[1:], rel=1e-8, abs=1e-30)


def test_generated_driver_builds_and_runs(generated):
    if shutil.which("make") is None:
        pytest.skip("needs make")
    build = subprocess.run(["make"], cwd=generated.folder, capture_output=True, text=True)
    assert build.returncode == 0, build.stdout + build.stderr
    # The driver's placeholder inputs (every parameter 1, 100 K) are not physical:
    # GOW needs its tables, and COthin overflows below ~200 K (kl11 ~ 1e-321).
    if generated.net.label != "demo2":
        pytest.skip("the generated driver's placeholder inputs suit demo2 only")

    proc = subprocess.run(["./test"], cwd=generated.folder, capture_output=True, text=True)

    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "evolution saved to solution.dat" in proc.stdout
    lines = (generated.folder / "solution.dat").read_text().splitlines()
    assert lines[0].split() == ["#", "t[s]", "T[K]", *generated.net.species.names()]
    assert len(lines) == 1 + 101
