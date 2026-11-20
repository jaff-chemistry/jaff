# ABOUTME: The kokkos_ode template (header-only-integrators VODE) as a standalone solver:
# ABOUTME: generated for real networks it builds, its Jacobian matches and VODE integrates
"""End-to-end tests of the generated C++ (``kokkos_ode``).

Each network is generated and configured with CMake, which fetches
header-only-integrators (or uses the copy named by the environment variable
``FETCHCONTENT_SOURCE_DIR_HEADER_ONLY_INTEGRATORS``). Small drivers compiled with
clang++ against the generated header and the library call ``ChemistryODE`` and
VODE. The state is the species number densities [cm^-3] followed by the gas
temperature [K]. Unknown functions (GOW's tables) return 1, their partial
derivatives 0, by rewriting the generated stubs.
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
from tests import kokkos_build
from tests._paths import FIXTURES_DIR, REPO_ROOT
from tests.python_solver import generate

pytestmark = pytest.mark.skipif(
    not sys.platform.startswith("linux")
    or shutil.which("cmake") is None
    or shutil.which("clang++") is None,
    reason="needs Linux with cmake and clang++ (and network access, or a local "
    "header-only-integrators in FETCHCONTENT_SOURCE_DIR_HEADER_ONLY_INTEGRATORS)",
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
    """C++ statements setting every parameter of ``ChemistryODE::params``."""
    names = net.ode_system(thermal="dtdt").parameters
    return "\n".join(
        f"    ChemistryODE::params.{name} = {_value(str(name))!r};" for name in names
    )


def _neutral_start(net: Network) -> str:
    """Every species at 1 cm^-3, electrons equal to the total ion charge."""
    electron = net.species.e_idx()
    charge = sum(s.charge for s in net.species if s.charge > 0)
    values = [float(charge) if i == electron else 1.0 for i in range(net.species.count)]

    return "\n".join(f"    y[{i}] = {v!r};" for i, v in enumerate(values))


def _run(built: SimpleNamespace, name: str, body: str) -> subprocess.CompletedProcess:
    """Compile a driver whose main() is *body*; run it in the generated folder."""
    source = built.folder / f"{name}.cpp"
    source.write_text(
        "#include <cstdio>\n#include <string>\n#include <vector>\n"
        '#include "chemistry_ode.hpp"\nusing namespace integrators;\n\n'
        "int main() {\n" + body + "\n}\n",
        encoding="utf-8",
    )
    exe = built.folder / name
    subprocess.run(
        ["clang++", "-std=c++20", "-O2", f"-I{built.folder}", f"-I{built.include}"]
        + [str(source), "-o", str(exe)],
        check=True,
        capture_output=True,
        text=True,
    )

    return subprocess.run([str(exe)], cwd=built.folder, capture_output=True, text=True)


def _numbers(text: str) -> np.ndarray:
    return np.array([float(v) for v in text.split()])


@pytest.fixture(scope="module", params=list(NETWORKS))
def built(request, tmp_path_factory) -> SimpleNamespace:
    network = NETWORKS[request.param]
    net = Network(network)
    folder = generate(network, tmp_path_factory.mktemp(request.param), "kokkos_ode")
    _stub_unknown_functions(folder)
    configure = kokkos_build.configure(folder)
    build = subprocess.run(
        ["cmake", "--build", "build", "-j"], cwd=folder, capture_output=True, text=True
    )
    include = kokkos_build.include_dir(folder)

    return SimpleNamespace(
        net=net, folder=folder, configure=configure, build=build, include=include
    )


def test_writes_the_solver_files(built):
    for name in ("chemistry_ode.hpp", "user_functions.hpp", "chemistry_ode.cpp"):
        assert (built.folder / name).is_file()


def test_builds_with_clang_by_default(built):
    assert built.configure.returncode == 0, (
        built.configure.stdout + built.configure.stderr
    )
    assert built.build.returncode == 0, built.build.stdout + built.build.stderr
    compiler = next(
        (built.folder / "build" / "CMakeFiles").glob("*/CMakeCXXCompiler.cmake")
    )
    assert 'set(CMAKE_CXX_COMPILER_ID "Clang")' in compiler.read_text(encoding="utf-8")


def test_species_and_parameters_by_name(built):
    names = built.net.species.names()
    body = """
    std::printf("%d %d\\n", ChemistryODE::nspecs, static_cast<int>(ChemistryODE::neqs));
    for (int i = 0; i < ChemistryODE::nspecs; ++i)
        std::printf("%s\\n", std::string(ChemistryODE::species[i]).c_str());
    for (auto name : ChemistryODE::unset_parameters())
        std::printf("unset %s\\n", std::string(name).c_str());
"""
    lines = _run(built, "names", body).stdout.splitlines()

    assert lines[0] == f"{len(names)} {len(names) + 1}"
    assert lines[1 : 1 + len(names)] == names
    params = built.net.ode_system(thermal="dtdt").parameters
    assert lines[1 + len(names) :] == [f"unset {p}" for p in params]


def test_jacobian_matches_finite_differences(built):
    nspec = built.net.species.count
    body = f"""
{_parameters(built.net)}
    ChemistryODE::evolve_temperature = true;
    ChemistryODE::state_type y{{}}, up{{}}, down{{}};
    ChemistryODE::rhs_type fu{{}}, fd{{}};
    for (int i = 0; i < ChemistryODE::nspecs; ++i) y[i] = 1e2 * (1.0 + 0.37 * i);
    y[ChemistryODE::nspecs] = {TEMPERATURE!r};
    ChemistryODE::jacobian_type jac{{}};
    ChemistryODE::jacobian(0.0, y, jac);
    for (int i = 0; i < ChemistryODE::neqs; ++i) {{
        for (int j = 0; j < ChemistryODE::neqs; ++j) std::printf("%.17e ", jac[i][j]);
        std::printf("\\n");
    }}
    for (int i = 0; i < ChemistryODE::neqs; ++i) {{
        for (int j = 0; j < ChemistryODE::neqs; ++j) {{
            const double h = 1e-6 * std::abs(y[j]);
            up = y;
            down = y;
            up[j] += h;
            down[j] -= h;
            ChemistryODE::rhs(0.0, up, fu);
            ChemistryODE::rhs(0.0, down, fd);
            std::printf("%.17e ", (fu[i] - fd[i]) / (2 * h));
        }}
        std::printf("\\n");
    }}
"""
    jac, fd = _numbers(_run(built, "check_jacobian", body).stdout).reshape(
        2, nspec + 1, nspec + 1
    )

    for row in range(nspec + 1):
        scale = max(np.abs(fd[row]).max(), np.abs(jac[row]).max(), 1e-300)
        assert jac[row] == pytest.approx(fd[row], rel=1e-5, abs=1e-7 * scale)


def _has_unknown_functions(net: Network) -> bool:
    return bool(net.symbols.applied_functions)


@pytest.mark.parametrize("evolve", [False, True])
def test_vode_integrates(built, evolve):
    if evolve and _has_unknown_functions(built.net):
        # Stubbed tables (= 1) make GOW cool at ~1e14 K/s: T -> 0 within 1e-11 s.
        pytest.skip("temperature evolution needs the real tables (benchmarks)")
    flag = "true" if evolve else "false"
    body = f"""
{_parameters(built.net)}
    ChemistryODE::evolve_temperature = {flag};
    ChemistryODE::state_type y{{}};
{_neutral_start(built.net)}
    y[ChemistryODE::nspecs] = {TEMPERATURE!r};
    VODE<ChemistryODE> integrator{{}};
    VODEState<ChemistryODE::neqs> state{{}};
    state.t = 0.0;
    state.tout = 1e8;
    state.y = y;
    state.rtol = 1e-8;  // this VODE stalls below rtol ~ 1e-10 on stiff networks
    state.atol = 1e-20;
    state.max_steps = 1000000;
    state.jacobian_analytic = true;
    ChemistryODE::state_type problem{{}};
    const auto result = integrator.integrate(problem, state);
    std::printf("%d", static_cast<int>(result));
    for (int i = 0; i < ChemistryODE::neqs; ++i) std::printf(" %.17e", state.y[i]);
"""
    status, *values = _numbers(_run(built, "integrate", body).stdout)

    assert status == 1  # IntegratorResult::SUCCESS
    assert np.all(np.isfinite(values))
    if not evolve:
        assert values[-1] == pytest.approx(TEMPERATURE, rel=1e-12)


def test_generated_driver_runs(built):
    # The driver's placeholder inputs (every parameter 1, 100 K) are not physical:
    # GOW needs its tables, and COthin overflows below ~200 K (kl11 ~ 1e-321).
    if built.net.label != "demo2":
        pytest.skip("the generated driver's placeholder inputs suit demo2 only")

    proc = subprocess.run(
        ["./build/chemistry_ode"], cwd=built.folder, capture_output=True, text=True
    )

    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "net_charge:" in proc.stdout
    assert "evolution saved to solution.dat" in proc.stdout
    lines = (built.folder / "solution.dat").read_text().splitlines()
    assert lines[0].split() == ["#", "t[s]", "T[K]", *built.net.species.names()]
    assert len(lines) == 1 + 101
