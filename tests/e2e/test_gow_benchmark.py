# ABOUTME: Benchmark: the generated Python, Fortran, C++, Microphysics and FLASH solvers: GOW
# ABOUTME: to equilibrium at n_H = 100, chi = 0.1, 1, 10 and match the reference to 1e-6
"""GOW one-zone equilibrium against reference runs.

The network is ``tests/fixtures/GOW_benchmark``: ``networks/GOW`` without the
He + CR heating term (``deltaE2``), as in the reference runs. The references
(``reference/state_over_time_n100_00_chi*_met1_00.txt``, Microphysics runs of the
same network) start from the Gong+17 abundances at n_H = 100 cm^-3 and 100 K and
reach equilibrium by ~2e13 s; their last row (7e20 s) is the equilibrium state.

The RHS and Jacobian are generated from the network by the current code, for
the ``python_solve_ivp``, ``fortran_dlsodes``, ``rosenbrock_cpp``, ``kokkos_ode``,
``microphysics`` and ``flash`` templates; the test supplies only the table functions
(:mod:`tests.gow_tables`, and the same tables in
``tests/fixtures/GOW_benchmark/user_functions.f90`` and ``.hpp``), the start state
(the reference's first row) and the parameters below. The Rosenbrock integrator
evolves the specific internal energy instead of the temperature; VODE
(``kokkos_ode``) evolves the temperature, at rtol 1e-8 (it stalls below ~1e-10);
Microphysics evolves the specific energy with the generated EOS (rtol 1e-8, one
burn in the one-zone test, ``tests/fixtures/microphysics_cell``); FLASH burns the
generated network with its DLSODES once per FLASH step (rtol 1e-8), in a temporary
overlay of the checkout named by ``FLASH_HOME`` (:mod:`tests.flash_overlay`).
The gas sits at the cloud face: no extinction, all column densities 0, so every
shielding factor is 1.
"""

import re
import shutil
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from tests import flash_overlay, gow_tables, kokkos_build, microphysics_cell
from tests._paths import FIXTURES_DIR
from tests.python_solver import generate, load

BENCHMARK = FIXTURES_DIR / "GOW_benchmark"

#: Reference file name tag -> FUV field chi (Draine units)
CHI = {"0_10": 0.1, "1_00": 1.0, "10_00": 10.0}

#: Reference column -> species name in the network
SPECIES = {
    "H": "H", "Hp": "H+", "E": "e-", "H2": "H2", "H2p": "H2+", "HE": "He",
    "HEp": "He+", "Carbon": "C", "Cp": "C+", "CO": "CO", "HCOp": "HCO+", "O": "O",
    "SI": "Si", "SIp": "Si+", "CH": "CHx", "OH": "OHx", "H3p": "H3+", "Op": "O+",
}  # fmt: skip

#: End of the integration [s]: 5x the slowest approach to equilibrium (chi = 10)
T_END = 1.0e14

RTOL = 1e-12
ATOL = 1e-20
#: Required agreement with the reference equilibrium (relative)
AGREEMENT = 1e-6


def _parameters(chi: float) -> dict[str, float]:
    return {
        "av": 0.0,
        "chi": chi,
        "crate": 2e-16,
        "d2g": 0.01,
        "gradv": 3.0 * 3.0e-14,
        "vdisp": 3.0 * 1.0e5 / 2**0.5,
        "ncol_C": 0.0,
        "ncol_CO": 0.0,
        "ncol_H": 0.0,
        "ncol_H2": 0.0,
    }


def _reference(tag: str) -> dict[str, np.ndarray]:
    """The reference run for *tag* as ``{column: values over time}``."""
    path = BENCHMARK / "reference" / f"state_over_time_n100_00_chi{tag}_met1_00.txt"
    header = path.read_text(encoding="utf-8").splitlines()[0].replace("#", "").split()
    data = np.loadtxt(path)

    return {name: data[:, i] for i, name in enumerate(header)}


@pytest.fixture(scope="module")
def solver(tmp_path_factory) -> SimpleNamespace:
    folder = generate(str(BENCHMARK / "GOW.jet"), tmp_path_factory.mktemp("gow"))
    modules = load(folder)
    gow_tables.install(modules.chemistry, folder)

    return modules


@pytest.mark.parametrize("tag", list(CHI))
def test_equilibrium_matches_the_reference(solver, tag):
    ref = _reference(tag)
    n0 = {name: ref[column][0] for column, name in SPECIES.items()}

    result = solver.solver.solve(
        n0,
        ref["Tgas"][0],
        T_END,
        _parameters(CHI[tag]),
        evolve_temperature=True,
        t_eval=[T_END],
        rtol=RTOL,
        atol=ATOL,
    )

    assert result.temperature[-1] == pytest.approx(ref["Tgas"][-1], rel=AGREEMENT)
    for column, name in SPECIES.items():
        assert result.n[name][-1] == pytest.approx(ref[column][-1], rel=AGREEMENT), name


# --------------------------------------------------------------------------- #
# Fortran (fortran_dlsodes)                                                    #
# --------------------------------------------------------------------------- #

needs_gfortran = pytest.mark.skipif(
    not sys.platform.startswith("linux")
    or shutil.which("gfortran") is None
    or shutil.which("make") is None,
    reason="needs Linux with gfortran and make",
)

FFLAGS = ["-ffree-line-length-none", "-O2", "-fallow-argument-mismatch"]


def _real(value: float) -> str:
    """*value* as a Fortran double-precision literal."""
    return f"{value:.17e}".replace("e", "d")


def _fortran_driver(tag: str, species: list[str]) -> str:
    """A program solving the benchmark for *tag* and printing T then the densities."""
    ref = _reference(tag)
    start = {name: ref[column][0] for column, name in SPECIES.items()}
    densities = "\n".join(
        f"    n({i + 1}) = {_real(start[name])}  ! {name}"
        for i, name in enumerate(species)
    )
    parameters = "\n".join(
        f"    {name} = {_real(value)}" for name, value in _parameters(CHI[tag]).items()
    )

    return f"""
program benchmark
    use chemistry
    use solver
    implicit none
    real*8 :: n(nspecs), temperature
{densities}
{parameters}
    temperature = {_real(ref["Tgas"][0])}
    evolve_temperature = .true.
    call solve(n, temperature, {_real(T_END)}, "benchmark.dat", noutput=2, &
               rtol={_real(RTOL)}, atol={_real(ATOL)})
    print '(*(es25.16e3))', temperature, n
end program benchmark
"""


@pytest.fixture(scope="module")
def fortran(tmp_path_factory) -> SimpleNamespace:
    """The generated Fortran, built with the benchmark's table functions."""
    folder = generate(
        str(BENCHMARK / "GOW.jet"), tmp_path_factory.mktemp("gow_f90"), "fortran_dlsodes"
    )
    shutil.copy(BENCHMARK / "user_functions.f90", folder / "user_functions.f90")
    subprocess.run(["make"], cwd=folder, check=True, capture_output=True, text=True)
    species = [name for _, name in sorted((i, n) for n, i in _species_index(folder))]

    return SimpleNamespace(folder=folder, species=species)


def _species_index(folder) -> list[tuple[str, int]]:
    """(name, 1-based index) of every species, from the generated species_name()."""
    source = (folder / "chemistry.f90").read_text(encoding="utf-8")
    return [
        (name, int(i))
        for i, name in re.findall(r'case \((\d+)\); name = "([^"]+)"', source)
    ]


@needs_gfortran
@pytest.mark.parametrize("tag", list(CHI))
def test_fortran_equilibrium_matches_the_reference(fortran, tag):
    ref = _reference(tag)
    source = fortran.folder / f"benchmark_{tag}.f90"
    source.write_text(_fortran_driver(tag, fortran.species), encoding="utf-8")
    objects = [str(p) for p in fortran.folder.glob("*.o") if p.name != "main.o"]
    exe = fortran.folder / f"benchmark_{tag}"
    subprocess.run(
        ["gfortran", *FFLAGS, source.name, *objects, "-o", str(exe)],
        cwd=fortran.folder,
        check=True,
        capture_output=True,
        text=True,
    )

    proc = subprocess.run([str(exe)], cwd=fortran.folder, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    final = [float(v) for v in proc.stdout.split()]
    temperature, densities = final[0], dict(zip(fortran.species, final[1:]))

    assert temperature == pytest.approx(ref["Tgas"][-1], rel=AGREEMENT)
    for column, name in SPECIES.items():
        assert densities[name] == pytest.approx(ref[column][-1], rel=AGREEMENT), name


# --------------------------------------------------------------------------- #
# C++ (rosenbrock_cpp)                                                         #
# --------------------------------------------------------------------------- #

needs_gxx = pytest.mark.skipif(
    not sys.platform.startswith("linux") or shutil.which("g++") is None,
    reason="needs Linux with g++",
)


def _cpp_driver(tag: str, lu: str) -> str:
    """A program solving the benchmark for *tag* with linear solver *lu* and
    printing the status, T, then the densities."""
    ref = _reference(tag)
    densities = "\n".join(
        f'    n[chem::species_index("{name}")] = {float(ref[column][0])!r};'
        for column, name in SPECIES.items()
    )
    parameters = "\n".join(
        f"    p.{name} = {value!r};" for name, value in _parameters(CHI[tag]).items()
    )

    return f"""#include <cstdio>
#include "solver.hpp"

int main() {{
    chem::Parameters p;
{parameters}
    double n[chem::nspecs] = {{}};
{densities}
    double e = chem::energy(n, {float(ref["Tgas"][0])!r});
    chem::Options options;
    options.rtol = {RTOL!r};
    options.atol = {ATOL!r};
    options.lu = chem::LinearSolver::{lu};
    const chem::Status status = chem::solve(n, e, {T_END!r}, p, options);
    std::printf("%d %.17e", static_cast<int>(status), chem::temperature(n, e));
    for (int i = 0; i < chem::nspecs; ++i) std::printf(" %.17e", n[i]);
    std::printf("\\n");
}}
"""


@pytest.fixture(scope="module")
def cpp(tmp_path_factory) -> SimpleNamespace:
    """The generated C++, with the benchmark's table functions."""
    folder = generate(
        str(BENCHMARK / "GOW.jet"), tmp_path_factory.mktemp("gow_cpp"), "rosenbrock_cpp"
    )
    shutil.copy(BENCHMARK / "user_functions.hpp", folder / "user_functions.hpp")
    names = re.findall(r'^\s*"([^"]+)",$', (folder / "chemistry.hpp").read_text(), re.M)

    return SimpleNamespace(folder=folder, species=names)


@needs_gxx
@pytest.mark.parametrize("lu", ["sparse", "dense"])
@pytest.mark.parametrize("tag", list(CHI))
def test_cpp_equilibrium_matches_the_reference(cpp, tag, lu):
    ref = _reference(tag)
    source = cpp.folder / f"benchmark_{tag}_{lu}.cpp"
    source.write_text(_cpp_driver(tag, lu), encoding="utf-8")
    exe = cpp.folder / f"benchmark_{tag}_{lu}"
    subprocess.run(
        ["g++", "-std=c++17", "-O2", f"-I{cpp.folder}", str(source), "-o", str(exe)],
        check=True,
        capture_output=True,
        text=True,
    )

    proc = subprocess.run([str(exe)], cwd=cpp.folder, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    status, temperature, *values = [float(v) for v in proc.stdout.split()]
    densities = dict(zip(cpp.species, values))

    assert status == 0, proc.stdout
    assert temperature == pytest.approx(ref["Tgas"][-1], rel=AGREEMENT)
    for column, name in SPECIES.items():
        assert densities[name] == pytest.approx(ref[column][-1], rel=AGREEMENT), name


# --------------------------------------------------------------------------- #
# C++ (kokkos_ode: header-only-integrators VODE)                               #
# --------------------------------------------------------------------------- #

needs_clang = pytest.mark.skipif(
    not sys.platform.startswith("linux")
    or shutil.which("cmake") is None
    or shutil.which("clang++") is None,
    reason="needs Linux with cmake and clang++ (and network access, or a local "
    "header-only-integrators in FETCHCONTENT_SOURCE_DIR_HEADER_ONLY_INTEGRATORS)",
)

#: This VODE stalls at order 1 below rtol ~ 1e-10, so it runs at the template's 1e-8
KOKKOS_RTOL = 1e-8


def _kokkos_driver(tag: str, parameters: list[str]) -> str:
    """A program solving the benchmark for *tag* with VODE and printing the result
    code, T, then the densities; *parameters* are the generated ones."""
    ref = _reference(tag)
    values = _parameters(CHI[tag])
    densities = "\n".join(
        f'    y[ChemistryODE::species_index("{name}")] = {float(ref[column][0])!r};'
        for column, name in SPECIES.items()
    )
    assignments = "\n".join(
        f"    ChemistryODE::params.{name} = {values[name]!r};" for name in parameters
    )

    return f"""#include <cstdio>
#include "chemistry_ode.hpp"
using namespace integrators;

int main() {{
{assignments}
    ChemistryODE::evolve_temperature = true;
    ChemistryODE::state_type y{{}};
{densities}
    y[ChemistryODE::nspecs] = {float(ref["Tgas"][0])!r};
    VODE<ChemistryODE> integrator{{}};
    VODEState<ChemistryODE::neqs> state{{}};
    state.t = 0.0;
    state.tout = {T_END!r};
    state.y = y;
    state.rtol = {KOKKOS_RTOL!r};
    state.atol = {ATOL!r};
    state.max_steps = 1000000;
    state.jacobian_analytic = true;
    ChemistryODE::state_type problem{{}};
    const auto result = integrator.integrate(problem, state);
    std::printf("%d %.17e", static_cast<int>(result), state.y[ChemistryODE::nspecs]);
    for (int i = 0; i < ChemistryODE::nspecs; ++i) std::printf(" %.17e", state.y[i]);
    std::printf("\\n");
}}
"""


@pytest.fixture(scope="module")
def kokkos(tmp_path_factory) -> SimpleNamespace:
    """The generated kokkos_ode C++, with the benchmark's table functions."""
    network = str(BENCHMARK / "GOW.jet")
    folder = generate(network, tmp_path_factory.mktemp("gow_kokkos"), "kokkos_ode")
    shutil.copy(BENCHMARK / "user_functions.hpp", folder / "user_functions.hpp")
    configure = kokkos_build.configure(folder)
    assert configure.returncode == 0, configure.stdout + configure.stderr
    source = (folder / "chemistry_ode.hpp").read_text(encoding="utf-8")
    block = source[source.index("struct ChemistryParameters") : source.index("};")]

    return SimpleNamespace(
        folder=folder,
        include=kokkos_build.include_dir(folder),
        species=re.findall(r'^\s*"([^"]+)",$', source, re.M),
        parameters=re.findall(r"double (\w+) = std::numeric_limits", block),
    )


@needs_clang
@pytest.mark.parametrize("tag", list(CHI))
def test_kokkos_equilibrium_matches_the_reference(kokkos, tag):
    ref = _reference(tag)
    source = kokkos.folder / f"benchmark_{tag}.cpp"
    source.write_text(_kokkos_driver(tag, kokkos.parameters), encoding="utf-8")
    exe = kokkos.folder / f"benchmark_{tag}"
    subprocess.run(
        ["clang++", "-std=c++20", "-O2", f"-I{kokkos.folder}", f"-I{kokkos.include}"]
        + [str(source), "-o", str(exe)],
        check=True,
        capture_output=True,
        text=True,
    )

    proc = subprocess.run([str(exe)], cwd=kokkos.folder, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    status, temperature, *values = [float(v) for v in proc.stdout.split()]
    densities = dict(zip(kokkos.species, values))

    assert status == 1, proc.stdout  # IntegratorResult::SUCCESS
    assert temperature == pytest.approx(ref["Tgas"][-1], rel=AGREEMENT)
    for column, name in SPECIES.items():
        assert densities[name] == pytest.approx(ref[column][-1], rel=AGREEMENT), name


# --------------------------------------------------------------------------- #
# AMReX Microphysics (microphysics: network + its EOS, one-zone VODE)          #
# --------------------------------------------------------------------------- #

needs_microphysics = pytest.mark.skipif(
    not sys.platform.startswith("linux") or not microphysics_cell.available(),
    reason="needs Linux, make, and MICROPHYSICS_HOME and AMREX_HOME set",
)


@pytest.fixture(scope="module")
def microphysics(tmp_path_factory):
    """The one-zone Microphysics test built with the generated network and tables."""
    network = str(BENCHMARK / "GOW.jet")
    folder = generate(network, tmp_path_factory.mktemp("gow_mp"), "microphysics")
    shutil.copy(BENCHMARK / "user_functions.hpp", folder / "user_functions.H")
    workdir = tmp_path_factory.mktemp("gow_mp_cell")
    build = microphysics_cell.build(folder, workdir)
    assert build.returncode == 0, build.stdout[-3000:] + build.stderr[-3000:]

    return workdir


@needs_microphysics
@pytest.mark.parametrize("tag", list(CHI))
def test_microphysics_equilibrium_matches_the_reference(microphysics, tag):
    ref = _reference(tag)
    start = {name: float(ref[column][0]) for column, name in SPECIES.items()}

    proc = microphysics_cell.run(
        microphysics,
        start,
        _parameters(CHI[tag]),
        f"unit_test.temperature={float(ref['Tgas'][0])!r}",
        f"unit_test.tmax={T_END!r}",
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    success, temperature, densities = microphysics_cell.final_state(proc.stdout)

    assert success
    assert temperature == pytest.approx(ref["Tgas"][-1], rel=AGREEMENT)
    for column, name in SPECIES.items():
        assert densities[name] == pytest.approx(ref[column][-1], rel=AGREEMENT), name


# --------------------------------------------------------------------------- #
# FLASH (flash: KromeChemistry network + one-zone simulation, in an overlay)   #
# --------------------------------------------------------------------------- #

needs_flash = pytest.mark.skipif(
    not sys.platform.startswith("linux") or not flash_overlay.available(),
    reason="needs Linux, make, mpif90 and FLASH_HOME naming a FLASH checkout",
)


@pytest.fixture(scope="module")
def flash(tmp_path_factory):
    """FLASH's one-zone simulation built with the generated network and tables."""
    folder = generate(
        str(BENCHMARK / "GOW.jet"), tmp_path_factory.mktemp("gow_fl"), "flash"
    )
    unit = folder / "source/physics/sourceTerms/KromeChemistry/KromeChemistryMain/GOW"
    shutil.copy(BENCHMARK / "user_functions.f90", unit / "user_functions.f90")
    home = flash_overlay.overlay(
        folder, tmp_path_factory.mktemp("gow_fl_home") / "f", "GOW"
    )
    build = flash_overlay.build(home, "GOW")
    assert build.returncode == 0, build.stdout[-3000:] + build.stderr[-3000:]

    return home / "obj"


@needs_flash
@pytest.mark.parametrize("tag", list(CHI))
def test_flash_equilibrium_matches_the_reference(flash, tag, tmp_path):
    ref = _reference(tag)
    start = {name: float(ref[column][0]) for column, name in SPECIES.items()}
    settings = {f"jaff_{k}": v for k, v in _parameters(CHI[tag]).items()}
    settings |= {"sim_initTemp": float(ref["Tgas"][0]), "tmax": T_END}

    proc = flash_overlay.run(flash, start, settings, tmp_path)
    assert proc.returncode == 0, proc.stdout[-3000:] + proc.stderr[-3000:]
    t, temperature, densities = flash_overlay.final_state(tmp_path)

    assert t == pytest.approx(T_END)
    assert temperature == pytest.approx(ref["Tgas"][-1], rel=AGREEMENT)
    for column, name in SPECIES.items():
        assert densities[name] == pytest.approx(ref[column][-1], rel=AGREEMENT), name
