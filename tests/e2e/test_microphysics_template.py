# ABOUTME: The microphysics template (an AMReX Microphysics network plus its EOS): the
# ABOUTME: generated files, and a one-zone Microphysics build and run of real networks
"""End-to-end tests of the generated Microphysics network (``microphysics``).

The file checks need only JAFF. The build tests compile the one-zone test in
``tests/fixtures/microphysics_cell`` against the generated network and its
``eos/`` with Microphysics' GNU make; they run only when ``MICROPHYSICS_HOME``
and ``AMREX_HOME`` name the checkouts.
"""

import re
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from jaff import Network
from tests import microphysics_cell
from tests._paths import FIXTURES_DIR, REPO_ROOT
from tests.python_solver import generate, reference

DEMO2 = str(REPO_ROOT / "networks" / "demos" / "demo2.jet")
GOW = str(FIXTURES_DIR / "GOW_benchmark" / "GOW.jet")

needs_microphysics = pytest.mark.skipif(
    not sys.platform.startswith("linux") or not microphysics_cell.available(),
    reason="needs Linux, make, and MICROPHYSICS_HOME and AMREX_HOME set",
)


@pytest.fixture(scope="module")
def gow(tmp_path_factory) -> SimpleNamespace:
    folder = generate(GOW, tmp_path_factory.mktemp("gow"), "microphysics")
    return SimpleNamespace(net=Network(GOW), folder=folder)


def test_writes_the_network_and_its_eos(gow):
    names = {str(p.relative_to(gow.folder)) for p in gow.folder.rglob("*")}

    assert {
        "pynucastro.net",
        "_parameters",
        "actual_network.H",
        "actual_network_data.cpp",
        "actual_rhs.H",
        "user_functions.H",
        "Make.package",
        "eos/actual_eos.H",
        "eos/actual_eos_data.H",
        "eos/_parameters",
        "eos/Make.package",
    } <= names
    assert "jaffgen_template_functions.py" not in names


def test_species_file_lists_every_species_with_an_identifier_short_name(gow):
    rows = [
        line.split()
        for line in (gow.folder / "pynucastro.net").read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]

    assert [row[0] for row in rows] == gow.net.species.names()
    assert all(re.fullmatch(r"[A-Z]\w*", row[1]) for row in rows)
    assert {row[0]: row[1] for row in rows}["C"] == "Carbon"  # C: constants namespace
    masses = {row[0]: float(row[2]) for row in rows}
    assert masses["H2"] == pytest.approx(2.0159, rel=1e-4)
    assert {row[0]: int(row[3]) for row in rows}["e-"] == -1


def test_parameters_are_runtime_parameters_unset_by_default(gow):
    text = (gow.folder / "_parameters").read_text()
    params = [str(p) for p in gow.net.ode_system(thermal="dedt").parameters]

    assert "@namespace: network" in text
    for name in params:
        assert re.search(rf"^{name}\s+real\s+-1\.e300$", text, re.M), name


def test_eos_takes_the_species_gammas_and_masses(gow):
    text = (gow.folder / "eos" / "actual_eos.H").read_text()
    gammas = [float(v) for v in re.findall(r"gammas\[\d+\] = (\S+);", text)]
    masses = [float(v) for v in re.findall(r"spmasses\[\d+\] = (\S+);", text)]

    assert gammas == [
        pytest.approx(float(s.thermodynamics.gamma)) for s in gow.net.species
    ]
    assert masses == [pytest.approx(s.mass, rel=1e-6, abs=0) for s in gow.net.species]


def test_unknown_functions_get_device_stubs(gow):
    text = (gow.folder / "user_functions.H").read_text()

    assert "AMREX_GPU_HOST_DEVICE AMREX_INLINE amrex::Real l0_co_interp1d(" in text
    assert "l0_co_interp1d_partial_0(" in text


def test_species_sharing_a_short_name_are_an_error(tmp_path):
    network = tmp_path / "clash.dat"
    network.write_text("CO + Co -> C + O + Co [10,1e4] 1e-10\n", encoding="utf-8")

    with pytest.raises(Exception, match="short name"):
        generate(str(network), tmp_path / "out", "microphysics")


# --------------------------------------------------------------------------- #
# One-zone Microphysics build and run                                          #
# --------------------------------------------------------------------------- #


@pytest.fixture(scope="module")
def demo2_cell(tmp_path_factory) -> Path:
    folder = generate(DEMO2, tmp_path_factory.mktemp("demo2"), "microphysics")
    workdir = tmp_path_factory.mktemp("demo2_cell")
    build = microphysics_cell.build(folder, workdir)
    assert build.returncode == 0, build.stdout[-3000:] + build.stderr[-3000:]

    return workdir


@needs_microphysics
def test_demo2_matches_the_python_template(demo2_cell, tmp_path):
    start = {"e-": 1.0, "CO+": 1.0, "CO": 10.0}
    proc = microphysics_cell.run(demo2_cell, start, {"av": 1.0}, "unit_test.tmax=1e10")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    success, temperature, n = microphysics_cell.final_state(proc.stdout)
    ref_t, ref_n = reference(DEMO2, tmp_path, start, 100.0, 1e10, {"av": 1.0})

    assert success
    # recombination heats the gas by ~0.5 K: T checks the energy bookkeeping
    assert temperature == pytest.approx(ref_t, rel=1e-8)
    for name, value in ref_n.items():
        assert n[name] == pytest.approx(value, rel=1e-6), name


@needs_microphysics
def test_unset_parameters_stop_the_run_named(demo2_cell):
    proc = microphysics_cell.run(demo2_cell, {"e-": 1.0, "CO+": 1.0}, {})

    assert proc.returncode != 0
    assert "network.av" in proc.stdout + proc.stderr


@needs_microphysics
def test_gow_builds_with_its_stubs(gow, tmp_path):
    build = microphysics_cell.build(gow.folder, tmp_path)

    assert build.returncode == 0, build.stdout[-3000:] + build.stderr[-3000:]


@needs_microphysics
def test_cothin_burns_and_conserves_elements_and_charge(tmp_path_factory):
    net = Network("COthin")
    folder = generate("COthin", tmp_path_factory.mktemp("cothin"), "microphysics")
    workdir = tmp_path_factory.mktemp("cothin_cell")
    build = microphysics_cell.build(folder, workdir)
    assert build.returncode == 0, build.stdout[-3000:] + build.stderr[-3000:]
    # every species at 1 cm^-3, electrons balancing the other charges: neutral
    charge = sum(s.charge for s in net.species if s.name != "e-")
    start = {s.name: float(charge) if s.name == "e-" else 1.0 for s in net.species}

    proc = microphysics_cell.run(
        workdir,
        start,
        {"av": 1.0, "crate": 1e-16, "tdust": 20.0},
        "unit_test.temperature=1000.0",  # COthin rate fits overflow below ~200 K
        "unit_test.tmax=1e10",
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    success, _, final = microphysics_cell.final_state(proc.stdout)

    assert success
    names = net.species.names()
    for row in net.elements.density_matrix():
        before = sum(c * start[n] for c, n in zip(row, names))
        after = sum(c * final[n] for c, n in zip(row, names))
        assert after == pytest.approx(before, rel=1e-6)
    assert sum(s.charge * start[s.name] for s in net.species) == 0.0
    net_charge = sum(s.charge * final[s.name] for s in net.species)
    assert net_charge == pytest.approx(0.0, abs=1e-6 * charge)
