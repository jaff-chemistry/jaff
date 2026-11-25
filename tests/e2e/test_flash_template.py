# ABOUTME: The flash template (a network for FLASH's KromeChemistry unit plus a one-zone
# ABOUTME: simulation): the generated files, and a build and run in a FLASH overlay
"""End-to-end tests of the generated FLASH network (``flash``).

The file checks need only JAFF. The build tests set up and compile the one-zone
simulation in a temporary overlay of the FLASH checkout named by ``FLASH_HOME``
(:mod:`tests.helpers.flash_overlay`), which they leave untouched.
"""

import os
import re
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from jaff import Network
from tests.helpers import flash_overlay
from tests._paths import REPO_ROOT
from tests.helpers.python_solver import generate, reference

DEMO2 = str(REPO_ROOT / "networks" / "demos" / "demo2.jet")
UNIT = Path("source/physics/sourceTerms/KromeChemistry/KromeChemistryMain")
COMPOSITION = Path("source/Simulation/SimulationComposition/KromeChemistry")
SIMULATION = Path("source/Simulation/SimulationMain")

needs_flash = pytest.mark.skipif(
    not sys.platform.startswith("linux") or not flash_overlay.available(),
    reason="needs Linux, make, mpif90 and FLASH_HOME naming a FLASH checkout",
)


@pytest.fixture(scope="module")
def cothin(tmp_path_factory) -> SimpleNamespace:
    folder = generate("COthin", tmp_path_factory.mktemp("cothin"), "flash")
    return SimpleNamespace(net=Network("COthin"), folder=folder)


def test_writes_the_units_under_flash_paths_named_after_the_network(cothin):
    # Forward slashes on every OS, as the expected names are written.
    files = {p.relative_to(cothin.folder).as_posix() for p in cothin.folder.rglob("*.*")}
    label = cothin.net.label  # the network file's stem, react_COthin
    unit, comp, sim = ((d / label).as_posix() for d in (UNIT, COMPOSITION, SIMULATION))

    assert {
        f"{unit}/chemistry.f90",
        f"{unit}/chemistry_dlsodes.f90",
        f"{unit}/user_functions.f90",
        f"{unit}/KromeChemistry.F90",
        f"{unit}/KromeChemistry_init.F90",
        f"{unit}/KromeChemistry_data.F90",
        f"{comp}/Simulation_initSpecies.F90",
        f"{sim}/Simulation_init.F90",
        f"{sim}/Simulation_initBlock.F90",
        f"{sim}/IO_writeIntegralQuantities.F90",
        f"{sim}/flash.par",
    } <= files
    assert (cothin.folder / unit / "Config").is_file()
    assert not any("$" in f or "jaffgen_template_functions" in f for f in files)


def test_config_declares_every_parameter_unset_and_every_species(cothin):
    text = (cothin.folder / UNIT / cothin.net.label / "Config").read_text()
    params = [str(p) for p in cothin.net.ode_system(thermal="dtdt").parameters]

    for name in params:
        assert re.search(rf"^PARAMETER jaff_{name} REAL -1\.0e300$", text, re.M), name
    species = re.findall(r"^SPECIES (\w+)$", text, re.M)
    assert species == list(cothin.net.species.normalized_names())


def test_species_get_the_eos_gammas_and_molar_masses(cothin):
    text = (
        cothin.folder / COMPOSITION / cothin.net.label / "Simulation_initSpecies.F90"
    ).read_text()
    gammas = [
        float(v.replace("d", "e")) for v in re.findall(r"gammas\(\d+\) = (\S+)", text)
    ]
    masses = [
        float(v.replace("d", "e")) for v in re.findall(r"masses\(\d+\) = (\S+)", text)
    ]

    assert gammas == [
        pytest.approx(float(s.thermodynamics.gamma)) for s in cothin.net.species
    ]
    assert masses == [pytest.approx(s.mass, rel=1e-6, abs=0) for s in cothin.net.species]
    assert "masses(n) * avogadro" in text


def test_species_sharing_a_flash_name_are_an_error(tmp_path):
    network = tmp_path / "clash.dat"
    network.write_text("CO + Co -> C + O + Co [10,1e4] 1e-10\n", encoding="utf-8")

    with pytest.raises(Exception, match="FLASH name"):
        generate(str(network), tmp_path / "out", "flash")


# --------------------------------------------------------------------------- #
# Build and run in a FLASH overlay                                             #
# --------------------------------------------------------------------------- #


def _checkout_status() -> str:
    """``git status`` of the FLASH checkout (to show the tests leave it alone)."""
    home = os.environ["FLASH_HOME"]
    return subprocess.run(
        ["git", "-C", home, "status", "--porcelain", "--untracked-files=all"],
        capture_output=True,
        text=True,
    ).stdout


@pytest.fixture(scope="module")
def demo2_flash(tmp_path_factory) -> SimpleNamespace:
    before = _checkout_status()
    folder = generate(DEMO2, tmp_path_factory.mktemp("demo2"), "flash")
    home = flash_overlay.overlay(folder, tmp_path_factory.mktemp("home") / "f", "demo2")
    build = flash_overlay.build(home, "demo2")
    assert build.returncode == 0, build.stdout[-3000:] + build.stderr[-3000:]

    return SimpleNamespace(obj=home / "obj", before=before)


@needs_flash
def test_demo2_matches_the_python_template(demo2_flash, tmp_path):
    start = {"e-": 1.0, "CO+": 1.0, "CO": 10.0}
    proc = flash_overlay.run(
        demo2_flash.obj, start, {"jaff_av": 1.0, "tmax": 1.0e10}, tmp_path / "run"
    )
    assert proc.returncode == 0, proc.stdout[-3000:] + proc.stderr[-3000:]
    t, temperature, n = flash_overlay.final_state(tmp_path / "run")
    ref_t, ref_n = reference(DEMO2, tmp_path / "ref", start, 100.0, 1e10, {"av": 1.0})

    assert t == pytest.approx(1.0e10)
    # recombination heats the gas by ~0.5 K: T checks the energy bookkeeping
    assert temperature == pytest.approx(ref_t, rel=1e-7)
    for name, value in ref_n.items():
        assert n[name] == pytest.approx(value, rel=1e-6), name


@needs_flash
def test_unset_parameters_stop_the_run_named(demo2_flash, tmp_path):
    proc = flash_overlay.run(
        demo2_flash.obj, {"e-": 1.0, "CO+": 1.0}, {"jaff_av": -1.0e300}, tmp_path
    )

    assert proc.returncode != 0
    assert "set in flash.par: jaff_av" in proc.stdout + proc.stderr


@needs_flash
def test_the_flash_checkout_is_left_untouched(demo2_flash):
    assert _checkout_status() == demo2_flash.before


@needs_flash
def test_cothin_burns_and_conserves_elements_and_charge(cothin, tmp_path):
    net, label = cothin.net, cothin.net.label
    home = flash_overlay.overlay(cothin.folder, tmp_path / "home", label)
    build = flash_overlay.build(home, label)
    assert build.returncode == 0, build.stdout[-3000:] + build.stderr[-3000:]
    # every species at 1 cm^-3, electrons balancing the other charges: neutral
    charge = sum(s.charge for s in net.species if s.name != "e-")
    start = {s.name: float(charge) if s.name == "e-" else 1.0 for s in net.species}
    settings = {"jaff_av": 1.0, "jaff_crate": 1e-16, "jaff_tdust": 20.0}
    # COthin rate fits overflow below ~200 K
    settings |= {"sim_initTemp": 1000.0, "tmax": 1.0e10}

    proc = flash_overlay.run(home / "obj", start, settings, tmp_path / "run")
    assert proc.returncode == 0, proc.stdout[-3000:] + proc.stderr[-3000:]
    _, _, final = flash_overlay.final_state(tmp_path / "run")

    names = net.species.names()
    for row in net.elements.density_matrix():
        before = sum(c * start[n] for c, n in zip(row, names))
        after = sum(c * final[n] for c, n in zip(row, names))
        assert after == pytest.approx(before, rel=1e-6)
    net_charge = sum(s.charge * final[s.name] for s in net.species)
    assert net_charge == pytest.approx(0.0, abs=1e-6 * charge)
