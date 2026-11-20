# ABOUTME: A built-in template is one folder, <repo>/templates/<name>: jaffgen
# ABOUTME: collects every file there (static solver files too); no preprocessor/Builder

import importlib
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

import jaff
from jaff.cli.jaffgen._engine import JaffGen
from jaff.config import JAFF_DIR, TEMPLATES_DIR
from tests._paths import REPO_ROOT
from tests.jaffgen_run import run_folder

#: Static solver files each template ships besides its $$JAFF templates.
STATIC = {
    "fortran_dlsodes": {"main.f90", "Makefile", "opkda1.f", "opkda2.f", "opkdmain.f"},
    "python_solve_ivp": {"main.py", "solver.py"},
    "microphysics": {"actual_network.H", "Make.package"},
}


@pytest.mark.parametrize("template", STATIC)
def test_static_solver_files_are_collected(template):
    names = {f.relpath.as_posix() for f in JaffGen._template_files(template)}

    assert STATIC[template] <= names


@pytest.mark.parametrize("template", STATIC)
def test_every_file_comes_from_the_template_folder(template):
    folder = (TEMPLATES_DIR / template).resolve()

    assert all(
        f.abspath.is_relative_to(folder) for f in JaffGen._template_files(template)
    )


def test_templates_live_at_the_repository_root():
    assert TEMPLATES_DIR.resolve() == (REPO_ROOT / "templates").resolve()
    assert not (JAFF_DIR / "templates").exists()


def test_one_level_of_folders():
    assert not (TEMPLATES_DIR / "generator").exists()
    assert not (TEMPLATES_DIR / "preprocessor").exists()


def test_builder_and_preprocessor_are_gone():
    assert not hasattr(jaff, "Builder") and not hasattr(jaff, "Preprocessor")
    assert not hasattr(jaff.codegen, "Builder")
    assert not hasattr(jaff.codegen, "Preprocessor")
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("jaff.plugins")


def _generate_from(indir, outdir, network):
    """Run the jaffgen CLI on the template folder *indir*; return the run folder."""
    jaffgen = shutil.which("jaffgen", path=str(Path(sys.executable).parent))
    assert jaffgen is not None
    proc = subprocess.run(
        [jaffgen, "--network", str(network), "--indir", str(indir)]
        + ["--outdir", str(outdir)],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr

    return run_folder(outdir)


def test_label_in_a_template_path_becomes_the_network_label(tmp_path):
    network = REPO_ROOT / "networks" / "demos" / "demo2.jet"
    indir = tmp_path / "templates"
    (indir / "units" / "$label$").mkdir(parents=True)
    (indir / "units" / "$label$" / "notes_$label$.txt").write_text("plain\n")
    (indir / "top.txt").write_text("top\n")

    out = _generate_from(indir, tmp_path / "out", network)

    assert (out / "units" / "demo2" / "notes_demo2.txt").read_text() == "plain\n"
    assert (out / "top.txt").is_file()
    assert not any("$" in str(p) for p in out.rglob("*"))
