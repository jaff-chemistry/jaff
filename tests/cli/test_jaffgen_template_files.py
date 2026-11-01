# ABOUTME: A built-in template is one folder, templates/generator/<name>: jaffgen
# ABOUTME: collects every file there (static solver files too); no preprocessor/Builder

import importlib

import pytest

import jaff
from jaff.cli.jaffgen._engine import JaffGen
from jaff.config import TEMPLATES_DIR

#: Static solver files each template ships besides its $$JAFF templates.
STATIC = {
    "fortran_dlsodes": {"main.f90", "Makefile", "opkda1.f", "opkda2.f", "opkdmain.f"},
    "python_solve_ivp": {"main.py"},
    "microphysics": {"actual_network.H", "Make.package"},
}


@pytest.mark.parametrize("template", STATIC)
def test_static_solver_files_are_collected(template):
    names = {f.relpath.as_posix() for f in JaffGen._template_files(template)}

    assert STATIC[template] <= names


@pytest.mark.parametrize("template", STATIC)
def test_every_file_comes_from_the_generator_folder(template):
    folder = (TEMPLATES_DIR / "generator" / template).resolve()

    assert all(
        f.abspath.is_relative_to(folder) for f in JaffGen._template_files(template)
    )


def test_no_preprocessor_templates():
    assert not (TEMPLATES_DIR / "preprocessor").exists()


def test_builder_and_preprocessor_are_gone():
    assert not hasattr(jaff, "Builder") and not hasattr(jaff, "Preprocessor")
    assert not hasattr(jaff.codegen, "Builder")
    assert not hasattr(jaff.codegen, "Preprocessor")
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("jaff.plugins")
