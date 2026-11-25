# ABOUTME: Tests for the Language registry: file-extension lookup per target language
# ABOUTME: Each language declares its own template extensions; duplicates are rejected

import pytest

from jaff.codegen._languages import Language

EXTENSION_MAP = [
    ("cpp", "cxx"),
    ("cxx", "cxx"),
    ("cc", "cxx"),
    ("hpp", "cxx"),
    ("hxx", "cxx"),
    ("hh", "cxx"),
    ("h", "cxx"),
    ("c", "c"),
    ("f", "fortran"),
    ("for", "fortran"),
    ("f90", "fortran"),
    ("f95", "fortran"),
    ("f03", "fortran"),
    ("f08", "fortran"),
    ("py", "python"),
    ("rs", "rust"),
    ("jl", "julia"),
    ("r", "r"),
]


@pytest.fixture
def clean_registry():
    yield
    Language.LOOKUP.pop("dup_test_lang", None)
    Language._register.pop("dup_test_lang", None)


@pytest.mark.parametrize("ext,name", EXTENSION_MAP)
def test_extension_resolves_to_language(ext, name):
    lang = Language.from_extension(ext)
    assert lang is not None
    assert lang.name == name


@pytest.mark.parametrize("ext,name", [(".F90", "fortran"), (".HPP", "cxx")])
def test_leading_dot_and_case_are_accepted(ext, name):
    assert Language.from_extension(ext).name == name


@pytest.mark.parametrize("ext", ["txt", "", "."])
def test_unknown_extension_returns_none(ext):
    assert Language.from_extension(ext) is None


def test_from_extension_returns_the_singleton():
    assert Language.from_extension("cpp") is Language("cxx")


def test_duplicate_extension_is_rejected_without_registering(clean_registry):
    with pytest.raises(ValueError, match="already belongs to"):

        class Dup(Language):
            name = "dup_test_lang"
            extensions = ("py",)
            brac = "[]"
            matrix_sep = "]["

    assert "dup_test_lang" not in Language.LOOKUP
    assert "dup_test_lang" not in Language._register
    assert Language.from_extension("py").name == "python"
