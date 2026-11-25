# ABOUTME: Codegen token overrides: None keeps the language's token, "" is a real value
# ABOUTME: (drops C++'s ;); var_prefix / def_prefix replace the CSE declaration's type

import re

import pytest

from jaff import Network
from jaff.codegen import Codegen, Language


@pytest.fixture(scope="module")
def cg():
    net = Network("h_photoionization", radiation={"bands": [1, 13.6, 100, "inf"]})

    return Codegen(net, lang="cxx")


class TestDerive:
    def test_no_override_is_the_language_itself(self):
        cxx = Language("c++")

        assert cxx.derive() is cxx

    def test_an_empty_string_is_a_real_override(self):
        view = Language("c++").derive(line_end="", assignment_op="")

        assert (view.line_end, view.assignment_op) == ("", "")
        assert Language("c++").line_end == ";"


class TestEmptyLineEnd:
    def test_default_keeps_the_terminator(self, cg):
        assert cg.get_ode_str().splitlines()[-1].endswith(";")

    def test_empty_line_end_drops_it(self, cg):
        lines = cg.get_ode_str(line_end="").splitlines()

        assert lines
        assert not any(line.endswith(";") for line in lines)


@pytest.mark.parametrize(
    ("method", "prefix_arg", "options"),
    [
        ("get_ode_str", "def_prefix", {}),
        ("get_rhs_str", "def_prefix", {"thermal": "dedt"}),
        ("get_radode_str", "def_prefix", {}),
        ("get_jacobian_str", "var_prefix", {"thermal": "dtdt"}),
    ],
)
def test_prefix_replaces_the_declaration_type(cg, method, prefix_arg, options):
    lines = getattr(cg, method)(**{prefix_arg: "amrex::Real "}, **options).splitlines()
    declared = [line for line in lines if re.search(r"\br?cse\d+ = ", line)]

    assert declared, f"{method} declared no temporaries"
    assert all(re.match(r"amrex::Real r?cse\d+ = ", line) for line in declared), declared


def test_rates_var_prefix():
    cg = Codegen(Network("GOW"), lang="cxx")
    lines = cg.get_rates_str(var_prefix="amrex::Real ").splitlines()
    declared = [line for line in lines if re.match(r"\S+ x\d+ = ", line)]

    assert declared
    assert all(line.startswith("amrex::Real x") for line in declared), declared[:3]


@pytest.mark.parametrize("prefix", ["", "amrex::Real "])
def test_declarations_keep_the_language_assignment(cg, prefix):
    """An accumulating Jacobian (+=) still declares its temporaries with =."""
    lines = cg.get_jacobian_str(
        thermal="dtdt", assignment_op="+=", var_prefix=prefix
    ).splitlines()
    entries = [line for line in lines if line.startswith("J[")]
    declared = [line for line in lines if line and not line.startswith("J[")]

    assert declared and entries
    assert all(re.search(r"\bcse\d+ = ", line) for line in declared), declared[:2]
    assert all(" += " in line for line in entries), entries[:2]
