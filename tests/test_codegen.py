# ABOUTME: Tests for multi-language code generation, CSE, and ODE/Jacobian output
# ABOUTME: Language dialect table is data-driven; numeric ODE/Jac checks use fixed fixtures

import re
from pathlib import Path
from typing import List

import pytest
import sympy as sp

from jaff import Network
from jaff.codegen import Codegen

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def jac_network():
    """A single-reaction network used for language and ODE/Jacobian checks."""
    return Network(str(FIXTURES / "test_jac.dat"))


@pytest.fixture
def jac_codegen(jac_network):
    return Codegen(jac_network, lang="c++")


@pytest.fixture
def dedt_network():
    """A single-reaction network that also carries an internal-energy rate."""
    return Network(str(FIXTURES / "test_jac_dedt.dat"))


@pytest.fixture
def dedt_codegen(dedt_network):
    return Codegen(dedt_network, lang="c++")


@pytest.fixture
def cse_network():
    """A network with shared subexpressions across reaction rates."""
    return Network(str(FIXTURES / "test_cse.dat"))


@pytest.fixture
def cse_codegen(cse_network):
    return Codegen(cse_network, lang="c++")


# --------------------------------------------------------------------------- #
# Language dialects                                                            #
# --------------------------------------------------------------------------- #

# name -> dialect properties.  ``sep=None`` means the separator is unspecified
# for that language and is not asserted.
LANGS = {
    "c": dict(
        lb="[",
        rb="]",
        offset=0,
        line_end=";",
        comment="//",
        sep="][",
        types={"int": "int ", "float": "float ", "double": "double ", "bool": "_Bool "},
    ),
    "cxx": dict(
        lb="[",
        rb="]",
        offset=0,
        line_end=";",
        comment="//",
        sep="][",
        types={"int": "int ", "float": "float ", "double": "double ", "bool": "bool "},
    ),
    "fortran": dict(
        lb="(",
        rb=")",
        offset=1,
        line_end="",
        comment="!",
        sep=", ",
        types={"int": None, "float": None, "double": None, "bool": None},
    ),
    "python": dict(
        lb="[",
        rb="]",
        offset=0,
        line_end="",
        comment="#",
        sep="][",
        types={"int": None, "float": None, "double": None, "bool": None},
    ),
    "rust": dict(
        lb="[",
        rb="]",
        offset=0,
        line_end=";",
        comment="//",
        sep=None,
        types={"int": "i32 ", "float": "f32 ", "double": "f64 ", "bool": "bool "},
    ),
    "julia": dict(
        lb="[",
        rb="]",
        offset=1,
        line_end="",
        comment="#",
        sep=", ",
        types={
            "int": "Int64 ",
            "float": "Float32 ",
            "double": "Float64 ",
            "bool": "Bool ",
        },
    ),
    "r": dict(
        lb="[",
        rb="]",
        offset=1,
        line_end="",
        comment="#",
        sep=", ",
        types={"int": None, "float": None, "double": None, "bool": None},
    ),
}

# Canonical aliases accepted by Codegen(lang=...).
ALIASES = {
    "c++": "cxx",
    "cpp": "cxx",
    "cxx": "cxx",
    "c": "c",
    "fortran": "fortran",
    "f90": "fortran",
    "python": "python",
    "py": "python",
    "rust": "rust",
    "rs": "rust",
    "julia": "julia",
    "jl": "julia",
    "r": "r",
}

SEMICOLON_LANGS = {"c", "cxx", "rust"}
ONE_BASED_LANGS = {"fortran", "julia", "r"}


@pytest.mark.parametrize("name", list(LANGS))
def test_language_dialect(jac_network, name):
    """Each language reports the expected brackets, indexing, and separators."""
    spec = LANGS[name]
    lang = Codegen(jac_network, lang=name).lang
    assert lang.name == name
    assert (lang.lb, lang.rb) == (spec["lb"], spec["rb"])
    assert lang.idx_offset == spec["offset"]
    assert lang.line_end == spec["line_end"]
    assert lang.comment == spec["comment"]
    if spec["sep"] is not None:
        assert lang.sep == spec["sep"]


@pytest.mark.parametrize("name", list(LANGS))
def test_language_types(jac_network, name):
    """Type declaration strings match each language's convention."""
    types = Codegen(jac_network, lang=name).lang.types
    for key, expected in LANGS[name]["types"].items():
        assert types.get(key) == expected


@pytest.mark.parametrize("alias,canonical", list(ALIASES.items()))
def test_language_aliases(jac_network, alias, canonical):
    assert Codegen(jac_network, lang=alias).lang.name == canonical


@pytest.mark.parametrize("name", list(LANGS))
def test_rate_generation_basic_syntax(jac_network, name):
    """Rate strings use the language's array syntax and terminator."""
    cg = Codegen(jac_network, lang=name)
    rates = cg.get_rates_str(rate_variable="k", use_cse=False)
    index_token = "k(" if name == "fortran" else "k["
    assert index_token in rates
    # R assigns with "<-"; every other language uses "=".
    assert ("<-" if name == "r" else "=") in rates
    if name in SEMICOLON_LANGS:
        assert rates.count(";") > 0


@pytest.mark.parametrize("name", sorted(ONE_BASED_LANGS))
def test_one_based_indexing(jac_network, name):
    cg = Codegen(jac_network, lang=name)
    rates = cg.get_rates_str(idx_offset=-1, rate_variable="k", use_cse=False)
    open_b = "(" if name == "fortran" else "["
    close_b = ")" if name == "fortran" else "]"
    assert f"k{open_b}1{close_b}" in rates or f"{open_b}1{close_b}" in rates
    assert f"k{open_b}0{close_b}" not in rates


@pytest.mark.parametrize("name", list(LANGS))
def test_flux_generation(jac_network, name):
    fluxes = Codegen(jac_network, lang=name).get_flux_expressions_str(flux_var="flux")
    assert len(fluxes) > 0
    assert ("flux(" if name == "fortran" else "flux[") in fluxes


@pytest.mark.parametrize("name", list(LANGS))
def test_ode_generation(jac_network, name):
    odes = Codegen(jac_network, lang=name).get_ode_str(ode_var="f", use_cse=False)
    assert len(odes) > 0


# --------------------------------------------------------------------------- #
# Common subexpression elimination (C++)                                       #
# --------------------------------------------------------------------------- #


class TestCSE:
    """CSE extraction in C++ rate generation (test_cse.dat has shared terms)."""

    def test_network_has_expected_reactions(self, cse_network):
        assert cse_network.reactions.count == 8

    def test_extracts_common_subexpressions(self, cse_codegen):
        rates = cse_codegen.get_rates_str(use_cse=True)
        assert rates.count("const double x") > 0

    def test_extracts_temperature_dependent_terms(self, cse_codegen):
        rates = cse_codegen.get_rates_str(use_cse=True)
        assert any(
            "tgas" in line and "const double x" in line for line in rates.splitlines()
        )

    def test_reduces_redundant_exp_calls(self, cse_codegen):
        no_cse = cse_codegen.get_rates_str(use_cse=False)
        with_cse = cse_codegen.get_rates_str(use_cse=True)
        assert with_cse.count("exp(") <= no_cse.count("exp(")

    def test_same_number_of_rate_assignments(self, cse_codegen):
        no_cse = cse_codegen.get_rates_str(use_cse=False)
        with_cse = cse_codegen.get_rates_str(use_cse=True)
        assert no_cse.count("k[") == with_cse.count("k[") > 0

    def test_output_statements_terminated(self, cse_codegen):
        rates = cse_codegen.get_rates_str(use_cse=True)
        assert "const double" in rates
        for line in rates.strip().splitlines():
            stripped = line.strip()
            if stripped and not stripped.startswith("//"):
                assert stripped.endswith(";"), f"unterminated: {line}"


def _cse_names(code: str, prefix: str):
    """Return (declared, used) CSE temporaries named ``<prefix><digits>``."""
    pattern = re.compile(rf"\b{re.escape(prefix)}\d+\b")
    declared: List[str] = []
    used = set()
    for line in code.splitlines():
        if "=" not in line:
            continue
        lhs, rhs = line.split("=", 1)
        target = lhs.strip().split()[-1] if lhs.strip() else ""
        if pattern.fullmatch(target):
            declared.append(target)
        used.update(pattern.findall(rhs))
    return declared, used


_CSE_STR_APIS = {
    "rates": lambda cg, p: cg.get_rates_str(use_cse=True, cse_var=p),
    "odes": lambda cg, p: cg.get_ode_str(use_cse=True, cse_var=p),
    "rhs": lambda cg, p: cg.get_rhs_str(use_cse=True, cse_var=p),
    "jacobian": lambda cg, p: cg.get_jacobian_str(use_cse=True, cse_var=p),
}


class TestCSEPrefixWithDigits:
    """CSE temporaries must be declared under the same name they are used by,
    even when the user-supplied prefix itself contains digits."""

    @pytest.mark.parametrize("prefix", ["tmp2", "stage2_cse", "x"])
    @pytest.mark.parametrize("api", sorted(_CSE_STR_APIS))
    def test_declared_temporaries_match_used(self, cse_network, api, prefix):
        code = _CSE_STR_APIS[api](Codegen(cse_network, lang="python"), prefix)
        declared, used = _cse_names(code, prefix)

        assert declared, f"no CSE temporaries emitted:\n{code}"
        assert len(declared) == len(set(declared)), f"duplicate temps:\n{code}"
        assert used <= set(declared), (
            f"undeclared temps {sorted(used - set(declared))}:\n{code}"
        )


# --------------------------------------------------------------------------- #
# Unknown-function derivative conversion                                       #
# --------------------------------------------------------------------------- #

_convert_derivs = Codegen._Codegen__convert_unknown_derivatives
_x, _y, _z = sp.symbols("x y z")
_fun = sp.Function("fun")


def _eval_partials(expr: sp.Expr) -> sp.Expr:
    """Evaluate ``fun_partial_k`` calls for the asymmetric ``fun(a, b) = a + 2*b``."""
    return expr.replace(sp.Function("fun_partial_0"), sp.Lambda((_x, _y), 1)).replace(
        sp.Function("fun_partial_1"), sp.Lambda((_x, _y), 2)
    )


class TestUnknownDerivatives:
    """Partial positions must come from the original derivative, not from the
    arguments after a ``Subs`` evaluation point made them coincide."""

    def test_subs_onto_repeated_argument_keeps_position(self):
        expr = sp.Subs(sp.Derivative(_fun(_x, _z), _z), _z, _x)
        assert _convert_derivs(expr) == sp.Function("fun_partial_1")(_x, _x)

    def test_diff_of_repeated_argument_sums_both_partials(self):
        converted = _convert_derivs(sp.diff(_fun(_y, _y), _y))
        assert _eval_partials(converted) == 3

    def test_repeated_argument_after_cse(self):
        jac = [sp.diff(_fun(_x * _y, _x * _y), _y)]
        replacements, reduced = sp.cse(jac, symbols=sp.numbered_symbols("cse"))
        defs = {str(k): v for k, v in replacements}
        converted = _convert_derivs(reduced[0], defs)
        subs_back = {k: v for k, v in replacements}
        assert _eval_partials(converted).xreplace(subs_back) == 3 * _x


# --------------------------------------------------------------------------- #
# ODE / Jacobian numeric output (C++)                                          #
# --------------------------------------------------------------------------- #


def _rhs_terms(text: str) -> List[str]:
    """Strip ``lhs =`` and trailing ``;`` from each generated statement."""
    return [line.split("=")[-1].strip().strip(";") for line in text.strip().splitlines()]


class TestOdeJacobian:
    """Exact generated RHS/Jacobian strings for the fake-rate test network."""

    def test_single_reaction_loaded(self, jac_network):
        assert jac_network.reactions.count == 1

    def test_rate_string(self, jac_codegen):
        rates = jac_codegen.get_rates_str().strip().splitlines()
        assert len(rates) == 1
        assert rates[-1].split("=")[-1].strip().rstrip(";") == "nden[0]"

    def test_ode_and_jacobian(self, jac_codegen):
        ode = _rhs_terms(jac_codegen.get_ode_str(use_cse=False))
        jac = _rhs_terms(jac_codegen.get_jacobian_str(use_cse=False))
        assert ode == [
            "-std::pow(nden[0], 2)*nden[1]",
            "-std::pow(nden[0], 2)*nden[1]",
            "std::pow(nden[0], 2)*nden[1]",
        ]
        assert jac == [
            "-2*nden[0]*nden[1]",
            "-std::pow(nden[0], 2)",
            "-2*nden[0]*nden[1]",
            "-std::pow(nden[0], 2)",
            "2*nden[0]*nden[1]",
            "std::pow(nden[0], 2)",
        ]


class TestOdeJacobianWithInternalEnergy:
    """Same, for a network carrying an internal-energy (dEdt) contribution."""

    def test_single_reaction_loaded(self, dedt_network):
        assert dedt_network.reactions.count == 1

    def test_rate_string(self, dedt_codegen):
        rates = dedt_codegen.get_rates_str().strip().splitlines()
        assert len(rates) == 1
        assert rates[-1].split("=")[-1].strip().rstrip(";") == "nden[0]"

    def test_dedt_expression(self, dedt_network):
        assert str(dedt_network.dEdt_chem) == "nden[0, 0]**3*nden[1, 0]"

    def test_rhs_and_jacobian(self, dedt_codegen):
        rhs = _rhs_terms(dedt_codegen.get_rhs_str(use_cse=False))
        jac = _rhs_terms(dedt_codegen.get_jacobian_str(use_cse=False, use_dedt=True))
        assert rhs == [
            "-std::pow(nden[0], 2)*nden[1]",
            "-std::pow(nden[0], 2)*nden[1]",
            "std::pow(nden[0], 2)*nden[1]",
            "std::pow(nden[0], 3)*nden[1]",
        ]
        assert jac == [
            "-2*nden[0]*nden[1]",
            "-std::pow(nden[0], 2)",
            "-2*nden[0]*nden[1]",
            "-std::pow(nden[0], 2)",
            "2*nden[0]*nden[1]",
            "std::pow(nden[0], 2)",
            "3*std::pow(nden[0], 2)*nden[1]",
            "std::pow(nden[0], 3)",
        ]
