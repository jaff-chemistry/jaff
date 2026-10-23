# ABOUTME: Tests for SUB, GET and HAS validation and rendering on a small network
# ABOUTME: Bodies are rendered line by line with indentation kept and REPLACE applied

from pathlib import Path

import pytest

from jaff import Network
from jaff.codegen.codegen import Codegen
from jaff.codegen.templating.commands import Get, Has, Sub
from jaff.codegen.templating.errors import TemplateError
from jaff.codegen.templating.context import RenderContext
from jaff.codegen.templating.options import Options
from jaff.codegen.templating.resolve import Resolved, compile_replacements
from jaff.codegen.templating.syntax import parse_header
from tests._paths import FIXTURES_DIR

FIXTURES = FIXTURES_DIR


@pytest.fixture(scope="module")
def ctx() -> RenderContext:
    net = Network(str(FIXTURES / "react_cie_hepp.jet"))
    return RenderContext(net, Codegen(network=net, lang="python"), Path("t.py"))


def resolved(opts: Options = None) -> Resolved:
    """Resolved block options for *opts* (none by default)."""
    opts = opts or Options()
    return Resolved(opts, compile_replacements(opts))


def replace(text: str) -> Options:
    """Options holding the REPLACE pairs written in *text*."""
    return parse_header(f"$$JAFF SUB x {{ REPLACE {{ {text} }} }}", 1, Path("t")).options


def test_sub_substitutes_listed_tokens_with_arithmetic(ctx):
    body = ["  n = $nspec$;\n", "    m = $nspec+1$ $nreact$\n"]
    out = Sub(("nspec", "nreact")).render(body, resolved(), ctx)
    assert out == "  n = 5;\n    m = 6 3\n"


def test_sub_leaves_unlisted_tokens(ctx):
    assert Sub(("nspec",)).render(["$nreact$\n"], resolved(), ctx) == "$nreact$\n"


def test_sub_keeps_tab_indentation_and_the_rest_of_the_line(ctx):
    body = ["\tx = $nspec$  \n", "\t\tif x:\r\n", "\n"]
    assert Sub(("nspec",)).render(body, resolved(), ctx) == "\tx = 5  \n\t\tif x:\n\n"


def test_sub_applies_replace_after_substitution(ctx):
    out = Sub(("nspec",)).render(
        ["n_s = $nspec$\n"], resolved(replace(r'"n_" -> "N_" "(\d+)" -> "<\1>"')), ctx
    )
    assert out == "N_s = <5>\n"


def test_sub_filename_and_label(ctx):
    out = Sub(("filename", "label")).render(['"$filename$" "$label$"\n'], resolved(), ctx)
    assert out == '"t.py" "react_cie_hepp"\n'


def test_sub_rejects_unknown_token(ctx):
    with pytest.raises(TemplateError, match="SUB has no token nspecs"):
        Sub(("nspec", "nspecs")).validate([], resolved(), ctx)


@pytest.mark.parametrize(
    "props, entity, line, expected",
    [
        (("specie_idx",), "H+", "i = $specie_idx+1$", "i = 3"),
        (("specie_charge",), "He++", "q = $specie_charge$", "q = 2"),
        (("element_idx",), "He", "e = $element_idx$", "e = 1"),
        (("reaction_verbatim",), "H+ + e- -> H", "$reaction_verbatim$", "H+ + e- -> H"),
    ],
)
def test_get_substitutes_entity_properties(ctx, props, entity, line, expected):
    assert Get(props, entity).render([line + "\n"], resolved(), ctx) == expected + "\n"


def test_get_rejects_unknown_property(ctx):
    with pytest.raises(TemplateError, match="GET has no property mass"):
        Get(("mass",), "H").validate([], resolved(), ctx)


@pytest.mark.parametrize(
    "kind, entity, expected",
    [("specie", "He+", "1"), ("specie", "CO", "0"), ("element", "C", "0")],
)
def test_has_substitutes_one_or_zero(ctx, kind, entity, expected):
    out = Has(kind, entity).render([f"x = ${kind}$\n"], resolved(), ctx)
    assert out == f"x = {expected}\n"


def test_has_rejects_unknown_kind(ctx):
    with pytest.raises(TemplateError, match="HAS has no kind species"):
        Has("species", "H").validate([], resolved(), ctx)


@pytest.mark.parametrize(
    "prop, entity, message",
    [
        ("specie_idx", "Xyz", "GET specie_idx needs a species; Xyz is not one"),
        ("element_idx", "Xyz", "GET element_idx needs an element; Xyz is not one"),
        (
            "reaction_tmin",
            "H -> H2",
            "GET reaction_tmin needs a reaction; H -> H2 is not one",
        ),
        ("aux_func", "nope", "GET aux_func needs an auxiliary function; nope is not one"),
        ("reaction_tmin", "H", "GET reaction_tmin needs a reaction; H is not one"),
        ("specie_idx", "H + e- -> H+ + e- + e-", "needs a species"),
        ("element_idx", "H+", "needs an element"),
    ],
)
def test_get_rejects_entity_of_the_wrong_kind_or_unknown(ctx, prop, entity, message):
    with pytest.raises(TemplateError, match=message):
        Get((prop,), entity).validate([], resolved(), ctx)


@pytest.mark.parametrize(
    "prop, entity",
    [
        ("specie_idx", "H+"),
        ("element_idx", "He"),
        ("reaction_tmin", "H+ + e- -> H"),
        ("aux_func", "chemrate2"),
    ],
)
def test_get_accepts_existing_entities(ctx, prop, entity):
    Get((prop,), entity).validate([], resolved(), ctx)
