# ABOUTME: Tests for the REDUCE command of the new template engine
# ABOUTME: $( ... )$ expands to a parenthesised sum; other lines pass through unchanged

from pathlib import Path

import pytest

from jaff import Network
from jaff.codegen.codegen import Codegen
from jaff.codegen.templating.commands import Reduce
from jaff.codegen.templating.context import RenderContext
from jaff.codegen.templating.errors import TemplateError
from jaff.codegen.templating.options import Options
from jaff.codegen.templating.resolve import Resolved, compile_replacements
from jaff.codegen.templating.syntax import parse_header
from tests._paths import FIXTURES_DIR

FIXTURES = FIXTURES_DIR


@pytest.fixture(scope="module")
def ctx() -> RenderContext:
    # species H, e-, H+, He++, He+ with charges 0, -1, 1, 2, 1
    net = Network(str(FIXTURES / "react_cie_hepp.jet"))
    return RenderContext(net, Codegen(network=net, lang="python"), Path("t.py"))


def resolved(replace: str = "") -> Resolved:
    """Resolved options, with the REPLACE pairs written in *replace* if any."""
    opts = Options()
    if replace:
        header = f"$$JAFF SUB x {{ REPLACE {{ {replace} }} }}"
        opts = parse_header(header, 1, Path("t")).options
    return Resolved(opts, compile_replacements(opts))


def charges() -> Reduce:
    return Reduce(("specie_charge",), ("specie_charges",))


def test_reduction_sums_each_value_in_parentheses(ctx):
    out = charges().render(["  q = $($specie_charge$)$\n"], resolved(), ctx)
    assert out == "  q = ((0) + (-1) + (1) + (2) + (1))\n"


def test_summand_template_is_repeated_per_value(ctx):
    out = charges().render(["q = 2 * $($specie_charge$ * 3)$\n"], resolved(), ctx)
    assert out == "q = 2 * ((0 * 3) + (-1 * 3) + (1 * 3) + (2 * 3) + (1 * 3))\n"


def test_lines_without_a_reduction_are_copied_with_replace(ctx):
    body = ["  x = $specie_charge$  \n", "\ty\n"]
    out = charges().render(body, resolved('"x" -> "X"'), ctx)
    assert out == "  X = $specie_charge$  \n\ty\n"


def test_several_properties_are_zipped(ctx):
    command = Reduce(("specie_charge",), ("specie_charges", "charge_truths"))
    out = command.render(["$($specie_charge$*$charge_truth$)$\n"], resolved(), ctx)
    assert out == "((0*0) + (-1*1) + (1*1) + (2*1) + (1*1))\n"


def test_replace_runs_on_expanded_lines(ctx):
    out = charges().render(
        ["q = $($specie_charge$)$\n"], resolved('"\\+" -> "plus"'), ctx
    )
    assert out == "q = ((0) plus (-1) plus (1) plus (2) plus (1))\n"


def test_empty_property_reduces_to_zero(ctx):
    command = Reduce(("photo_reaction_index",), ("photo_reaction_indices",))
    assert command.render(["n = $($photo_reaction_index$)$\n"], resolved(), ctx) == (
        "n = 0\n"
    )


def test_properties_of_different_length_are_an_error_at_the_line(ctx):
    command = Reduce(("specie_charge",), ("specie_charges", "neutral_specie_masses"))
    body = ["a\n", "q = $($specie_charge$)$\n"]
    with pytest.raises(TemplateError, match="differ in length") as err:
        command.render(body, resolved(), ctx)
    assert err.value.body_index == 1


def test_unknown_property_is_rejected(ctx):
    with pytest.raises(TemplateError, match="REDUCE has no property charges"):
        Reduce(("specie_charge",), ("charges",)).validate([], resolved(), ctx)


def test_variable_must_belong_to_a_listed_property(ctx):
    with pytest.raises(TemplateError, match="REDUCE has no variable specie_mass"):
        Reduce(("specie_mass",), ("specie_charges",)).validate([], resolved(), ctx)


def test_every_reduction_on_a_line_is_expanded(ctx):
    body = ["a = $($specie_charge$)$ + $($specie_charge$*2)$ + $($other$)$\n"]
    total = "((0) + (-1) + (1) + (2) + (1))"
    doubled = "((0*2) + (-1*2) + (1*2) + (2*2) + (1*2))"
    assert charges().render(body, resolved(), ctx) == (
        f"a = {total} + {doubled} + $($other$)$\n"
    )


def test_reduction_without_a_listed_variable_is_copied_with_replace(ctx):
    body = ["  m = $($specie_mass$)$\n"]
    out = charges().render(body, resolved('"m =" -> "M ="'), ctx)
    assert out == "  M = $($specie_mass$)$\n"
