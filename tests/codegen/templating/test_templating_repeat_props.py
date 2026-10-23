# ABOUTME: Tests for REPEAT's property table and its per-property option groups
# ABOUTME: Pure checks on parsed options; no network or codegen is involved

from pathlib import Path

import pytest

from jaff.codegen.templating.commands.repeat.props import (
    GROUPS,
    PROPERTIES,
    option_kwargs,
)
from jaff.codegen.templating.errors import TemplateError
from jaff.codegen.templating.syntax import parse_header

# Every property the old engine's REPEAT offered.
OLD_PROPS = {
    "rates", "flux_expressions", "ode_expressions", "odes", "radodes", "rhses",
    "jacobian", "reactions", "species", "species_with_normalized_sign", "elements",
    "specie_masses", "reactants", "products", "photo_reactions",
    "photo_reaction_truths", "photo_reaction_indices", "specie_charges",
    "neutral_species", "charged_species", "tmins", "tmaxes", "element_density_matrix",
    "element_truth_matrix", "charged_indices", "neutral_indices", "charge_truths",
    "specie_masses_ne", "specie_charges_ne", "charge_truths_ne",
    "neutral_specie_indices", "charged_specie_indices", "charged_specie_indices_ne",
    "charged_specie_masses_ne", "neutral_specie_masses", "charged_specie_masses",
    "radiation_boundaries",
}  # fmt: skip


def options(text: str, line: int = 1):
    """Options of ``$$JAFF REPEAT x IN p { <text> }`` starting at template *line*."""
    return parse_header(f"$$JAFF REPEAT x IN p {{{text}}}", line, Path("t")).options


def test_properties_match_the_old_engine():
    assert set(PROPERTIES) == OLD_PROPS


@pytest.mark.parametrize("name", ["rates", "odes", "radodes", "rhses", "jacobian"])
def test_expression_properties_offer_cse(name):
    assert PROPERTIES[name].special == ("cse",)


def test_item_variable():
    assert PROPERTIES["element_density_matrix"].item == "element"
    assert PROPERTIES["species"].vars == ("idx", "specie")


def test_groups_are_the_union_over_properties():
    assert GROUPS == {
        "THERMAL": {"MODE", "ENERGY"},
        "RADIATION": {"ENABLED", "ORDER"},
        "SIGN": {"POS", "NEG"},
    }


def test_thermal_and_radiation_become_rhs_kwargs():
    opts = options(
        " THERMAL { MODE -> DTDT  ENERGY -> specific } "
        "RADIATION { ENABLED -> true  ORDER -> 1 } "
    )
    assert option_kwargs("rhses", opts) == {
        "thermal": "dtdt",
        "energy": "specific",
        "radiation": True,
        "rad_order": 1,
    }


def test_radodes_order_is_passed_as_order():
    assert option_kwargs("radodes", options(" RADIATION { ORDER -> 2 } ")) == {"order": 2}


def test_sign_becomes_normalized_name_kwargs():
    opts = options(" SIGN { POS -> p  NEG -> m } ")
    assert option_kwargs("species_with_normalized_sign", opts) == {
        "pos": "p",
        "neg": "m",
    }


def test_plain_lists_take_no_option_groups():
    with pytest.raises(
        TemplateError, match="REPEAT IN species takes no option group SORT"
    ):
        option_kwargs("species", options(" SORT { ENABLED -> True } "))


def test_group_a_property_does_not_take_is_an_error_at_its_line():
    with pytest.raises(TemplateError, match="REPEAT IN rates takes no option group") as e:
        option_kwargs("rates", options("\n THERMAL { MODE -> dtdt }\n", line=3))
    assert e.value.nline == 4


def test_key_a_property_does_not_take_is_an_error():
    with pytest.raises(TemplateError, match="RADIATION has no setting ENABLED"):
        option_kwargs("radodes", options(" RADIATION { ENABLED -> True } "))


def test_bad_value_is_an_error_at_the_group_line():
    with pytest.raises(TemplateError, match="RADIATION ENABLED expects True or False"):
        option_kwargs("rhses", options(" RADIATION { ENABLED -> maybe } "))

    with pytest.raises(TemplateError, match="THERMAL ENERGY expects one of") as e:
        option_kwargs("jacobian", options("\n\n THERMAL { ENERGY -> joules }\n", 7))
    assert e.value.nline == 9
