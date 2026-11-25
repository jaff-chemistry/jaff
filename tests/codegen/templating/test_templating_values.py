# ABOUTME: Tests for the option-value parsers used by commands' option groups
# ABOUTME: Booleans and choices are case-insensitive; bad values raise ValueError

import pytest

from jaff.codegen.templating.values import (
    parse_bool,
    parse_choice,
    parse_int,
    parse_text,
)


@pytest.mark.parametrize(
    "text, expected",
    [("True", True), ("TRUE", True), ("true", True), ("False", False), ("fAlSe", False)],
)
def test_bool_is_case_insensitive(text, expected):
    assert parse_bool(text) is expected


@pytest.mark.parametrize("text", ["yes", "1", "", " true"])
def test_bool_rejects_other_text(text):
    with pytest.raises(ValueError, match="expects True or False"):
        parse_bool(text)


def test_int_parses_decimal_integers():
    assert parse_int("2") == 2
    assert parse_int("-1") == -1


@pytest.mark.parametrize("text", ["one", "", "1.5"])
def test_int_rejects_other_text(text):
    with pytest.raises(ValueError, match=f"expects an integer, got {text!r}"):
        parse_int(text)


def test_text_is_unchanged():
    assert parse_text("Ab_c") == "Ab_c"


def test_choice_is_case_insensitive_and_lowercase():
    parse = parse_choice("none", "dedt", "dtdt")
    assert parse("DtDt") == "dtdt"


def test_choice_lists_the_choices_on_error():
    with pytest.raises(ValueError, match="expects one of none, dedt, dtdt, got 'x'"):
        parse_choice("none", "dedt", "dtdt")("x")


def test_choice_list_in_any_case_is_matched_lowercase():
    assert parse_choice("None", "DEDT")("dedt") == "dedt"
