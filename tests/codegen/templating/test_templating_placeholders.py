# ABOUTME: Tests for $name$ / $name op N$ substitution and body-line splitting
# ABOUTME: Pure string inputs: no network or codegen is involved

import pytest

from jaff.codegen.templating.placeholders import split_line, substitute


def const(value):
    return lambda: value


@pytest.mark.parametrize(
    "text, expected",
    [
        ("n = $nspec$", "n = 5"),
        ("n = $nspec+1$", "n = 6"),
        ("n = $nspec - 1$", "n = 4"),
        ("n = $nspec*2$", "n = 10"),
        ("n = $nspec/2$", "n = 2"),
        ("n = $nspec $", "n = 5"),
        ("$nspec$ and $nspec+1$", "5 and 6"),
        ("$nspecx$ $other$", "$nspecx$ $other$"),
        ("$nspec,1$", "$nspec,1$"),
    ],
)
def test_substitute_integer_token(text, expected):
    assert substitute(text, {"nspec": const(5)}) == expected


def test_arithmetic_is_ignored_for_non_integers():
    assert substitute("$label+1$", {"label": const("net")}) == "net"


def test_values_are_computed_only_when_used():
    def boom():
        raise AssertionError("evaluated")

    assert substitute("x = $a$", {"a": const(1), "b": boom}) == "x = 1"


def test_longer_name_sharing_a_prefix_still_matches():
    values = {"specie": const("H"), "specie_idx": const(0)}
    assert substitute("$specie$ $specie_idx$", values) == "H 0"


def test_no_names_leaves_text_unchanged():
    assert substitute("$a$", {}) == "$a$"


@pytest.mark.parametrize(
    "line, expected",
    [
        ("    x = 1;\n", ("    ", "x = 1;")),
        ("x  \n", ("", "x  ")),
        ("\tx\n", ("\t", "x")),
        ("  \t x\r\n", ("  \t ", "x")),
        ("\n", ("", "")),
        ("   \n", ("   ", "")),
    ],
)
def test_split_line(line, expected):
    assert split_line(line) == expected
