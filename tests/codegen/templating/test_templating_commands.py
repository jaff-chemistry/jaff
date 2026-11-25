# ABOUTME: Tests for the $$JAFF command classes and the registry the syntax layer uses
# ABOUTME: Checks registry keys match the grammar's aliases and each command's shape

import re
from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest

import jaff.codegen.templating as templating
from jaff.codegen.templating.commands import (
    COMMANDS,
    Comment,
    Define,
    End,
    Get,
    Has,
    Lang,
    Reduce,
    Repeat,
    Sub,
)

GRAMMAR = Path(templating.__file__).with_name("grammar.lark")


def _command_aliases() -> set:
    """Aliases of the grammar's ``command`` rule alternatives."""
    rule = re.search(r"^command:(.*?)^\s*$", GRAMMAR.read_text(), re.M | re.S)
    return set(re.findall(r"->\s*(\w+)", rule.group(1)))


def test_registry_keys_match_grammar_command_aliases():
    assert set(COMMANDS) == _command_aliases()


def test_command_aliases_do_not_shadow_transformer_rules():
    # Explicit transformer methods win over __default__, so a command alias with
    # one of these names would never reach the registry
    rules = {"start", "names", "replacement", "setting", "replace", "use", "group"}
    assert not rules & set(COMMANDS)


@pytest.mark.parametrize(
    "alias, items, command",
    [
        ("repeat", [("idx", "rhs"), "rhses"], Repeat(("idx", "rhs"), "rhses")),
        ("sub", [("nspec",)], Sub(("nspec",))),
        ("get", [("specie_idx",), "H+"], Get(("specie_idx",), "H+")),
        ("has", ["specie", "e-"], Has("specie", "e-")),
        ("reduce", [("a",), ("p",)], Reduce(("a",), ("p",))),
        ("define", ["micro"], Define("micro")),
        ("end", [], End()),
        ("lang", ["fortran"], Lang("fortran")),
        ("comment", ["#"], Comment("#")),
    ],
)
def test_from_items(alias, items, command):
    assert COMMANDS[alias].from_items(items) == command


@pytest.mark.parametrize(
    "cls, has_body",
    [
        (Repeat, True),
        (Sub, True),
        (Get, True),
        (Has, True),
        (Reduce, True),
        (Define, False),
        (End, False),
    ],
)
def test_has_body(cls, has_body):
    assert cls.has_body is has_body


def test_commands_are_immutable():
    with pytest.raises(FrozenInstanceError):
        Repeat(("idx",), "x").prop = "y"
