# ABOUTME: Tests for the $$JAFF directive header grammar and its transformer
# ABOUTME: Pure string inputs: no network or codegen is involved

from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest

from jaff.codegen.templating.commands import (
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
from jaff.codegen.templating.syntax import (
    DirectiveSyntaxError,
    IncompleteHeaderError,
    Options,
    SplitOptionError,
    parse_header,
)
from jaff.errors import ParserError

FILE = Path("template.cpp")


def parse(text: str, line: int = 1):
    """Parse *text* as the header following a ``$$JAFF`` marker."""
    return parse_header(f"$$JAFF {text}", line, FILE)


@pytest.mark.parametrize(
    "text, command",
    [
        ("REPEAT idx, rhs, cse IN rhses", Repeat(("idx", "rhs", "cse"), "rhses")),
        ("SUB nspec, label", Sub(("nspec", "label"))),
        (
            "GET specie_idx, specie_mass IN H+",
            Get(("specie_idx", "specie_mass"), "H+"),
        ),
        (
            'GET reaction_idx IN "H + H -> H2"',
            Get(("reaction_idx",), "H + H -> H2"),
        ),
        ("HAS specie e-", Has("specie", "e-")),
        ("REDUCE a, b IN p, q", Reduce(("a", "b"), ("p", "q"))),
        ("DEFINE micro { }", Define("micro")),
        ("END", End()),
        ("LANG c++", Lang("c++")),
        ('COMMENT "//"', Comment("//")),
    ],
)
def test_commands(text, command):
    assert parse(text).command == command


def test_no_block_means_no_options():
    assert parse("SUB nspec").options is None


def test_empty_block_gives_empty_options():
    assert parse("SUB nspec { }").options == Options()


def test_prop_named_like_a_keyword_suffix():
    # Old engine split on the substring "IN"; INF must stay one name (bug 5)
    assert parse("REPEAT idx IN INF").command == Repeat(("idx",), "INF")


def test_replace_pairs_keep_regex_escapes():
    options = parse(
        r'REPEAT idx IN x { REPLACE { "nden\[\s*(\d+)\s*\]" -> "state.xn[\1]" } }'
    ).options
    assert options.replace == [(r"nden\[\s*(\d+)\s*\]", r"state.xn[\1]")]


def test_escaped_quote_in_replacement_becomes_quote():
    options = parse(r'SUB x { REPLACE { "a" -> "say \"hi\"" } }').options
    assert options.replace == [("a", 'say "hi"')]


def test_replace_blocks_concatenate_in_order():
    options = parse('SUB x { REPLACE { "a" -> "b" } REPLACE { "c" -> "d" } }').options
    assert options.replace == [("a", "b"), ("c", "d")]


def test_groups_and_uses():
    options = parse(
        "REPEAT idx IN x {\n"
        "    USE micro\n"
        "    THERMAL { MODE -> dtdt  ENERGY -> specific }\n"
        "    USE other\n"
        "}"
    ).options
    assert options.uses == ["micro", "other"]
    assert options.groups == {"THERMAL": {"MODE": "dtdt", "ENERGY": "specific"}}


def test_quoted_setting_value_loses_quotes():
    options = parse('SUB x { SIGN { POS -> "p" } }').options
    assert options.groups == {"SIGN": {"POS": "p"}}


def test_comments_are_ignored_but_not_inside_strings():
    options = parse(
        "SUB x {  # trailing\n"
        "    REPLACE {\n"
        '        # "tgas" -> "T"\n'
        '        "a#b" -> "c"   # note\n'
        "    }\n"
        "}"
    ).options
    assert options.replace == [("a#b", "c")]


@pytest.mark.parametrize(
    "text",
    [
        'SUB x { REPLACE { tgas -> "T" } }',  # unquoted pattern
        "SUB x { REPLACE { 'a' -> 'b' } }",  # single quotes
        "SUB x { RADIATION -> True }",  # top-level setting
        'SUB x { DEDT { REPLACE { "a" -> "b" } } }',  # nested group
        "SUB x { DEDT { ENABLED } }",  # setting without value
        "REPEAT idx rhses",  # missing IN
        "FOO bar",  # unknown command
        'SUB x { REPLACE { "a -> "b" } }',  # unterminated string
        "SUB x { SORT { ENABLED -> True }",  # missing closing brace
        "",  # empty text
    ],
)
def test_syntax_errors_raise_parser_error(text):
    with pytest.raises(ParserError, match="Invalid directive"):
        parse(text)


def test_error_line_is_offset_into_template():
    with pytest.raises(ParserError) as err:
        parse("SUB x {\n    DEDT { ENABLED }\n}", line=10)
    assert err.value.nline == 11


def test_repeated_group_is_an_error():
    with pytest.raises(ParserError, match="SORT given twice"):
        parse("SUB x { SORT { ENABLED -> True } SORT { ENABLED -> False } }")


def test_repeated_key_is_an_error():
    with pytest.raises(ParserError, match="sets ENABLED more than once"):
        parse("SUB x { SORT { ENABLED -> True ENABLED -> False } }")


def test_repeated_group_reports_its_own_line():
    text = "SUB x {\n  SORT { ENABLED -> True }\n  SORT { ENABLED -> False }\n}"
    with pytest.raises(ParserError, match="SORT given twice") as err:
        parse(text, line=10)
    assert err.value.nline == 12
    assert err.value.line.strip() == "SORT { ENABLED -> False }"


def test_repeated_key_reports_its_own_line():
    text = "SUB x {\n  SORT {\n    ENABLED -> True\n    ENABLED -> False\n  }\n}"
    with pytest.raises(ParserError, match="more than once") as err:
        parse(text, line=10)
    assert err.value.nline == 13
    assert err.value.line.strip() == "ENABLED -> False"


def test_lark_error_message_is_clean():
    with pytest.raises(ParserError) as err:
        parse("SUB x {\n    DEDT { ENABLED }\n}", line=10)
    message = str(err.value)
    assert "column" in message
    assert "Token(" not in message
    assert "at line" not in message


def test_unexpected_end_of_header_is_named():
    with pytest.raises(ParserError, match="end of header"):
        parse("SUB x {")


def test_quoted_setting_value_unescapes_quotes():
    options = parse(r'SUB x { SIGN { POS -> "v\"w" } }').options
    assert options.groups == {"SIGN": {"POS": 'v"w'}}


def test_directive_and_headers_are_immutable():
    directive = parse("SUB a, b")
    with pytest.raises(FrozenInstanceError):
        directive.options = Options()
    assert isinstance(directive.command.tokens, tuple)


@pytest.mark.parametrize("text", ["SUB x {", 'SUB x { REPLACE { "a" -> "b"'])
def test_end_inside_open_block_is_incomplete(text):
    with pytest.raises(IncompleteHeaderError):
        parse(text)


@pytest.mark.parametrize(
    "text", ["REPEAT idx", "SUB x { SORT", 'SUB x { REPLACE { "a" ->', "SUB x { USE", ""]
)
def test_other_end_of_input_is_final_error(text):
    with pytest.raises(ParserError) as err:
        parse(text)
    assert not isinstance(err.value, IncompleteHeaderError)


def test_error_shows_full_source_line_and_absolute_column():
    text = "    $$JAFF SUB x { SORT { K } }"
    with pytest.raises(ParserError) as err:
        parse_header(text, 7, FILE)
    assert err.value.line == text
    assert err.value.nline == 7
    assert f"column {text.index('}') + 1}" in err.value.message


def test_marker_needs_a_space_after_it():
    with pytest.raises(ParserError, match="Expected a space after"):
        parse_header("$$JAFFREPEAT idx IN x", 1, FILE)


def test_crlf_header_error_line_has_no_carriage_return():
    with pytest.raises(ParserError) as err:
        parse_header("$$JAFF SUB x {\r\n  SORT { K }\r\n}", 1, FILE)
    assert err.value.line == "  SORT { K }"


def test_grammar_errors_are_directive_syntax_errors():
    with pytest.raises(DirectiveSyntaxError):
        parse("REPEAT idx rhses")


def test_repeated_group_is_not_a_syntax_error():
    with pytest.raises(ParserError) as err:
        parse("SUB x { SORT { A -> 1 } SORT { A -> 2 } }")
    assert not isinstance(err.value, DirectiveSyntaxError)


def test_options_record_template_lines_of_uses_and_groups():
    options = parse(
        "SUB x {\n    USE micro\n    SORT { ENABLED -> True }\n}", line=10
    ).options
    assert options.use_lines == {"micro": 11}
    assert options.group_lines == {"SORT": 12}


def test_option_lines_do_not_affect_equality():
    assert (
        parse("SUB x { USE m }", line=1).options
        == parse("SUB x {\n USE m }", line=5).options
    )


@pytest.mark.parametrize(
    "text, command",
    [
        (r'GET reaction_idx IN "a \"b\""', Get(("reaction_idx",), 'a "b"')),
        ('HAS reaction "H + H -> H2"', Has("reaction", "H + H -> H2")),
    ],
)
def test_quoted_command_values_unescape_quotes(text, command):
    assert parse(text).command == command


def test_bare_marker_reports_missing_command():
    with pytest.raises(DirectiveSyntaxError, match="missing command after"):
        parse_header("  $$JAFF   ", 1, FILE)


@pytest.mark.parametrize(
    "text",
    [
        "SUB x { SORT { K ->",
        "SUB x { SORT",
        "SUB x { USE",
        'SUB x { REPLACE { "a"',
    ],
)
def test_line_ending_inside_an_option_is_a_split_option_error(text):
    with pytest.raises(SplitOptionError, match="keep each pair"):
        parse(text)


def test_line_ending_outside_a_block_is_a_plain_syntax_error():
    with pytest.raises(DirectiveSyntaxError) as err:
        parse("REPEAT idx IN")
    assert not isinstance(err.value, SplitOptionError)


@pytest.mark.parametrize(
    "text",
    ['$$JAFF SUB x { REPLACE { "a -> "b" } }', "$$JAFFSUB x"],
)
def test_lexer_and_marker_errors_are_directive_syntax_errors(text):
    with pytest.raises(DirectiveSyntaxError):
        parse_header(text, 1, FILE)


def test_tab_indented_marker_keeps_absolute_column():
    text = "\t$$JAFF SUB x { SORT { K } }"
    with pytest.raises(ParserError) as err:
        parse_header(text, 1, FILE)
    assert f"column {text.index('}') + 1}" in err.value.message


def test_crlf_header_parses():
    directive = parse_header("$$JAFF SUB x {\r\n  USE m\r\n}\r\n", 1, FILE)
    assert directive.options.uses == ["m"]


def test_options_record_template_lines_of_replace_pairs():
    options = parse(
        'SUB x {\n  REPLACE {\n    "a" -> "b"\n    "c" -> "d"\n  }\n}', line=20
    ).options
    assert options.replace_lines == [22, 23]


def test_repeated_use_reports_its_own_line():
    text = "SUB x {\n  USE m\n  USE m\n}"
    with pytest.raises(ParserError, match="USE m given twice") as err:
        parse(text, line=10)
    assert err.value.nline == 12
