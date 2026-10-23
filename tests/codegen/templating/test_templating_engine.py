# ABOUTME: Tests for the template engine's single pass splitting text and directive blocks
# ABOUTME: Pure string inputs: no network or codegen is involved

from pathlib import Path

import pytest

from jaff.codegen.templating.blocks import CommandBlock, TextBlock
from jaff.codegen.templating.commands import Comment, Define, Lang, Repeat, Sub
from jaff.codegen.templating.engine import split_blocks
from jaff.errors import ParserError

FILE = Path("template.cpp")


def split(text: str):
    return split_blocks(text, FILE)


def test_plain_text_is_one_block():
    blocks = split("int a;\n// comment $nspec$\n")
    assert blocks == [TextBlock(["int a;\n", "// comment $nspec$\n"], 1)]


def test_single_line_directive_with_body():
    blocks = split("before\n$$JAFF SUB nspec\nn = $nspec$;\n$$JAFF END\nafter\n")
    assert [type(b) for b in blocks] == [TextBlock, CommandBlock, TextBlock]
    block = blocks[1]

    assert block.command == Sub(("nspec",))
    assert block.header_lines == ["$$JAFF SUB nspec\n"]
    assert block.body_lines == ["n = $nspec$;\n"]
    assert block.line == 2
    assert blocks[2] == TextBlock(["after\n"], 5)


def test_multiline_header_and_indented_body():
    text = (
        "  $$JAFF REPEAT idx, rhs IN rhses {\n"
        '      REPLACE { "tgas" -> "T" }\n'
        "  }\n"
        "  f[$idx$] = $rhs$;\n"
        "\n"
        "  $$JAFF END\n"
    )
    (block,) = split(text)
    assert block.command == Repeat(("idx", "rhs"), "rhses")
    assert block.options.replace == [("tgas", "T")]
    assert len(block.header_lines) == 3
    assert block.body_lines == ["  f[$idx$] = $rhs$;\n", "\n"]


def test_brace_on_header_line_only():
    # A { on the next line is body text, not an option block
    (block,) = split("$$JAFF SUB nspec\n{ $nspec$ }\n$$JAFF END\n")
    assert block.options is None
    assert block.body_lines == ["{ $nspec$ }\n"]


def test_define_has_no_body():
    blocks = split("$$JAFF DEFINE micro {\n    RADIATION { ENABLED -> True }\n}\nx\n")
    assert blocks[0].command == Define("micro")
    assert blocks[0].body_lines == []
    assert blocks[1] == TextBlock(["x\n"], 4)


def test_end_may_carry_a_comment():
    (block,) = split("$$JAFF SUB nspec\nx\n$$JAFF END  # sub\n")
    assert block.body_lines == ["x\n"]


@pytest.mark.parametrize(
    "text, message, nline",
    [
        ("$$JAFF SUB nspec\nx\n", "Missing \\$\\$JAFF END", 1),
        ("$$JAFF SUB x {\n  SORT { ENABLED -> True }\n", "Unclosed", 1),
        ("$$JAFF END\n", "END without a directive", 1),
        ("$$JAFF SUB x\n$$JAFF SUB y\n$$JAFF END\n", "Nested", 2),
        ("$$JAFF SUB x\n$$JAFF END { USE m }\n", "END takes no options", 2),
        ("a\n$$JAFF SUB x {\n  BAD { K }\n}\n", "Invalid directive", 3),
    ],
)
def test_structure_errors(text, message, nline):
    with pytest.raises(ParserError, match=message) as err:
        split(text)
    assert err.value.nline == nline


def test_unclosed_header_stops_at_next_directive():
    text = "$$JAFF SUB x {\n  SORT { ENABLED -> True }\n$$JAFF END\n"
    with pytest.raises(ParserError, match="Unclosed") as err:
        split(text)
    assert err.value.nline == 1


def test_unclosed_header_swallowing_body_line_reports_it_as_invalid():
    text = "$$JAFF SUB x {\n  SORT { ENABLED -> True }\ncode;\n$$JAFF END\n"
    with pytest.raises(ParserError, match="Invalid directive") as err:
        split(text)
    assert err.value.nline == 3
    assert "opened at line 1" in str(err.value)


def test_unclosed_header_stops_at_directive_line():
    text = "$$JAFF SUB x {\n$$JAFF END\n"
    with pytest.raises(ParserError, match="Unclosed") as err:
        split(text)
    assert err.value.nline == 1


def test_form_feed_is_not_a_line_break():
    with pytest.raises(ParserError) as err:
        split("a\x0cb\n$$JAFF END\n")
    assert err.value.nline == 2
    assert split("a\x0cb\n") == [TextBlock(["a\x0cb\n"], 1)]


@pytest.mark.parametrize(
    "text",
    [
        "$$JAFF SUB x\nx\n$$JAFF SUB y\n$$JAFF END\n",
        "$$JAFF END\r\n",
        "$$JAFF SUB x\r\n",
    ],
)
def test_error_line_has_no_trailing_newline(text):
    with pytest.raises(ParserError) as err:
        split(text)
    assert not err.value.line.endswith(("\n", "\r"))


def test_empty_template():
    assert split("") == []


def test_crlf_is_preserved():
    (block,) = split("$$JAFF SUB nspec\r\nx\r\n$$JAFF END\r\n")
    assert block.command == Sub(("nspec",))
    assert block.header_lines == ["$$JAFF SUB nspec\r\n"]
    assert block.body_lines == ["x\r\n"]


def test_no_trailing_newline_after_end_and_define():
    (block,) = split("$$JAFF SUB x\nb\n$$JAFF END")
    assert block.body_lines == ["b\n"]
    blocks = split("$$JAFF DEFINE m {\n  RADIATION { ENABLED -> True }\n}")
    assert blocks[0].command == Define("m")


def test_missing_end_reports_header_line_after_text():
    with pytest.raises(ParserError, match="Missing") as err:
        split("a\nb\n$$JAFF SUB x\nbody\n")
    assert err.value.nline == 3


def test_lone_carriage_return_is_not_a_line_break():
    assert split("a\rb\n") == [TextBlock(["a\rb\n"], 1)]


def test_unclosed_header_not_balanced_by_later_brace():
    text = "$$JAFF SUB x {\n  SORT { ENABLED -> True }\n$$JAFF END\n}\n"
    with pytest.raises(ParserError, match="Unclosed") as err:
        split(text)
    assert err.value.nline == 1


def test_header_without_brace_does_not_swallow_body():
    with pytest.raises(ParserError, match="Invalid directive") as err:
        split("$$JAFF REPEAT idx\nf[$idx$] = 1;\n$$JAFF END\n")
    assert err.value.nline == 1
    assert "f[" not in err.value.line
    assert "opened at line" not in str(err.value)


def test_group_brace_on_next_line_is_invalid_not_unclosed():
    text = "$$JAFF SUB x {\n  SORT\n  { ENABLED -> True }\n}\n$$JAFF END\n"
    with pytest.raises(ParserError, match="Invalid directive"):
        split(text)


def test_braces_in_strings_and_comments_across_lines():
    text = '$$JAFF SUB x {\n  REPLACE { "a{" -> "}b" }  # } {\n}\nbody\n$$JAFF END\n'
    (block,) = split(text)
    assert isinstance(block, CommandBlock)
    assert list(block.options.replace) == [("a{", "}b")]
    assert block.body_lines == ["body\n"]


def test_duplicate_group_across_lines_hint_has_no_missing_brace():
    text = (
        "$$JAFF SUB x {\n  SORT { ENABLED -> True }\n"
        "  SORT { ENABLED -> False }\n}\nb\n$$JAFF END\n"
    )
    with pytest.raises(ParserError) as err:
        split(text)
    assert "opened at line 1" in str(err.value)
    assert "missing" not in str(err.value)


@pytest.mark.parametrize(
    "text, message",
    [
        ("$$JAFF DEFINE micro\n", "DEFINE needs an option block"),
        ("$$JAFF DEFINE micro { USE other }\n", "DEFINE cannot contain USE"),
    ],
)
def test_define_rules(text, message):
    with pytest.raises(ParserError, match=message) as err:
        split(text)
    assert err.value.nline == 1


def test_empty_define_block_is_allowed():
    (block,) = split("$$JAFF DEFINE micro { }\n")
    assert block.command == Define("micro")
    assert block.body_lines == []


def test_malformed_nested_directive_reports_its_own_error():
    # A nested $$JAFF line is parsed before the nesting check, so its own syntax
    # error wins over "Nested $$JAFF directive"
    with pytest.raises(ParserError, match="Invalid directive") as err:
        split("$$JAFF SUB x\n$$JAFF SUB {\n$$JAFF END\n")
    assert err.value.nline == 2


def test_end_line_is_kept_with_its_indentation():
    (block,) = split("  $$JAFF SUB x\nbody\n    $$JAFF END  # done\n")
    assert block.end_lines == ["    $$JAFF END  # done\n"]


def test_define_has_no_end_lines():
    (block,) = split("$$JAFF DEFINE m { }\n")
    assert block.end_lines == []


def test_error_on_continuation_line_reports_that_line():
    with pytest.raises(ParserError) as err:
        split("$$JAFF SUB x {\n  SORT { ENABLED -> True }\n  BAD { K }\n}\n")
    assert err.value.nline == 3
    assert err.value.line == "  BAD { K }"


def test_split_option_on_continuation_line_has_no_brace_hint():
    with pytest.raises(ParserError) as err:
        split("$$JAFF SUB x {\n  SORT { K ->\n  V }\n}\n")
    assert "keep each pair" in err.value.message
    assert "missing" not in err.value.message


def test_define_use_error_points_at_the_use_line():
    with pytest.raises(ParserError, match="DEFINE cannot contain USE") as err:
        split("$$JAFF DEFINE m {\n  SORT { ENABLED -> True }\n  USE other\n}\n")
    assert err.value.nline == 3
    assert err.value.line == "  USE other"


def test_define_brace_on_next_line_needs_an_option_block():
    with pytest.raises(ParserError, match="DEFINE needs an option block"):
        split("$$JAFF DEFINE m\n{ SORT { ENABLED -> True } }\n")


def test_lang_and_comment_are_header_only_blocks():
    text = '#!/bin/sh\n$$JAFF LANG fortran\n$$JAFF COMMENT "#"\nx\n'
    blocks = split(text)
    assert [type(b) for b in blocks] == [TextBlock, CommandBlock, CommandBlock, TextBlock]
    assert blocks[1].command == Lang("fortran")
    assert blocks[2].command == Comment("#")
    assert blocks[2].body_lines == []


def test_file_setting_after_another_directive_is_an_error():
    text = "$$JAFF SUB nspec\nx\n$$JAFF END\n$$JAFF LANG python\n"
    with pytest.raises(ParserError, match="must come before the other") as err:
        split(text)
    assert err.value.nline == 4


def test_file_setting_given_twice_is_an_error():
    with pytest.raises(ParserError, match="COMMENT given twice") as err:
        split('$$JAFF COMMENT "#"\n$$JAFF COMMENT "!"\n')
    assert err.value.nline == 2


def test_file_setting_takes_no_options():
    with pytest.raises(ParserError, match="LANG takes no options"):
        split("$$JAFF LANG python { USE m }\n")


@pytest.mark.parametrize("token", ['""', '" "'])
def test_comment_needs_a_token(token):
    with pytest.raises(ParserError, match="COMMENT needs a comment token"):
        split(f"$$JAFF COMMENT {token}\n")
