# ABOUTME: Tests for rendering whole templates with the new engine's TemplateParser
# ABOUTME: Echo of directives, DEFINE/USE layering, validation errors and the facade


import pytest

from jaff import Network
from jaff.codegen.templating import TemplateParser
from jaff.codegen.templating.commands import Sub
from jaff.codegen.templating.errors import TemplateError
from jaff.errors import ParserError
from tests._paths import FIXTURES_DIR

FIXTURES = FIXTURES_DIR


@pytest.fixture(scope="module")
def net() -> Network:
    return Network(str(FIXTURES / "react_cie_hepp.jet"))


def render(net, tmp_path, text, name="t.py", echo=False, **kw) -> str:
    template = tmp_path / name
    template.write_text(text)
    return TemplateParser(net, template, echo_directives=echo, **kw).parse_file()


def test_text_outside_directives_is_copied_verbatim(net, tmp_path):
    text = "a = $nspec$\n\t# comment  \n"
    assert render(net, tmp_path, text) == text


def test_sub_block_without_echo(net, tmp_path):
    text = "x\n$$JAFF SUB nspec\nn = $nspec$\n$$JAFF END\ny\n"
    assert render(net, tmp_path, text) == "x\nn = 5\ny\n"


@pytest.mark.parametrize(
    "name, comment", [("t.py", "#"), ("t.cpp", "//"), ("t.f90", "!")]
)
def test_echo_comments_header_and_end_in_the_language(net, tmp_path, name, comment):
    text = (
        '  $$JAFF SUB nspec {\n      REPLACE { "5" -> "five" }\n  }\n'
        "n = $nspec$\n  $$JAFF END\n"
    )
    out = render(net, tmp_path, text, name=name, echo=True)
    assert out == (
        f"  {comment} $$JAFF SUB nspec {{\n"
        f'      {comment} REPLACE {{ "5" -> "five" }}\n'
        f"  {comment} }}\n"
        "n = five\n"
        f"  {comment} $$JAFF END\n"
    )


def test_echo_adds_missing_final_newline(net, tmp_path):
    out = render(net, tmp_path, "$$JAFF DEFINE m { }", echo=True)
    assert out == "# $$JAFF DEFINE m { }\n"


def test_use_pulls_in_a_define(net, tmp_path):
    text = (
        '$$JAFF DEFINE up {\n    REPLACE { "n" -> "N" }\n}\n'
        "$$JAFF SUB nspec { USE up }\nn = $nspec$\n$$JAFF END\n"
    )
    assert render(net, tmp_path, text) == "N = 5\n"


def test_use_before_define_is_an_error_at_the_use_line(net, tmp_path):
    text = (
        "$$JAFF SUB nspec {\n    USE up\n}\nx\n$$JAFF END\n"
        '$$JAFF DEFINE up { REPLACE { "n" -> "N" } }\n'
    )
    with pytest.raises(ParserError, match="undefined option set up") as err:
        render(net, tmp_path, text)
    assert err.value.nline == 2
    assert err.value.line == "    USE up"


def test_define_twice_is_an_error(net, tmp_path):
    with pytest.raises(ParserError, match="DEFINE m given twice") as err:
        render(net, tmp_path, "$$JAFF DEFINE m { }\n$$JAFF DEFINE m { }\n")
    assert err.value.nline == 2


def test_group_a_command_does_not_take_is_an_error_at_its_line(net, tmp_path):
    text = "$$JAFF SUB nspec {\n    SORT { ENABLED -> True }\n}\nx\n$$JAFF END\n"
    with pytest.raises(ParserError, match="SUB takes no option group SORT") as err:
        render(net, tmp_path, text)
    assert err.value.nline == 2


def test_group_from_use_is_checked_too(net, tmp_path):
    text = (
        "$$JAFF DEFINE m { SORT { ENABLED -> True } }\n"
        "$$JAFF SUB nspec { USE m }\nx\n$$JAFF END\n"
    )
    with pytest.raises(ParserError, match="takes no option group SORT"):
        render(net, tmp_path, text)


def test_unknown_token_is_an_error_at_the_header(net, tmp_path):
    with pytest.raises(ParserError, match="SUB has no token nspecs") as err:
        render(net, tmp_path, "a\n$$JAFF SUB nspecs\nx\n$$JAFF END\n")
    assert err.value.nline == 2


def test_every_block_is_validated_before_rendering(net, tmp_path, monkeypatch):
    from jaff.codegen.templating.commands import Sub

    calls = []
    monkeypatch.setattr(Sub, "render", lambda *a: calls.append(1) or "")
    text = "$$JAFF SUB nspec\nx\n$$JAFF END\n$$JAFF SUB bogus\ny\n$$JAFF END\n"
    with pytest.raises(ParserError, match="bogus"):
        render(net, tmp_path, text)
    assert calls == []


def test_get_and_has_render_through_the_facade(net, tmp_path):
    text = (
        "$$JAFF GET specie_idx IN H+\ni = $specie_idx$\n$$JAFF END\n"
        "$$JAFF HAS specie CO\nh = $specie$\n$$JAFF END\n"
    )
    assert render(net, tmp_path, text) == "i = 2\nh = 0\n"


SUB_BLOCK = "$$JAFF SUB nspec\nn = $nspec$\n$$JAFF END\n"


def test_unknown_extension_needs_a_lang_directive(net, tmp_path):
    with pytest.raises(ParserError, match="add \\$\\$JAFF LANG") as err:
        render(net, tmp_path, "x\n" + SUB_BLOCK, name="t.xyz")

    assert err.value.nline == 1
    text = "$$JAFF LANG python\n" + SUB_BLOCK
    assert render(net, tmp_path, text, name="t.xyz") == "n = 5\n"


def test_file_without_extension_needs_a_lang_directive(net, tmp_path):
    with pytest.raises(ParserError, match="No language for template"):
        render(net, tmp_path, SUB_BLOCK, name="template")


def test_lang_overrides_the_extension(net, tmp_path):
    text = "$$JAFF LANG fortran\n$$JAFF SUB nspec\nn = $nspec$\n$$JAFF END\n"
    assert render(net, tmp_path, text, name="t.py", echo=True) == (
        "! $$JAFF LANG fortran\n! $$JAFF SUB nspec\nn = 5\n! $$JAFF END\n"
    )


def test_lang_accepts_aliases(net, tmp_path):
    out = render(net, tmp_path, "$$JAFF LANG c++\nx\n", name="t.txt", echo=True)
    assert out == "// $$JAFF LANG c++\nx\n"


def test_unknown_lang_is_an_error_at_its_line(net, tmp_path):
    with pytest.raises(ParserError, match="not a supported language") as err:
        render(net, tmp_path, "x\n$$JAFF LANG cobol\n", name="t.txt")
    assert err.value.nline == 2


def test_comment_sets_the_echo_token(net, tmp_path):
    text = '$$JAFF COMMENT "#"\n$$JAFF SUB nspec\nn = $nspec$\n$$JAFF END\n'
    assert render(net, tmp_path, text, name="Config.cpp", echo=True) == (
        '# $$JAFF COMMENT "#"\n# $$JAFF SUB nspec\nn = 5\n# $$JAFF END\n'
    )


def test_get_of_unknown_entity_is_a_parser_error_at_the_header(net, tmp_path):
    text = "x\n$$JAFF GET specie_idx IN Xyz\ni = $specie_idx$\n$$JAFF END\n"
    with pytest.raises(ParserError, match="Xyz is not one") as err:
        render(net, tmp_path, text)
    assert err.value.nline == 2


def test_replace_patterns_are_compiled_once_per_block(net, tmp_path, monkeypatch):
    import jaff.codegen.templating.commands.base as base
    import jaff.codegen.templating.render as render_mod

    calls = []
    original = render_mod.compile_replacements

    def counting(options):
        calls.append(options)
        return original(options)

    monkeypatch.setattr(render_mod, "compile_replacements", counting)
    monkeypatch.setattr(base, "compile_replacements", counting, raising=False)

    text = '$$JAFF SUB nspec { REPLACE { "5" -> "five" } }\nn = $nspec$\n$$JAFF END\n'
    assert render(net, tmp_path, text) == "n = five\n"
    assert len(calls) == 1


def test_validation_error_on_a_body_line_reports_that_line(net, tmp_path, monkeypatch):
    def validate(self, body, resolved, ctx):
        raise TemplateError("bad body line", body_index=1)

    monkeypatch.setattr(Sub, "validate", validate)
    text = "x\n$$JAFF SUB nspec {\n}\na\nb = 1\n$$JAFF END\n"
    with pytest.raises(ParserError, match="bad body line") as err:
        render(net, tmp_path, text)
    assert err.value.nline == 5
    assert "b = 1" in str(err.value)


def test_validate_receives_the_body(net, tmp_path, monkeypatch):
    seen = []
    monkeypatch.setattr(Sub, "validate", lambda self, body, r, c: seen.append(body))
    render(net, tmp_path, "$$JAFF SUB nspec\na\nb\n$$JAFF END\n")
    assert seen == [["a\n", "b\n"]]


def test_render_error_becomes_a_parser_error(net, tmp_path, monkeypatch):
    def render_fails(self, body, resolved, ctx):
        raise TemplateError("cannot expand", body_index=0)

    monkeypatch.setattr(Sub, "render", render_fails)
    with pytest.raises(ParserError, match="cannot expand") as err:
        render(net, tmp_path, "$$JAFF SUB nspec\nn = $nspec$\n$$JAFF END\n")
    assert err.value.nline == 2


def test_failure_inside_codegen_is_a_parser_error_at_the_directive(net, tmp_path):
    text = "x\n$$JAFF REPEAT idx, rhs IN rhses { RADIATION { ENABLED -> True } }\n"
    text += "f[$idx$] = $rhs$\n$$JAFF END\n"
    with pytest.raises(ParserError, match="REPEAT failed: ValueError") as err:
        render(net, tmp_path, text)

    assert err.value.nline == 2
    assert isinstance(err.value.__cause__, ValueError)


def test_arithmetic_error_is_a_parser_error(net, tmp_path):
    text = "$$JAFF SUB nspec\nn = $nspec/0$\n$$JAFF END\n"
    with pytest.raises(ParserError, match="SUB failed: ZeroDivisionError"):
        render(net, tmp_path, text)


def test_lang_is_case_insensitive(net, tmp_path):
    out = render(net, tmp_path, "$$JAFF LANG Fortran\nx\n", name="t.txt", echo=True)
    assert out == "! $$JAFF LANG Fortran\nx\n"


def test_fortran_wrapping_leaves_echoed_comment_lines_alone(net, tmp_path):
    pattern = "x" * 150
    header = f'$$JAFF SUB nspec {{ REPLACE {{ "{pattern}" -> "y" }} }}'
    text = f'$$JAFF LANG fortran\n$$JAFF COMMENT "#"\n{header}\nn = $nspec$\n$$JAFF END\n'
    out = render(net, tmp_path, text, name="Config", echo=True)
    assert f"# {header}\n" in out


@pytest.mark.parametrize("name", ["Makefile", "CMakeLists.txt", "jaffgen.toml"])
def test_file_without_directives_or_language_is_copied(net, tmp_path, name):
    text = "all:\n\t$(CC) -o $nspec$ main.c\n"
    assert render(net, tmp_path, text, name=name, echo=True) == text
