# ABOUTME: Parses one $$JAFF directive header with Lark into a command and its options
# ABOUTME: Builds commands from the COMMANDS registry; no network or codegen knowledge

from dataclasses import dataclass
from functools import cache
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Type

from lark import Lark, Token, Transformer, UnexpectedInput, UnexpectedToken
from lark.tree import Meta

from ...errors import ParserError
from .commands import COMMANDS, Command
from .errors import TemplateError
from .options import Options

GRAMMAR = Path(__file__).with_name("grammar.lark")
MARKER = "$$JAFF"


@dataclass(frozen=True)
class Directive:
    """A parsed directive header: its command and its option block, if any."""

    command: Command
    options: Optional[Options]


def _text(token: Token) -> str:
    """Token text; double-quoted strings lose their quotes, escapes stay raw."""
    if token.type == "STRING":
        return token[1:-1]

    return str(token)


def _unescape(token: Token) -> str:
    """Token text with quotes removed and ``\\"`` turned into ``"``."""
    return _text(token).replace('\\"', '"')


# One raw option: (kind, value); REPLACE carries [(pattern, replacement, line)],
# USE (name, line), groups (name, line, [(key, value, line)]); lines are relative
# to the header text
_Part = Tuple[str, Any]


class _ToDirective(Transformer):
    """Lark transformer turning the parse tree into a command and raw option parts."""

    def __default__(self, data: str, children: List[Any], meta: Meta) -> Any:
        """Build a command for command aliases; leave other rules as trees."""
        if data in COMMANDS:
            items = [_unescape(c) if isinstance(c, Token) else c for c in children]
            return COMMANDS[data].from_items(items)

        return super().__default__(data, children, meta)

    def start(self, items: List[Any]) -> Tuple[Command, Optional[List[_Part]]]:
        """Pair the command with its raw option parts (``None`` without a block)."""
        return items[0], (items[1] if len(items) > 1 else None)

    def names(self, items: List[Any]) -> Tuple[str, ...]:
        """Comma-separated names as a tuple of strings."""
        return tuple(str(t) for t in items)

    def replacement(self, items: List[Any]) -> Tuple[str, str, int]:
        """A REPLACE pair: raw regex pattern, unescaped replacement, relative line."""
        return (_text(items[0]), _unescape(items[1]), items[0].line)

    def setting(self, items: List[Any]) -> Tuple[str, str, int]:
        """A group setting: key, unescaped value and the key's relative line."""
        return (str(items[0]), _unescape(items[1]), items[0].line)

    def replace(self, items: List[Any]) -> _Part:
        """A REPLACE block's pairs."""
        return ("replace", items)

    def use(self, items: List[Any]) -> _Part:
        """A USE of a named option set."""
        return ("use", (str(items[0]), items[0].line))

    def group(self, items: List[Any]) -> _Part:
        """A named group: name, its relative line and its settings."""
        return ("group", (str(items[0]), items[0].line, items[1:]))

    def options(self, items: List[Any]) -> List[_Part]:
        """The block's option parts in written order."""
        return items


@cache
def _parser() -> Lark:
    """The header parser, built once per process."""
    return Lark(
        GRAMMAR.read_text(encoding="utf-8"), parser="lalr", propagate_positions=True
    )


def _add_group(
    options: Options, first: int, name: str, rel_line: int, settings: List[Any]
) -> None:
    """Add one group to *options*; a repeated group or key is an error.

    *first* is the template line of the header's first line.
    """
    if name in options.groups:
        raise TemplateError(f"option group {name} given twice", first + rel_line - 1)
    seen: Dict[str, int] = {}
    for key, _, key_line in settings:
        if key in seen:
            raise TemplateError(f"{name} sets {key} more than once", first + key_line - 1)

        seen[key] = key_line
    options.groups[name] = {key: value for key, value, _ in settings}
    options.group_lines[name] = first + rel_line - 1


def _build_options(parts: List[_Part], first: int) -> Options:
    """Collect option parts; a repeated group or key in one block is an error.

    *first* is the template line of the header's first line.
    """
    options = Options()
    for kind, value in parts:
        if kind == "replace":
            for pattern, replacement, rel_line in value:
                options.replace.append((pattern, replacement))
                options.replace_lines.append(first + rel_line - 1)
        elif kind == "use":
            name, rel_line = value
            if name in options.use_lines:
                raise TemplateError(f"USE {name} given twice", first + rel_line - 1)
            options.uses.append(name)
            options.use_lines[name] = first + rel_line - 1
        else:
            _add_group(options, first, *value)

    return options


def _terminal_label(name: str) -> str:
    """Terminal name; anonymous terminals (e.g. ``->``) show their literal text."""
    if not name.startswith("__ANON"):
        return name
    for terminal in _parser().terminals:
        if terminal.name == name:
            return f"'{terminal.pattern.value}'"

    return name


def _describe(e: UnexpectedInput) -> str:
    """Readable message for a Lark failure: what was found and what was expected."""
    if isinstance(e, UnexpectedToken):
        found = "end of header" if e.token.type == "$END" else f"'{e.token}'"
        names = e.expected
    else:
        found = f"'{getattr(e, 'char', '')}'"
        names = getattr(e, "allowed", None)

    message = f"Invalid directive: unexpected {found} at column {max(e.column, 1)}"
    if names:
        shown = sorted(_terminal_label(name) for name in names)[:8]
        message += f"; expected one of: {', '.join(shown)}"

    return message


class DirectiveSyntaxError(ParserError):
    """The header text does not match the directive grammar."""


class IncompleteHeaderError(DirectiveSyntaxError):
    """The header ended inside an open option block; more lines may complete it."""


class SplitOptionError(DirectiveSyntaxError):
    """A line ended in the middle of an option (a pair, ``NAME {`` or ``USE name``)."""


def _inside_block(e: UnexpectedToken) -> bool:
    """Whether the parser had an option block ``{`` open when it failed."""
    state = getattr(getattr(e, "interactive_parser", None), "parser_state", None)
    stack = getattr(state, "value_stack", [])
    return any(getattr(item, "type", None) == "LBRACE" for item in stack)


def _end_of_header_error(e: UnexpectedToken) -> Tuple[str, Type[ParserError]]:
    """Message and error class for a header that ran out of text."""
    if "RBRACE" in e.expected:
        return _describe(e), IncompleteHeaderError

    if _inside_block(e):
        message = (
            f"Invalid directive: line ends inside an option at column "
            f"{max(e.column, 1)}; keep each pair, 'NAME {{' and 'USE name' on one line"
        )
        return message, SplitOptionError

    return _describe(e), DirectiveSyntaxError


def _blank_marker(text: str) -> Optional[str]:
    """Replace the leading ``$$JAFF`` with spaces so Lark columns match the source.

    *text* must start with optional whitespace and ``$$JAFF`` (the engine only
    passes lines for which ``is_directive`` holds). Returns ``None`` when the
    marker is not followed by whitespace.
    """
    body = text.lstrip()
    indent = len(text) - len(body)
    rest = body[len(MARKER) :]
    if rest and not rest[0].isspace():
        return None

    return text[:indent] + " " * len(MARKER) + rest


def parse_header(text: str, line: int, file: Path) -> Directive:
    """Parse a directive header, including its ``$$JAFF`` marker.

    Args:
        text: Raw header source, starting with the (indented) ``$$JAFF`` line;
            may span lines.
        line: 1-based template line on which the header starts.
        file: Template path, for error messages.

    Returns:
        The parsed directive.

    Raises:
        IncompleteHeaderError: If the text ends inside an open option block.
        SplitOptionError: If a line ends in the middle of an option.
        DirectiveSyntaxError: If the text does not match the grammar.
        ParserError: On a repeated option group, key or ``USE`` name.
    """
    lines = [part.rstrip("\r") for part in text.split("\n")]

    def error(
        message: str, rel_line: int, cls: Type[ParserError] = ParserError
    ) -> ParserError:
        rel = rel_line if 0 < rel_line <= len(lines) else 1
        return cls(message, lines[rel - 1], line + rel - 1, file)

    source = _blank_marker(text)
    if source is None:
        raise error(f"Expected a space after {MARKER}", 1, DirectiveSyntaxError)
    if not source.strip():
        raise error(
            f"Invalid directive: missing command after {MARKER}", 1, DirectiveSyntaxError
        )

    try:
        tree = _parser().parse(source)
    except UnexpectedToken as e:
        if e.token.type == "$END":
            message, cls = _end_of_header_error(e)
            raise error(message, e.line, cls) from None

        raise error(_describe(e), e.line, DirectiveSyntaxError) from None
    except UnexpectedInput as e:
        raise error(_describe(e), e.line, DirectiveSyntaxError) from None

    command, parts = _ToDirective().transform(tree)
    if parts is None:
        return Directive(command, None)
    try:
        return Directive(command, _build_options(parts, line))
    except TemplateError as e:
        # Option errors always carry their absolute line; fall back to the header.
        nline = e.nline if e.nline is not None else line
        raise error(str(e), nline - line + 1) from None
