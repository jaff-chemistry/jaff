# ABOUTME: Base class for $$JAFF commands, built from a parsed directive header
# ABOUTME: Subclasses add the header's fields plus validate and render for their body

from __future__ import annotations

from dataclasses import dataclass
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    ClassVar,
    Dict,
    List,
    Pattern,
    Self,
    Sequence,
    Set,
    Tuple,
)

from ..errors import TemplateError
from ..placeholders import split_line
from ..resolve import apply_replacements

if TYPE_CHECKING:
    from ..context import RenderContext
    from ..resolve import Resolved


@dataclass(frozen=True)
class Command:
    """A parsed ``$$JAFF`` command; subclasses add the header's fields.

    Attributes:
        has_body: Whether the command owns the lines up to ``$$JAFF END``.
        groups: Option groups the command accepts, with their keys.
    """

    has_body: ClassVar[bool] = True
    groups: ClassVar[Dict[str, Set[str]]] = {}

    @classmethod
    def from_items(cls, items: Sequence[Any]) -> Self:
        """Build the command from its header parts, in grammar order.

        Args:
            items: Header parts as plain strings or tuples of names.

        Returns:
            The command instance.
        """
        return cls(*items)

    @property
    def keyword(self) -> str:
        """The command's keyword as written in templates, e.g. ``SUB``."""
        return type(self).__name__.upper()

    def validate(self, body: List[str], resolved: Resolved, ctx: RenderContext) -> None:
        """Check the block against the network before anything is rendered.

        Args:
            body: Raw body lines, with their line endings.
            resolved: Layered options and compiled REPLACE of the block.
            ctx: Render context.

        Raises:
            TemplateError: If the block names something the command cannot use;
                ``body_index`` set when a body line is at fault.
        """

    def render(self, body: List[str], resolved: Resolved, ctx: RenderContext) -> str:
        """Render the body lines into output text.

        Args:
            body: Raw body lines, with their line endings.
            resolved: Layered options and compiled REPLACE of the block.
            ctx: Render context.

        Returns:
            The generated text, one ``\\n``-terminated line per output line.

        Raises:
            TemplateError: If a body line cannot be expanded.
        """
        raise NotImplementedError(f"{self.keyword} rendering is not implemented")


@dataclass(frozen=True)
class FileSetting(Command):
    """A directive that configures the whole template, such as ``LANG``.

    It has no body and renders nothing; the engine requires it before every
    other directive, at most once per keyword.
    """

    has_body: ClassVar[bool] = False

    def render(self, body: List[str], resolved: Resolved, ctx: RenderContext) -> str:
        """A file setting produces no output of its own."""
        return ""


def finish_line(
    indent: str, content: str, replacements: List[Tuple[Pattern[str], str]]
) -> str:
    """An output line: *indent*, then *content* after the REPLACE pairs, then ``\\n``.

    REPLACE never sees the indentation, so every command treats lines alike. An
    empty line gets no indentation.
    """
    text = apply_replacements(content, replacements)
    return (indent + text if text else "") + "\n"


def render_each_line(
    body: List[str],
    transform: Callable[[str], str],
    replacements: List[Tuple[Pattern[str], str]],
) -> str:
    """Render every body line: indentation kept, content through *transform*.

    Used by SUB, GET and HAS; see :func:`finish_line` for REPLACE.
    """
    out = []
    for line in body:
        indent, content = split_line(line)
        out.append(finish_line(indent, transform(content), replacements))

    return "".join(out)


def check_names(
    keyword: str, kind: str, names: Sequence[str], known: Sequence[str]
) -> None:
    """Every name in *names* must be one of *known*.

    Raises:
        TemplateError: Listing the unknown names and the known ones.
    """
    unknown = [name for name in names if name not in known]
    if unknown:
        raise TemplateError(
            f"{keyword} has no {kind} {', '.join(unknown)} "
            f"(known: {', '.join(sorted(known))})"
        )
