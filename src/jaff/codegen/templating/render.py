# ABOUTME: Renders split template blocks: validates every block first, then renders
# ABOUTME: Resolves USE layering, registers DEFINEs and CSE groups, echoes directives

from typing import Dict, List, Set

from ...errors import ParserError
from .blocks import Block, CommandBlock, TextBlock
from .commands import COMMANDS, Define
from .commands.repeat.groups import problem as cse_group_problem
from .context import RenderContext
from .placeholders import split_line
from .errors import TemplateError
from .resolve import Resolved, check_groups, compile_replacements, layer


def _all_groups() -> Dict[str, Set[str]]:
    """Every option group any command accepts, with the union of its keys."""
    groups: Dict[str, Set[str]] = {}
    for cls in COMMANDS.values():
        for name, keys in cls.groups.items():
            groups.setdefault(name, set()).update(keys)

    return groups


def _error(block: CommandBlock, error: TemplateError, ctx: RenderContext) -> ParserError:
    """ParserError at the body line, the line or the header *error* points at."""
    if error.body_index is not None:
        nline = block.line + len(block.header_lines) + error.body_index
        text = block.body_lines[error.body_index]
    else:
        nline = error.nline or block.line
        index = nline - block.line
        lines = block.header_lines
        text = lines[index] if 0 <= index < len(lines) else lines[0]

    return ParserError(str(error), text.rstrip("\r\n"), nline, ctx.file)


def _prepare(block: CommandBlock, ctx: RenderContext) -> Resolved:
    """Validate one block and return its resolved options; register DEFINEs."""
    command = block.command
    try:
        options = layer(block.options, ctx.defines)
        accepted = _all_groups() if isinstance(command, Define) else command.groups
        check_groups(command.keyword, options, accepted)
        indent, _ = split_line(block.header_lines[0])
        resolved = Resolved(options, compile_replacements(options), indent)
        command.validate(block.body_lines, resolved, ctx)
    except TemplateError as e:
        raise _error(block, e, ctx) from e.__cause__

    if isinstance(command, Define):
        if command.name in ctx.defines:
            error = TemplateError(f"DEFINE {command.name} given twice")
            raise _error(block, error, ctx)

        ctx.defines[command.name] = options

    return resolved


def echo(lines: List[str], comment: str) -> str:
    """Directive lines as comments, each keeping its own indentation."""
    out = []
    for line in lines:
        text = line.rstrip("\r\n")
        content = text.lstrip()
        indent = text[: len(text) - len(content)]
        out.append(f"{indent}{comment} {content}".rstrip() + "\n")

    return "".join(out)


def _check_cse_groups(blocks: List[Block], ctx: RenderContext) -> None:
    """Every CSE group has its ``REPEAT IN cses`` block and at least one member."""
    found = cse_group_problem(ctx)
    if found is None:
        return

    command, message = found
    block = next(
        b for b in blocks if isinstance(b, CommandBlock) and id(b.command) == command
    )
    raise _error(block, TemplateError(message), ctx)


def render_blocks(blocks: List[Block], ctx: RenderContext) -> str:
    """Validate every block, then render them all in order.

    Args:
        blocks: Blocks from ``split_blocks``.
        ctx: Render context; its ``defines`` fill up as DEFINEs are met.

    Returns:
        The generated text.

    Raises:
        ParserError: On the first invalid block, before anything is rendered; on
            a body line that cannot be expanded; or at the directive whose
            rendering failed in the network or codegen (the original exception
            is the cause).
    """
    prepared = [
        _prepare(block, ctx) if isinstance(block, CommandBlock) else None
        for block in blocks
    ]
    _check_cse_groups(blocks, ctx)
    comment = ctx.comment or ctx.cg.lang.comment
    out: List[str] = []
    for block, resolved in zip(blocks, prepared):
        if isinstance(block, TextBlock):
            out.extend(block.lines)
            continue

        if ctx.echo_directives:
            out.append(echo(block.header_lines, comment))
        try:
            out.append(block.command.render(block.body_lines, resolved, ctx))
        except TemplateError as e:
            raise _error(block, e, ctx) from e.__cause__
        except Exception as e:
            # Network or codegen failures (e.g. no radiation bands): report them at
            # the directive, keeping the original exception as the cause.
            keyword = block.command.keyword
            error = TemplateError(f"{keyword} failed: {type(e).__name__}: {e}")
            raise _error(block, error, ctx) from e

        if ctx.echo_directives:
            out.append(echo(block.end_lines, comment))

    return "".join(out)
