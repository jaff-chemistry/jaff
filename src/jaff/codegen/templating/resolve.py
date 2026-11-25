# ABOUTME: Resolves a block's options: USE layering, per-command group checks, REPLACE
# ABOUTME: Pure functions over Options; errors carry the template line they belong to

import re
from dataclasses import dataclass
from typing import List, Mapping, Optional, Pattern, Set, Tuple

from .errors import TemplateError
from .options import Options


@dataclass(frozen=True)
class Resolved:
    """A block's options after USE layering, with REPLACE compiled once.

    Attributes:
        options: Layered options of the block.
        replacements: Compiled REPLACE pairs of *options*, in order.
        indent: Indentation of the directive's first header line, for commands
            without a body whose output stands in for the directive (``CALL``).
    """

    options: Options
    replacements: List[Tuple[Pattern[str], str]]
    indent: str = ""


def _merge(target: Options, source: Options, line: Optional[int]) -> None:
    """Lay *source* over *target*; group lines become *line* when given."""
    target.replace.extend(source.replace)
    target.replace_lines.extend(
        [line] * len(source.replace) if line is not None else source.replace_lines
    )

    for name, settings in source.groups.items():
        target.groups.setdefault(name, {}).update(settings)
        target.group_lines[name] = line if line is not None else source.group_lines[name]


def layer(local: Optional[Options], defines: Mapping[str, Options]) -> Options:
    """Options in effect for a block: its USE sets in order, then its own options.

    REPLACE lists concatenate; groups merge key by key and later sets win. Options
    taken from a USE report the USE's line in later errors.

    Args:
        local: The block's own options, or ``None`` without an option block.
        defines: Option sets registered so far by ``DEFINE``.

    Returns:
        A new ``Options`` without ``uses``.

    Raises:
        TemplateError: If a USE names a set that has not been defined yet.
    """
    result = Options()
    if local is None:
        return result

    for name in local.uses:
        line = local.use_lines[name]
        if name not in defines:
            raise TemplateError(
                f"USE of undefined option set {name} (DEFINE it above this block)", line
            )

        _merge(result, defines[name], line)
    _merge(result, local, None)

    return result


def check_groups(
    command: str, options: Options, accepted: Mapping[str, Set[str]]
) -> None:
    """Every group must be one *command* takes, with keys it knows.

    Args:
        command: Command name for messages, e.g. ``"SUB"``.
        options: Layered options of the block.
        accepted: Group name to the set of keys the command accepts.

    Raises:
        TemplateError: On an unknown group or key.
    """
    for name, settings in options.groups.items():
        line = options.group_lines[name]
        if name not in accepted:
            allowed = ", ".join(sorted(accepted)) or "none"
            raise TemplateError(
                f"{command} takes no option group {name} (accepted: {allowed})", line
            )
        unknown = sorted(set(settings) - accepted[name])
        if unknown:
            raise TemplateError(
                f"{name} has no setting {', '.join(unknown)} "
                f"(accepted: {', '.join(sorted(accepted[name]))})",
                line,
            )


def compile_replacements(options: Options) -> List[Tuple[Pattern[str], str]]:
    """Compile the REPLACE patterns in order and check their replacements.

    Raises:
        TemplateError: If a pattern is not a valid regular expression, or a
            replacement refers to a group the pattern does not have.
    """
    compiled: List[Tuple[Pattern[str], str]] = []
    for (pattern, replacement), line in zip(options.replace, options.replace_lines):
        try:
            compiled_pattern = re.compile(pattern)
        except re.error as e:
            raise TemplateError(
                f"Invalid REPLACE pattern {pattern!r}: {e}", line
            ) from None

        try:
            compiled_pattern.sub(replacement, "")  # validates group references
        except re.error as e:
            raise TemplateError(
                f"Invalid REPLACE replacement {replacement!r}: {e}", line
            ) from None

        compiled.append((compiled_pattern, replacement))

    return compiled


def apply_replacements(text: str, compiled: List[Tuple[Pattern[str], str]]) -> str:
    """Apply compiled REPLACE pairs to *text*, one after another."""
    for pattern, replacement in compiled:
        text = pattern.sub(replacement, text)

    return text
