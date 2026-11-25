# ABOUTME: Block types produced by the template engine's single pass over a template
# ABOUTME: TextBlock holds plain lines; CommandBlock a command, its options and body

from dataclasses import dataclass
from typing import List, Optional, Union

from .commands import Command
from .options import Options


@dataclass
class TextBlock:
    """Template lines copied to the output unchanged."""

    lines: List[str]
    line: int


@dataclass
class CommandBlock:
    """One ``$$JAFF`` directive: command, options, raw header lines and body lines.

    ``body_lines`` and ``end_lines`` are empty for commands without a body
    (``DEFINE`` and the file settings ``LANG`` and ``COMMENT``); ``end_lines``
    holds the raw ``$$JAFF END`` line otherwise.
    """

    command: Command
    options: Optional[Options]
    header_lines: List[str]
    body_lines: List[str]
    end_lines: List[str]
    line: int


Block = Union[TextBlock, CommandBlock]
