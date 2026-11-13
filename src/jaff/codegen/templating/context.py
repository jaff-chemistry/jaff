# ABOUTME: Render context shared by every command while one template is rendered
# ABOUTME: Network, codegen, template path, echo settings, DEFINEs, CALL functions, CSE groups

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Optional, Tuple

from .options import Options

if TYPE_CHECKING:
    from ...core.network.network import Network
    from ..codegen import Codegen
    from .commands.repeat.groups import CseGroupSpec
    from .functions import TemplateFunctions


@dataclass
class RenderContext:
    """State shared by the commands of one template.

    Attributes:
        net: Network the template is rendered for.
        cg: Code generator for the template's language.
        file: Template path, for messages and the ``filename``/``filepath``
            tokens.
        echo_directives: Whether directives are echoed as comments in the output.
        comment: Comment token for the echo, from ``$$JAFF COMMENT``; ``None``
            means the language's own comment.
        defines: Option sets registered by ``DEFINE``, by name.
        functions: Functions for ``$$JAFF CALL``, or ``None`` without a
            ``jaffgen_template_functions.py``.
        cse_groups: CSE groups by name, filled while the blocks are validated.
        cse_members: ``id`` of each block joining a CSE group, to its group's name
            and its position among the members.
    """

    net: Network
    cg: Codegen
    file: Path
    echo_directives: bool = True
    comment: Optional[str] = None
    defines: Dict[str, Options] = field(default_factory=dict)
    functions: Optional[TemplateFunctions] = None
    cse_groups: Dict[str, CseGroupSpec] = field(default_factory=dict)
    cse_members: Dict[int, Tuple[str, int]] = field(default_factory=dict)
