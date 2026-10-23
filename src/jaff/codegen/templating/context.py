# ABOUTME: Render context shared by every command while one template is rendered
# ABOUTME: Network, code generator, template path, echo settings and DEFINE registry

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Optional

from .options import Options

if TYPE_CHECKING:
    from ...core.network.network import Network
    from ..codegen import Codegen


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
    """

    net: Network
    cg: Codegen
    file: Path
    echo_directives: bool = True
    comment: Optional[str] = None
    defines: Dict[str, Options] = field(default_factory=dict)
