# ABOUTME: CALL command: insert the text a @jaffgen_function returns for the network
# ABOUTME: No body; every returned line gets the directive's indentation and REPLACE

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Callable, ClassVar, List

from ...errors import TemplateError
from ...functions import FUNCTIONS_FILE
from ..base import Command, finish_line

if TYPE_CHECKING:
    from ...context import RenderContext
    from ...resolve import Resolved


@dataclass(frozen=True)
class Call(Command):
    """``CALL name``."""

    has_body: ClassVar[bool] = False
    name: str

    def _function(self, ctx: RenderContext) -> Callable[..., str]:
        """The function called *name*; loads the functions file on first use."""
        if ctx.functions is None:
            raise TemplateError(f"CALL {self.name} needs a {FUNCTIONS_FILE}")

        try:
            function = ctx.functions.get(self.name)
            if function is None:
                known = ", ".join(ctx.functions.names()) or "none"
                raise TemplateError(f"CALL has no function {self.name} (known: {known})")
        except RuntimeError as e:
            raise TemplateError(str(e)) from e.__cause__ or e

        return function

    def validate(self, body: List[str], resolved: Resolved, ctx: RenderContext) -> None:
        """The functions file must load and register *name*."""
        self._function(ctx)

    def render(self, body: List[str], resolved: Resolved, ctx: RenderContext) -> str:
        """The function's text, indented like the directive, REPLACE applied."""
        text = self._function(ctx)(ctx.net)
        if not isinstance(text, str):
            raise TemplateError(
                f"CALL {self.name} returned {type(text).__name__}, expected str"
            )

        lines = text.splitlines()

        return "".join(
            finish_line(resolved.indent, line, resolved.replacements) for line in lines
        )
