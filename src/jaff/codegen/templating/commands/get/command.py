# ABOUTME: GET command: substitute properties of a species, reaction, element or aux func
# ABOUTME: Each body line keeps its indentation; REPLACE runs on every rendered line

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, List, Tuple

from ...placeholders import substitute
from ...errors import TemplateError
from ..base import Command, check_names, render_each_line
from ..has.props import KINDS
from .props import PROPERTIES, PROPS

if TYPE_CHECKING:
    from ...context import RenderContext
    from ...resolve import Resolved

_KIND_NAMES = {
    "element": "an element",
    "specie": "a species",
    "reaction": "a reaction",
    "aux": "an auxiliary function",
}


def _has_entity(ctx: RenderContext, kind: str, entity: str) -> bool:
    """Whether *entity* exists as *kind* and resolves to a single item for GET."""
    if kind == "aux":
        return entity in ctx.net.spec.aux_funcs

    if entity not in KINDS[kind](ctx):
        return False

    # A reaction name shared by several mechanisms gives a list, which has no tmin etc.
    return kind != "reaction" or not isinstance(ctx.net.reactions[entity], list)


@dataclass(frozen=True)
class Get(Command):
    """``GET props IN entity``."""

    props: Tuple[str, ...]
    entity: str

    def validate(self, body: List[str], resolved: Resolved, ctx: RenderContext) -> None:
        """Every property must be known and the entity must exist as its kind."""
        check_names(self.keyword, "property", self.props, PROPS)
        for name in self.props:
            kind = PROPERTIES[name][0]
            if not _has_entity(ctx, kind, self.entity):
                raise TemplateError(
                    f"{self.keyword} {name} needs {_KIND_NAMES[kind]}; "
                    f"{self.entity} is not one in this network"
                )

    def render(self, body: List[str], resolved: Resolved, ctx: RenderContext) -> str:
        """Substitute the listed properties of the entity on every body line."""
        values = {
            name: (lambda get=PROPERTIES[name][1]: get(ctx, self.entity))
            for name in self.props
        }

        return render_each_line(
            body, lambda text: substitute(text, values), resolved.replacements
        )
