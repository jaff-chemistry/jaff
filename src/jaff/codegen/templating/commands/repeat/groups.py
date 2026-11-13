# ABOUTME: CSE groups of a template: REPEAT IN cses blocks declare them, expression
# ABOUTME: blocks join them with CSE { GROUP }; one Codegen.cse_group per group

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

from .....types import IndexedList
from ...errors import TemplateError

if TYPE_CHECKING:
    from ...context import RenderContext


@dataclass
class CseGroupSpec:
    """One CSE group of a template, filled while its blocks are validated.

    Attributes:
        name: The group's name.
        declared_by: ``id`` of the ``REPEAT IN cses`` command declaring it.
        cse_var, cse_suffix: The temporaries' names, from that block's ``$cse$``
            line.
        members: ``(kind, options)`` of each joining block, in template order.
        joined_by: ``id`` of each joining command, in the same order.
    """

    name: str
    declared_by: Optional[int] = None
    cse_var: str = "cse"
    cse_suffix: str = ""
    members: List[Tuple[str, Dict[str, Any]]] = field(default_factory=list)
    joined_by: List[int] = field(default_factory=list)


def declare(
    ctx: RenderContext, command: object, name: str, cse_var: str, cse_suffix: str
) -> None:
    """Record *command* as the ``REPEAT IN cses`` block of group *name*.

    Raises:
        TemplateError: If the group already has one.
    """
    spec = ctx.cse_groups.setdefault(name, CseGroupSpec(name))
    if spec.declared_by is not None:
        raise TemplateError(f"CSE group {name} is declared twice")

    spec.declared_by = id(command)
    spec.cse_var, spec.cse_suffix = cse_var, cse_suffix


def join(
    ctx: RenderContext, command: object, name: str, kind: str, options: Dict[str, Any]
) -> None:
    """Record *command*, printing *kind* with *options*, as a member of group *name*."""
    spec = ctx.cse_groups.setdefault(name, CseGroupSpec(name))
    ctx.cse_members[id(command)] = (name, len(spec.members))
    spec.members.append((kind, dict(options)))
    spec.joined_by.append(id(command))


def problem(ctx: RenderContext) -> Optional[Tuple[int, str]]:
    """The first incomplete group: ``(id of the command to blame, message)``."""
    for spec in ctx.cse_groups.values():
        if spec.declared_by is None:
            return spec.joined_by[0], (
                f"CSE group {spec.name} has no REPEAT idx, cse IN cses block "
                f"(add one with GROUP {{ NAME -> {spec.name} }})"
            )

        if not spec.members:
            return spec.declared_by, (
                f"no block joins CSE group {spec.name} "
                f"(add CSE {{ GROUP -> {spec.name} }} to the blocks it serves)"
            )

    return None


def temporaries(ctx: RenderContext, name: str) -> IndexedList:
    """The shared temporaries of group *name*."""
    return _group(ctx, ctx.cse_groups[name]).temporaries


def expressions(ctx: RenderContext, command: object) -> IndexedList:
    """The expressions of the joining *command*, reading its group's temporaries."""
    name, index = ctx.cse_members[id(command)]

    return _group(ctx, ctx.cse_groups[name]).expressions[index]


def _group(ctx: RenderContext, spec: CseGroupSpec) -> Any:
    """``Codegen.cse_group`` of *spec* (cached there)."""
    return ctx.cg.cse_group(spec.members, spec.cse_var, spec.cse_suffix)
