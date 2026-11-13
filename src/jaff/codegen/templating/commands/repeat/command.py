# ABOUTME: REPEAT command: expand body lines once per item of a network property
# ABOUTME: Vertical with idx (one line per item), else horizontal into a bracketed list

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, ClassVar, Dict, List, Set, Tuple

from ...errors import TemplateError
from ...placeholders import split_line
from ..base import Command, check_names, finish_line
from . import groups
from .props import GROUPS, PROPERTIES, RepeatProp, option_kwargs
from .expand import as_flattened, expand_horizontal, expand_vertical, find_idx

if TYPE_CHECKING:
    from ...context import RenderContext
    from ...resolve import Resolved


def _active(line: str, names: Tuple[str, ...]) -> bool:
    """Whether a body line holds a placeholder of this REPEAT.

    Placeholders are ``$name$`` for *names* (the item variable), ``$idx$`` with
    or without arithmetic, and ``$IDX$``. Other lines are copied, with REPLACE
    applied.
    """
    pattern = "|".join(re.escape(name) for name in names)
    placeholder = re.search(rf"\$(?:{pattern})\$", line)
    return bool(placeholder or find_idx(line) or "$IDX$" in line)


@dataclass(frozen=True)
class Repeat(Command):
    """``REPEAT vars IN prop``."""

    groups: ClassVar[Dict[str, Set[str]]] = GROUPS
    vars: Tuple[str, ...]
    prop: str

    def validate(self, body: List[str], resolved: Resolved, ctx: RenderContext) -> None:
        """Check the property, its variables, its options and the CSE line.

        A ``REPEAT IN cses`` block declares its CSE group, and an expression block
        with ``CSE { GROUP -> name }`` joins one (see :mod:`.groups`).
        """
        check_names(self.keyword, "property", (self.prop,), tuple(PROPERTIES))
        prop = PROPERTIES[self.prop]
        if prop.kind and "cse" in self.vars:
            raise TemplateError(
                f"REPEAT IN {self.prop} has no cse variable: declare the temporaries "
                "in their own block, REPEAT idx, cse IN cses { GROUP { NAME -> g } }, "
                "and join it with CSE { GROUP -> g }"
            )

        check_names(f"{self.keyword} IN {self.prop}", "variable", self.vars, prop.vars)
        kwargs = self._kwargs(body, resolved)
        if self.prop == "cses":
            if "name" not in kwargs:
                raise TemplateError("REPEAT IN cses needs its group: GROUP { NAME -> g }")
            groups.declare(
                ctx, self, kwargs["name"], kwargs["cse_var"], kwargs["cse_suffix"]
            )
        elif "cse_group" in kwargs:
            name = kwargs.pop("cse_group")
            assert prop.kind is not None  # only expression properties take CSE
            groups.join(ctx, self, name, prop.kind, kwargs)

    def render(self, body: List[str], resolved: Resolved, ctx: RenderContext) -> str:
        """Expand every line holding a placeholder; copy the others once.

        REPLACE runs on the content of every output line, expanded or copied.
        """
        prop = PROPERTIES[self.prop]
        items: Any = None
        out: List[str] = []
        for index, raw in enumerate(body):
            indent, line = split_line(raw)
            contents = [line]
            if _active(line, (prop.item,)):
                if items is None:
                    items = self._items(body, resolved, ctx)
                contents = self._expand(prop, items, line, index)

            out.extend(finish_line(indent, c, resolved.replacements) for c in contents)

        return "".join(out)

    def _items(self, body: List[str], resolved: Resolved, ctx: RenderContext) -> Any:
        """The items: from the CSE group, or from the property's source."""
        if self.prop == "cses":
            return groups.temporaries(ctx, self._kwargs(body, resolved)["name"])

        if id(self) in ctx.cse_members:
            return groups.expressions(ctx, self)

        # Recomputed rather than kept from validate: commands are immutable, and
        # the options are cheap next to codegen.
        return PROPERTIES[self.prop].source(ctx, self._kwargs(body, resolved))

    def _expand(self, prop: RepeatProp, items: Any, line: str, index: int) -> List[str]:
        """The contents of the output lines for one line holding a placeholder."""
        if "idx" not in self.vars:
            expanded = expand_horizontal(items, prop.item, line)
            return [line if expanded is None else expanded]

        return expand_vertical(as_flattened(items), line, f"${prop.item}$", index)

    def _kwargs(self, body: List[str], resolved: Resolved) -> Dict[str, Any]:
        """Keyword arguments from the options; for ``cses``, the temporaries' names."""
        kwargs = option_kwargs(self.prop, resolved.options)
        if self.prop == "cses":
            kwargs.update(self._cse_names(body))

        return kwargs

    def _cse_names(self, body: List[str]) -> Dict[str, Any]:
        """The temporaries' names, taken from the ``$cse$`` line of a ``cses`` block.

        The ``$cse$`` line must be the first line that expands, and names its
        temporaries with an identifier around a plain ``$idx$``, e.g. ``cse$idx$``.

        Raises:
            TemplateError: At the offending body line.
        """
        index = next(
            (i for i, raw in enumerate(body) if _active(raw.strip(), ("cse",))), None
        )
        if index is None:
            return {"cse_var": "cse", "cse_suffix": ""}

        line = body[index].strip()

        def error(message: str) -> TemplateError:
            return TemplateError(message, body_index=index)

        if "$cse$" not in line:
            raise error(
                "The $cse$ line must be the first line of the REPEAT block "
                "so CSE temporaries are declared before they are used"
            )
        spans = find_idx(line)
        if not spans:
            raise error("No valid idx variable detected")
        (begin, end), op, _ = spans[0]
        if op:
            raise error(
                "CSE temporaries cannot use arithmetic on $idx$; use $idx$ in the $cse$ name"
            )

        before = re.search(r"\w*$", line[:begin])
        after = re.match(r"\w*", line[end:])
        prefix = before.group() if before else ""
        suffix = after.group() if after else ""
        if not re.fullmatch(r"[A-Za-z_]\w*", prefix):
            raise error(
                "CSE temporary must be an identifier like cse$idx$ or tmp$idx$_value; "
                f"got '{prefix}$idx${suffix}'"
            )

        return {"cse_var": prefix, "cse_suffix": suffix}
