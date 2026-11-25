# ABOUTME: REDUCE command: expand $( ... )$ into a sum over network properties
# ABOUTME: The i-th summand takes the i-th value of every listed property

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

from ...errors import TemplateError
from ...placeholders import split_line
from ..base import Command, check_names, finish_line
from .props import PROPERTIES, PROPS

if TYPE_CHECKING:
    from ...context import RenderContext
    from ...resolve import Resolved

# $( ... )$: group 1 is the whole reduction, group 2 the summand template
_REDUCTION = re.compile(r"(\$\((.*?)\)\$)")


@dataclass(frozen=True)
class Reduce(Command):
    """``REDUCE vars IN props``."""

    vars: Tuple[str, ...]
    props: Tuple[str, ...]

    def validate(self, body: List[str], resolved: Resolved, ctx: RenderContext) -> None:
        """Properties must be known and every variable must belong to one of them."""
        check_names(self.keyword, "property", self.props, PROPS)
        known = [PROPERTIES[prop][0] for prop in self.props]
        check_names(self.keyword, "variable", self.vars, known)

    def render(self, body: List[str], resolved: Resolved, ctx: RenderContext) -> str:
        """Expand every ``$( ... )$`` that uses a listed variable, on every line.

        Other text is copied. REPLACE runs on the content of every line.
        """
        columns: Optional[Dict[str, List[Any]]] = None

        def expand(match: re.Match[str], index: int) -> str:
            nonlocal columns
            summand = match.group(2)
            if not any(f"${var}$" in summand for var in self.vars):
                return match.group(1)
            if columns is None:
                columns = self._columns(ctx)
            return self._sum(summand, columns, index)

        out = []
        for index, raw in enumerate(body):
            indent, line = split_line(raw)
            line = _REDUCTION.sub(lambda m: expand(m, index), line)
            out.append(finish_line(indent, line, resolved.replacements))

        return "".join(out)

    def _columns(self, ctx: RenderContext) -> Dict[str, List[Any]]:
        """Placeholder variable to its list of values, in the order of the props."""
        columns = {}
        for prop in self.props:
            var, getter = PROPERTIES[prop]
            columns[var] = getter(ctx)

        return columns

    def _sum(self, summand: str, columns: Dict[str, List[Any]], index: int) -> str:
        """``((s0) + (s1) + ...)`` over the values, or ``0`` without any.

        The first property sets the number of terms; extra values in later ones
        are ignored, and a later property with fewer values raises a
        :class:`~jaff.codegen.templating.errors.TemplateError`.
        """
        count = len(next(iter(columns.values())))
        terms = []
        for i in range(count):
            term = summand
            for var, values in columns.items():
                if i >= len(values):
                    raise TemplateError(
                        f"{self.keyword} properties differ in length: "
                        f"{', '.join(self.props)}",
                        body_index=index,
                    )

                term = term.replace(f"${var}$", str(values[i]))
            terms.append(f"({term})")

        return f"({' + '.join(terms)})" if terms else "0"
