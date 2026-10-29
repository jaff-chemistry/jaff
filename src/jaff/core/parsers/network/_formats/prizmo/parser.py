# ABOUTME: PRIZMO parser: arrow-notation reactions and VARIABLES { } directives
# ABOUTME: Seeds temperature shorthand globals (t32, invtgas, ...)

"""PRIZMO parser: arrow-notation reactions + ``VARIABLES { }`` directives."""

from sympy import Basic, parse_expr
from sympy.core.function import AppliedUndef

from ......errors import ParserError

from .._parser import Parser, register
from .._record import ParsedRecord, ParseResult
from .reaction import PrizmoReaction
from .vars import PrizmoVars


@register
class PrizmoParser(Parser):
    """Parses the PRIZMO arrow-notation format and its variables block."""

    name = "prizmo"
    priority = 30

    #: Temperature shorthands common to the KROME/PRIZMO conventions. The
    #: density accessors and ``user_*`` aliases are KROME-specific and live
    #: only on the KROME parser. A file's own ``VARIABLES`` block overrides
    #: these (applied after seeding).
    BASE_GLOBALS = {
        "invt32": "1e0 / t32",
        "invte": "1e0 / te",
        "t32": "tgas/3e2",
        "te": "tgas*8.617343e-5",
        "invtgas": "1e0 / tgas",
        "sqrtgas": "sqrt(tgas)",
    }

    def __init__(self):
        self.handlers = [PrizmoVars(), PrizmoReaction()]

    def process(self, records) -> ParseResult:
        """Walk the PRIZMO bucket in file order.

        ``VARIABLES`` directives toggle the local ``parse_vars`` state and
        populate ``globals``; reaction lines are parsed into
        :class:`ParsedRecord` objects.
        """
        state = self._initial_state()
        globals_: dict = {k: parse_expr(v) for k, v in self.BASE_GLOBALS.items()}
        by_name = {h.name: h for h in self.handlers}
        reactions: list[ParsedRecord] = []

        for rec in records:
            handler = by_name[rec.format]
            if handler.is_reaction:
                fields = handler.parse(rec.line, rec.nline, state, self.file)
                reactions.append(
                    ParsedRecord(
                        **fields,
                        source_index=rec.source_index,
                        nline=rec.nline,
                        sub_order=rec.sub_order,
                    )
                )
            else:
                handler.apply(
                    rec.line, rec.nline, state, globals_, self.file, self.logger
                )

        return ParseResult(reactions, globals_)

    def post_process(self, network) -> None:
        """Number photo reactions: ``prizmo_photo(lo, hi)`` → ``photorates(n, lo, hi)``.

        ``n`` is the reaction's 1-based position in
        ``network.reactions.with_type("photo")``, so every format numbers photo
        reactions the same way.  A photo reaction without the placeholder (e.g.
        one whose rate comes from a ``.jfunc`` ``chemRate``) leaves its number
        unused.

        Raises
        ------
        ParserError
            If a placeholder is left in any of this parser's reactions.
        """
        for n, rea in enumerate(network.reactions.with_type("photo"), start=1):
            rea.rate_segments.map_rates(lambda rate, n=n: self._number_photo(rate, n))

        for rea in network.reactions.with_format(self.name):
            exprs = [rea.rate_coefficient, *(seg.rate for seg in rea.rate_segments)]
            if any(self._placeholders(expr) for expr in exprs):
                raise ParserError(f"Unresolved photo placeholder in {rea}")

    @staticmethod
    def _placeholders(expr) -> list:
        """The ``prizmo_photo(lo, hi)`` calls in *expr*."""
        if not isinstance(expr, Basic):
            return []

        return [f for f in expr.atoms(AppliedUndef) if f.func.__name__ == "prizmo_photo"]

    @classmethod
    def _number_photo(cls, expr, n: int):
        """Replace every placeholder in *expr* with ``photorates(n, lo, hi)``."""
        # Local import: a module-level one would cycle via jaff.core.network.
        from .....network import NetworkSymbols

        placeholders = cls._placeholders(expr)
        if not placeholders:
            return expr

        return expr.xreplace(
            {f: NetworkSymbols.photorates(n, *f.args) for f in placeholders}
        )
