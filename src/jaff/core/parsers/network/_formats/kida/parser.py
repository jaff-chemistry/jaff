"""KIDA parser: a single fixed-width reaction line-type."""

from .._parser import Parser, register
from .._record import ParsedRecord, ParseResult
from .reaction import KIDA_OPEN_TMAX, KidaReaction


@register
class KidaParser(Parser):
    """Parses the KIDA fixed-width reaction format."""

    name = "kida"
    priority = 80

    def __init__(self):
        self.handlers = [KidaReaction()]

    def process(self, records) -> ParseResult:
        """Parse each KIDA reaction record into a :class:`ParsedRecord`."""
        state = self._initial_state()
        by_name = {h.name: h for h in self.handlers}
        parsed = []

        for rec in records:
            handler = by_name[rec.format]
            parsed.append((rec, handler.parse(rec.line, rec.nline, state, self.file)))

        _resolve_open_tmax([fields for _, fields in parsed])

        reactions = [
            ParsedRecord(**fields, source_index=rec.source_index, sub_order=rec.sub_order)
            for rec, fields in parsed
        ]

        return ParseResult(reactions, {})


def _resolve_open_tmax(rows: list[dict]) -> None:
    """Map KIDA's ``Tmax = 9999`` "no upper limit" to ``None`` in place.

    9999 is a real bound when another row of the same reaction starts at a
    higher temperature (e.g. ``3000-9999`` followed by ``10000-20000``);
    otherwise the row is open-ended.
    """
    starts: dict[tuple, list[float]] = {}
    for f in rows:
        if f["tmin"] is not None:
            starts.setdefault((tuple(f["r"]), tuple(f["p"])), []).append(f["tmin"])

    for f in rows:
        if f["tmax"] != KIDA_OPEN_TMAX:
            continue

        lo = float("-inf") if f["tmin"] is None else f["tmin"]
        later = starts.get((tuple(f["r"]), tuple(f["p"])), [])
        if not any(t > lo for t in later):
            f["tmax"] = None
