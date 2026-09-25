"""UDFA parser: a single colon-delimited reaction line-type."""

from .._parser import Parser, register
from .._record import ParsedRecord, ParseResult
from .reaction import UdfaReaction


@register
class UdfaParser(Parser):
    """Parses the UDFA (UMIST) colon-delimited reaction format."""

    name = "udfa"
    priority = 50

    def __init__(self):
        self.handlers = [UdfaReaction()]

    def process(self, records) -> ParseResult:
        """Parse each UDFA reaction record into a :class:`ParsedRecord`."""
        state = self._initial_state()
        by_name = {h.name: h for h in self.handlers}
        reactions: list[ParsedRecord] = []

        for rec in records:
            handler = by_name[rec.format]
            segments = handler.parse(rec.line, rec.nline, state, self.file)
            for i, fields in enumerate(segments):
                reactions.append(
                    ParsedRecord(
                        **fields,
                        source_index=rec.source_index,
                        sub_order=rec.sub_order + i,
                    )
                )

        return ParseResult(reactions, {})
