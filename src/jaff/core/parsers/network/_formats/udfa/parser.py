# ABOUTME: UDFA parser: processes the bucket of colon-delimited reaction lines
# ABOUTME: Registered with the format-parser registry

"""UDFA parser: a single colon-delimited reaction line-type."""

from .._parser import Parser, register
from .._record import ParsedRecord, ParseResult
from .reaction import UdfaReaction, UdfaReactionFallback


@register
class UdfaParser(Parser):
    """Parses the UDFA (UMIST) colon-delimited reaction format."""

    name = "udfa"
    priority = 50

    def __init__(self):
        self.handlers = [UdfaReaction(), UdfaReactionFallback()]

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
                        nline=rec.nline,
                        sub_order=rec.sub_order + i,
                    )
                )

        return ParseResult(reactions, {})

    def post_process(self, network) -> None:
        """No post-processing needed for this format."""
        return
