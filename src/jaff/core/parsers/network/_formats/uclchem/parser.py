# ABOUTME: UCLCHEM parser: processes the bucket of NAN-sentinel reaction lines
# ABOUTME: Registered with the format-parser registry

"""UCLCHEM parser: a single ``NAN``-sentinel reaction line-type."""

from .._parser import Parser, register
from .._record import ParsedRecord, ParseResult
from .reaction import UclchemReaction


@register
class UclchemParser(Parser):
    """Parses the UCLCHEM comma-delimited reaction format."""

    name = "uclchem"
    priority = 70

    def __init__(self):
        self.handlers = [UclchemReaction()]

    def process(self, records) -> ParseResult:
        """Parse each UCLCHEM reaction record into a :class:`ParsedRecord`."""
        state = self._initial_state()
        by_name = {h.name: h for h in self.handlers}
        reactions: list[ParsedRecord] = []

        for rec in records:
            handler = by_name[rec.format]
            fields = handler.parse(rec.line, rec.nline, state, self.file)
            reactions.append(
                ParsedRecord(
                    **fields,
                    source_index=rec.source_index,
                    sub_order=rec.sub_order,
                )
            )

        return ParseResult(reactions, {})

    def post_process(self, network) -> None:
        """No post-processing needed for this format."""
        return
