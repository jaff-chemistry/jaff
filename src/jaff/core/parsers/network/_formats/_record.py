# ABOUTME: Record types passed between the network-parser engine and format parsers
# ABOUTME: Record (detected line), ParsedRecord (parsed reaction), ParseResult


class Record:
    """A detected-but-not-yet-processed line: a reaction or a directive.

    Captured in the engine's detection phase and handed to the owning format's
    ``process`` method.  The engine leaves ``metadata`` empty and no bundled
    parser reads it.

    Attributes
    ----------
    source_index : int
        File-side reaction number (0-based, stable under dedup); ``-1`` for a
        directive line.
    line : str
        Raw file line.
    nline : int
        1-based file line number (for error messages).
    format : str
        Name of the handler that detected the line; the owning parser uses it
        to dispatch the record.
    metadata : dict
        Optional per-line data; ``{}`` unless passed in (the engine never
        does).
    sub_order : int
        Tie-breaker for a single line expanding into several reactions.
    """

    def __init__(
        self,
        source_index: int,
        line: str,
        nline: int,
        format: str,
        metadata: dict | None = None,
        sub_order: int = 0,
    ):
        self.source_index: int = source_index
        self.line: str = line
        self.nline: int = nline
        self.format: str = format
        self.metadata: dict = metadata if metadata is not None else {}
        self.sub_order: int = sub_order


class ParsedRecord:
    """A fully parsed reaction produced by a format family's ``process``.

    The engine sorts these by ``(source_index, sub_order)`` and converts each to
    the ``parsedListProps`` dict consumed by ``Network``.  ``metadata`` is
    parser-private data that ends up on the reaction as ``Reaction._metadata``
    for the parser's post-processing hook (``{}`` for every bundled parser).
    """

    def __init__(
        self,
        r: list[str],
        p: list[str],
        tmin: float | None,
        tmax: float | None,
        rate: str,
        type: str,
        string: str,
        source_index: int,
        sub_order: int = 0,
        metadata: dict | None = None,
    ):
        self.r = r
        self.p = p
        self.tmin = tmin
        self.tmax = tmax
        self.rate = rate
        self.type = type
        self.string = string
        self.source_index = source_index
        self.sub_order = sub_order
        self.metadata: dict = metadata if metadata is not None else {}

    def as_parsed_props(self, format: str) -> dict:
        """Convert to the engine's ``parsedListProps`` dict, tagged with *format*.

        Parameters
        ----------
        format : str
            Name of the parser that produced this record (``Parser.name``).
        """
        return {
            "r": self.r,
            "p": self.p,
            "tmin": self.tmin,
            "tmax": self.tmax,
            "rate": self.rate,
            "type": self.type,
            "string": self.string,
            "source_index": self.source_index,
            "format": format,
            "metadata": self.metadata,
        }


class ParseResult:
    """The output of a :class:`~._parser.Parser`'s ``process``.

    Attributes
    ----------
    reactions : list[ParsedRecord]
        Fully parsed reactions produced from the parser's bucket.
    globals : dict
        Symbolic globals (e.g. ``@var`` / ``VARIABLES`` entries) the parser
        discovered; merged into the engine's global symbol map.
    """

    def __init__(self, reactions: list["ParsedRecord"], globals: dict):
        self.reactions: list["ParsedRecord"] = reactions
        self.globals: dict = globals
