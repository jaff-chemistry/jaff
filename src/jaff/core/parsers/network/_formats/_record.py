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

    The engine sets :attr:`format`, lower-cases :attr:`rate`, sorts the records
    by ``(source_index, sub_order)`` and hands them to ``Network``, which builds
    the species and reactions from them.

    Attributes
    ----------
    reactants : list[str]
        Reactant names, as written (pseudo-species such as ``_PHOTON``
        included).
    products : list[str]
        Product names, as written.
    tmin : float or None
        Lower temperature bound (K), or ``None``.
    tmax : float or None
        Upper temperature bound (K), or ``None``.
    rate : str
        The rate expression, a SymPy-parsable string.
    type : str
        Reaction type concluded by the parser (e.g. ``"photo"``).
    line : str
        The network-file line the reaction came from.
    source_index : int
        File-side reaction number (0-based, in file order).
    sub_order : int
        Tie-breaker for a single line expanding into several reactions.
    metadata : dict
        Parser-private data that ends up on the reaction as
        ``Reaction._metadata`` for the parser's post-processing hook (``{}``
        for every bundled parser).
    format : str or None
        Name of the parser that produced the record; set by the engine
        (``None`` until then).
    """

    def __init__(
        self,
        reactants: list[str],
        products: list[str],
        tmin: float | None,
        tmax: float | None,
        rate: str,
        type: str,
        line: str,
        source_index: int,
        sub_order: int = 0,
        metadata: dict | None = None,
        format: str | None = None,
    ):
        """Store the parsed fields; see the class attributes."""
        self.reactants: list[str] = reactants
        self.products: list[str] = products
        self.tmin: float | None = tmin
        self.tmax: float | None = tmax
        self.rate: str = rate
        self.type: str = type
        self.line: str = line
        self.source_index: int = source_index
        self.sub_order: int = sub_order
        self.metadata: dict = metadata if metadata is not None else {}
        self.format: str | None = format


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
