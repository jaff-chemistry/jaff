# ABOUTME: NetworkParser: detects each network-file line's format and dispatches it
# ABOUTME: Buckets lines per format parser and collects parsed reactions + globals

"""Low-level reaction network file parser for multiple astrochemical formats.

``NetworkParser`` reads a single network file and turns each reaction line
into a format-independent
:class:`~jaff.core.parsers.network._formats.ParsedRecord` (reactants,
products, temperature bounds, rate string, type, source line, file-side
number, the parser's name in ``format`` and parser-private ``metadata``).

Supported file formats
----------------------
The parser auto-detects the format from line patterns.  Each format is a
self-contained parser under ``parsers.network._formats``; the engine discovers
them through :func:`~jaff.core.parsers.network._formats.all_parsers`.  Detection is
driven by each parser's handler descriptors, ordered by their declared
``priority`` (lower is matched first):

1. **PRIZMO** — arrow-notation (``->``) followed by a ``[tmin, tmax]``
   bracket, both before any comma (the brackets are required; either bound may be
   left empty).  Variables in a ``VARIABLES { }`` block.
2. **UDFA** — colon-delimited, fixed-column database from the UMIST project; a
   line starts with its index and a colon.
3. **KROME** — comma-separated, declared via ``@format:`` header.
   Variable aliases via ``@var:``.
4. **UCLCHEM** — comma-separated with a ``NAN`` sentinel, includes grain
   surface reactions.
5. **KIDA** — fixed-width column format from the KIDA database.

Last, any other line with ``->`` or ``:`` goes to a PRIZMO or UDFA fallback handler,
which reports it as an invalid reaction of that format instead of skipping it.

Rate normalization
------------------
Each format parser normalizes its own rate strings: format-specific symbols
(e.g. ``user_crflux``, ``user_av``) are rewritten by the parser, and shorthand
variables such as ``t32`` or ``invtgas`` are supplied as globals by the KROME
and PRIZMO parsers (their ``BASE_GLOBALS``).  The engine itself only
lower-cases every rate string after parsing, before ``Network`` passes it to
SymPy.
"""

import logging
from pathlib import Path

from sympy import Basic

from ....common import resolve_symbolic_dependencies
from ....io import JaffLogger, jaff_progress
from ._formats import ParsedRecord, Parser, Record, all_parsers


class NetworkParser:
    """Auto-detecting parser for astrochemical reaction network files.

    On construction the file is read, reactions are extracted, and rate
    strings are normalised.  Use as a context manager to ensure internal
    pattern state is freed after use.

    Parameters
    ----------
    file : str | Path
        Path to the network file.
    logger : logging.Logger | None, optional
        Logger instance.  A new JAFF logger is created if ``None``.

    Raises
    ------
    ValueError
        If *file* is not a ``str`` or ``Path``.
    FileNotFoundError
        If *file* does not exist on disk.
    ParserError
        On syntax errors encountered while parsing the file.
    """

    def __init__(self, file: str | Path, logger: logging.Logger | None = None):
        """Parse *file* and prepare the internal parsed-reaction list.

        Parameters
        ----------
        file : str | Path
            Path to the network file.
        logger : logging.Logger | None, optional
            External logger.  Defaults to a new JAFF logger.
        """
        if isinstance(file, str):
            file = Path(file)
        if not isinstance(file, (str, Path)):
            raise ValueError(f"Invalid file type detected for {file}: {type(file)}")

        file = file.resolve()
        if not file.exists():
            raise FileNotFoundError(file)

        self.__file: Path = file
        self.__logger: logging.Logger = logger or JaffLogger().get_logger()
        self.__globals: dict[str, Basic] = {}

        self.__parsed_list: list[ParsedRecord] = []

        self.__parsers = all_parsers()
        for parser in self.__parsers:
            parser.file = self.__file
            parser.logger = self.__logger

        self.__descriptors = sorted(
            (
                (
                    handler.priority,
                    handler.global_re,
                    handler.is_reaction,
                    handler.name,
                    parser,
                )
                for parser in self.__parsers
                for handler in parser.handlers
            ),
            key=lambda d: d[0],
        )

        self.__source_counter: int = 0
        self.__buckets: dict[str, list[Record]] = {}

        self.__parse_file()

        # Each parsed record is kept with the name of the parser that made it.
        reactions = []
        self.__used_parsers: list[Parser] = []
        for parser in self.__parsers:
            bucket = self.__buckets.get(parser.name)
            if not bucket:
                continue

            self.__used_parsers.append(parser)
            result = parser.process(bucket)
            reactions.extend((pr, parser.name) for pr in result.reactions)
            self.__globals.update(result.globals)

        reactions.sort(key=lambda item: (item[0].source_index, item[0].sub_order))
        for record, name in reactions:
            record.format = name
        self.__parsed_list[:] = [record for record, _ in reactions]

        self.__normalize_rates()
        self.__globals = resolve_symbolic_dependencies(self.__globals, fname=self.__file)

    def __enter__(self) -> "NetworkParser":
        """Return self when entering a ``with`` block.

        Returns
        -------
        NetworkParser
        """
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Free the registered parsers on context manager exit."""
        self.__parsers.clear()

        return

    def get_parsers(self) -> list[Parser]:
        """Return the parsers that produced records, in priority order.

        ``Network`` calls each one's ``post_process`` once the network is built.
        """
        return list(self.__used_parsers)

    def get_parsed(self) -> tuple[list[ParsedRecord], dict[str, Basic]]:
        """Return the parsed reactions and the resolved global variable map.

        Returns
        -------
        tuple[list[ParsedRecord], dict[str, Basic]]
            - ``list[ParsedRecord]``: one record per reaction, in
              ``(source_index, sub_order)`` order, its ``format`` set and its
              ``rate`` lower-cased.
            - ``dict[str, Basic]``: global symbolic constants defined in the
              file (e.g. ``@var`` entries), with all inter-dependencies
              resolved via SymPy substitution.
        """
        return self.__parsed_list, resolve_symbolic_dependencies(
            dep_map=self.__globals, fname=self.__file
        )

    def __parse_file(self) -> None:
        """Read the network file line-by-line and detect/bucket each line.

        Iterates over every line of :attr:`__file`, advancing the line counter
        and calling :meth:`__parse_line` for each.
        """
        with open(self.__file, "r", encoding="utf-8") as f:
            lines = f.readlines()
            for i, line in enumerate(
                jaff_progress.track(lines, description=f"Parsing {self.__file.name}")
            ):
                self.__parse_line(line, i + 1)

    def __parse_line(self, line: str, nline: int) -> None:
        """Detect the owning parser for *line* and bucket a raw :class:`Record`.

        Iterates through the handler descriptors in priority order.  The first
        handler whose global regex matches owns the line: a :class:`Record` is
        appended to its parser's bucket, tagged with the matched handler name
        and — for a reaction handler — the next ``source_index``.  Directives
        get ``source_index = -1``.  If no handler matches the line is silently
        skipped.
        """
        if not line.strip():
            return

        for _priority, global_re, is_reaction, name, parser in self.__descriptors:
            if not global_re.match(line):
                continue

            if is_reaction:
                source_index = self.__source_counter
                self.__source_counter += 1
            else:
                source_index = -1

            self.__buckets.setdefault(parser.name, []).append(
                Record(
                    source_index=source_index,
                    line=line,
                    nline=nline,
                    format=name,
                )
            )

            break

    def __normalize_rates(self):
        """Lower-case all rate strings so SymPy ``parse_expr`` is case-insensitive."""
        for record in self.__parsed_list:
            assert isinstance(record.rate, str)
            record.rate = record.rate.lower()
