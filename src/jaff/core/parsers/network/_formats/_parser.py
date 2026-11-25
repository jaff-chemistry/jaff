# ABOUTME: Parser ABC for one network-file format, plus the parser registry
# ABOUTME: all_parsers() discovers the format subpackages and sorts them by priority

"""The single ABC for a per-subpackage network parser plus its registry.

Each format subpackage (``kida``, ``krome``, ``prizmo``, ``udfa``,
``uclchem``) defines exactly one concrete :class:`Parser` in its ``parser.py``
and registers it with :func:`register`.  A parser owns a list of plain
**handler** objects (one per line-type) and walks its own bucket of
:class:`~._record.Record` in file order, dispatching each record to the
handler that detected it.

The engine discovers parsers through :func:`all_parsers`, which imports the
subpackages (triggering ``@register``) and returns one instance per parser,
sorted by ``priority``.  Handlers carry no base class; only :class:`Parser` is
abstract.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ....network import Network
    from ._record import ParseResult, Record

_PARSER_REGISTRY: list[type["Parser"]] = []


def register(cls: type["Parser"]) -> type["Parser"]:
    """Register *cls* as an available network parser.

    Parameters
    ----------
    cls : type[Parser]
        The parser class to register.

    Returns
    -------
    type[Parser]
        *cls* unchanged, so the decorator is transparent.
    """
    _PARSER_REGISTRY.append(cls)

    return cls


class Parser(ABC):
    """Owns and processes one subpackage's bucket of detected records.

    Class attributes
    ----------------
    name : str
        Unique parser identifier (the subpackage folder name); the engine
        buckets records under this name.
    priority : int
        Ordering hint (unused for line detection, which is driven by each
        handler's own ``priority``).
    handlers : list
        The parser's handler objects.  Each handler exposes ``name``,
        ``priority``, ``is_reaction``, and a static ``global_re`` for detection,
        plus a ``parse`` (reaction) or ``apply`` (directive) method.

    The engine sets :attr:`file` and :attr:`logger` before calling
    :meth:`process`, so handlers can raise contextual
    :class:`~jaff.errors.ParserError` and emit warnings.
    """

    name: str
    priority: int
    handlers: list

    #: Set by the engine before :meth:`process`.
    file = None
    logger = None

    def _initial_state(self) -> dict:
        """Seed the parser's local state from any handler ``default_state``.

        Handlers that need mutable, file-order state (e.g. the KROME
        ``@format`` column layout, PRIZMO's ``parse_vars`` toggle) expose a
        ``default_state()``; their defaults are merged so a reaction line
        that precedes any directive is parsed with those defaults.
        """
        state: dict = {}
        for handler in self.handlers:
            default_state = getattr(handler, "default_state", None)
            if default_state is not None:
                state.update(default_state())

        return state

    @abstractmethod
    def process(self, records: list["Record"]) -> "ParseResult":
        """Walk *records* in file order and return the parsed reactions + globals."""
        ...

    @abstractmethod
    def post_process(self, network: "Network") -> None:
        """Adjust this parser's reactions once the network has been built.

        Called once, only for parsers that produced records, after every
        reaction is built from its standardised rate and before ``rc_<N>``
        references are resolved: a hook writing a shorthand symbol (``n_X``,
        ...) standardises it itself (``network.symbols.standardize``).  The
        reaction and
        species catalogues are frozen: a hook may change reactions (e.g. their
        rate segments, props or ``_metadata``) but must not add or remove
        reactions or species.  Use ``network.reactions.with_format(self.name)``
        to select this parser's reactions.  Parsers with nothing to do
        implement it as a no-op.
        """
        ...


def all_parsers() -> list[Parser]:
    """Instantiate every registered parser, sorted by priority.

    Importing the format subpackages here triggers their ``@register``
    decorators (via each subpackage ``__init__`` importing its ``parser``),
    populating :data:`_PARSER_REGISTRY`.

    Returns
    -------
    list[Parser]
        One instance per registered parser, in ascending priority order.
    """
    from jaff.common._helper import import_subpackages

    parent = __name__.rsplit(".", 1)[0]
    import_subpackages(parent)

    return sorted((cls() for cls in _PARSER_REGISTRY), key=lambda p: p.priority)
