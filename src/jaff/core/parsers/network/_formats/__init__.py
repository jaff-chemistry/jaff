# ABOUTME: Registry of network-file format parsers, one subpackage per format
# ABOUTME: Exposes Parser, register, all_parsers and the Record types

"""Registry for network-file parsers.

Each supported format lives in its own subpackage (``krome``, ``prizmo``,
``udfa``, ``uclchem``, ``kida``) and registers one concrete
:class:`~._parser.Parser` subclass (in its ``parser.py``) via the
:func:`~._parser.register` decorator.  A parser owns plain **handler** objects
(one per line-type, kept in files like ``krome/header.py``); handlers carry no
base class.

The engine discovers parsers through :func:`~._parser.all_parsers`, which
returns them sorted by ``priority``.  Adding a new format requires only a new
subpackage with a ``@register``-ed ``Parser``; no engine edits.
"""

from ._parser import Parser, all_parsers, register
from ._record import ParsedRecord, ParseResult, Record

__all__ = [
    "Parser",
    "ParseResult",
    "ParsedRecord",
    "Record",
    "register",
    "all_parsers",
]
