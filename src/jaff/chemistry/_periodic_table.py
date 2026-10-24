# ABOUTME: periodic_table(): the bundled element table as Element records, read once
# ABOUTME: per process; with a symbol, returns that one Element

"""The :func:`periodic_table` of :class:`~jaff.core.Element` records."""

from __future__ import annotations

from functools import cache
from types import MappingProxyType
from typing import TYPE_CHECKING, Mapping, overload

if TYPE_CHECKING:
    from ..core.elements import Element


@cache
def _load() -> Mapping[str, Element]:
    """Read the bundled ``atomic_masses`` table once; a read-only symbol → Element map."""
    # Imported here: jaff.core imports this package (species parsing).
    from ..core.elements.element import Element
    from ..drivers.sqlite import JaffDb

    with JaffDb() as jdb:
        rows = jdb.table("atomic_masses").all_rows()

    return MappingProxyType({row["element"]: Element._from_row(row) for row in rows})


@overload
def periodic_table() -> Mapping[str, Element]: ...


@overload
def periodic_table(symbol: str) -> Element: ...


def periodic_table(symbol: str | None = None) -> Mapping[str, Element] | Element:
    """The element table, or one element of it.

    Read from the bundled ``atomic_masses`` table on first use and shared by
    every caller for the rest of the process.  Besides the chemical elements it
    holds pseudo-elements such as ``"e-"`` and ``"GRAIN"``.

    Parameters
    ----------
    symbol : str or None, optional
        An element symbol (case-sensitive, e.g. ``"He"``); ``None`` (default)
        for the whole table.

    Returns
    -------
    Mapping[str, Element] or Element
        Without *symbol*, the whole table: a read-only mapping from symbol to
        :class:`Element`.  With *symbol*, that element.

    Raises
    ------
    KeyError
        If *symbol* is not in the table.

    Examples
    --------
    >>> from jaff.chemistry import periodic_table
    >>> periodic_table("H").mass
    1.673773e-24
    >>> periodic_table()["H"] is periodic_table("H")
    True
    """
    table = _load()
    if symbol is None:
        return table

    try:
        return table[symbol]
    except KeyError:
        raise KeyError(f"{symbol!r} is not in the periodic table") from None
