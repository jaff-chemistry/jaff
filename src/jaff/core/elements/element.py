# ABOUTME: Element: one record of the periodic table (mass, atomic mass and
# ABOUTME: proton/neutron/electron counts); Element("H") looks it up in the table

"""Single chemical :class:`Element`: one record of :func:`~jaff.chemistry.periodic_table`."""

from __future__ import annotations

from typing import Any, Mapping

from ...chemistry import periodic_table


class Element:
    """A chemical element: one record of :func:`~jaff.chemistry.periodic_table`.

    Every element is built once, when the table is read, and shared:
    ``Element("H")`` returns the table's record for ``"H"`` (the same object as
    ``periodic_table("H")``), so constructing it twice gives the same object.

    Elements compare, sort and hash by :attr:`symbol`; comparing an
    ``Element`` with any other type raises :exc:`TypeError`.  ``str(element)``
    returns the symbol.

    Attributes
    ----------
    symbol : str
        Periodic-table symbol (e.g. ``"H"``, ``"He"``).
    name : str
        Full element name (e.g. ``"Hydrogen"``).
    mass : float
        Mass of the element in grams (CGS), i.e. the standard atomic weight
        expressed in grams.
    atomic_mass : float or None
        Standard atomic weight in atomic mass units.
    protons : int or None
        Number of protons (atomic number).
    neutrons : int or None
        Number of neutrons in the most common isotope.
    electrons : int or None
        Number of electrons in the neutral atom.

    Notes
    -----
    The table also holds pseudo-elements such as ``"e-"`` and ``"GRAIN"``; for
    these only ``name`` and ``mass`` are set and the other numeric attributes
    are ``None``.

    Examples
    --------
    >>> from jaff.chemistry import periodic_table
    >>> from jaff.core import Element
    >>> h = Element("H")
    >>> h.name, h.protons, h.atomic_mass
    ('Hydrogen', 1, 1.008)
    >>> h is periodic_table("H")
    True
    >>> str(Element("He"))
    'He'
    >>> sorted([Element("O"), Element("C")])
    [ElementObject(symbol='C'), ElementObject(symbol='O')]
    """

    symbol: str
    name: str
    mass: float
    atomic_mass: float | None
    protons: int | None
    neutrons: int | None
    electrons: int | None

    def __new__(cls, symbol: str) -> Element:
        """Return the periodic table's record for *symbol*.

        Parameters
        ----------
        symbol : str
            Periodic-table symbol (case-sensitive, e.g. ``"He"``).

        Returns
        -------
        Element
            The shared record.

        Raises
        ------
        KeyError
            If *symbol* is not in the periodic table.
        """
        return periodic_table(symbol)

    @classmethod
    def _from_row(cls, row: Mapping[str, Any]) -> Element:
        """Build the record of one row of the bundled ``atomic_masses`` table."""
        element = object.__new__(cls)
        element.symbol = row["element"]
        element.name = row["name"]
        element.mass = row["mass"]
        element.atomic_mass = row["atomic_mass"]
        element.protons = row["protons"]
        element.neutrons = row["neutrons"]
        element.electrons = row["electrons"]

        return element

    def __repr__(self) -> str:
        """Return detailed string representation of this element.

        Returns
        -------
        str
            String of the form ``"ElementObject(symbol='C')"``.
        """
        return f"ElementObject(symbol={self.symbol!r})"

    def __str__(self) -> str:
        """Return the periodic-table symbol.

        Returns
        -------
        str
            Element symbol (e.g. ``"He"``).
        """
        return self.symbol

    def __eq__(self, other) -> bool:
        """Check equality by comparing element symbols.

        Parameters
        ----------
        other : Element
            Element to compare against.

        Returns
        -------
        bool

        Raises
        ------
        TypeError
            If *other* is not an ``Element`` instance.
        """
        if not isinstance(other, Element):
            raise TypeError(
                f"'==' not supported between instances of 'Element' and '{other}'"
            )

        return self.symbol == other.symbol

    def __lt__(self, other) -> bool:
        """Compare elements lexicographically by symbol.

        Parameters
        ----------
        other : Element
            Element to compare against.

        Returns
        -------
        bool

        Raises
        ------
        TypeError
            If *other* is not an ``Element`` instance.
        """
        if not isinstance(other, Element):
            raise TypeError(
                f"'<' not supported between instances of 'Element' and '{other}'"
            )

        return self.symbol < other.symbol

    def __hash__(self):
        """Return hash based on the element symbol.

        Returns
        -------
        int
        """
        return hash(self.symbol)
