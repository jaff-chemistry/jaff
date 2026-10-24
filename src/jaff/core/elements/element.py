# ABOUTME: Element: flyweight for one chemical element loaded from the JAFF mass table
# ABOUTME: Exposes mass, atomic mass and proton/neutron/electron counts

"""Single chemical :class:`Element` flyweight for JAFF networks.

An ``Element`` represents a single chemical element, loaded from the JAFF mass
dictionary (CGS units throughout).  Instances are flyweights: constructing
``Element("H")`` twice returns the same object.
"""

from __future__ import annotations

from ..chemistry import element_table


class Element:
    """A chemical element loaded from the JAFF mass dictionary.

    Instances are flyweights: constructing ``Element("H")`` twice returns the
    same object.  The first construction populates all attributes from the
    mass dictionary; subsequent constructions are no-ops.

    Elements compare, sort and hash by :attr:`symbol`; comparing an
    ``Element`` with any other type raises :exc:`TypeError`.  ``str(element)``
    returns the symbol.

    Attributes
    ----------
    symbol : str
        Periodic-table symbol (e.g. ``"H"``, ``"He"``).
    name : str
        Full element name as stored in the mass dictionary (e.g.
        ``"Hydrogen"``).
    mass : float
        Mass of the element in grams (CGS), i.e. the standard atomic weight
        expressed in grams.
    atomic_mass : float
        Standard atomic weight in atomic mass units.
    protons : int
        Number of protons (atomic number).
    neutrons : int
        Number of neutrons in the most common isotope.
    electrons : int
        Number of electrons in the neutral atom.

    Notes
    -----
    The mass dictionary also holds pseudo-elements such as ``"e-"`` and
    ``"GRAIN"``; for these only ``name`` and ``mass`` are set and the other
    numeric attributes are ``None``.

    Examples
    --------
    >>> from jaff.core import Element
    >>> h = Element("H")
    >>> h.name, h.protons, h.atomic_mass
    ('Hydrogen', 1, 1.008)
    >>> Element("H") is h
    True
    >>> str(Element("He"))
    'He'
    >>> sorted([Element("O"), Element("C")])
    [ElementObject(symbol='C'), ElementObject(symbol='O')]
    """

    _register: dict = {}

    def __new__(cls, symbol: str):
        """Return the flyweight instance for *symbol*, creating it if absent.

        Parameters
        ----------
        symbol : str
            Periodic-table symbol (case-sensitive).

        Returns
        -------
        Element
            Existing cached instance, or a newly allocated one registered for
            future calls.
        """
        # Return the cached instance if this element has already been built;
        # otherwise create a fresh one and register it for future look-ups.
        if symbol in cls._register:
            return cls._register[symbol]

        instance = super().__new__(cls)
        cls._register[symbol] = instance

        return instance

    def __init__(self, symbol: str):
        """Initialise an Element from the mass dictionary.

        Parameters
        ----------
        symbol : str
            Periodic-table symbol of the element (case-sensitive, e.g. ``"He"``).

        Raises
        ------
        KeyError
            If *symbol* is not present in the mass dictionary.
        """
        if getattr(self, "_initialized", False):
            return

        mass_dict = element_table()

        if symbol not in mass_dict:
            raise KeyError(f"No specie found in mass dictionary: {symbol}")

        self.symbol: str = symbol
        self.name: str = mass_dict[symbol]["name"]
        self.mass: float = mass_dict[symbol]["mass"]
        self.atomic_mass: float = mass_dict[symbol]["atomic_mass"]
        self.protons: int = mass_dict[symbol]["protons"]
        self.neutrons: int = mass_dict[symbol]["neutrons"]
        self.electrons: int = mass_dict[symbol]["electrons"]
        self._initialized = True

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
