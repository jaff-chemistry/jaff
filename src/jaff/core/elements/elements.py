# ABOUTME: Elements: ordered, unique collection of the elements in a network's species
# ABOUTME: Derived from a Specie list, with lookup and per-element helpers

"""Ordered, unique :class:`Elements` collection for JAFF networks.

``Elements`` is derived from a list of ``Specie`` objects, with helpers for
building truth and density matrices used in stoichiometry calculations.  It is
a flyweight: instances built from the same ordered species sequence are reused.
"""

from __future__ import annotations

from functools import cache
from typing import TYPE_CHECKING, Iterable, Union

from ...types import Catalogue, Vector
from .element import Element

if TYPE_CHECKING:
    from ..species import Specie

SpeciesInput = Union["Specie", str, Iterable[Union["Specie", str]]]
"""One species or name, or any iterable of them (e.g. a ``Species`` catalogue)."""


class Elements(Catalogue):
    """Sorted, deduplicated collection of elements derived from a species list.

    ``Elements`` is also a flyweight: instances built from the same ordered
    species sequence are reused.  The internal element order follows
    alphabetical sort on element symbol, which fixes the row order of the
    composition matrices (:meth:`truth_matrix`, :meth:`density_matrix`) and
    of every per-element accessor (:meth:`symbols`, :meth:`masses`, ...); the
    column order follows the caller's species order.

    Elements are extracted from the exploded atom list of each species; charge
    tokens (``+``, ``-``) and the electron are ignored.  Items can be looked up
    by position or by symbol (``elements[0]``, ``elements["H"]``), and by full
    name via :meth:`from_name`.  A :class:`~jaff.core.Network` exposes the
    collection for its species as ``net.elements``.

    Attributes
    ----------
    species : list[Specie]
        The input species whose atoms were used to populate this collection,
        stored as a list of ``Specie`` objects.
    count : int
        Number of unique elements present (inherited from ``Catalogue``).

    Examples
    --------
    >>> from jaff.core import Elements
    >>> elements = Elements(["CO", "H2O", "CH4"])
    >>> elements.symbols()
    ['C', 'H', 'O']
    >>> elements.count
    3
    >>> elements["H"].name
    'Hydrogen'
    >>> Elements(["CO", "H2O", "CH4"]) is elements
    True
    """

    _register: dict = {}

    @staticmethod
    def __get_species_list(species: SpeciesInput) -> list[Specie]:
        """Normalise the *species* argument to a ``list[Specie]``.

        A name or a ``Specie`` is one species; any other iterable (``list``,
        ``tuple``, a ``Species`` catalogue, a ``Vector``) holds names and/or
        ``Specie`` objects, names becoming ``Specie`` objects at their position.
        """
        from ..species import Specie as _Specie

        if isinstance(species, (str, _Specie)):
            species = [species]

        return [
            _Specie(item, idx) if isinstance(item, str) else item
            for idx, item in enumerate(species)
        ]

    def __new__(cls, species: SpeciesInput):
        """Return the flyweight ``Elements`` instance for the given species set.

        Parameters
        ----------
        species : Specie | str | Iterable[Specie | str]
            Species whose atoms define the element set: one species or name,
            or any iterable of them (e.g. a ``Species`` catalogue).

        Returns
        -------
        Elements
            Existing cached instance keyed by the ordered species sequence,
            or a newly allocated one.
        """
        _species = cls.__get_species_list(species)
        _serialized: tuple[str, ...] = tuple(str(s) for s in _species)
        if _serialized in cls._register:
            return cls._register[_serialized]

        instance = super().__new__(cls)
        cls._register[_serialized] = instance

        return instance

    def __init__(self, species: SpeciesInput) -> None:
        """Build the element collection from *species*.

        Parameters
        ----------
        species : Specie | str | Iterable[Specie | str]
            One species or name, or any iterable of them (``list``, ``tuple``, a
            ``Species`` catalogue, ...), whose constituent atoms define the
            element set.  Plain strings are converted to ``Specie`` objects on
            the fly.

        Raises
        ------
        KeyError
            If a species (or one of its atoms) is not in the mass dictionary.
        """
        if getattr(self, "_initialized", False):
            return

        self.species: list[Specie] = self.__get_species_list(species)

        self.__set_elements()
        self._initialized = True

    def __repr__(self):
        return f"Catalogue({self.symbols()!r})"

    def __set_elements(self) -> None:
        """Collect unique alphabetic atoms across all species and build indices."""
        elements: set[str] = set()

        for specie in self.species:
            elements |= set(specie.exploded)  # type: ignore[arg-type]

        # Filter out charge tokens ('+', '-') — only real element symbols remain.
        _list = sorted(list({Element(e) for e in elements if e.isalpha()}))

        _by_name = {e.name: e for e in _list}
        _by_symbol = {e.symbol: e for e in _list}

        super().__init__(_list, _by_symbol)
        # _by_serialized is reused here to store the name→Element mapping.
        self._by_serialized = _by_name

    @cache
    def truth_matrix(self) -> list[list[int]]:
        """Binary element-presence matrix.

        Returns
        -------
        list[list[int]]
            A 2-D integer matrix of shape ``(n_elements, n_species)``.
            Entry ``[i][j]`` is ``1`` if element *i* appears at least once in
            species *j*, otherwise ``0``.

        Notes
        -----
        The result is cached after the first call (via ``functools.cache``).
        Row order matches the sorted element list; column order matches the
        order of *species* passed to ``__init__``.

        Examples
        --------
        >>> elements = Elements(["CO", "H2O", "CH4"])
        >>> elements.symbols()
        ['C', 'H', 'O']
        >>> elements.truth_matrix()
        [[1, 0, 1], [0, 1, 1], [1, 1, 0]]
        """
        element_truth_matrix: list[list[int]] = [
            [0] * len(self.species) for _ in range(self.count)
        ]

        for i, element in enumerate(self._list):
            for j, specie in enumerate(self.species):
                element_truth_matrix[i][j] = int(str(element) in specie.exploded)

        return element_truth_matrix

    @cache
    def density_matrix(self) -> list[list[int]]:
        """Atom-count matrix (stoichiometric composition matrix).

        Returns
        -------
        list[list[int]]
            A 2-D integer matrix of shape ``(n_elements, n_species)``.
            Entry ``[i][j]`` is the number of atoms of element *i* contained
            in species *j* (e.g. entry for H in H2O is ``2``).

        Notes
        -----
        The result is cached after the first call (via ``functools.cache``).
        Row order matches the sorted element list; column order matches the
        order of *species* passed to ``__init__``.

        Examples
        --------
        >>> elements = Elements(["CO", "H2O", "CH4"])
        >>> elements.symbols()
        ['C', 'H', 'O']
        >>> elements.density_matrix()
        [[1, 0, 1], [0, 2, 4], [1, 1, 0]]
        """
        element_density_matrix: list[list[int]] = [
            [0] * len(self.species) for _ in range(self.count)
        ]

        for i, element in enumerate(self._list):
            for j, specie in enumerate(self.species):
                element_density_matrix[i][j] = specie.exploded.count(element.symbol)

        return element_density_matrix

    def from_name(self, name: str) -> Element:
        """Return the ``Element`` with the given full name (e.g. ``"Hydrogen"``).

        Parameters
        ----------
        name : str
            Full element name as stored in the mass dictionary (case-sensitive,
            capitalised, e.g. ``"Hydrogen"``, ``"Carbon"``).

        Returns
        -------
        Element
            The ``Element`` flyweight for the requested element.

        Raises
        ------
        KeyError
            If no element of this collection has that name.

        Examples
        --------
        >>> Elements(["CO", "H2O"]).from_name("Oxygen")
        ElementObject(symbol='O')
        """
        return self._by_serialized[name]

    def from_symbol(self, symbol: str) -> Element:
        """Return the ``Element`` with the given periodic-table symbol.

        Parameters
        ----------
        symbol : str
            Element symbol (case-sensitive, e.g. ``"H"``, ``"He"``).

        Returns
        -------
        Element
            The ``Element`` flyweight for the requested element.

        Raises
        ------
        KeyError
            If no element of this collection has that symbol.
        """
        return self._by_prop[symbol]

    def get_list(self) -> list[Element]:
        """Return the underlying ordered list of ``Element`` objects.

        The order is alphabetical by symbol and matches the row order of the
        composition matrices (:meth:`truth_matrix`, :meth:`density_matrix`).

        Returns
        -------
        list[Element]
            All unique elements, sorted alphabetically by symbol.
        """
        return self._list

    def symbols(self) -> Vector[str]:
        """Return a ``Vector`` of element symbols in sorted order.

        This alphabetical order is the canonical order shared by all other
        per-element accessors and by the composition matrices.

        Returns
        -------
        Vector[str]
            Periodic-table symbols (e.g. ``["C", "H", "O"]``).
        """
        return Vector([e.symbol for e in self])

    def names(self) -> Vector[str]:
        """Return a ``Vector`` of full element names in sorted order.

        Returns
        -------
        Vector[str]
            Full element names in the order of :meth:`symbols` (e.g.
            ``["Carbon", "Hydrogen", "Oxygen"]``).
        """
        return Vector([e.name for e in self])

    def masses(self) -> Vector[float]:
        """Return a ``Vector`` of element masses in grams (CGS).

        Unlike :meth:`atomic_masses`, values are in grams rather than atomic
        mass units.

        Returns
        -------
        Vector[float]
            Mass per element in grams, in the order of :meth:`symbols`.
        """
        return Vector([e.mass for e in self])

    def atomic_masses(self) -> Vector[float]:
        """Return a ``Vector`` of standard atomic weights in atomic mass units.

        Unlike :meth:`masses`, values are in atomic mass units rather than
        grams.

        Returns
        -------
        Vector[float]
            Standard atomic weight per element, in the order of :meth:`symbols`.
        """
        return Vector([e.atomic_mass for e in self])

    def protons(self) -> Vector[int]:
        """Return a ``Vector`` of proton counts (atomic numbers).

        Returns
        -------
        Vector[int]
            Atomic number per element, in the order of :meth:`symbols`.
        """
        return Vector([e.protons for e in self])

    def neutrons(self) -> Vector[int]:
        """Return a ``Vector`` of neutron counts for the most common isotope.

        Returns
        -------
        Vector[int]
            Neutron count per element, in the order of :meth:`symbols`.
        """
        return Vector([e.neutrons for e in self])

    def electrons(self) -> Vector[int]:
        """Return a ``Vector`` of electron counts for the neutral atom.

        Returns
        -------
        Vector[int]
            Electron count per element (equal to the atomic number), in the order of :meth:`symbols`.
        """
        return Vector([e.electrons for e in self])
