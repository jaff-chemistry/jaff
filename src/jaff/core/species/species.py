# ABOUTME: Species: ordered, name-indexed catalogue of a network's Specie objects
# ABOUTME: Lookup by name, serialized form or index

"""Ordered, name-indexed :class:`Species` catalogue for JAFF networks.

``Species`` supports look-up by name (``species["H2O"]``), by serialized form
(``species["+/H/H/O"]``), and by integer index (``species[0]``).
"""

from __future__ import annotations

from functools import cached_property

from ...types import Catalogue, Vector
from ..elements import Elements
from .specie import Specie


class Species(Catalogue[Specie]):
    """Ordered, name-indexed collection of ``Specie`` objects.

    Species supports look-up by name (``species["H2O"]``), by serialized form
    (``species["+/H/H/O"]``), and by integer index (``species[0]``).  The
    ``ne`` parameter on many accessor methods excludes the electron species
    (``"e-"``), which is often treated separately in network solvers.

    Vector accessors (:meth:`names`, :meth:`masses`, :meth:`charges`, ...)
    return one entry per species in catalogue order.

    Attributes
    ----------
    count : int
        Total number of species in the collection.

    Examples
    --------
    >>> from jaff.core.species import Species
    >>> species = Species(["H2", "H+", "e-"])
    >>> species["H+"].name
    'H+'
    >>> species["H/H"].name
    'H2'
    >>> species[0].name
    'H2'
    >>> species.names(ne=True)
    ['H2', 'H+']
    >>> species.count
    3
    """

    def __init__(
        self,
        species: list[Specie] | list[str] | None = None,
        check_length: bool = True,
    ):
        """Initialise the species catalogue.

        Parameters
        ----------
        species : list[Specie] | list[str] | None, optional
            Initial species.  Plain strings are converted to ``Specie``
            objects with indices assigned in list order.  If ``None``, an
            empty catalogue is created (items can be added with ``add``).
        check_length : bool, optional
            If ``True`` (default), the base ``Catalogue`` verifies that the
            list and name-dict have the same length.  Set to ``False`` when
            constructing a catalogue from reactants/products that may contain
            duplicate species.
        """
        _by_name: dict[str, Specie] | None = None
        _by_serialized: dict[str, Specie] = {}

        if species is not None:
            if species and isinstance(species[0], str):
                species = [Specie(name, idx) for idx, name in enumerate(species)]  # type: ignore[arg-type]
            _by_name = {sp.name: sp for sp in species}  # type: ignore
            _by_serialized = {sp.serialized: sp for sp in species}  # type: ignore

        _species: list[Specie] = species  # type: ignore

        super().__init__(_species, _by_name, check_length)
        self._by_serialized = _by_serialized
        self.__frozen: bool = False

    def _freeze(self) -> None:
        """Disallow further :meth:`add` calls.

        ``Network`` freezes its catalogue before running the parsers'
        post-processing hooks, which may not add species.
        """
        self.__frozen = True

    def __repr__(self):
        return f"Catalogue({self.names()!r})"

    def add(self, specie: Specie) -> None:
        """Append a new species to the catalogue if not already present.

        Duplicates (by name) are silently ignored.

        Parameters
        ----------
        specie : Specie
            The species to add.

        Raises
        ------
        ValueError
            If *specie* is not a ``Specie`` instance.
        RuntimeError
            If the catalogue has been frozen (``Network`` freezes it before
            running the parsers' post-processing hooks).
        """
        if not isinstance(specie, Specie):
            raise ValueError(f"'{specie}' must be an instance of 'Specie'")

        if self.__frozen:
            raise RuntimeError(
                "Species catalogue is frozen: Species cannot be further added to the network"
            )

        if specie.name not in self._by_prop:
            self._by_prop[specie.name] = specie
            self._by_serialized[specie.serialized] = specie
            self._list.append(specie)
            self.count = len(self._list)
            # Invalidate the cached core/special sub-catalogues on mutation.
            self.__dict__.pop("core", None)
            self.__dict__.pop("special", None)

    @cached_property
    def special(self) -> "Species":
        """Sub-catalogue of the special pseudo-species (``is_special``).

        Returns
        -------
        Species
        """
        return Species([s for s in self._list if s.is_special], check_length=False)

    @cached_property
    def core(self) -> "Species":
        """Sub-catalogue of the core (real) species (``is_core``).

        Used wherever only physically integrated species should participate —
        e.g. the mass-action density product and the ODE assembly — so the
        special pseudo-species are excluded from the kinetics.

        Returns
        -------
        Species
        """
        return Species([s for s in self._list if s.is_core], check_length=False)

    def from_serialized(self, serialized: str) -> Specie:
        """Return the species matching the given serialized form.

        Parameters
        ----------
        serialized : str
            Canonical serialized string (e.g. ``"+/H/H/O"``).

        Returns
        -------
        Specie

        Raises
        ------
        KeyError
            If no species has this serialized form.
        """
        return self._by_serialized[serialized]

    def from_name(self, name: str) -> Specie:
        """Return the species with the given chemical name.

        Parameters
        ----------
        name : str
            Species name (e.g. ``"H2O+"``).

        Returns
        -------
        Specie

        Raises
        ------
        KeyError
            If no species has this name.
        """
        return self._by_prop[name]

    def get_list(self) -> list[Specie]:
        """Return the underlying ordered list of ``Specie`` objects.

        The list itself is returned, not a copy, so mutating it mutates the
        catalogue.

        Returns
        -------
        list[Specie]
        """
        return self._list

    def names(self, ne: bool = False) -> Vector[str]:
        """Return a ``Vector`` of species names.

        Names are returned in catalogue (insertion) order.  ``ne=True`` is
        useful when building solver state vectors that treat electrons
        separately.

        Parameters
        ----------
        ne : bool, optional
            If ``True``, exclude the electron species (``"e-"``),
            by default ``False``.

        Returns
        -------
        Vector[str]
        """
        return Vector([s.name for s in self if not (ne and s.name == "e-")])

    def masses(self, ne: bool = False) -> Vector[float | None]:
        """Return a ``Vector`` of species masses in grams (CGS).

        Each mass is the sum over the species' constituent atoms (see
        :attr:`Specie.mass`).

        Parameters
        ----------
        ne : bool, optional
            If ``True``, exclude the electron species, by default ``False``.

        Returns
        -------
        Vector[float | None]
        """
        return Vector([s.mass for s in self if not (ne and s.name == "e-")])

    def exploded(self, ne: bool = False) -> Vector[list[str]]:
        """Return a ``Vector`` of atom lists (one per species).

        Each entry is the sorted atomic-symbol list with repetition (e.g.
        ``["H", "H", "O"]`` for H2O); charge tokens (``"+"``, ``"-"``) are
        included.

        Parameters
        ----------
        ne : bool, optional
            If ``True``, exclude the electron species, by default ``False``.

        Returns
        -------
        Vector[list[str]]
        """
        return Vector([s.exploded for s in self if not (ne and s.name == "e-")])

    def latex(self, dollars: bool = True, ne: bool = False) -> Vector[str]:
        """Return a ``Vector`` of LaTeX species strings.

        See :meth:`Specie.latex` for the formatting rules.  Unlike
        :meth:`Specie.latex`, *dollars* defaults to ``True`` here so the
        strings render inline.

        Parameters
        ----------
        dollars : bool, optional
            Wrap each string in ``$...$`` math delimiters, by default ``True``.
        ne : bool, optional
            If ``True``, exclude the electron species, by default ``False``.

        Returns
        -------
        Vector[str]
        """
        return Vector([s.latex(dollars) for s in self if not (ne and s.name == "e-")])

    def charges(self, ne: bool = False) -> Vector[int]:
        """Return a ``Vector`` of net charges in units of elementary charge.

        The charge is derived by counting trailing ``"+"`` and ``"-"``
        characters in the species name; ``"e-"`` is always ``-1``.

        Parameters
        ----------
        ne : bool, optional
            If ``True``, exclude the electron species, by default ``False``.

        Returns
        -------
        Vector[int]
        """
        return Vector([s.charge for s in self if not (ne and s.name == "e-")])

    def serialized(self, ne: bool = False) -> Vector[str]:
        """Return a ``Vector`` of canonical serialized species strings.

        Each entry is ``"/".join(sorted(exploded))``.  Isomers share the same
        serialized key, making these strings suitable for isomer-insensitive
        comparison.

        Parameters
        ----------
        ne : bool, optional
            If ``True``, exclude the electron species, by default ``False``.

        Returns
        -------
        Vector[str]
        """
        return Vector([s.serialized for s in self if not (ne and s.name == "e-")])

    def elements(self, ne: bool = False) -> Vector[Elements]:
        """Return a ``Vector`` of ``Elements`` collections (one per species).

        Each ``Elements`` instance lists the unique chemical elements present
        in the corresponding species and is used to build composition
        matrices.

        Parameters
        ----------
        ne : bool, optional
            If ``True``, exclude the electron species, by default ``False``.

        Returns
        -------
        Vector[Elements]
        """
        return Vector([s.elements for s in self if not (ne and s.name == "e-")])

    def e_idx(self) -> int | None:
        """Return the integer index of the electron species, or ``None``.

        Needed when building solver state vectors that place the electron
        density at a fixed position.  The value is the electron's
        :attr:`Specie.index`.

        Returns
        -------
        int or None
            ``None`` when ``"e-"`` is not in this catalogue.

        Examples
        --------
        >>> Species(["H2", "H+", "e-"]).e_idx()
        2
        >>> Species(["H2"]).e_idx() is None
        True
        """
        if "e-" in self:
            return self["e-"].index

    def normalized_names(self, pos: str = "j", neg: str = "k") -> Vector[str]:
        """Return species names normalized for use as code identifiers.

        Each name is stripped and lower-cased; ``"+"`` is replaced with *pos*
        and ``"-"`` with *neg*.  The defaults ``j``/``k`` mirror the scheme in
        :meth:`Specie.get_fidx` (``.strip().lower()`` with the same
        replacements), so the two agree for ordinary species -- though only
        partially: ``get_fidx`` special-cases the electron to ``"idx_e"``,
        whereas here ``e-`` normalizes to ``"ek"``.

        Choosing ``j``/``k`` avoids the charge-marker-versus-element-letter
        collision that ``p``/``n`` suffered (e.g. neutral tin ``Sn`` versus the
        sulfur anion ``S-`` -> ``sk``).  It is *not*, however, injective:
        because names are lower-cased, case-distinct species collapse onto the
        same key (``CO``/``Co``, ``CS``/``Cs``, ``NO``/``No`` -> ``co``/``cs``/
        ``no``).  See :meth:`charge_reverse_map`, which raises on such a clash.

        Parameters
        ----------
        pos : str, optional
            Replacement for ``"+"``, by default ``"j"``.
        neg : str, optional
            Replacement for ``"-"``, by default ``"k"``.

        Returns
        -------
        Vector[str]

        Examples
        --------
        >>> Species(["HCO+", "e-"]).normalized_names()
        ['hcoj', 'ek']
        """
        return Vector(
            [s.name.replace("+", pos).replace("-", neg).strip().lower() for s in self]
        )

    def charge_reverse_map(self) -> dict[str, Specie]:
        """Map each species' j/k-normalized identifier back to the ``Specie``.

        Built from :meth:`normalized_names` with the default ``j``/``k`` scheme.
        The mapping avoids charge-marker/element-letter collisions but is not
        injective in general: because :meth:`normalized_names` lower-cases,
        case-distinct species (e.g. ``CO`` and ``Co``) collapse onto the same
        key.  Any such duplicate key raises :exc:`ValueError`.

        Returns
        -------
        dict[str, Specie]

        Raises
        ------
        ValueError
            If two species normalize to the same identifier.
        """
        out: dict[str, Specie] = {}
        for name, sp in zip(self.normalized_names(), self):
            if name in out:
                raise ValueError(
                    f"Charge-normalized name collision: {name!r} maps to both "
                    f"{out[name].name!r} and {sp.name!r}"
                )
            out[name] = sp

        return out

    def neutral(self, attr: str = "") -> Vector[Specie | int]:
        """Return neutral (charge == 0) species or one of their attributes.

        Parameters
        ----------
        attr : str, optional
            If given, return the named attribute of each neutral species
            instead of the ``Specie`` object itself.  Must be one of the
            values in ``Specie._ATTRS``: ``charge``, ``elements``,
            ``exploded``, ``fidx``, ``index``, ``mass``, ``name``,
            ``serialized``.

        Returns
        -------
        Vector[Specie | int]

        Raises
        ------
        ValueError
            If *attr* is not a valid ``Specie`` attribute name.
        """
        if attr:
            if attr not in Specie._ATTRS:
                raise ValueError(f"Invalid attribute passed: {attr}")

            return Vector([getattr(s, attr) for s in self if s.charge == 0])

        return Vector([s for s in self if s.charge == 0])

    def charged(self, attr: str = "", ne: bool = False) -> Vector[Specie]:
        """Return charged (charge != 0) species or one of their attributes.

        Parameters
        ----------
        attr : str, optional
            If given, return the named attribute of each charged species.
            Must be one of the values in ``Specie._ATTRS`` (see
            :meth:`neutral`).
        ne : bool, optional
            If ``True``, exclude the electron species (``"e-"``),
            by default ``False``.

        Returns
        -------
        Vector[Specie]

        Raises
        ------
        ValueError
            If *attr* is not a valid ``Specie`` attribute name.

        Examples
        --------
        >>> species = Species(["H2", "H+", "e-"])
        >>> species.charged("name")
        ['H+', 'e-']
        >>> species.charged("name", ne=True)
        ['H+']
        """
        if attr:
            if attr not in Specie._ATTRS:
                raise ValueError(f"Invalid attribute passed: {attr}")

            return Vector(
                [
                    getattr(s, attr)
                    for s in self
                    if s.charge != 0 and not (ne and s.name == "e-")
                ]
            )

        return Vector([s for s in self if s.charge != 0 and not (ne and s.name == "e-")])

    def charge_truths(self, ne: bool = False) -> Vector[int]:
        """Return a binary ``Vector`` indicating whether each species is charged.

        Parameters
        ----------
        ne : bool, optional
            If ``True``, exclude the electron species, by default ``False``.

        Returns
        -------
        Vector[int]
            ``1`` for charged species, ``0`` for neutral.

        Examples
        --------
        >>> species = Species(["H2", "H+", "e-"])
        >>> species.charge_truths()
        [0, 1, 1]
        >>> species.charge_truths(ne=True)
        [0, 1]
        """
        return Vector([int(bool(s.charge)) for s in self if not (ne and s.name == "e-")])
