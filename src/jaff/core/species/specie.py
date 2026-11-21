# ABOUTME: Specie class: one chemical species' name, mass, composition and charge
# ABOUTME: Parses formulas into exploded atom lists and serialized identity forms

"""Single chemical species (:class:`Specie`) for JAFF networks.

A *specie* is a single molecular or atomic entity present in a reaction
network.  Each ``Specie`` instance carries its chemical name, mass (CGS),
elemental composition (the ``exploded`` atom list), charge, and a canonical
serialized form used for identity comparisons.

The serialized form sorts atom tokens alphabetically and joins them with
``"/"``.  For example ``H2O+`` becomes ``"+/H/H/O"`` — so structural isomers
(e.g. HCO+ and HOC+) share the same serialized form.

The ``fidx`` attribute encodes the species name as a C/Fortran-compatible
identifier (``"idx_X"`` for most species; ``"idx_e"`` for electrons).
"""

from __future__ import annotations

import re
import sys
from functools import cached_property
from itertools import product
from typing import TYPE_CHECKING, Mapping

from ...chemistry import periodic_table
from ...io import JaffLogger
from ..elements import Elements

if TYPE_CHECKING:
    import logging

    from ...physics.thermodynamics import SpecieThermodynamics
    from ..elements import Element


class Specie:
    """A single chemical species in a JAFF reaction network.

    A ``Specie`` represents an atom, molecule, ion, grain or special
    pseudo-species.  On construction the name is parsed to determine the
    elemental composition, total mass, net charge and LaTeX formatting.  An
    independent ``Specie`` can be created outside of a network; the *index*
    then only matters when it belongs to a :class:`Species` catalogue.

    Attributes
    ----------
    name : str
        Species name as it appears in the network file (e.g. ``"H2O+"``,
        ``"e-"``).
    mass : float or None
        Total mass in grams (CGS), summed over constituent atoms.  Populated
        by the constructor, so it is always a ``float`` on a constructed
        instance.
    exploded : list[str]
        Sorted list of atomic symbols with repetition
        (e.g. ``["H", "H", "O"]`` for H2O).  Charge tokens (``"+"``, ``"-"``)
        are included so that the serialized form encodes charge.
    charge : int
        Net charge in units of elementary charge.  Trailing ``"+"`` or
        ``"-"`` characters at the end of *name* are counted; ``e-`` always
        has charge ``-1``.
    index : int
        Position of this species in the parent ``Species`` catalogue.
    fidx : str
        Flat index identifier, safe for use in generated C/Fortran/Python
        source.  ``"+"`` maps to ``"j"`` and ``"-"`` maps to ``"k"``, so
        H2O+ becomes ``"idx_h2oj"``.  Electrons are special-cased to
        ``"idx_e"``.
    serialized : str
        Canonical form built by sorting *exploded* and joining with ``"/"``.
        Isomers share the same serialized form (e.g. HCO+ and HOC+ both
        become ``"+/C/H/O"``).
    elements : Elements
        Lazy-loaded ``Elements`` collection for this single species.
    thermodynamics : SpecieThermodynamics or None
        Thermodynamic properties under the network's EOS (adiabatic index),
        attached by :class:`~jaff.physics.thermodynamics.Thermodynamics`; ``None``
        for a species outside a network or an EOS without per-species properties.

    Raises
    ------
    KeyError
        If *name* contains a token that is not a known element symbol (or
        stoichiometric number) in the mass dictionary.

    Examples
    --------
    >>> from jaff.core.species import Specie
    >>> sp = Specie("H2O+")
    >>> sp.exploded
    ['+', 'H', 'H', 'O']
    >>> sp.serialized
    '+/H/H/O'
    >>> sp.charge
    1
    >>> sp.fidx
    'idx_h2oj'
    >>> Specie("HCO+").serialized == Specie("HOC+").serialized
    True
    """

    _ATTRS: frozenset[str] = frozenset(
        {"name", "mass", "exploded", "charge", "index", "fidx", "serialized", "elements"}
    )

    def __init__(self, name: str, index: int = 0):
        """Construct a ``Specie`` and immediately parse its composition.

        Parameters
        ----------
        name : str
            Chemical formula as used in the network file.  The electron must
            be given as ``"e-"``; the names ``"E"``, ``"E-"`` and (in any case)
            ``"e"``, ``"eletron"``, ``"electrons"``, ``"el"`` and ``"els"``
            raise a fatal error.  ``"electron"`` itself is not in that list.
        index : int, optional
            Position in the parent ``Species`` catalogue, by default ``0``.

        Raises
        ------
        SystemExit
            If *name* is one of the electron misspellings listed under *name*
            (the correct form is ``"e-"``).
        """
        self.logger: logging.Logger = JaffLogger().get_logger()
        if name.lower() in ["e", "eletron", "electrons", "el", "els"] or name in [
            "E",
            "E-",
        ]:
            self.logger.error(f"Electrons found with name: {name}. Use e- instead")
            sys.exit(1)

        self.name: str = name
        self.mass: float | None = None
        self.exploded: list[str] = []
        self.charge: int = 0
        self.index: int = index
        self.__latex: str = ""
        self.fidx: str = self.get_fidx()
        self.serialized: str = ""
        self.thermodynamics: SpecieThermodynamics | None = None

        self._parse(periodic_table())
        self.serialize()

    def __repr__(self):
        """Return detailed string representation of this species.

        Returns
        -------
        str
            String of the form ``"SpecieObject(<name>)"``.
        """
        return f"SpecieObject({self.name!r})"

    def __str__(self):
        """Return the species chemical name.

        Returns
        -------
        str
            Chemical name as used in the network file.
        """
        return self.name

    def __eq__(self, other):
        """Check equality by comparing species names.

        Parameters
        ----------
        other : Specie
            Species to compare against.

        Returns
        -------
        bool

        Raises
        ------
        TypeError
            If *other* is not a ``Specie`` instance.
        """
        if not isinstance(other, Specie):
            raise TypeError(
                f"'==' not supported between instances of 'Specie' and '{other}'"
            )

        return self.name == other.name

    def __hash__(self):
        """Return hash based on the species name.

        Returns
        -------
        int
        """
        return hash(self.name)

    def __lt__(self, other):
        """Compare species lexicographically by name.

        Parameters
        ----------
        other : Specie
            Species to compare against.

        Returns
        -------
        bool

        Raises
        ------
        TypeError
            If *other* is not a ``Specie`` instance.
        """
        if not isinstance(other, Specie):
            raise TypeError(
                f"'<' not supported between instances of 'Specie' and '{other}'"
            )

        return self.name < other.name

    @cached_property
    def elements(self) -> Elements:
        """Lazy ``Elements`` collection derived from this species' atoms.

        Returns
        -------
        Elements
        """
        return Elements(self)

    def get_fidx(self) -> str:
        """Compute the flat index identifier string for generated code.

        Since ``+`` and ``-`` are not valid in identifiers, they are replaced
        with ``j`` and ``k`` respectively and the result is lower-cased, giving
        a name usable in generated C, Fortran and Python source.  Electrons
        (``e-``) are special-cased to ``"idx_e"``.

        Returns
        -------
        str
            ``"idx_e"`` for electrons; otherwise ``"idx_<name>"`` with
            ``"+"`` replaced by ``"j"`` and ``"-"`` replaced by ``"k"``.

        Examples
        --------
        >>> Specie("HCO+").get_fidx()
        'idx_hcoj'
        >>> Specie("H2O+").get_fidx()
        'idx_h2oj'
        >>> Specie("e-").get_fidx()
        'idx_e'
        """
        return (
            "idx_e"
            if self.name == "e-"
            else f"idx_{self.name.replace('+', 'j').replace('-', 'k').strip().lower()}"
        )

    @property
    def is_special(self) -> bool:
        """Whether this is a special pseudo-species.

        Special pseudo-species (``_PHOTON``, ``_CR``, ``_GRAIN``, ``_DUMMY``,
        ...) are radiation/cosmic-ray/grain agents and markers; they carry the
        reaction's identity but do not participate in the mass-action kinetics
        or the integrated ODE state.  They are identified by a leading
        underscore — real species never start with ``_`` (underscore-suffixed
        grain/ice species such as ``H2O_DUST`` have the underscore mid-name).

        Returns
        -------
        bool
        """
        return self.name.startswith("_")

    @property
    def is_core(self) -> bool:
        """Whether this is a core (real) species, i.e. not :attr:`is_special`.

        Returns
        -------
        bool
        """
        return not self.is_special

    def serialize(self) -> str:
        """Build and store the canonical serialized form of this species.

        The serialized form is ``"/".join(sorted(self.exploded))``, which
        places isomers into the same equivalence class.  It is used for
        lookups and comparisons that must not depend on how the species was
        originally written.  Called automatically by the constructor.

        Returns
        -------
        str
            The serialized string (also stored in ``self.serialized``).

        Examples
        --------
        >>> Specie("H2O").serialize()
        'H/H/O'
        >>> Specie("H2O+").serialize()
        '+/H/H/O'
        """
        self.serialized = "/".join(sorted(self.exploded))

        return self.serialized

    def latex(self, dollars: bool = False):
        """Return a LaTeX representation of this species.

        Parameters
        ----------
        dollars : bool, optional
            When ``True``, wrap the string in ``$...$`` math delimiters,
            by default ``False``.

        Returns
        -------
        str
            LaTeX string with subscript digits, superscript charges, and
            ``\\rm`` roman font for element names.  Suffixes ``_ORTHO``,
            ``_PARA``, ``_META`` become ``o``, ``p``, ``m`` prefixes;
            ``_DUST`` becomes an ``ice`` suffix; ``_GRAIN`` becomes ``g``.
            Repeated charges are collapsed (``Fe++`` gives ``Fe^{2+}``) and
            ``_PHOTON`` is rendered as ``\\gamma``.

        Examples
        --------
        >>> print(Specie("H2O").latex())
        {\\rm H_{2}O}
        >>> print(Specie("H2O+").latex(dollars=True))
        ${\\rm H_{2}O^{+}}$
        >>> print(Specie("H2_ORTHO").latex())
        {\\rm oH_{2}}
        """
        return f"${self.__latex}$" if dollars else self.__latex

    def _parse(self, table: Mapping[str, Element]) -> None:
        """Parse *name* into elemental composition, mass, and charge.

        This is called automatically during ``__init__``.  It uses a proxy
        substitution strategy to avoid ambiguous greedy matches when element
        symbols are substrings of each other (e.g. ``"C"`` inside ``"Ca"``):
        each element symbol is temporarily replaced by a unique 4-character
        proxy drawn from the alphabet ``{q, z, x, j}``, the formula is
        tokenized, then proxies are reversed back to real symbols.  Symbols
        are substituted longest-first so that multi-character symbols (e.g.
        ``"He"``) match before single-character ones (e.g. ``"H"``).

        Parameters
        ----------
        table : Mapping[str, Element]
            The periodic table (:func:`~jaff.chemistry.periodic_table`).

        Notes
        -----
        Numeric stoichiometry coefficients immediately following an element
        token are expanded into repeated atom entries (e.g. ``H2`` → ``["H",
        "H"]``).  The special coefficient ``"x"`` is treated as 1 so that
        wildcard species names are tolerated.

        Charge is determined by counting trailing ``"+"`` and ``"-"``
        characters at the *end* of the name only, which avoids misreading
        embedded signs (e.g. the lone minus in ``"e-"`` which is
        special-cased separately).
        """
        # Sort atoms longest-first so that multi-character symbols (e.g. "He")
        # are matched before single-character ones (e.g. "H").
        atoms = sorted(table.keys(), key=lambda x: len(x), reverse=True)
        ps = ["".join(x) for x in product("qzxj", repeat=4)][: len(atoms)]
        proxy = {a: p for a, p in zip(atoms, ps)}
        proxy_rev = {p: a for a, p in proxy.items()}

        pname = self.name.strip()
        for a in atoms:
            pname = pname.replace(a, f"${proxy[a]}$")

        def is_number(s: str) -> bool:
            if s == "x":
                return True
            try:
                float(s)
                return True
            except ValueError:
                return False

        alist = [x for x in pname.split("$") if x != ""]
        expl = []
        pold = None
        for p in alist:
            if not is_number(p):
                expl += [p]
            else:
                if p != "x":
                    # Repeat the previous atom (n-1) extra times to account
                    # for the first copy already appended in the previous step
                    # (an explicit count of 1, e.g. "C1", adds none).
                    expl += [pold] * (int(p) - 1)
            pold = p
        self.exploded = sorted([proxy_rev[x] for x in expl])
        self.mass = sum([table[x].mass for x in self.exploded])

        # --- Build LaTeX representation ---
        latex = self.name.strip()
        latex = re.sub(r"\d+", lambda m: "_{" + m.group() + "}", latex)
        # Trailing proxy marker "x" (CHx, OHx; see
        # Reaction.normalized_proxy_reaction_str) is a subscript, placed
        # before any charge; after a count it joins that subscript (C2x -> C_{2x}).
        latex = re.sub(r"_\{(\d+)\}x(?=[+-]*$)", r"_{\1x}", latex)
        latex = re.sub(r"(?<=[A-Za-z])x(?=[+-]*$)", "_{x}", latex)
        latex = re.sub(
            r"\++",
            lambda m: f"^{{{len(m.group()) if len(m.group()) > 1 else ''}+}}",
            latex,
        )
        latex = re.sub(
            r"-+",
            lambda m: f"^{{{len(m.group()) if len(m.group()) > 1 else ''}-}}",
            latex,
        )
        if "_ORTHO" in latex:
            latex = "o" + latex.replace("_ORTHO", "")
        if "_PARA" in latex:
            latex = "p" + latex.replace("_PARA", "")
        if "_META" in latex:
            latex = "m" + latex.replace("_META", "")
        if "_DUST" in latex:
            latex = latex.replace("_DUST", "") + "ice"

        latex = latex.replace("_GRAIN", "g")
        # Any remaining underscore (pseudo-species such as _PHOTON, suffixes
        # such as _BULK) is literal text, not a mathtext subscript.
        latex = re.sub(r"_(?!\{)", r"\\_", latex)

        if self.name.strip() == "_PHOTON":
            latex = "\\gamma"

        self.__latex = f"{{\\rm {latex}}}"

        # --- Determine net charge ---
        if self.name == "e-":
            self.charge = -1
            return

        self.charge = 0
        name = self.name
        # Count charge only at the end of the name so that mid-name hyphens
        # (e.g. compound names) are not misidentified as negative charges.
        while name.endswith("+") or name.endswith("-"):
            if name.endswith("+"):
                self.charge += 1
            elif name.endswith("-"):
                self.charge -= 1
            name = name[:-1]
        # self.charge = self.name.count("+") - self.name.count("-")
