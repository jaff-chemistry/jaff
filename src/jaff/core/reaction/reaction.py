# ABOUTME: Reaction class: one reaction's species, rate segments and photo props
# ABOUTME: Serialization, conservation checks, flux strings and rate plotting

"""Single astrochemical :class:`Reaction` for JAFF chemical networks.

A ``Reaction`` holds the symbolic rate expression, reactants, products, and
optional temperature bounds for a single astrochemical reaction.

Serialized form
---------------
The serialized form of a reaction is::

    "<sorted_reactant_names>__<sorted_product_names>"

where species names are sorted alphabetically and joined with ``"."``, and
the two sides are separated by ``"__"``.  For example
``H + H2O+ -> H2O + H+`` serializes as ``"H.H2O+__H+.H2O"``.  This canonical
form is used for hashing and duplicate detection.  Equality compares the
``(serialized, type)`` pair, so the same reaction under two mechanisms is two
distinct reactions; the hash uses ``serialized`` alone.

The ``"."`` species joiner (rather than ``"_"``) is required because special
pseudo-species names start with an underscore (e.g. ``_PHOTON``, ``_GRAIN``)
and underscore-suffixed grain/ice species exist (``X_DUST``); a ``"_"`` joiner
would collide with those and make the form ambiguous.

Reaction types
--------------
The reaction type is concluded by the network-format parser and passed to the
``Reaction`` constructor; the ``type`` attribute holds that stored value (it no
longer inspects the rate expression).  The generic gas-phase types are:

- ``"photo"``       — radiation-driven (photodissociation/ionisation)
- ``"cosmic_ray"``  — cosmic-ray driven
- ``"3_body"``      — three-body reaction
- ``"unknown"``     — unclassified

KIDA networks map their ``itype`` column to ``"cosmic_ray"``, ``"photo"``,
``"two_body"``, ``"radiative_association"``, ``"associative_detachment"``,
``"recombination"``, ``"grain"`` and ``"special"``.

Grain/ice networks (e.g. UCLCHEM) additionally classify by surface mechanism:
``"freeze"``, ``"desorption_thermal"``, ``"desorption_cr"``,
``"desorption_uvcr"``, ``"desorption_h2"``, ``"eley_rideal"``,
``"eley_rideal_desorption"``, ``"langmuir_hinshelwood"``,
``"langmuir_hinshelwood_desorption"``, ``"h2_formation"``, ``"bulk_swap"``,
``"surface_swap"``.  The set is open-ended: any string a parser supplies is
stored verbatim.  See the UCLCHEM parser README for that format's mapping.
"""

from __future__ import annotations

import sys
from functools import cache, cached_property
from typing import TYPE_CHECKING, Any, Mapping, cast

import numpy as np
from sympy import (
    Basic,
    Expr,
    Float,
    sympify,
)

from ...errors import ParserError
from ...io import JaffLogger
from ..elements import Elements
from ..species import Specie, Species
from .config import PhotoConfig, ReactionConfig
from .photo_reaction import PhotoReaction
from .rate import Rate
from .reaction_radiation import ReactionRadiation
from .types import RateSegment, RateSegments, ReactionSource

if TYPE_CHECKING:
    import matplotlib.pyplot as plt

    from ...core import Network


class Reaction:
    """A single astrochemical reaction with a symbolic rate expression.

    Holds the reaction's reactants, products, symbolic rate coefficient, valid
    temperature bounds, energy change and radiation terms.  Reactions are
    normally built by :class:`~jaff.core.Network` while it parses a network
    file and are reached through ``net.reactions`` (a :class:`Reactions`
    catalogue).  Its settings (``t_cutoff``, photo-reaction options,
    shielding) are :attr:`config`: the global ``[network.reactions]`` table
    with any ``[network.reactions.reaction."<serialized>"]`` entry applied over
    it (see :meth:`~jaff.core.reaction.config.ReactionsConfig.resolve`).

    Attributes
    ----------
    net : Network
        Network the reaction belongs to.  Supplies the ``.jfunc`` auxiliary
        functions, the radiation field and the cross-section lookup.
    source : ReactionSource
        The network-file line that produced this reaction: raw ``line``, its
        line number ``nline``, ``index`` and the parser ``format``.  A
        reaction merged from several lines keeps each line on its rate
        segments (``rate_segments[i].source``); this is the first.
    reactants : Species
        Ordered species catalogue of reactant species.
    products : Species
        Ordered species catalogue of product species.
    rate_coefficient : Expr
        SymPy expression for the reaction rate coefficient (units depend on
        the reaction order; typically cm³ s⁻¹ for two-body reactions).
        Read-only: the photo band rate when a radiation field sets it
        (``photo.band_rate``), else the :attr:`rate_segments` collapsed into
        one expression.
    rate : Rate
        The reaction rate per unit volume, ``rate_coefficient`` times the core
        reactants' number densities, used by every equation builder:
        ``rate.expanded`` writes the coefficient out, ``rate.reduced`` uses its
        symbol ``krate[r]``.
    tmin : float or None
        Minimum gas temperature at which the rate is valid (Kelvin), from the
        lowest rate segment.  ``None`` means no lower bound.  Read-only.
    tmax : float or None
        Maximum gas temperature at which the rate is valid (Kelvin), from the
        highest rate segment.  ``None`` means no upper bound.  Read-only.
    config : ReactionConfig
        The reaction's resolved settings, every key set: ``t_cutoff``
        (out-of-range behaviour of the rate segments, see
        :class:`~jaff.core.reaction.types.RateSegments`), ``duplicate_policy``
        and ``photo`` (``pi_database``, ``override_custom_rate``,
        ``override_custom_cross_section``, ``use_proxy_photoreaction``,
        ``shielding``).  Read-only.
    rate_segments : RateSegments
        Per-temperature-range rate pieces.  A single-range reaction holds one
        segment; reactions merged across disjoint ranges accumulate several,
        collapsed into :attr:`rate_coefficient` (``rate_segments.equivalent_rate``).
        Change the rates through it (``map_rates``, ``Reactions.add_rate``).
    dE : Basic
        SymPy expression for the energy released per reaction event (erg),
        from a ``.jfunc`` ``deltaE<N>`` function if there is one; zero
        otherwise.
    verbatim : str
        Human-readable string ``"R1 + R2 -> P1 + P2"``.
    index : int
        File-side source number (``source_index``) of the reaction, i.e. its
        file-order row number.  Stable per parsed row and gapped when rows
        merge (``duplicate_policy`` or a temperature-range merge), so it is
        NOT a valid array subscript.  Used for ``.jfunc`` lookups (``chemRateCoeffN``,
        ``deltaEN``, ...; see :meth:`function`) and
        :meth:`Reactions.by_source_index`.
    catalogue_index : int
        Dense 0..count-1 position of this reaction in the parent ``Reactions``
        catalogue, assigned by :meth:`Reactions.add`.  Safe to use as an array
        subscript in generated code (what the ``reaction_idx`` template
        property emits).  ``-1`` until catalogued.
    type : str
        Reaction type concluded by the parser and stored verbatim.  Generic
        gas-phase values are ``"photo"``, ``"cosmic_ray"``, ``"3_body"``,
        ``"unknown"``; KIDA networks add ``"two_body"``,
        ``"radiative_association"``, ``"associative_detachment"``,
        ``"recombination"``, ``"grain"`` and ``"special"``; grain/ice networks
        (e.g. UCLCHEM) add surface-mechanism
        types: ``"freeze"``, ``"desorption_thermal"``, ``"desorption_cr"``,
        ``"desorption_uvcr"``, ``"desorption_h2"``, ``"eley_rideal"``,
        ``"eley_rideal_desorption"``, ``"langmuir_hinshelwood"``,
        ``"langmuir_hinshelwood_desorption"``, ``"h2_formation"``,
        ``"bulk_swap"``, ``"surface_swap"``.  The set is open-ended.
    serialized : str
        Canonical name-level (isomer-sensitive) form
        ``"<sorted_reactants>__<sorted_products>"``; see :meth:`serialize`.
    serialized_exploded : str
        Like ``serialized`` but built from the atom-level serialized forms of
        each species (isomer-insensitive comparison); see
        :meth:`serialize_exploded`.
    photo : PhotoReaction or None
        The photo-reaction component, ``None`` for non-photo reactions: its
        resolved photo settings (``photo.config``, the same object as
        ``config.photo``), cross sections (``photo.xsec``), per-band data
        (``photo.bands``, ``photo.band_xsecs``), ``custom_rad_rate`` and
        ``shielding_tables``.  See :class:`PhotoReaction`.
    rad : ReactionRadiation
        The reaction's radiation component: the radiation it emits
        (``rad.spectral_yield``, its ``.jfunc`` ``spectralYield<N>``
        function, erg eV⁻¹ per event as a function of ``ph_E``, or ``None``)
        and its per-group integrals (``rad.group_yields``, set when a
        radiation field is configured).  Any reaction type can emit, not only
        photo reactions.  See :class:`ReactionRadiation`.

    Notes
    -----
    For a photo reaction the constructor also looks up its cross sections
    (from a ``.jfunc`` ``photoCrossSection<N>`` function or the
    cross-section databases) and, when the network has a radiation field,
    replaces :attr:`rate_coefficient` with the sum of its per-band rate
    coefficients.  A photo reaction whose rate comes from a ``.jfunc``
    ``chemRateCoeff<N>`` keeps that rate (``photo.custom_rad_rate`` is ``True``)
    and is not coupled to the bands, unless ``override_custom_rate`` is set.

    Two reactions compare equal when both their :attr:`serialized` form and
    :attr:`type` match; ordering and hashing use :attr:`serialized` alone.

    Examples
    --------
    >>> from jaff import Network
    >>> net = Network("h_photoionization")
    >>> rxn = net.reactions["H+.e-__H"]
    >>> rxn.verbatim
    'H+ + e- -> H'
    >>> rxn.type, rxn.index, rxn.catalogue_index
    ('unknown', 1, 1)
    >>> rxn.rate_coefficient
    2.63e-13/(tgas/10000.0)**0.7
    >>> net.reactions[0].type
    'photo'
    """

    def __init__(
        self,
        net: "Network",
        reactants: list[Specie],
        products: list[Specie],
        rate_coefficient: Expr,
        tmin: float | None,
        tmax: float | None,
        source: ReactionSource,
        dE: Expr | None = None,
        type: str = "unknown",
        metadata: dict | None = None,
    ):
        """Construct a ``Reaction`` from its species and rate.

        Parameters
        ----------
        net : Network
            Network the reaction belongs to; supplies the auxiliary functions,
            radiation field and cross-section lookup.
        reactants : list[Specie]
            Reactant species (may contain duplicates for three-body reactions).
        products : list[Specie]
            Product species (may contain duplicates).
        rate_coefficient : Basic
            SymPy expression for the rate coefficient.
        tmin : float or None
            Lower temperature bound in Kelvin, or ``None``.
        tmax : float or None
            Upper temperature bound in Kelvin, or ``None``.
        source : ReactionSource
            The network-file line that produced this reaction.  Its ``index``
            becomes the reaction's ``index``.
        dE : Expr or None, optional
            Energy change per event (erg).  A ``.jfunc`` ``deltaE<N>``
            function takes precedence; ``None`` (default) without one means
            zero.
        type : str, optional
            Reaction type as concluded by the network-format parser (e.g.
            ``"photo"``, ``"cosmic_ray"``, ``"3_body"``, ``"unknown"``, or a
            grain surface-mechanism type).  Stored verbatim on the ``type``
            attribute; defaults to ``"unknown"``.
        metadata : dict or None, optional
            Parser-private data attached at parse time; ``None`` (default)
            starts empty.

        Raises
        ------
        ParserError
            For invalid photo-reaction or shielding configuration (e.g. a
            ``pi_database`` set on a reaction whose cross section comes from
            a ``photoCrossSection<N>`` function).
        """
        self.logger = JaffLogger().get_logger()
        self.net: "Network" = net
        self._metadata: dict = metadata if metadata is not None else {}
        self.source: ReactionSource = source
        self.index: int = source.index
        self.__check_functions()
        self.reactants: Species = Species(reactants, check_length=False)
        self.products: Species = Species(products, check_length=False)
        self.serialized: str = self.serialize(self.reactants, self.products)
        self.config: ReactionConfig = self.net.reactions.config.resolve(self.serialized)
        self.rate_segments: RateSegments = RateSegments(
            self, [RateSegment(rate_coefficient, tmin, tmax, source)]
        )
        function_dE = self.function("deltaE")
        self.dE: Expr = (
            function_dE
            if function_dE is not None
            else (dE if dE is not None else Float(0.0))
        )
        self.rad: ReactionRadiation = ReactionRadiation(self)
        # verbatim is the "R1 + R2 -> P1 + P2" form built from the species
        # names; source.line keeps the raw source line
        self.verbatim: str = self.get_verbatim()
        self.catalogue_index: int = -1

        self.serialized_exploded: str = self.serialize_exploded()
        # The reaction type is concluded by the parser and supplied here, not
        # inferred from the rate expression.
        self.type: str = type

        self.photo: PhotoReaction | None = (
            PhotoReaction(self, self.__photo_config()) if type == "photo" else None
        )
        if self.photo is not None:
            self.photo.compute()

        self.rad.compute()

    def __repr__(self):
        """Return detailed string representation of this reaction.

        Returns
        -------
        str
            String of the form ``"ReactionObject(<verbatim>)"``.
        """
        return f"ReactionObject({self.verbatim})"

    def __str__(self):
        """Return the human-readable verbatim reaction string.

        Returns
        -------
        str
            Verbatim form ``"R1 + R2 -> P1 + P2"``.
        """
        return self.verbatim

    def __eq__(self, other):
        """Check equality by comparing serialized (name-level) forms and types.

        Two reactions are equal when both :attr:`serialized` and :attr:`type`
        match.

        Parameters
        ----------
        other : Reaction
            Reaction to compare against.

        Returns
        -------
        bool

        Raises
        ------
        TypeError
            If *other* is not a ``Reaction`` instance.
        """
        if not isinstance(other, Reaction):
            raise TypeError(
                f"'==' not supported between instances of 'Reaction' and '{other}'"
            )

        return self.serialized == other.serialized and self.type == other.type

    def __hash__(self):
        """Return hash based on the serialized (name-level) form.

        The type is not hashed, so reactions that differ only in type share a
        hash (but compare unequal).

        Returns
        -------
        int
        """
        return hash(self.serialized)

    def __lt__(self, other):
        """Compare reactions lexicographically by serialized form.

        Parameters
        ----------
        other : Reaction
            Reaction to compare against.

        Returns
        -------
        bool

        Raises
        ------
        TypeError
            If *other* is not a ``Reaction`` instance.
        """
        if not isinstance(other, Reaction):
            raise TypeError(
                f"'<' not supported between instances of 'Reaction' and '{other}'"
            )

        return self.serialized < other.serialized

    def function(self, kind: str) -> Expr | None:
        """The reaction's ``.jfunc`` ``<kind><N>`` function, read as a network input.

        ``<N>`` is :attr:`index`; e.g. ``function("deltaE")`` is the reaction's
        ``deltaE<N>``, with its parameters read as network quantities (see
        :meth:`~jaff.core.network.NetworkFunctions.value`).

        Parameters
        ----------
        kind : str
            The function kind (case-insensitive): ``"chemRateCoeff"``, ``"deltaE"``,
            ``"photoCrossSection"``, ``"spectralYield"``, ...

        Returns
        -------
        sympy.Expr or None
            The standardized function, or ``None`` when the ``.jfunc`` file
            defines none for this reaction.

        Examples
        --------
        >>> net = Network("GOW")
        >>> rxn = net.reactions.by_source_index(0)
        >>> rxn.verbatim
        'H + _CR -> H+ + e-'
        >>> rxn.function("deltaE") is None, rxn.function("photoCrossSection") is None
        (False, True)
        """
        name = f"{kind}{self.index}"
        functions = self.net.functions

        return functions.value(name) if name in functions else None

    def expressions(self) -> list[Expr]:
        """Every expression the reaction and its components hold.

        Its rate segments, :attr:`dE`, then those of :attr:`rad` and
        :attr:`photo` (see their ``expressions``).  :meth:`substitute` changes
        exactly these.

        Returns
        -------
        list[Expr]
            The expressions.
        """
        exprs = [seg.rate for seg in self.rate_segments] + [self.dE]
        exprs += self.rad.expressions()
        if self.photo is not None:
            exprs += self.photo.expressions()

        return exprs

    def substitute(self, mapping: Mapping[Basic, Basic]) -> None:
        """Apply *mapping* (``xreplace``) to every expression in :meth:`expressions`.

        Parameters
        ----------
        mapping : Mapping[Basic, Basic]
            Old symbol or expression to its replacement.
        """
        self.rate_segments.map_rates(lambda rate: cast(Expr, rate.xreplace(mapping)))
        self.dE = cast(Expr, self.dE.xreplace(mapping))
        self.rad.substitute(mapping)
        if self.photo is not None:
            self.photo.substitute(mapping)

    def __check_functions(self) -> None:
        """Reject a reaction with both a custom rate and a custom cross section."""
        functions = self.net.functions
        if (
            f"chemratecoeff{self.index}" in functions
            and f"photocrosssection{self.index}" in functions
        ):
            raise ParserError(
                "A reaction cannot have both a chemRateCoeff and a photoCrossSection "
                "auxiliary function"
            )

    def __photo_config(self) -> PhotoConfig:
        """The resolved photo settings (always set once resolved)."""
        photo = self.config.photo
        assert photo is not None, f"{self}: unresolved config"

        return photo

    @property
    def rate_coefficient(self) -> Expr:
        """SymPy rate coefficient; see the class attributes.

        The photo band rate when a radiation field sets it, else the rate
        segments collapsed into one expression.  Units depend on the reaction
        order (e.g. cm³ s⁻¹ for a two-body reaction).  Change it through
        :attr:`rate_segments`.

        Returns
        -------
        Expr
        """
        if self.photo is not None and self.photo.band_rate is not None:
            return self.photo.band_rate

        return self.rate_segments.equivalent_rate

    @property
    def tmin(self) -> float | None:
        """Lower temperature bound of the lowest rate segment (K); see the attributes."""
        return self.rate_segments.sort()[0].tmin

    @property
    def tmax(self) -> float | None:
        """Upper temperature bound of the highest rate segment (K); see the attributes."""
        return self.rate_segments.sort()[-1].tmax

    @cached_property
    def rate(self) -> Rate:
        """Reaction rate per unit volume, ``rate_coefficient · Π nden[reactants]``.

        A :class:`~jaff.core.reaction.Rate`: ``rate.expanded`` writes the rate
        coefficient out, ``rate.reduced`` uses its symbol ``krate[r]``.  The densities
        are those (``net.symbols.ndens``) of the core reactants; pseudo-species such as
        ``_PHOTON`` or ``_CR`` carry no density.

        Returns
        -------
        Rate

        Examples
        --------
        >>> from jaff import Network
        >>> net = Network("h_photoionization")
        >>> net.reactions["H+.e-__H"].rate.expanded
        1.65941781598291e-10*nden[1]*nden[2]/tgas**0.7
        """
        return Rate(self)

    @cached_property
    def species(self) -> Species:
        """All unique species involved in this reaction (reactants ∪ products).

        Merges the reactant and product catalogues without double-counting a
        species that appears on both sides.  Cached after the first access.

        Returns
        -------
        Species
        """
        return Species(list(set(self.reactants._list) | set(self.products._list)))

    @cached_property
    def elements(self) -> Elements:
        """All elements present across reactants and products.

        Cached after the first access.

        Returns
        -------
        Elements
            Catalogue of every unique chemical element in the reaction's
            reactants and products.
        """
        return Elements(self.reactants._list + self.products._list)

    def is_isomer_version(self, other: "Reaction") -> bool:
        """Check whether *other* is an isomer variant of this reaction.

        Two reactions are considered isomer versions of each other when they
        involve the same set of atoms on each side (same ``serialized_exploded``
        form) but differ in at least one species name.

        Parameters
        ----------
        other : Reaction
            The reaction to compare against.

        Returns
        -------
        bool
            ``True`` if *other* is an isomeric version of this reaction (e.g.
            the same reaction with HOC+ in place of HCO+).
        """
        # Atom-level comparison ignores isomer distinctions (e.g. HCO+ ≡ HOC+).
        is_same_serialized = self.serialized_exploded == other.serialized_exploded

        # Name-level comparison detects the isomer distinction.
        rp1 = sorted([x.name for x in self.reactants._list + self.products._list])
        rp2 = sorted([x.name for x in other.reactants._list + other.products._list])
        has_different_species_names = rp1 != rp2

        return is_same_serialized and has_different_species_names

    def serialize_exploded(self) -> str:
        """Build the atom-level serialized form (isomer-insensitive).

        Each species is replaced by its ``Specie.serialized`` form (e.g.
        H2O+ → ``"+/H/H/O"``), then species tokens are sorted and joined
        with ``"."``.  Reactants and products are separated by ``"__"``.
        Two reactions that differ only in isomers (e.g. HCO+ vs HOC+) share
        this form.

        Returns
        -------
        str
            Atom-level canonical key, e.g. ``"+/H/H/O.H__+/H.H/H/O"`` for
            ``H + H2O+ -> H2O + H+``; also stored as
            :attr:`serialized_exploded`.
        """
        sr = ".".join(sorted([x.serialized for x in self.reactants]))
        sp = ".".join(sorted([x.serialized for x in self.products]))

        return f"{sr}__{sp}"

    @staticmethod
    def serialize(
        reactants: Species | list[Specie],
        products: Species | list[Specie],
    ) -> str:
        """Build the name-level serialized form (isomer-sensitive).

        Species names are sorted alphabetically and joined with ``"."``.
        Reactants and products are separated by ``"__"``.  The ``"."``
        joiner (not ``"_"``) is used because special pseudo-species names
        start with ``_`` (e.g. ``_PHOTON``) and grain/ice species end in
        ``_DUST``-style suffixes.  A static method: the instance's form is
        stored as :attr:`serialized`.

        Parameters
        ----------
        reactants : Species or list[Specie]
            Reactant species.
        products : Species or list[Specie]
            Product species.

        Returns
        -------
        str
            Name-level canonical key, e.g. ``"H.H2O+__H+.H2O"`` for
            ``H + H2O+ -> H2O + H+``.

        Examples
        --------
        >>> net = Network("h_photoionization")
        >>> rxn = net.reactions[0]
        >>> Reaction.serialize(rxn.reactants, rxn.products)
        'H._PHOTON__H+.e-'
        """
        sr = ".".join(sorted([r.name for r in reactants]))
        sp = ".".join(sorted([p.name for p in products]))

        return f"{sr}__{sp}"

    @cache
    def normalized_proxy_reaction_str(self):
        """Serialize the reaction with proxy markers stripped from species names.

        Any reactant or product whose name ends in the proxy marker ``"x"`` has
        that trailing character removed, and the resulting species are
        re-serialized.  The returned string is used as the cross-section lookup
        key for proxy photo-reactions.

        Returns
        -------
        str
            The serialized reaction string with proxy markers removed.
        """
        rs = [Specie(r.name[:-1]) if r.name.endswith("x") else r for r in self.reactants]
        ps = [Specie(p.name[:-1]) if p.name.endswith("x") else p for p in self.products]

        return self.serialize(rs, ps)

    def check(self, strict: bool) -> None:
        """Validate mass and charge conservation for this reaction.

        Compares the summed mass and net charge of the reactants against
        those of the products (see :meth:`check_mass` and
        :meth:`check_charge`).  Violations are logged; with *strict* the
        first one is fatal.  :meth:`~jaff.core.Network.check_conservation` runs
        it for every reaction of a network.

        Parameters
        ----------
        strict : bool
            When ``True``, terminate the process on any conservation failure;
            when ``False``, only emit a warning.

        Raises
        ------
        SystemExit
            If *strict* is ``True`` and mass or charge is not conserved.
        """
        if not self.check_mass():
            message = f"Mass not conserved in reaction {self.index}: {self.verbatim}"
            if strict:
                self.logger.error(message)
                sys.exit(1)

            self.logger.warning(message)

        if not self.check_charge():
            message = f"Charge not conserved in reaction {self.index}: {self.verbatim}"
            if strict:
                self.logger.error(message)
                sys.exit(1)

            self.logger.warning(message)

    def check_mass(self) -> bool:
        """Return ``True`` if mass is conserved within one electron mass.

        The tolerance (9.1093837e-28 g) is chosen to accommodate reactions
        that appear to gain or lose a single electron mass due to ionisation
        (the electron mass is negligible for chemistry purposes).

        Returns
        -------
        bool
            ``True`` if the reactant–product mass difference is below
            9.1093837e-28 g.
        """
        return (
            abs(
                np.sum([r.mass for r in self.reactants])
                - np.sum([p.mass for p in self.products])
            )
            < 9.1093837e-28
        )

    def check_charge(self) -> bool:
        """Return ``True`` if the net charge is conserved.

        Returns
        -------
        bool
            ``True`` if the summed reactant charge equals the summed product
            charge.
        """
        return (
            np.sum([x.charge for x in self.reactants])
            - np.sum([x.charge for x in self.products])
        ) == 0

    def get_verbatim(self) -> str:
        """Return a human-readable reaction string ``"R1 + R2 -> P1 + P2"``.

        Built from the species names (pseudo-species included); the same
        string is stored as :attr:`verbatim`.

        Returns
        -------
        str
            E.g. ``"H + _PHOTON -> H+ + e-"``.
        """
        return (
            f"{' + '.join([x.name for x in self.reactants])}"
            " -> "
            f"{' + '.join([x.name for x in self.products])}"
        )

    def get_latex(self) -> str:
        """Return a LaTeX-formatted reaction equation wrapped in ``$...$``.

        Reactants and products are joined with ``+`` and separated by
        ``\\,\\to\\,``; each species is formatted by its own
        :meth:`Specie.latex`.  Used as the panel title by
        :meth:`plot_rate_coefficient`.

        Returns
        -------
        str
            Inline-math LaTeX string, e.g.
            ``"${\\rm H^{+}} + {\\rm e^{-}}\\,\\to\\,{\\rm H}$"``.
        """
        latex = (
            f"{' + '.join([r.latex() for r in self.reactants])}"
            "\\,\\to\\,"
            f"{' + '.join([x.latex() for x in self.products])}"
        )
        return f"${latex}$"

    def get_flux_expression(
        self,
        idx: int = 0,
        rate_variable: str = "k",
        species_variable: str = "y",
        brackets: str = "[]",
        idx_prefix: str = "",
    ) -> str:
        """Return a source-code string for the reaction flux.

        The flux has the form ``k[idx] * y[idx_R1] * y[idx_R2] * ...``,
        where ``idx_Ri`` is derived from each core reactant's ``fidx``
        attribute.  Special pseudo-species (``_PHOTON``, ``_CR``, ...) are
        excluded — they carry the reaction's identity but not its kinetics.

        Parameters
        ----------
        idx : int, optional
            Index into the rate-coefficient array, by default ``0``.
        rate_variable : str, optional
            Name of the rate array variable, by default ``"k"``.
        species_variable : str, optional
            Name of the species density array variable, by default ``"y"``.
        brackets : str, optional
            Two-character string whose first character is the left bracket and
            second is the right bracket (e.g. ``"[]"`` or ``"()"``),
            by default ``"[]"``.
        idx_prefix : str, optional
            Optional prefix prepended to each species index token, by default
            ``""``.

        Returns
        -------
        str
            E.g. ``"k[0] * y[idx_hj] * y[idx_e]"`` for ``H+ + e- -> H``.

        Raises
        ------
        SystemExit
            If *brackets* is not exactly 2 characters.
        """
        if len(brackets) != 2:
            self.logger.error("Brackets must be a string of length 2, e.g. '[]'")
            sys.exit(1)

        lb, rb = brackets[0], brackets[1]
        flux = f"{rate_variable}{lb}{idx}{rb} * " + " * ".join(
            [
                f"{species_variable}{lb}{idx_prefix + x.fidx}{rb}"
                for x in self.reactants.core
            ]
        )

        return flux

    def has_any_species(self, species: list[Specie | str] | str | Specie) -> bool:
        """Return ``True`` if *any* of *species* appear in reactants or products.

        Parameters
        ----------
        species : list[Specie | str] | str | Specie
            One or more species (objects or names) to test.

        Returns
        -------
        bool
            ``True`` if at least one of *species* is in this reaction.
        """
        sp_list: list[str] = []
        if isinstance(species, Specie):
            sp_list.append(species.name)
        elif isinstance(species, str):
            sp_list.append(species)
        elif isinstance(species, list):
            sp_list = [sp.name if isinstance(sp, Specie) else sp for sp in species]

        return any(
            [x.name in sp_list for x in self.reactants._list + self.products._list]
        )

    def has_reactant(self, species: list[Specie | str] | str | Specie) -> bool:
        """Return ``True`` if *all* of *species* appear in the reactants.

        Parameters
        ----------
        species : list[Specie | str] | str | Specie
            One or more species (objects or names) to test.

        Returns
        -------
        bool
            ``True`` if every one of *species* is among the reactants.
        """
        sp_list: list[str] = []
        if isinstance(species, Specie):
            sp_list.append(species.name)
        elif isinstance(species, str):
            sp_list.append(species)
        elif isinstance(species, list):
            sp_list = [sp.name if isinstance(sp, Specie) else sp for sp in species]

        return all([s in self.reactants for s in sp_list])

    def has_product(self, species: list[Specie | str] | str | Specie) -> bool:
        """Return ``True`` if *all* of *species* appear in the products.

        Parameters
        ----------
        species : list[Specie | str] | str | Specie
            One or more species (objects or names) to test.

        Returns
        -------
        bool
            ``True`` if every one of *species* is among the products.
        """
        sp_list: list[str] = []
        if isinstance(species, Specie):
            sp_list.append(species.name)
        elif isinstance(species, str):
            sp_list.append(species)
        elif isinstance(species, list):
            sp_list = [sp.name if isinstance(sp, Specie) else sp for sp in species]

        return all([s in self.products for s in sp_list])

    def get_code(self, lang="cpp") -> str:
        """Generate source code for the reaction rate expression.

        Parameters
        ----------
        lang : str, optional
            Target programming language, by default ``"cpp"``.
            Supported values: ``"python"``, ``"c"``, ``"cxx"``,
            ``"fortran"``, ``"rust"``, ``"julia"``, ``"r"``, and the aliases
            ``"py"``, ``"cpp"``, ``"c++"``, ``"f90"``, ``"rs"``, ``"jl"``.

        Returns
        -------
        str
            Source code string for the rate expression (of
            :meth:`get_sympy`).  Unknown functions are printed as calls.

        Raises
        ------
        InvalidLanguageError
            If *lang* is not one of the supported language keys.
        """
        from ...codegen import Language

        language = Language(lang)

        return language.code_gen(
            self.get_sympy(), strict=False, allow_unknown_functions=True
        )

    def get_sympy(self) -> Basic:
        """Return the rate as a canonical SymPy expression.

        The :attr:`rate_coefficient` passed through ``sympify``; it can go
        straight to SymPy's code printers or differentiation routines.

        Returns
        -------
        Basic
        """
        return sympify(self.rate_coefficient)

    def plot_rate_coefficient(self, **kwargs: Any) -> tuple[plt.Figure, Any] | None:
        """Plot the rate coefficient as a function of gas temperature.

        Forwards to :func:`jaff.plotting.plot_rates` with this reaction as the
        only item; every ``plot_rates`` keyword is accepted (``grid``, ``tmin``,
        ``tmax``, ``xscale``, ``yscale``, ``labels``, ``style``, ``palette``,
        ``save``, ``show``, ``close`` ...).
        The panel is titled with the LaTeX reaction equation (see
        :meth:`get_latex`).

        Returns
        -------
        tuple[matplotlib.figure.Figure, numpy.ndarray] or None
            ``(fig, axes)``, or ``None`` if the rate cannot be evaluated
            numerically (e.g. a photo reaction, whose rate carries the
            symbolic radiation-density variable).

        Notes
        -----
        The temperature axis spans [``tmin``, ``tmax``] on a log scale.
        When the ``tmin`` or ``tmax`` keyword is not given, the reaction's own
        bound is used, and when that is ``None``, 2.73 K and 1e6 K
        respectively.

        Examples
        --------
        >>> fig, axes = net.reactions["H+.e-__H"].plot_rate_coefficient(
        ...     tmin=10, tmax=1e4, save="rate.pdf", show=False
        ... )
        """
        from ...plotting import plot_rates

        return plot_rates(self, **kwargs)
