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
from typing import TYPE_CHECKING, Any, cast

import numpy as np
from astropy import units as u
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
from ._helper import to_float_or_none
from .reaction_props import (
    AuxFunctions,
    PhotoReactionProps,
    ReactionRadiationProps,
)
from .types import RateSegment, RateSegments

if TYPE_CHECKING:
    import matplotlib.pyplot as plt
    import pandas as pd

    from ...core import Network


class Reaction:
    """A single astrochemical reaction with a symbolic rate expression.

    Holds the reaction's reactants, products, symbolic rate coefficient, valid
    temperature bounds, energy change and radiation terms.  Reactions are
    normally built by :class:`~jaff.core.Network` while it parses a network
    file and are reached through ``net.reactions`` (a :class:`Reactions`
    catalogue).  Per-reaction settings (``t_cutoff``, photo-reaction options,
    shielding) come from the network config: the global
    ``[network.reactions]`` table with any
    ``[network.reactions.reaction."<serialized>"]`` entry applied over it
    (see :meth:`~jaff.core.reaction.config.ReactionsConfig.resolve`).

    Attributes
    ----------
    net : Network
        Network the reaction belongs to.  Supplies the ``.jfunc`` auxiliary
        functions, the radiation field, the cross-section lookup and the
        conservation-check mode (``net.spec.strict``).
    original_string : str
        The raw network-file line that produced this reaction.
    reactants : Species
        Ordered species catalogue of reactant species.
    products : Species
        Ordered species catalogue of product species.
    rate_coefficient : Basic
        SymPy expression for the reaction rate coefficient (units depend on
        the reaction order; typically cm³ s⁻¹ for two-body reactions).
        Reassigning it refreshes :attr:`rate`.
    rate : Expr
        ``rate_coefficient`` times the core reactants' number densities, the
        reaction rate per unit volume used by every equation builder.
    tmin : float or None
        Minimum gas temperature at which the rate is valid (Kelvin).
        ``None`` means no lower bound.
    tmax : float or None
        Maximum gas temperature at which the rate is valid (Kelvin).
        ``None`` means no upper bound.
    t_cutoff : str
        Out-of-range behaviour for the rate segment(s): ``"clip"`` holds the
        boundary value, ``"extrapolate"`` lets the outermost segment extend
        unbounded.  Set from the ``t_cutoff`` reaction config key.  See
        :class:`~jaff.core.reaction.types.RateSegments`.
    rate_segments : RateSegments
        Per-temperature-range rate pieces.  A single-range reaction holds one
        segment; reactions merged across disjoint ranges accumulate several,
        collapsed into :attr:`rate_coefficient` via ``evaluate_equivalent_rate``.
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
        NOT a valid array subscript.  Used for aux lookups (``chemRateN``,
        ``deltaEN``, ...) and :meth:`Reactions.by_source_index`.
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
    photo : PhotoReactionProps or None
        Photo-reaction state, ``None`` for non-photo reactions: ``xsec``
        (cross sections of the single decay channel, an
        :class:`~jaff.physics.XsecProps`), ``pi_database`` (per-reaction
        photoionization database override), ``override_custom_rate`` /
        ``override_custom_cross_section`` (per-reaction overrides of the
        global flags), ``custom_rad_rate`` (rate supplied by a ``.jfunc`` aux
        function), ``rgrouped`` (per-band data; see :attr:`band_xsecs`) and
        ``shielding``; the spectrum-averaged cross section is ``xsec.total``
        (cm²).  See :class:`PhotoReactionProps`.
    rad : ReactionRadiationProps
        Radiation the reaction emits: ``spectral_yield`` (from a ``.jfunc``
        ``spectralYield<N>`` function, erg eV⁻¹ per event as a function of
        ``ph_E``, or ``None``) and its per-group integrals ``group_yields``
        (set when a radiation field is configured).  Any reaction type can
        emit, not only photo reactions.

    Notes
    -----
    For a photo reaction the constructor also looks up its cross sections
    (from a ``.jfunc`` ``photoCrossSection<N>`` function or the
    cross-section databases) and, when the network has a radiation field,
    replaces :attr:`rate_coefficient` with the sum of its per-band rate
    coefficients.  A photo reaction whose rate comes from a ``.jfunc``
    ``chemRate<N>`` keeps that rate (``photo.custom_rad_rate`` is ``True``)
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
        original_string: str,
        index: int,
        dE: Expr | None = None,
        t_cutoff: str = "clip",
        type: str = "unknown",
        photo: PhotoReactionProps | None = None,
        rad: ReactionRadiationProps | None = None,
        format: str | None = None,
        metadata: dict | None = None,
    ):
        """Construct a ``Reaction`` and validate mass/charge conservation.

        Parameters
        ----------
        net : Network
            Network the reaction belongs to; supplies the auxiliary functions,
            radiation field, cross-section lookup and conservation-check mode.
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
        original_string : str
            The raw network-file line that produced this reaction.
        index : int
            File-side source number (``source_index``): the file-order row
            number of this reaction, gapped when rows merge.  Not an array
            subscript; see the ``index`` / ``catalogue_index`` attributes.
        dE : Expr or None, optional
            Energy change per event (erg).  A ``.jfunc`` ``deltaE<N>``
            function takes precedence; ``None`` (default) without one means
            zero.
        t_cutoff : str, optional
            Out-of-range behaviour for the initial rate segment (``"clip"`` by
            default); see the :attr:`t_cutoff` attribute.
        type : str, optional
            Reaction type as concluded by the network-format parser (e.g.
            ``"photo"``, ``"cosmic_ray"``, ``"3_body"``, ``"unknown"``, or a
            grain surface-mechanism type).  Stored verbatim on the ``type``
            attribute; defaults to ``"unknown"``.
        photo : PhotoReactionProps or None, optional
            Photo-reaction state; pass one for photo reactions, ``None``
            (default) otherwise.  :class:`~jaff.core.Network` passes one for
            every photo reaction, filled from the reaction's photo config.
        rad : ReactionRadiationProps or None, optional
            Emitted radiation.  ``None`` (default) builds one from the
            ``.jfunc`` ``spectralYield<N>`` function, if any.
        format : str or None, optional
            Name of the network-format parser that produced the reaction
            (e.g. ``"prizmo"``); ``None`` (default) when unknown.
        metadata : dict or None, optional
            Parser-private data attached at parse time; ``None`` (default)
            starts empty.

        Raises
        ------
        SystemExit
            When ``net.spec.strict`` is set and mass or charge is not conserved.
        ParserError
            For invalid photo-reaction or shielding configuration (e.g. a
            ``pi_database`` set on a reaction whose cross section comes from
            a ``photoCrossSection<N>`` function).
        """
        self.logger = JaffLogger().get_logger()
        self.net: "Network" = net
        self._format: str | None = format
        self._metadata: dict = metadata if metadata is not None else {}
        self.index: int = index
        self.aux_funcs: AuxFunctions = self.__get_aux_funcs()
        self.reactants: Species = Species(reactants, check_length=False)
        self.products: Species = Species(products, check_length=False)
        self.rate_coefficient = rate_coefficient
        self.tmin: float | None = tmin
        self.tmax: float | None = tmax
        self.t_cutoff: str = t_cutoff
        self.rate_segments: RateSegments = RateSegments(
            [RateSegment(rate_coefficient, tmin, tmax)], t_cutoff
        )
        self.dE: Expr = (
            cast(Expr, self.aux_funcs.deltae.expr)
            if self.aux_funcs.deltae is not None
            else (dE if dE is not None else Float(0.0))
        )
        self.photo: PhotoReactionProps | None = photo
        self.rad: ReactionRadiationProps = (
            rad
            if rad is not None
            else ReactionRadiationProps(
                spectral_yield=cast(Expr, self.aux_funcs.spectralyield.expr)
                if self.aux_funcs.spectralyield is not None
                else None
            )
        )
        self.original_string = original_string
        # verbatim is the "R1 + R2 -> P1 + P2" form built from the species
        # names; original_string keeps the raw source line
        self.verbatim: str = self.get_verbatim()
        self.catalogue_index: int = -1

        self.check(self.net.spec.strict)
        self.serialized_exploded: str = self.serialize_exploded()
        self.serialized: str = self.serialize(self.reactants, self.products)
        # The reaction type is concluded by the parser and supplied here, not
        # inferred from the rate expression.
        self.type: str = type

        if self.photo is not None:
            if self.aux_funcs.chemrate is not None:
                self.photo.custom_rad_rate = not self.__overrides_custom_rate()

            self.__set_xsec_props()

            if self.net.radiation is not None:
                self.__set_photo_reaction_props()

        if self.net.radiation is not None and self.rad.spectral_yield is not None:
            self.rad.group_yields = self.net.radiation.integrate_over_groups(
                self.rad.spectral_yield
            )

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

    def __get_aux_funcs(self) -> AuxFunctions:
        af = self.net.spec.aux_funcs

        return AuxFunctions(
            chemrate=af.get(f"chemrate{self.index}", None),
            deltae=af.get(f"deltae{self.index}", None),
            photocrosssection=af.get(f"photocrosssection{self.index}", None),
            spectralyield=af.get(f"spectralyield{self.index}", None),
        )

    def __overrides_custom_rate(self) -> bool:
        """Whether radiation recomputes this photo reaction's ``chemRate`` rate.

        The reaction's own ``override_custom_rate`` wins over the global
        ``network.reactions.photo.override_custom_rate``; without radiation
        there is no cross-section rate to use, so the ``chemRate`` is always
        kept.
        """
        assert self.photo is not None

        rad = self.net.radiation
        if rad is None:
            return False

        if self.photo.override_custom_rate is not None:
            return self.photo.override_custom_rate

        return bool(self.net.config.reactions.photo.override_custom_rate)

    def __overrides_custom_cross_section(self) -> bool:
        """Whether this photo reaction ignores its ``photoCrossSection``.

        The reaction's own ``override_custom_cross_section`` wins over the
        global ``network.reactions.photo.override_custom_cross_section``
        (``False`` without radiation).
        """
        assert self.photo is not None

        if self.photo.override_custom_cross_section is not None:
            return self.photo.override_custom_cross_section

        if self.net.radiation is None:
            return False

        return bool(self.net.config.reactions.photo.override_custom_cross_section)

    def __set_xsec_props(self):
        from ...physics import Photochemistry

        assert self.photo is not None

        if self.net._photochemistry is None:
            self.net._photochemistry = Photochemistry(self.net)

        custom = self.aux_funcs.photocrosssection
        if custom is not None and not self.__overrides_custom_cross_section():
            if self.photo.pi_database is not None:
                raise ParserError(
                    f"pi_database cannot be set for {self}: its cross section comes "
                    "from a photoCrossSection function"
                )

            self.photo.xsec = self.net._photochemistry.get_custom_xsec(
                self, cast(Expr, custom.expr)
            )
            return

        # A .jaff reload has no aux functions but carries the custom cross section.
        if self.photo.xsec is not None and self.photo.xsec.database == "custom":
            return

        # A cross section is mandatory only when radiation sets the rate.
        self.photo.xsec = (
            self.net._photochemistry.get_xsec(
                self,
                required=self.net.radiation is not None
                and not self.photo.custom_rad_rate,
            )
            or self.photo.xsec
        )

    def __set_photo_reaction_props(self):
        assert self.net.radiation is not None
        assert self.photo is not None

        # A custom rate is not driven by the transported field, so the reaction
        # is not coupled to the bands and gets no band data.
        if self.photo.custom_rad_rate:
            self.logger.warning(
                f"Photo reaction [cyan]{self}[/] (index {self.index}) uses the custom "
                f"rate chemRate{self.index}: it is not coupled to the radiation bands "
                "and removes no photons. Supply a photoCrossSection"
                f"{self.index} instead, or set override_custom_rate = true to "
                "compute its rate from the database cross sections."
            )
            return

        k_tot, xsec_tot, grouped = self.net.radiation.get_reaction_rate_props(self)
        self.rate_coefficient = k_tot

        assert self.photo.xsec is not None
        self.photo.xsec.total = xsec_tot
        self.photo.rgrouped = grouped

    @property
    def rate_coefficient(self) -> Expr:
        """SymPy rate coefficient; see the class attributes.

        Units depend on the reaction order (e.g. cm³ s⁻¹ for a two-body
        reaction).  Assigning a new value drops the cached :attr:`rate`.

        Returns
        -------
        Expr
        """
        return self._rate_coefficient

    @rate_coefficient.setter
    def rate_coefficient(self, value: Expr) -> None:
        self._rate_coefficient: Expr = value
        # Drop the cached rate so the next read recomputes it from the new value.
        self.__dict__.pop("rate", None)

    @cached_property
    def rate(self) -> Expr:
        """Reaction rate per unit volume, ``rate_coefficient · Π nden[reactants]``.

        The rate coefficient times the number density (``net.symbols.ndens``)
        of each core reactant; pseudo-species such as ``_PHOTON`` or ``_CR``
        carry no density.  Cached, and recomputed after
        :attr:`rate_coefficient` is reassigned.

        Returns
        -------
        Expr

        Examples
        --------
        >>> from jaff import Network
        >>> net = Network("h_photoionization")
        >>> net.reactions["H+.e-__H"].rate
        1.65941781598291e-10*nden[1]*nden[2]/tgas**0.7
        """
        nden = self.net.symbols.ndens
        rate = self.rate_coefficient
        for s in self.reactants.core:
            rate *= nden[self.net.species[s.name].index]

        return rate

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
        first one is fatal.  The constructor runs it with
        ``net.spec.strict``.

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

    @property
    def band_xsecs(self) -> pd.DataFrame:
        """Band-averaged cross sections for this reaction, one row per band.

        Pairs each per-band entry of ``photo.rgrouped`` (populated when a
        radiation field is configured, e.g. through the ``radiation``
        argument of :class:`~jaff.core.Network` or ``[network.radiation]``)
        with the :class:`~jaff.physics.photo_reactions._radiation.RadiationGroup`
        it was computed for.  A photo reaction whose rate comes from a
        ``.jfunc`` ``chemRate<N>`` is not coupled to the bands and has no band
        data, so its frame is empty.  Intended as the data source for band bar
        plots (the ``bands`` option of :meth:`plot_xsecs`).

        Returns
        -------
        pandas.DataFrame
            One row per radiation band this reaction contributes to, with
            columns:

            - ``lower`` : lower band edge in eV.
            - ``upper`` : upper band edge in eV (``inf`` for an open top band).
            - ``eavg``  : photon-number-weighted band-average energy in eV.
            - ``xsec``  : photon-number-weighted band-average cross section in
              cm² (``NaN`` when it is not a number, e.g. a
              ``photoCrossSection<N>`` that depends on ``tgas``).
            - ``xsec_frac`` : the band-average cross section divided by the
              spectrum-averaged one (``photo.xsec.total``); ``0`` when that is
              zero.  Not normalised to sum to 1 over the bands.

            Empty (with the columns above) when the reaction contributes to no
            band, e.g. no radiation field is configured.

        Notes
        -----
        The rows are ordered by ascending band index: ``photo.rgrouped`` is
        built by iterating ``net.radiation.groups`` in order, so the two are
        paired by position.

        Examples
        --------
        >>> net = Network(
        ...     "h_photoionization",
        ...     radiation={"bands": [1, 13.6, 100, "inf"]},
        ... )
        >>> rxn = net.reactions.photo_reactions()[0]
        >>> rxn.band_xsecs.columns.tolist()
        ['lower', 'upper', 'eavg', 'xsec', 'xsec_frac']
        >>> len(rxn.band_xsecs)
        3
        >>> rxn.plot_xsecs(bands=True)  # overlay the band bars on σ(E)
        """
        import pandas as pd

        columns = ["lower", "upper", "eavg", "xsec", "xsec_frac"]
        grouped = self.photo.rgrouped if self.photo is not None else None
        groups = self.net.radiation.groups if grouped and self.net.radiation else []
        rows = [
            {
                "lower": to_float_or_none(group.lower),
                "upper": to_float_or_none(group.upper),
                "eavg": (
                    None
                    if group.eavg is None
                    else to_float_or_none(group.eavg * u.erg.to(u.eV))
                ),
                "xsec": to_float_or_none(props.xsec),
                "xsec_frac": to_float_or_none(props.xsec_frac),
            }
            for group, props in zip(groups, grouped or [], strict=True)
        ]

        return pd.DataFrame(rows, columns=columns)

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

    def plot_xsecs(self, **kwargs: Any) -> tuple[plt.Figure, Any] | None:
        """Plot photo cross sections against photon energy or wavelength.

        Forwards to :func:`jaff.plotting.plot_xsecs` with this reaction as the
        only item; every ``plot_xsecs`` keyword is accepted (``grid``,
        ``processes``, ``energy_unit``, ``xsec_unit``, ``xscale``, ``yscale``,
        ``bands``, ``style``, ``palette``, ``save``, ``show``, ``close`` ...).
        The panel is titled with the reaction equation and the curves are
        labelled by process.

        Returns
        -------
        tuple[matplotlib.figure.Figure, numpy.ndarray] or None
            ``(fig, axes)``, or ``None`` (with a log message) if
            ``self.photo`` or its ``xsec`` is ``None`` or no requested process has data.

        Examples
        --------
        >>> photo = net.reactions.photo_reactions()[0]
        >>> photo.plot_xsecs()  # all processes, eV vs Mb
        >>> photo.plot_xsecs(processes="photodecay", yscale="linear")
        >>> photo.plot_xsecs(energy_unit="nm", xsec_unit="cm^2", bands=True)
        """
        from ...plotting import plot_xsecs

        return plot_xsecs(self, **kwargs)
