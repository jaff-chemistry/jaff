# ABOUTME: Reactions: ordered catalogue of a network's reactions
# ABOUTME: Lookup by verbatim string, serialized form, (serialized, type) or index

"""Ordered, doubly-indexed :class:`Reactions` catalogue for JAFF networks.

Reactions can be looked up by verbatim string
(``reactions["H + H2O+ -> H2 + OH+"]``) or by serialized form
(``reactions["H.H2O+__H2.OH+"]``).
"""

from __future__ import annotations

from functools import cached_property
from typing import TYPE_CHECKING, Any, Sequence, SupportsIndex, cast, overload

from sympy import Basic, Expr, Float

from ...config.base import dotted
from ...errors import ParserError
from ...io import JaffLogger
from ...types import Catalogue, Vector
from ..species import Species
from .config import ReactionsConfig
from .reaction import Reaction
from .types import RateSegment

if TYPE_CHECKING:
    import matplotlib.pyplot as plt

    from ...physics import Photochemistry


class Reactions(Catalogue[Reaction]):
    """Ordered, doubly-indexed catalogue of ``Reaction`` objects.

    Reactions can be looked up by verbatim string (``reactions["H + H2O+ -> H2 + OH+"]``)
    or by serialized form (``reactions["H.H2O+__H2.OH+"]``), and by position
    (``reactions[0]``).  The vector accessors (:meth:`verbatim`,
    :meth:`rate_coefficients`, :meth:`tmins`, ...) return one entry per
    reaction, in catalogue order.  A network's catalogue is ``net.reactions``.

    Both indexes map a key to a *list* of reactions: several reactions may
    share a verbatim string or serialized form when they differ only by
    mechanism/``type`` (e.g. thermal vs cosmic-ray desorption).  String
    indexing is therefore *scalar-or-list*: it returns a single
    :class:`Reaction` when the key is unique, or a ``list[Reaction]`` when
    several mechanisms share it (a missing key raises ``KeyError``).  Use the
    ``(id, type)`` tuple form to always pick exactly one, or
    :meth:`from_serialized` / :meth:`all` for a result that is always a
    collection.

    The catalogue holds the ``[network.reactions]`` settings of its
    reactions (:attr:`config`): each reaction's own settings are
    ``config.resolve(serialized)``, a second rate for a catalogued reaction
    follows its ``duplicate_policy`` (:meth:`add_rate`), and
    :meth:`check_config` checks that every per-reaction entry names one of its
    reactions.

    Attributes
    ----------
    config : ReactionsConfig
        The ``[network.reactions]`` settings (read-only).
    count : int
        Number of reactions in the catalogue.

    Examples
    --------
    >>> from jaff import Network
    >>> net = Network("h_photoionization")
    >>> net.reactions.verbatim()
    ['H + _PHOTON -> H+ + e-', 'H+ + e- -> H']
    >>> net.reactions["H+ + e- -> H"] is net.reactions["H+.e-__H"]
    True
    >>> net.reactions["H+.e-__H", "unknown"].index
    1
    >>> ("H+.e-__H", "photo") in net.reactions
    False
    """

    _by_prop: dict[str, list[Reaction]]  # type: ignore
    _by_serialized: dict[str, list[Reaction]]  # type: ignore

    def __init__(self, config: ReactionsConfig):
        """Initialise an empty reactions catalogue.

        Reactions are appended with :meth:`add`.

        Parameters
        ----------
        config : ReactionsConfig
            The ``[network.reactions]`` settings; :class:`~jaff.core.Network`
            passes ``net.config.reactions``.
        """
        self.config: ReactionsConfig = config
        super().__init__()
        self._by_serialized = {}
        self._by_source_index: dict[int, Reaction] = {}
        # Filter results keyed by (filter, argument); cleared by add().
        self.__filter_cache: dict[tuple[str, ...], Vector[Reaction]] = {}
        self.__frozen: bool = False

    def _freeze(self) -> None:
        """Disallow further :meth:`add` calls.

        ``Network`` freezes its catalogue before running the parsers'
        post-processing hooks, which may change reactions but not add them.
        Logs the :attr:`photochemistry` database fallbacks, if any.
        """
        self.__frozen = True
        if "photochemistry" in self.__dict__:
            self.photochemistry.log_fallback_summary()

    @cached_property
    def photochemistry(self) -> Photochemistry:
        """The cross-section lookup service of the photo reactions.

        Built on first use (it downloads the cross-section and shielding data on
        first use); its database fallbacks are logged once, when the catalogue
        is frozen.

        Returns
        -------
        Photochemistry
        """
        from ...physics import Photochemistry

        return Photochemistry(self.config)

    @overload
    def __getitem__(self, key: str) -> Reaction | list[Reaction]: ...

    @overload
    def __getitem__(self, key: SupportsIndex) -> Reaction: ...

    @overload
    def __getitem__(self, key: tuple[str, str]) -> Reaction: ...

    @overload
    def __getitem__(self, key: slice) -> list[Reaction]: ...

    def __getitem__(self, key):
        """Look up a reaction by serialized/verbatim string, ``(id, type)``, or index.

        Several reactions may share a serialized form or verbatim string when
        they differ only by mechanism/``type``, so string lookup is
        *scalar-or-list*:

        * ``str`` — the reactions with that serialized form (checked first) or
          verbatim string.  Returns the sole :class:`Reaction` when the key is
          unique, or a ``list`` of reactions when several mechanisms share it.
          Raise :exc:`KeyError` when the key is absent.  Use the ``(id, type)``
          form to always pick one, or :meth:`all` / :meth:`from_serialized` for
          an always-collection result.
        * ``tuple`` — ``(id, type)``: return the reaction whose serialized form
          or verbatim string is *id* and whose ``type`` is *type* (exactly one
          must match).
        * ``int`` / ``slice`` — positional access, delegated to the base list.

        Parameters
        ----------
        key : str, tuple[str, str], int, or slice

        Returns
        -------
        Reaction or list[Reaction]
            A single ``Reaction`` for a unique string, tuple, or int key; a
            list for an ambiguous string key or a slice.

        Raises
        ------
        KeyError
            If a string key is missing, or a tuple key does not match exactly
            one reaction.
        """
        if isinstance(key, tuple):
            ident, rtype = key
            group = self._by_serialized.get(ident) or self._by_prop.get(ident, [])
            matches = [r for r in group if r.type == rtype]

            if len(matches) != 1:
                raise KeyError(f"{key!r}: {len(matches)} matches")

            return matches[0]

        if isinstance(key, str):
            # Serialized form first, then verbatim string; both map to a list.
            # Unwrap to a scalar when unique, else hand back the whole group.
            for index in (self._by_serialized, self._by_prop):
                group = index.get(key)
                if group is not None:
                    return group[0] if len(group) == 1 else group

            raise KeyError(f"{key!r} not found in catalogue")

        # int / slice — delegate to base (positional list access)
        return super().__getitem__(key)

    def __repr__(self):
        return f"Catalogue({self.verbatim()!r})"

    def __contains__(self, item) -> bool:
        """Test membership by ``Reaction``, name/serialized string, or
        ``(serialized, type)`` tuple.

        Parameters
        ----------
        item : Reaction, str, or tuple[str, str]
            * ``str`` — checks verbatim names and serialized forms.
            * ``tuple`` — ``(serialized, type)``: true if a reaction with that
              serialized form *and* type exists.
            * ``Reaction`` — true if an equal reaction (same serialized form
              and type; see :meth:`Reaction.__eq__`) is in the catalogue.

        Returns
        -------
        bool
        """
        if isinstance(item, tuple):
            ident, rtype = item
            group = self._by_serialized.get(ident) or self._by_prop.get(ident, [])

            return any(r.type == rtype for r in group)

        if isinstance(item, str):
            return item in self._by_prop or item in self._by_serialized

        return item in self._list

    def add(self, reaction: Reaction) -> None:
        """Append a reaction to the catalogue (duplicates are not checked here).

        Updates the verbatim, serialized and source-number indexes, so the
        reaction can be looked up at once, and sets its ``catalogue_index``
        to its dense ``0..count-1`` position (distinct from
        ``Reaction.index``, the gapped file-side source number).  Duplicates
        over the same temperature range are resolved earlier, during network
        construction, by the ``duplicate_policy`` reaction setting.  Clears
        the cached filter results.

        Parameters
        ----------
        reaction : Reaction
            The reaction to append.

        Raises
        ------
        ValueError
            If *reaction* is not a ``Reaction`` instance.
        RuntimeError
            If the catalogue has been frozen; ``Network`` freezes it once
            parsing is done, before the parsers' post-processing hooks.
        ParserError
            If photo settings are set for *reaction* and it is not a photo
            reaction.
        """
        if not isinstance(reaction, Reaction):
            raise ValueError(f"'{reaction}' must be an instance of 'Reaction'")

        if self.__frozen:
            raise RuntimeError(
                "Reactions catalogue is frozen: Reactions cannot be further added to the network"
            )

        own = self.config.own(reaction.serialized)
        if own is not None and own.photo is not None and reaction.type != "photo":
            keys = ", ".join(sorted(own.photo.model_fields_set))
            raise ParserError(
                f"{keys} is only valid for photo reactions: {reaction.serialized}"
            )

        reaction.catalogue_index = len(self._list)
        self._by_prop.setdefault(reaction.verbatim, []).append(reaction)
        self._by_serialized.setdefault(reaction.serialized, []).append(reaction)
        self._by_source_index[reaction.index] = reaction
        self._list.append(reaction)
        self.count = len(self._list)
        self.__filter_cache.clear()

    def add_rate(self, reaction: Reaction, segment: RateSegment) -> None:
        """Add another rate coefficient to a catalogued reaction.

        A rate over a new temperature range becomes another of the reaction's
        rate segments.  A second rate over a range the reaction already has
        is resolved by the reaction's ``duplicate_policy``:
        ``"preserve-first"`` keeps the existing rate, ``"preserve-last"``
        replaces it and ``"error"`` raises.

        Parameters
        ----------
        reaction : Reaction
            A reaction of this catalogue (same reaction and mechanism).
        segment : RateSegment
            The new rate coefficient and its temperature range.

        Raises
        ------
        ParserError
            If the range is already covered and the policy is ``"error"``.
        """
        span = (segment.tmin, segment.tmax)
        existing = reaction.rate_segments._by_prop.get(span)
        if existing is None:
            reaction.rate_segments.add(segment)
            return

        policy = reaction.config.duplicate_policy
        name = f"{reaction.verbatim} [type={reaction.type}]"
        if policy == "error":
            raise ParserError(
                f"Duplicate reaction: {name} has two rate coefficients over the "
                f"same temperature range {span}:\n"
                f"  existing: {existing.rate}\n"
                f"  new:      {segment.rate}\n"
                "Same reaction, same mechanism and overlapping temperature range "
                "cannot be resolved automatically. Deduplicate the source network."
            )

        kept, dropped = (
            (existing.rate, segment.rate)
            if policy == "preserve-first"
            else (segment.rate, existing.rate)
        )
        JaffLogger().get_logger().warning(
            f"Duplicate rate for [cyan]{reaction.verbatim}[/] [type={reaction.type}] "
            f"over T={span} resolved by policy '{policy}': kept {kept}, "
            f"dropped {dropped}"
        )
        if policy == "preserve-last":
            reaction.rate_segments.add(segment)

    def check_config(self) -> None:
        """Check that every per-reaction config entry names one of the reactions.

        Raises
        ------
        ParserError
            Listing each ``[network.reactions.reaction."<serialized>"]`` entry
            whose reaction is not in the catalogue.
        """
        unknown = [key for key in self.config.reaction if key not in self._by_serialized]
        if not unknown:
            return

        problems = "\n  ".join(
            f"{dotted(('network', 'reactions', 'reaction', key))}: "
            "no such reaction in the network"
            for key in unknown
        )
        raise ParserError(f"Invalid config:\n  {problems}")

    def link(self) -> None:
        """Resolve every ``rc_<N>`` reference to reaction *N*'s final rate coefficient.

        Run once every reaction is built.  ``<N>`` is the reaction's file-side
        number (:meth:`by_source_index`).  Rate coefficients are resolved in
        dependency order, so ``rc_<N>`` may point at a reaction whose own rate
        uses another ``rc_``; every other expression a reaction holds (``dE``,
        radiation group yields, ...; see :meth:`Reaction.expressions`) then
        gets the final rates.

        Raises
        ------
        ParserError
            If an ``rc_`` is not followed by a number, names a reaction that is
            not in the network, or the rate coefficients reference each other in
            a cycle (the message shows it, e.g. ``rc_0 → rc_1 → rc_0``).
        """
        references = {id(r): self._references(r.expressions()) for r in self}
        done: set[int] = set()
        for reaction in self:
            self._link(reaction, [], done)

        for reaction in self:
            refs = references[id(reaction)]
            if refs:
                reaction.substitute({s: t.rate_coefficient for s, t in refs.items()})

    def _link(self, reaction: Reaction, path: list[Reaction], done: set[int]) -> None:
        """Resolve the ``rc_`` references of *reaction*'s rate, depth first."""
        if id(reaction) in done:
            return

        if any(r is reaction for r in path):
            start = next(i for i, r in enumerate(path) if r is reaction)
            cycle = " → ".join(f"rc_{r.index}" for r in [*path[start:], reaction])
            raise ParserError(
                f"Rate coefficients reference each other in a cycle: {cycle}"
            )

        refs = self._references([reaction.rate_coefficient])
        for target in refs.values():
            self._link(target, [*path, reaction], done)
        if refs:
            reaction.substitute({s: t.rate_coefficient for s, t in refs.items()})

        done.add(id(reaction))

    def resolve_rate_references(self, expr: Expr) -> Expr:
        """Replace each ``rc_<N>`` in *expr* with reaction *N*'s rate coefficient.

        For expressions outside the reactions (e.g. a ``.jfunc``
        ``heatingCoolingRate``), once :meth:`link` has run.

        Parameters
        ----------
        expr : Expr
            Expression to resolve.

        Returns
        -------
        Expr
            *expr* with every ``rc_<N>`` replaced.

        Raises
        ------
        ParserError
            If an ``rc_`` is not followed by a number or names a reaction that
            is not in the network.
        """
        refs = self._references([expr])

        return cast(Expr, expr.xreplace({s: t.rate_coefficient for s, t in refs.items()}))

    def _references(self, exprs: Sequence[Basic]) -> dict[Basic, Reaction]:
        """The ``rc_<N>`` symbols in *exprs*, each with its target reaction."""
        refs: dict[Basic, Reaction] = {}
        for expr in exprs:
            symbols = getattr(expr, "free_symbols", set())
            for symbol in sorted(symbols, key=str):
                name = str(symbol)
                if name.lower().startswith("rc_"):
                    refs[symbol] = self._reference_target(name)

        return refs

    def _reference_target(self, name: str) -> Reaction:
        """The reaction an ``rc_<N>`` symbol called *name* points at."""
        try:
            number = int(name[3:])
        except ValueError:
            raise ParserError(
                f"'{name}': 'rc_' must be followed by a reaction number, e.g. rc_3"
            ) from None

        target = self.by_source_index(number)
        if target is None:
            raise ParserError(
                f"'{name}' references reaction {number}, which is not in the network"
            )

        return target

    def by_source_index(self, n: int) -> Reaction | None:
        """Return the reaction whose file-side number (``source_index``) is *n*.

        This is the file-order reaction number (``Reaction.index``), the same
        one ``.jfunc`` ``chemRateN`` binds to — NOT the positional catalogue
        slot returned by ``self[n]``. Returns ``None`` if no reaction carries
        that number (e.g. a reaction merged away by ``duplicate_policy``).

        Parameters
        ----------
        n : int
            File-side reaction number.

        Returns
        -------
        Reaction | None
        """
        return self._by_source_index.get(n)

    def from_serialized(self, serialized: str) -> list[Reaction]:
        """Return all reactions sharing a serialized form.

        A serialized form (species-only) can be shared by several reactions
        that differ by mechanism/``type`` (e.g. thermal vs cosmic-ray
        desorption).  All matching reactions are returned as a list, whatever
        the count.  For scalar-or-list ergonomics use ``reactions[serialized]``;
        for a non-raising ``Vector`` use :meth:`all`.

        Parameters
        ----------
        serialized : str
            Canonical form ``"<sorted_reactants>__<sorted_products>"``, e.g.
            ``"H.H2O+__H2.OH+"``.

        Returns
        -------
        list[Reaction]
            Every reaction with that serialized form.

        Raises
        ------
        KeyError
            If no reaction has that serialized form.
        """
        return self._by_serialized[serialized]

    def all(self, serialized: str) -> Vector[Reaction]:
        """Return all reactions sharing *serialized* (empty if none).

        Unlike :meth:`from_serialized`, a missing key yields an empty
        ``Vector`` rather than raising.

        Parameters
        ----------
        serialized : str
            Canonical ``"<sorted_reactants>__<sorted_products>"`` form, e.g.
            ``"H.H2O+__H2.OH+"``.

        Returns
        -------
        Vector[Reaction]
            Every reaction with that serialized form, or an empty ``Vector``.
        """
        return Vector(self._by_serialized.get(serialized, []))

    def from_verbatim(self, verbatim: str, type: str | None = None) -> Reaction | None:
        """Look up a reaction by its verbatim string.

        Parameters
        ----------
        verbatim : str
            Human-readable string (e.g. ``"H + H2O+ -> H2 + OH+"``).
        type : str or None, optional
            If supplied, return ``None`` when the reaction type does not match.
            Also disambiguates when several reactions share a verbatim string.

        Returns
        -------
        Reaction or None
            The first matching reaction, or ``None`` if not found.
        """
        for rea in self._by_prop.get(verbatim, []):
            if type is None or rea.type == type:
                return rea

        return None

    def get_list(self) -> list[Reaction]:
        """Return the underlying ordered list of ``Reaction`` objects.

        Returns
        -------
        list[Reaction]
            The catalogue's own list (not a copy), in catalogue order.
        """
        return self._list

    def get(self, reaction: str, type: str | None = None) -> Reaction | None:
        """Look up a reaction by name or serialized form, with optional type filter.

        When *type* is supplied and *reaction* is a serialized form shared by
        several mechanisms, the exact ``(serialized, type)`` match is returned.
        Never raises: a missing lookup, or one with a *type* that does not
        match exactly one reaction, yields ``None``.

        Parameters
        ----------
        reaction : str
            Verbatim string (e.g. ``"H + H2O+ -> H2 + OH+"``) or serialized
            form (e.g. ``"H.H2O+__H2.OH+"``).
        type : str or None, optional
            If given, return ``None`` when the reaction type does not match.

        Returns
        -------
        Reaction, list of Reaction, or None
            The matching reaction, the list of matching reactions (without
            *type*, for a key shared by several mechanisms), or ``None`` if not
            found.

        Notes
        -----
        Without *type*, a key shared by several mechanisms returns the
        ``list`` of them, as ``reactions[key]`` does.
        """
        if type is not None and (
            reaction in self._by_serialized or reaction in self._by_prop
        ):
            try:
                return self[reaction, type]
            except KeyError:
                return None

        try:
            rea = self[reaction]
        except KeyError:
            return None

        if type is None or rea.type == type:
            return rea

    def with_type(self, type: str) -> Vector[Reaction]:
        """Return all reactions matching the given reaction type.

        Parameters
        ----------
        type : str
            A reaction type string, e.g. ``"photo"``, ``"cosmic_ray"``,
            ``"3_body"``, ``"unknown"``, or a grain surface-mechanism type
            (see :class:`~jaff.core.reaction.Reaction`).

        Returns
        -------
        Vector[Reaction]
            In catalogue order.  Cached: the same object is returned until a
            reaction is added, so treat it as read-only.
        """
        key = ("type", type)
        if key not in self.__filter_cache:
            self.__filter_cache[key] = Vector([r for r in self if r.type == type])

        return self.__filter_cache[key]

    def with_format(self, format: str) -> Vector[Reaction]:
        """Return all reactions produced by the named network-format parser.

        Parameters
        ----------
        format : str
            Parser name, e.g. ``"prizmo"``, ``"krome"``, ``"kida"``,
            ``"udfa"``, ``"uclchem"``; compared against the reaction's private
            ``_format`` field.

        Returns
        -------
        Vector[Reaction]
            In catalogue order.  Cached: the same object is returned until a
            reaction is added, so treat it as read-only.
        """
        key = ("format", format)
        if key not in self.__filter_cache:
            self.__filter_cache[key] = Vector([r for r in self if r._format == format])

        return self.__filter_cache[key]

    def verbatim(self) -> Vector[str]:
        """Return a ``Vector`` of verbatim reaction strings.

        Returns
        -------
        Vector[str]
            The ``"R1 + R2 -> P1 + P2"`` string of each reaction, in catalogue
            order.
        """
        return Vector([r.verbatim for r in self])

    def types(self) -> Vector[str]:
        """Return a ``Vector`` of reaction type strings.

        Returns
        -------
        Vector[str]
            The type label (e.g. ``"photo"``, ``"cosmic_ray"``) of each
            reaction, in catalogue order.
        """
        return Vector([r.type for r in self])

    def reactants(self) -> Vector[Species]:
        """Return a ``Vector`` of reactant ``Species`` catalogues.

        Returns
        -------
        Vector[Species]
        """
        return Vector([r.reactants for r in self])

    def products(self) -> Vector[Species]:
        """Return a ``Vector`` of product ``Species`` catalogues.

        Returns
        -------
        Vector[Species]
        """
        return Vector([r.products for r in self])

    def rate_coefficients(self) -> Vector[Basic]:
        """Return a ``Vector`` of rate-coefficient SymPy expressions.

        ``Reaction.rate_coefficient`` of each reaction, in catalogue order.
        Units depend on the reaction order (e.g. cm³ s⁻¹ for two-body
        reactions).  For the full rates, including the reactant densities,
        use :meth:`rates`.

        Returns
        -------
        Vector[Basic]
        """
        return Vector([r.rate_coefficient for r in self])

    def rates(self) -> Vector[Basic]:
        """Return a ``Vector`` of reaction rates, ``rate_coefficient · Π nden``.

        ``Reaction.rate`` of each reaction (the rate per unit volume: the rate
        coefficient times each core reactant's number density), in catalogue
        order.  For the coefficients alone use :meth:`rate_coefficients`.

        Returns
        -------
        Vector[Basic]
        """
        return Vector([r.rate for r in self])

    def odes(self, species: Species) -> list[Basic]:
        """Right-hand sides of the species ODEs, one per core species.

        ``d nden[s]/dt`` is the sum of the :meth:`rates` of the reactions
        producing *s* minus those destroying it; a species appearing several
        times on one side counts once per occurrence:

        ```
        d nden[s] / dt = sum_{i: s in products(i)} rate_i
                       - sum_{i: s in reactants(i)} rate_i
        ```

        Parameters
        ----------
        species : Species
            The network's species catalogue: the result has one entry per core
            species (pseudo-species such as ``_PHOTON`` are not integrated), in
            ``species.index`` order.

        Returns
        -------
        list of sympy.Basic
            The time derivative of each core species' number density
            (cm⁻³ s⁻¹).

        Examples
        --------
        >>> from jaff import Network
        >>> net = Network("h_photoionization")  # H + γ -> H+ + e-, H+ + e- -> H
        >>> [s.name for s in net.species.core]
        ['H', 'H+', 'e-']
        >>> odes = net.reactions.odes(net.species)
        >>> odes[1] == -odes[0]
        True
        """
        rates = self.rates()
        odes: list[Basic] = [Float(0.0) for _ in range(species.core.count)]
        for rate, reaction in zip(rates, self, strict=True):
            for reactant in reaction.reactants.core:
                odes[reactant.index] -= rate
            for product in reaction.products.core:
                odes[product.index] += rate

        return odes

    def tmins(self) -> Vector[float | None]:
        """Return a ``Vector`` of lower temperature bounds (K or ``None``).

        Returns
        -------
        Vector[float | None]
            Lower bound of each reaction in Kelvin, ``None`` where there is
            none, in catalogue order.
        """
        return Vector([r.tmin for r in self])

    def tmaxes(self) -> Vector[float | None]:
        """Return a ``Vector`` of upper temperature bounds (K or ``None``).

        Returns
        -------
        Vector[float | None]
            Upper bound of each reaction in Kelvin, ``None`` where there is
            none, in catalogue order.
        """
        return Vector([r.tmax for r in self])

    def dE(self) -> Vector[Basic]:
        """Return a ``Vector`` of energy-change SymPy expressions (erg).

        Returns
        -------
        Vector[Basic]
            ``Reaction.dE`` of each reaction (energy released per reaction
            event), in catalogue order.
        """
        return Vector([r.dE for r in self])

    def serialized(self) -> Vector[str]:
        """Return a ``Vector`` of name-level serialized reaction strings.

        These keys index the catalogue and are used for duplicate detection.

        Returns
        -------
        Vector[str]
            ``"<sorted_reactants>__<sorted_products>"`` of each reaction, in
            catalogue order.
        """
        return Vector([r.serialized for r in self])

    def serialized_exploded(self) -> Vector[str]:
        """Return a ``Vector`` of atom-level (isomer-insensitive) serialized strings.

        Built from ``Specie.serialized`` forms, so two reactions that differ
        only in isomers (e.g. HCO+ vs HOC+) share a key.

        Returns
        -------
        Vector[str]
            ``Reaction.serialized_exploded`` of each reaction, in catalogue
            order.
        """
        return Vector([r.serialized_exploded for r in self])

    def radiation_emitters(self) -> Vector[Reaction]:
        """Return all reactions that emit radiation (have a ``spectralYield``).

        A reaction emits when it has a ``.jfunc`` ``spectralYield<N>``
        function (``reaction.rad.spectral_yield`` is set); any reaction type
        can emit, not only photo reactions.

        Returns
        -------
        Vector[Reaction]
            Reactions whose ``rad.spectral_yield`` is set, in catalogue order.
            Cached: the same object is returned until a reaction is added, so
            treat it as read-only.
        """
        key = ("radiation_emitters",)
        if key not in self.__filter_cache:
            self.__filter_cache[key] = Vector(
                [r for r in self if r.rad.spectral_yield is not None]
            )

        return self.__filter_cache[key]

    def photo_reactions(self) -> Vector[Reaction]:
        """Return all photo-reactions (``type == "photo"``).

        Returns
        -------
        Vector[Reaction]
            In catalogue order.  Unlike :meth:`with_type`, not cached.
        """
        return Vector([r for r in self if r.type == "photo"])

    def photo_reaction_truths(self) -> Vector[int]:
        """Return a binary ``Vector`` marking photo-reactions with ``1``.

        Useful as a mask over the other vector accessors.

        Returns
        -------
        Vector[int]
            ``1`` for a photo-reaction, ``0`` otherwise; one entry per
            reaction, in catalogue order.
        """
        return Vector([int(reaction.type == "photo") for reaction in self])

    def photo_reaction_indices(self) -> Vector[int]:
        """Return the integer indices of photo-reactions within this catalogue.

        Use these to pick the photo-reaction entries out of the other vector
        accessors.

        Returns
        -------
        Vector[int]
            Zero-based catalogue position (``Reaction.catalogue_index``) of
            each photo-reaction, in ascending order.
        """
        return Vector([i for i, reaction in enumerate(self) if reaction.type == "photo"])

    def plot_rates(self, **kwargs: Any) -> tuple[plt.Figure, Any] | None:
        """Plot the rate coefficients of every reaction in the catalogue.

        Forwards to :func:`jaff.plotting.plot_rates` with all reactions as a
        flat list: one shared axes by default, one panel per reaction with
        ``grid=(rows, cols)``.  Accepts every ``plot_rates`` keyword
        (``tmin``, ``tmax``, ``xscale``, ``yscale``, ``labels``, ``style``,
        ``palette``, ``save``, ``show`` ...).

        Reactions whose rate cannot be evaluated numerically (e.g. photo
        reactions) are skipped with a warning.

        Returns
        -------
        tuple[matplotlib.figure.Figure, numpy.ndarray] or None
        """
        from ...plotting import plot_rates

        return plot_rates(list(self._list), **kwargs)

    def plot_xsecs(self, **kwargs: Any) -> tuple[plt.Figure, Any] | None:
        """Plot the photo cross sections of the catalogue's reactions.

        Forwards to :func:`jaff.plotting.plot_xsecs` with all reactions as a
        flat list.  Reactions without cross-section data are skipped.  Accepts
        every ``plot_xsecs`` keyword (``grid``, ``processes``, ``bands``,
        ``style``, ``palette``, ``save``, ``show`` ...).

        Returns
        -------
        tuple[matplotlib.figure.Figure, numpy.ndarray] or None
        """
        from ...plotting import plot_xsecs

        return plot_xsecs(list(self._list), **kwargs)
