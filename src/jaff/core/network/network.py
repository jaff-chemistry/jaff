"""Top-level ``Network`` class for loading and querying JAFF reaction networks.

A ``Network`` is the main entry point for users of the JAFF library.  It
reads a reaction network file in any supported format (KROME, PRIZMO, UDFA,
KIDA, UCLChem, or the binary ``.jaff`` serialisation), builds typed
``Species`` and ``Reactions`` catalogues, validates the network, and exposes
symbolic ODE and flux expressions for downstream code generation.

Auxiliary ``.jfunc`` files may accompany the network file.  They can supply
custom rate coefficients, chemical heating/cooling rates, and radiation
moment contributions using the ``@var`` / ``@function`` syntax parsed by
``AuxiliaryFunctionParser``.

Number densities are represented as ``nden[i]`` — a reference into a SymPy
``MatrixSymbol("nden", n_species, 1)``.  The symbol ``idx_X`` refers to the
integer position of species ``X`` within the network.
"""

from __future__ import annotations

import logging
import sys
from dataclasses import dataclass
from functools import cached_property, lru_cache
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, Mapping, Sequence, cast

import numpy as np
from sympy import (
    Basic,
    Expr,
    parse_expr,
    symbols,
)

from ...common import is_jaff_file, load_mass_dict, print_motd, resolve_dependencies
from ...config.loader import deep_merge
from ...errors import ParserError
from ...io import JaffLogger, jaff_progress
from ...io._io import JaffProps, from_jaff_file, to_jaff_file, write_data_table
from ...physics import (
    Dust,
    DustProps,
    EosProps,
    Photochemistry,
    Radiation,
    Thermodynamics,
    get_sfluxes,
    get_sodes,
    get_sradodes,
)
from ..elements import Elements
from ..parsers import NetworkParser
from ..reaction import (
    PhotoReactionProps,
    ReactionRadiationProps,
    RateSegment,
    RateSegments,
    Reaction,
    Reactions,
    ShieldingProps,
)
from ...physics.photo_reactions.shielding.config import ShieldingConfig
from ..species import Specie, Species
from ._spec import NetworkSpec
from ...config.base import ConfigModel
from ...physics.photo_reactions.config import RadiationConfig
from .config import ConfigLayer, NetworkConfig, load_network_config
from ._symbols import NetworkSymbols

if TYPE_CHECKING:
    from .._typing import ElementProps
    from ..parsers.auxiliary_func import AuxiliaryFunctionProps
    from ..parsers.network._formats import Parser


@lru_cache(maxsize=200000)
def _parse_rate_expr(rate: str) -> Expr:
    """Parse a rate string into a (non-evaluated) SymPy expression, memoized.

    Large networks contain many reactions with identical rate strings (≈40% of
    KIDA-2024 rates repeat), and the parsed expression depends only on the
    string — species substitution happens later in ``NetworkSymbols.standardize`` —
    so results are cached across reactions and networks.  SymPy expressions are
    immutable, making the shared objects safe to reuse.
    """
    return parse_expr(rate, evaluate=False)


def _keyword_layer(**kwargs: Any) -> Dict[str, Any]:
    """Return the config layer of :class:`Network`'s keyword arguments.

    Only the arguments given (not ``None``) are set, so an omitted one never
    overrides the ``jaff.toml`` or the ``config`` layers.
    """
    given = {key: value for key, value in kwargs.items() if value is not None}
    network = {
        k: given[k]
        for k in ("strict", "label", "funcfile", "expand_nuclei")
        if k in given
    }
    reactions: Dict[str, Any] = {}
    if "duplicate_policy" in given:
        reactions["duplicate_policy"] = given["duplicate_policy"]

    if "use_proxy_photoreaction" in given:
        reactions["photo"] = {"use_proxy_photoreaction": given["use_proxy_photoreaction"]}

    if reactions:
        network["reactions"] = reactions

    if "radiation" in given:
        network["radiation"] = _section(given["radiation"])

    return {"network": network} if network else {}


def _section(value: Mapping[str, Any] | ConfigModel) -> Dict[str, Any]:
    """Return a config section given as a dict or a model (its explicitly set keys)."""
    if isinstance(value, ConfigModel):
        return value.model_dump(exclude_unset=True)

    return dict(value)


def _jaffgen_metadata_layer(metadata: Mapping[str, Any]) -> Dict[str, Any]:
    """Return jaffgen's old-format reaction settings as a config layer.

    Transitional: jaffgen still passes its ``[network.reactions]`` /
    ``[network.rates]`` tables in the old format through ``_metadata``.
    """
    reactions: Dict[str, Any] = {}
    global_cutoff = metadata.get("rate_props", {}).get("T_cutoff")
    if global_cutoff is not None:
        reactions["t_cutoff"] = global_cutoff

    entries: Dict[str, Any] = {}
    for serialized, props in metadata.get("reaction_props", {}).items():
        entry: Dict[str, Any] = {}
        if "T_cutoff" in props:
            entry["t_cutoff"] = props["T_cutoff"]

        photo_keys = (
            "pi_database",
            "override_custom_rate",
            "override_custom_cross_section",
            "shielding",
        )
        photo = {key: props[key] for key in photo_keys if key in props}
        if photo:
            entry["photo"] = photo

        entries[serialized] = entry

    if entries:
        reactions["reaction"] = entries

    return {"network": {"reactions": reactions}} if reactions else {}


def _shielding_props(shielding: ShieldingConfig | None) -> ShieldingProps | None:
    """Return the reaction's :class:`ShieldingProps` from its shielding config."""
    if shielding is None:
        return None

    return ShieldingProps(type=shielding.type, props=dict(shielding.params))


@dataclass(frozen=True)
class SinkSourceReport:
    """Structured result of :meth:`Network.check_sink_sources`.

    Attributes
    ----------
    sinks : frozenset[str]
        Species names consumed but never produced.
    sources : frozenset[str]
        Species names produced but never consumed.
    """

    sinks: frozenset[str]
    sources: frozenset[str]

    @property
    def ok(self) -> bool:
        """``True`` when neither sinks nor sources were found."""
        return not self.sinks and not self.sources


@dataclass(frozen=True)
class RecombinationReport:
    """Structured result of :meth:`Network.check_recombinations`.

    Attributes
    ----------
    missing : frozenset[str]
        Positively charged species names lacking an electron recombination.
    """

    missing: frozenset[str]

    @property
    def ok(self) -> bool:
        """``True`` when every cation has an electron recombination."""
        return not self.missing


@dataclass(frozen=True)
class IsomerReport:
    """Structured result of :meth:`Network.check_isomers`.

    Attributes
    ----------
    groups : tuple[tuple[str, ...], ...]
        Groups of two or more species names sharing an atomic composition.
    """

    groups: tuple[tuple[str, ...], ...]

    @property
    def ok(self) -> bool:
        """``True`` when no isomer groups were found."""
        return not self.groups


class Network:
    """Astrochemical reaction network loaded from a file.

    After construction the following public attributes are available:

    Attributes
    ----------
    filename : Path
        Absolute path to the source network file.
    label : str
        Human-readable label for this network (defaults to the file stem).
    species : Species
        Ordered catalogue of the network's core (real) species.  Special
        pseudo-species (``_PHOTON``, ``_CR``, ``_GRAIN``, ...) are excluded;
        they live only inside each reaction's ``reactants``/``products`` and
        carry the reaction's identity without entering the integrated state.
    reactions : Reactions
        Ordered catalogue of all reactions in the network.
    elements : Elements
        Unique elements derived from all species.
    reactant_matrix : np.ndarray
        Integer stoichiometry matrix, shape ``(n_reactions, n_species)``.
        Entry ``[i, j]`` counts how many times species *j* appears as a
        reactant in reaction *i*.
    product_matrix : np.ndarray
        Integer stoichiometry matrix for products, same shape.
    thermodynamics : Thermodynamics
        Heating/cooling terms: ``dEdt_chemical`` (accumulated over all
        reactions) and ``dEdt_extra`` (from the ``heatingcoolingrate``
        auxiliary function, if present), both in erg cm⁻³ s⁻¹.
    symbols : NetworkSymbols
        Canonical symbols, density expressions, introspection sets and symbol
        standardization for this network.
    radiation : Radiation | None
        Radiation field object; ``None`` when no radiation bands are specified.
    mass_dict : dict[str, ElementProps]
        Element mass dictionary used during species construction.
    spec : NetworkSpec
        The normalized construction parameters this network was built from.
    config : NetworkConfig
        The validated config: the ``jaff.toml`` next to the network, the
        ``config`` layers, then the keyword arguments.
    """

    def __init__(
        self,
        fname: str | Path,
        config: ConfigLayer | Sequence[ConfigLayer] | None = None,
        strict: bool | None = None,
        label: str | None = None,
        funcfile: bool | str | Path | None = None,
        duplicate_policy: str | None = None,  # preserve-first, preserve-last, error
        expand_nuclei: bool | None = None,
        radiation_props: RadiationConfig | Mapping[str, Any] | None = None,
        dust_props: DustProps | None = None,
        use_proxy_photoreaction: bool | None = None,
        eos_props: EosProps | None = None,
        use_jaff_toml: bool = True,
        _from_cli: bool = False,
        _metadata: dict[str, Any] = {},
    ):
        """Load a reaction network from *fname*.

        Parameters
        ----------
        fname : str | Path
            Path to the network file, or the name of a predefined network (a
            sub-directory of :data:`~jaff.config.NETWORKS_DIR` containing a
            single ``.jet`` file).  A predefined network name wins over a
            same-named path on disk.  Supported extensions: any text format
            auto-detected by ``NetworkParser``, plus ``.jaff`` binary files.
        config : ConfigLayer or sequence of them, optional
            Config layers stacked over the ``jaff.toml`` next to the network
            file: a TOML file path, a ``dict`` shaped like a ``jaff.toml``
            (``{"network": {...}}``) or a ``JaffTomlConfig``; see
            :func:`~jaff.core.network.config.load_network_config`.  The keyword
            arguments below override them when given (not ``None``).
        strict : bool | None, optional
            If ``True``, treat conservation violations and validation problems
            as fatal errors (process exits).  Default (``network.strict``):
            ``False`` (warnings only).
        label : str | None, optional
            Human-readable name for this network.  Defaults to the file stem.
        funcfile : bool | str | Path | None, optional
            Path to a ``.jfunc`` auxiliary function file.  ``True`` (the
            default) scans the network's directory for ``<network>.jfunc``;
            ``False`` skips auxiliary-function loading entirely.
        duplicate_policy : str | None, optional
            How to resolve two rate coefficients sharing a reaction, mechanism,
            and temperature range: ``"preserve-first"`` (keep the first, drop
            later duplicates; the default), ``"preserve-last"`` (keep the last),
            or ``"error"`` (raise a ``ParserError``).  Sets
            ``network.reactions.duplicate_policy``; a reaction's own policy in
            the config still wins.
        expand_nuclei : bool | None, optional
            When ``True`` (default), an ``n_<element>_nuc`` symbol (e.g.
            ``n_H_nuc``, ``n_He_nuc``) in rate expressions is expanded to the
            element-nucleus sum — ``Σ atom-count · nden[i]`` over all species
            bearing that element.  Set to ``False`` to keep it as a free symbol
            ``n<element>_nuc`` (e.g. ``nh_nuc``).  (Plain ``n_X`` always resolves
            to species ``X``; only ``n_X_nuc`` is affected by this flag.)
        radiation_props : RadiationConfig | dict | None, optional
            Radiation field (``network.radiation``): enables the network's
            :class:`Radiation` and photochemistry.  Its keys override those of
            a ``[network.radiation]`` table in the config.  Without either, no
            radiation field is built and photochemistry is disabled.
        dust_props : DustProps | None, optional
            Dust configuration used to construct the network's :class:`Dust`
            object.  When ``None`` (default), the dust module is disabled.
        use_proxy_photoreaction : bool | None, optional
            When ``True`` (default), cross-sections for photo-reactions are looked up
            using the proxy photo-reaction string (via
            :meth:`Reaction.normalized_proxy_reaction_str`) rather than the
            reaction's standard serialized form, so proxy species such as
            ``CHx``/``OHx`` resolve to the ``CH``/``OH`` entries.  A custom
            ``.jfunc`` rate still takes precedence over the cross-sections.
            Sets ``network.reactions.photo.use_proxy_photoreaction``.
        eos_props : EosProps | None, optional
            Equation-of-state configuration used by
            :attr:`Thermodynamics.eos`.  When ``None`` (default), an ideal gas
            with the default adiabatic index (``EosProps("ideal")``) is used.
        use_jaff_toml : bool, optional
            Stack the ``jaff.toml`` next to the network file, if any, under the
            ``config`` layers (default ``True``).

        Raises
        ------
        FileNotFoundError
            If *fname* is neither an existing file nor a predefined network
            name, or a predefined network directory has no ``.jet`` file.
        ParserError
            If a predefined network directory has more than one ``.jet`` file,
            or the config is invalid.
        """
        self.logger: logging.Logger = JaffLogger().get_logger()

        overrides = deep_merge(
            _jaffgen_metadata_layer(_metadata),
            _keyword_layer(
                strict=strict,
                label=label,
                funcfile=funcfile,
                duplicate_policy=duplicate_policy,
                expand_nuclei=expand_nuclei,
                use_proxy_photoreaction=use_proxy_photoreaction,
                radiation=radiation_props,
            ),
        )
        self.config: NetworkConfig = load_network_config(
            Path(fname),
            config,
            use_jaff_toml=use_jaff_toml,
            overrides=overrides,
        )
        self.spec: NetworkSpec = NetworkSpec(fname, self.config, _from_cli, _metadata)

        jaff_props: JaffProps = {}  # type: ignore
        loaded_from_jaff_file = is_jaff_file(self.spec.fname)
        if loaded_from_jaff_file:
            jaff_props = from_jaff_file(self.spec.fname, self.spec.strict)
            self.spec.fname = jaff_props.get("file_name", self.spec.fname)
        # Finalize the label: a ``.jaff`` payload label wins, else the
        # configured label, else the file stem.
        self.spec.label = jaff_props.get(
            "label", self.config.label or self.spec.fname.stem
        )
        if not _from_cli:
            print_motd()

        self.mass_dict: dict[str, ElementProps] = {}
        self.species: Species = Species()
        self.reactions: Reactions = Reactions()
        self.symbols: NetworkSymbols = NetworkSymbols(self)
        self.reactant_matrix: np.ndarray | None = None
        self.product_matrix: np.ndarray | None = None
        self.radiation: Radiation | None = (
            Radiation(self, self.config.radiation)
            if self.config.radiation is not None
            else None
        )
        self.eos_props: EosProps = (
            eos_props if eos_props is not None else EosProps("ideal")
        )
        self._photochemistry: None | Photochemistry = None
        self.dust: Dust | None = (
            Dust(self, dust_props) if dust_props is not None else None
        )
        self.__tgas_clamp_cache: dict[tuple[float | None, float | None], Expr] = {}

        self.logger.info(f"Loading network from {self.spec.fname}")
        self.logger.info(f"Network label: [yellow]{self.label}[/]")

        self.mass_dict: dict[str, ElementProps] = load_mass_dict()
        Species.configure(self.mass_dict)

        parsers: list[Parser] = []
        if not loaded_from_jaff_file:
            parsers = self.__load_network()
        else:
            self.__load_network_from_jaff_file(jaff_props)

        self.elements: Elements = Elements(self.species._list)
        self.__post_process(parsers)
        self.__normalize_network_extras(loaded_from_jaff_file)
        self.__set_emissivity(jaff_props.get("emissivity"))
        self.__set_absorption_coefficient(jaff_props.get("absorption_coefficient"))
        self.thermodynamics: Thermodynamics = Thermodynamics(
            self, jaff_props.get("dEdt_extra")
        )
        self.symbols.log_summary()

        self.check_sink_sources(self.spec.strict)
        self.check_recombinations(self.spec.strict)
        self.check_isomers(self.spec.strict)

        self.__generate_reaction_matrices()

        self.logger.info("[green]Network loaded successfully![/]")

    @cached_property
    def filename(self) -> Path:
        """Absolute path to the source network file (from :attr:`spec`)."""
        return self.spec.fname

    @property
    def label(self) -> str:
        """Human-readable label for this network (from :attr:`spec`)."""
        # Finalized to a non-None value in __init__ (label or file stem).
        assert self.spec.label is not None
        return self.spec.label

    def __load_network(self) -> list[Parser]:
        """Parse the network file and build species, reactions, and auxiliary quantities.

        Returns
        -------
        list[Parser]
            The parsers that produced reactions, for :meth:`__post_process`.
        """
        specie_names = set()
        special_species: dict[str, Specie] = {}
        reactions_config = self.config.reactions

        with NetworkParser(self.spec.fname, self.logger) as netp:
            reactions_list, global_vars = netp.get_parsed()
            parsers = netp.get_parsers()

        aux_funcs = self.spec.aux_funcs

        global_vars = {
            var: resolve_dependencies(expr, {}, aux_funcs)
            for var, expr in global_vars.items()
        }
        subs_dict: dict[Basic, Basic] = {
            symbols(var.lower()): expr for var, expr in global_vars.items()
        }

        for i, reaction in enumerate(
            jaff_progress.track(
                reactions_list,
                description=f"Creating {self.label} network",
            )
        ):
            reactants: list[str] = reaction["r"]
            products: list[str] = reaction["p"]
            tmin: float | None = reaction["tmin"]
            tmax: float | None = reaction["tmax"]
            rate: str = reaction["rate"]
            si = reaction["source_index"]
            aux_chem_rate = f"chemrate{si}"

            for s in reactants + products:
                if s in specie_names:
                    continue
                specie_names.add(s)
                if s.startswith("_"):
                    special_species[s] = Specie(s, -1)
                else:
                    self.species.add(Specie(s, self.species.count))

            rr = [
                special_species[r] if r.startswith("_") else self.species[r]
                for r in reactants
            ]
            pp = [
                special_species[p] if p.startswith("_") else self.species[p]
                for p in products
            ]

            local_subs_dict = {**subs_dict}
            srxn = Reaction.serialize(rr, pp)
            # The global reaction settings with this reaction's own applied.
            settings = reactions_config.resolve(srxn)
            assert settings.t_cutoff is not None and settings.photo is not None

            rate_expr = self.__parse_rate(aux_chem_rate, rate, aux_funcs, global_vars)
            rate_expr = resolve_dependencies(rate_expr, local_subs_dict, aux_funcs)

            rtype = reaction.get("type", "unknown")
            if (srxn, rtype) in self.reactions:
                rea = self.reactions[srxn, rtype]
                if (tmin, tmax) in rea.rate_segments._by_prop:
                    # Same reaction, mechanism, and temperature range with a
                    # second rate coefficient — resolved per duplicate_policy.
                    existing = rea.rate_segments._by_prop[(tmin, tmax)]
                    policy = settings.duplicate_policy
                    if policy == "error":
                        raise ParserError(
                            f"Duplicate reaction: {rea.verbatim} [type={rtype}] "
                            f"has two rate coefficients over the same temperature "
                            f"range ({tmin}, {tmax}):\n"
                            f"  existing: {existing.rate}\n"
                            f"  new:      {rate_expr}\n"
                            "Same reaction, same mechanism and overlapping "
                            "temperature range cannot be resolved automatically. "
                            "Deduplicate the source network."
                        )

                    kept, dropped = (
                        (existing.rate, rate_expr)
                        if policy == "preserve-first"
                        else (rate_expr, existing.rate)
                    )
                    self.logger.warning(
                        f"Duplicate rate for [cyan]{rea.verbatim}[/] "
                        f"[type={rtype}] over T=({tmin}, {tmax}) resolved by "
                        f"policy '{policy}': kept {kept}, "
                        f"dropped {dropped}"
                    )
                    if policy == "preserve-first":
                        continue
                    # "preserve-last": fall through; add() replaces the segment.

                rea.rate_segments.add(RateSegment(rate_expr, tmin, tmax))
                continue

            # The reaction's own photo settings; unset ones fall back to the
            # global ones where they are used.
            entry = reactions_config.reaction.get(srxn)
            own_photo = entry.photo if entry is not None else None
            if own_photo is not None and rtype != "photo":
                keys = ", ".join(sorted(own_photo.model_fields_set))
                raise ParserError(f"{keys} is only valid for photo reactions: {srxn}")

            rea = Reaction(
                self,
                reactants=rr,
                products=pp,
                rate_coefficient=rate_expr,
                tmin=tmin,
                tmax=tmax,
                original_string=reaction["string"],
                index=si,
                type=rtype,
                t_cutoff=settings.t_cutoff,
                photo=PhotoReactionProps(
                    pi_database=own_photo.pi_database if own_photo else None,
                    override_custom_rate=(
                        own_photo.override_custom_rate if own_photo else None
                    ),
                    override_custom_cross_section=(
                        own_photo.override_custom_cross_section if own_photo else None
                    ),
                    shielding=_shielding_props(settings.photo.shielding),
                )
                if rtype == "photo"
                else None,
                format=reaction.get("format"),
                metadata=reaction.get("metadata"),
            )
            self.reactions.add(rea)

        if self._photochemistry is not None:
            self._photochemistry.log_fallback_summary()

        self.logger.info(f"Loaded {self.reactions.count} reactions")
        self.logger.info(
            f"Loaded {len(self.reactions.with_type('photo'))} photo-chemistry reactions"
        )

        return parsers

    def __load_network_from_jaff_file(self, jaff_props: JaffProps):
        """Restore species, reactions, and radiation state from a ``.jaff`` file payload.

        Parameters
        ----------
        jaff_props : JaffProps
            Deserialised property dictionary returned by
            :func:`~jaff.io._io.from_jaff_file`.
        """
        self.species = jaff_props["species"]

        for i, reaction in enumerate(jaff_props["reactions"]):
            rtype = reaction.get("reaction_type", "unknown")
            shielding = reaction.get("shielding")

            rea = Reaction(
                self,
                reactants=reaction["reactants"],
                products=reaction["products"],
                rate_coefficient=reaction["rate"],
                dE=reaction["dE"],
                tmin=reaction["tmin"],
                tmax=reaction["tmax"],
                t_cutoff=reaction.get("t_cutoff", "clip"),
                original_string=reaction["original_string"],
                index=reaction.get("index", i),
                type=rtype,
                format=reaction.get("format"),
                photo=PhotoReactionProps(
                    xsec=reaction.get("xsec_props"),
                    pi_database=reaction.get("pi_database"),
                    custom_rad_rate=bool(reaction.get("custom_rad_rate")),
                    shielding=ShieldingProps(**shielding) if shielding else None,
                )
                if rtype == "photo"
                else None,
                rad=ReactionRadiationProps(spectral_yield=reaction.get("spectral_yield")),
            )

            segments = reaction.get("rate_segments")
            if segments:
                rea.rate_segments = RateSegments(
                    [RateSegment(s["rate"], s["tmin"], s["tmax"]) for s in segments],
                    rea.t_cutoff,
                )
            self.reactions.add(rea)

        if self._photochemistry is not None:
            self._photochemistry.log_fallback_summary()

    def __post_process(self, parsers: list[Parser]) -> None:
        """Freeze the catalogues and run each parser's ``post_process`` hook.

        Species and reactions are final from here on: hooks may change
        reactions but not add any (``add()`` raises once frozen).  Runs before
        :meth:`__normalize_network_extras`, so whatever a hook writes is
        standardised like every other rate.

        Parameters
        ----------
        parsers : list[Parser]
            Parsers that produced reactions; empty for a ``.jaff`` network.
        """
        self.reactions._freeze()
        self.species._freeze()

        for parser in parsers:
            parser.post_process(self)

    def __normalize_network_extras(self, loaded_from_jaff: bool = False):
        """Standardize convenience symbols in all rate and auxiliary expressions.

        Replaces shorthand symbols (``n_X``, ``n_X_nuc``, ``n_e``, ``ntot``, …) with
        ``nden[i]`` references in every reaction rate.

        Parameters
        ----------
        loaded_from_jaff : bool, optional
            When ``True``, the network was restored from a ``.jaff`` file whose
            stored ``rate`` is already the final (piecewise-collapsed,
            standardized) expression.  Re-evaluating the rate segments would
            wrap the already-collapsed rate in a second piecewise, so the stored
            rate is used as-is instead.
        """
        for r in self.reactions:
            if loaded_from_jaff:
                r.rate_coefficient = self.symbols.standardize(r.rate_coefficient)
            elif r.type == "photo" and self.radiation is not None:
                r.rate_coefficient = self.symbols.standardize(r.rate_coefficient)
                r.rate_segments[0].rate = r.rate_coefficient
            else:
                for seg in r.rate_segments:
                    seg.rate = self.symbols.standardize(seg.rate)
                r.rate_coefficient = r.rate_segments.sort().evaluate_equivalent_rate(r)

            r.tmin, r.tmax = r.rate_segments[0].tmin, r.rate_segments[-1].tmax
            if r.rad.group_yields is not None:
                r.rad.group_yields = [
                    self.symbols.standardize(g) for g in r.rad.group_yields
                ]
            if r.photo is not None:
                self.__standardize_photo_band_data(r.photo)

    def __standardize_photo_band_data(self, photo: PhotoReactionProps) -> None:
        """Standardize a photo reaction's band data, as its rate already is.

        A symbolic ``photoCrossSection`` (e.g. one using ``n_X``) leaves those
        symbols in the band-averaged cross sections and rate coefficients that
        :func:`~jaff.physics._equations.get_sradodes` builds the equations from.
        """
        if photo.xsec is not None and isinstance(photo.xsec.total, Basic):
            photo.xsec.total = self.symbols.standardize(photo.xsec.total)

        for band in photo.rgrouped or []:
            band.k = self.symbols.standardize(band.k)
            if isinstance(band.xsec, Basic):
                band.xsec = self.symbols.standardize(band.xsec)

    def __set_emissivity(self, stored: Expr | None) -> None:
        """Integrate the ``globalEmissivity`` over the radiation groups.

        Stores the spectrum on ``radiation.emissivity`` and its standardized
        per-group integrals on each ``group.emissivity``.  The ``.jfunc``
        ``globalEmissivity`` function wins over the expression *stored* in a ``.jaff``
        file.  Does nothing without radiation or without an emissivity.
        """
        if self.radiation is None:
            return

        aux = self.spec.aux_funcs.get("globalemissivity")
        emissivity = cast(Expr, aux.expr) if aux is not None else stored
        if emissivity is None:
            return

        self.radiation.emissivity = emissivity
        groups = self.radiation.groups
        for group, value in zip(
            groups, self.radiation.integrate_over_groups(emissivity), strict=True
        ):
            group.emissivity = self.symbols.standardize(value)

    def __set_absorption_coefficient(self, stored: Expr | None) -> None:
        """Average the ``globalAbsorptionCoefficient`` over the radiation groups.

        Stores the spectrum on ``radiation.absorption_coefficient`` and its
        standardized per-group averages on each ``group.absorption_coefficient``.
        The ``.jfunc`` function wins over the expression *stored* in a ``.jaff``
        file.  Does nothing without radiation or without a coefficient.
        """
        if self.radiation is None:
            return

        aux = self.spec.aux_funcs.get("globalabsorptioncoefficient")
        coefficient = cast(Expr, aux.expr) if aux is not None else stored
        if coefficient is None:
            return

        self.radiation.absorption_coefficient = coefficient
        groups = self.radiation.groups
        for group, value in zip(
            groups, self.radiation.average_over_groups(coefficient), strict=True
        ):
            group.absorption_coefficient = self.symbols.standardize(value)

    @staticmethod
    def __parse_rate(
        aux_chem_rate: str,
        rate: str,
        aux_funcs: dict[str, AuxiliaryFunctionProps],
        global_vars: dict[str, Expr],
    ) -> Expr:
        """Convert a raw rate string to a SymPy expression.

        Checks, in priority order:

        1. Whether an auxiliary function named *aux_chem_rate* exists (custom rate).
        2. Whether *rate* is a global variable name.
        3. Falls back to ``sympy.parse_expr``.

        Format-specific rates (e.g. photo rates) arrive as placeholder calls that
        the producing parser resolves in its ``post_process`` hook.

        Parameters
        ----------
        aux_chem_rate : str
            Key for the optional custom-rate auxiliary function (e.g. ``"chemrate0"``).
        rate : str
            Raw rate string from the network file.
        aux_funcs : dict[str, AuxiliaryFunctionProps]
            Parsed auxiliary functions, keyed by name.
        global_vars : dict[str, Basic]
            Resolved global variable map from the network file.

        Returns
        -------
        Expr
            The rate expression.
        """
        if aux_chem_rate in aux_funcs:
            rate_expr = aux_funcs[aux_chem_rate].expr
        elif rate in global_vars:
            rate_expr = symbols(rate)
        else:
            rate_expr = _parse_rate_expr(rate)

        if not isinstance(rate_expr, Expr):
            raise ParserError(f"Rate expression is not an Expr: {rate_expr}")

        return rate_expr

    def to_jaff(self, filename: str | Path):
        """Serialise this network to a binary ``.jaff`` file.

        Parameters
        ----------
        filename : str | Path
            Output file path.  The ``.jaff`` extension is conventional but
            not enforced.
        """
        to_jaff_file(filename, self)

    def compare_reactions(self, other: Network, verbosity: int = 1):
        """Log reactions present in one network but not the other.

        Parameters
        ----------
        other : Network
            Network to compare against.
        verbosity : int, optional
            When ``1`` (default), print the reaction sets for common and
            unique reactions.  Other values suppress output.
        """
        self.logger.info(f'Comparing networks "{self.label}" and "{other.label}"...')

        self_reacts = {rea.serialized for rea in self.reactions}
        other_reacts = {rea.serialized for rea in other.reactions}

        common = self_reacts & other_reacts
        not_in_self = other_reacts - common
        not_in_other = self_reacts - common

        if verbosity == 1:
            self.logger.info(f"Reactions not present in {self.label}:")
            print(
                "\n".join(
                    str(r) for rea in not_in_self for r in other.reactions.all(rea)
                ),
                "\n",
            )

            self.logger.info(f"Reactions not present in {other.label}:")
            print(
                "\n".join(
                    str(r) for rea in not_in_other for r in self.reactions.all(rea)
                ),
                "\n",
            )

            self.logger.info(f"Reactions present in both {self.label} and {other.label}:")
            print(
                "\n".join(str(r) for rea in common for r in self.reactions.all(rea)),
                "\n",
            )

        self.logger.info(f"{len(common)} reactions are common in both networks")
        self.logger.info(f'{len(not_in_self)} reactions are missing in "{self.label}"')
        self.logger.info(f'{len(not_in_other)} reactions are missing in "{other.label}"')

    def compare_species(self, other: Network, verbosity: int = 1) -> None:
        """Log species present in one network but not the other.

        Parameters
        ----------
        other : Network
            Network to compare against.
        verbosity : int, optional
            When ``1`` (default), print the species sets.
        """
        self.logger.info(
            f'Comparing species in networks "{self.label}" and "{other.label}"...'
        )

        self_species = {sp.serialized for sp in self.species}
        other_species = {sp.serialized for sp in other.species}

        common = self_species & other_species
        not_in_self = other_species - common
        not_in_other = self_species - common

        if verbosity == 1:
            self.logger.info(f"Species not present in {self.label}:")
            print(
                ", ".join([str(other.species[sp]) for sp in not_in_self]),
                "\n",
            )

            self.logger.info(f"Species not present in {other.label}:")
            print(
                ", ".join([str(self.species[sp]) for sp in not_in_other]),
                "\n",
            )

            self.logger.info(f"Species present in both {self.label} and {other.label}:")
            print(
                ", ".join([str(self.species[sp]) for sp in common]),
                "\n",
            )

        self.logger.info(f"{len(common)} species are common in both networks")
        self.logger.info(f'{len(not_in_self)} species are missing in "{self.label}"')
        self.logger.info(f'{len(not_in_other)} species are missing in "{other.label}"')

    def check_sink_sources(self, errors: bool) -> SinkSourceReport:
        """Report (and warn, or abort) species never produced or never consumed.

        A *sink* species appears as a reactant in at least one reaction but
        is never produced.  A *source* species is produced but never consumed.
        The special species ``"_DUMMY"`` is excluded from the check.

        Parameters
        ----------
        errors : bool
            If ``True``, call ``sys.exit()`` when sinks or sources are found.

        Returns
        -------
        SinkSourceReport
            The sinks and sources found (empty when the network is balanced).
        """
        produced = {p.name for rea in self.reactions for p in rea.products}
        consumed = {r.name for rea in self.reactions for r in rea.reactants}
        species_names = {s.name for s in self.species if s.name != "_DUMMY"}

        report = SinkSourceReport(
            sinks=frozenset(species_names - produced),
            sources=frozenset(species_names - consumed),
        )

        for name in report.sinks:
            self.logger.info(f"Sink: [cyan]{name}[/]")

        for name in report.sources:
            self.logger.info(f"Source: [cyan]{name}[/]")

        if report.sinks:
            self.logger.warning("Sink detected")

        if report.sources:
            self.logger.warning("Source detected")

        if not report.ok and errors:
            self.logger.error("Exiting since errors are enabled")
            sys.exit()

        return report

    def check_recombinations(self, errors: bool) -> RecombinationReport:
        """Report (and warn, or abort) cations lacking an electron recombination.

        Parameters
        ----------
        errors : bool
            If ``True``, call ``sys.exit(1)`` when recombination reactions
            are missing.

        Returns
        -------
        RecombinationReport
            The positively charged species with no electron recombination.
        """
        electron_recomb_species = set()

        for rea in self.reactions:
            reactant_names = {r.name for r in rea.reactants}

            if "e-" in reactant_names:
                for r in rea.reactants:
                    if r.name != "e-":
                        electron_recomb_species.add(r.name)

        missing = {
            sp.name
            for sp in self.species
            if sp.charge > 0 and sp.name not in electron_recomb_species
        }

        for name in missing:
            self.logger.warning(f"Electron recombination not found for [cyan]{name}[/]")

        report = RecombinationReport(missing=frozenset(missing))

        if not report.ok and errors:
            self.logger.error("Recombination errors found")
            sys.exit(1)

        return report

    def check_isomers(self, errors: bool) -> IsomerReport:
        """Report (and warn, or abort) species sharing the same atomic composition.

        Isomers are detected by comparing ``Specie.exploded`` tuples.  For
        example, HCO+ and HOC+ both have ``["C", "H", "O", "+"]`` and would
        be reported here.

        Parameters
        ----------
        errors : bool
            If ``True``, call ``sys.exit(1)`` when isomers are found.

        Returns
        -------
        IsomerReport
            Groups of two or more species names sharing a composition.
        """
        buckets: dict[tuple, list[str]] = {}

        for sp in self.species:
            key = tuple(sp.exploded)
            buckets.setdefault(key, []).append(sp.name)

        groups = tuple(tuple(names) for names in buckets.values() if len(names) > 1)

        for names in groups:
            markup = ", ".join(f"[cyan]{name}[/]" for name in names)
            self.logger.warning(f"Isomers detected: {markup}")

        report = IsomerReport(groups=groups)

        if not report.ok and errors:
            self.logger.error("Isomer errors found")
            sys.exit(1)

        return report

    def __generate_reaction_matrices(self) -> None:
        """Build integer stoichiometry matrices: shape (n_reactions × n_species)."""
        self.reactant_matrix = np.zeros(
            (self.reactions.count, self.species.count), dtype=int
        )
        self.product_matrix = np.zeros(
            (self.reactions.count, self.species.count), dtype=int
        )

        for i, reaction in enumerate(self.reactions):
            for reactant in reaction.reactants.core:
                self.reactant_matrix[i, reactant.index] += 1

            for product in reaction.products.core:
                self.product_matrix[i, product.index] += 1

    def sfluxes(self) -> list[Expr]:
        """Return symbolic flux expressions for all reactions.

        The flux of reaction *i* is ``rate_i * nden[r1] * nden[r2] * ...``,
        where the product runs over all reactant species.

        Returns
        -------
        list[Expr]
            One SymPy expression per reaction, in reaction-index order.
        """
        return get_sfluxes(self.reactions, self.species, self.symbols.ndens)

    def sodes(self) -> list[Basic]:
        """Return symbolic ODE right-hand sides for all species.

        Each entry is the net rate of change of a species number density
        (cm⁻³ s⁻¹), formed by summing flux contributions over all reactions
        in which the species participates as a reactant or product.

        Returns
        -------
        list[Basic]
            One SymPy expression per species, in species-index order.
        """
        return get_sodes(self.reactions, self.species, self.symbols.ndens)

    def sradodes(self, order: int = 0) -> list[Expr]:
        """Return symbolic radiation moment ODE right-hand sides.

        Parameters
        ----------
        order : int, optional
            Radiation moment order (0 = number/energy density, 1 = flux),
            by default ``0``.

        Returns
        -------
        list[Expr]
            One SymPy expression per radiation band.
        """
        return get_sradodes(self, order)

    def to_hdf5(
        self,
        fname: str | Path,
        label: str | None = None,
        T_min=None,
        T_max=None,
        nT=64,
        err_tol=0.01,
        rate_min=1e-30,
        rate_max=1e100,
        fast_log=False,
        include_all=False,
        verbose=False,
    ) -> None:
        """Write a pre-tabulated rate coefficient table to an HDF5 file.

        Parameters
        ----------
        fname : str | Path
            Output file path.  The ``.hdf5`` extension is recommended.
        label : str | None, optional
            Dataset label stored inside the file.  Defaults to
            ``self.label``.
        T_min : float | None, optional
            Minimum temperature for the tabulation grid (K).
        T_max : float | None, optional
            Maximum temperature for the tabulation grid (K).
        nT : int, optional
            Number of temperature grid points, by default 64.
        err_tol : float, optional
            Maximum allowed relative error in rate interpolation, default
            0.01 (1 %).
        rate_min : float, optional
            Rates below this threshold are clamped to ``rate_min``,
            default ``1e-30``.
        rate_max : float, optional
            Rates above this threshold are clamped to ``rate_max``,
            default ``1e100``.
        fast_log : bool, optional
            Use a fast logarithm approximation, default ``False``.
        include_all : bool, optional
            Include all reactions, even those with temperature bounds,
            default ``False``.
        verbose : bool, optional
            Print per-reaction tabulation progress, default ``False``.
        """
        if isinstance(fname, str):
            fname = Path(fname)

        if fname.suffix not in [".hdf5", ".hdf"]:
            fname = fname.with_suffix(".hdf5")

        write_data_table(
            reactions=self.reactions,
            logger=self.logger,
            fname=fname,
            label=label or self.label,
            T_min=T_min,
            T_max=T_max,
            nT=nT,
            err_tol=err_tol,
            rate_min=rate_min,
            rate_max=rate_max,
            fast_log=fast_log,
            format="hdf5",
            include_all=include_all,
            verbose=verbose,
        )

    def to_txt(
        self,
        fname: str | Path,
        label: str | None = None,
        T_min=None,
        T_max=None,
        nT=64,
        err_tol=0.01,
        rate_min=1e-30,
        rate_max=1e100,
        fast_log=False,
        include_all=False,
        verbose=False,
    ) -> None:
        """Write a pre-tabulated rate coefficient table to a plain-text file.

        Parameters are identical to ``to_hdf5``, except the output format is
        a whitespace-separated text table.

        Parameters
        ----------
        fname : str | Path
            Output file path.  The ``.txt`` extension is recommended.
        label : str | None, optional
            Dataset label stored in the file header.
        T_min : float | None, optional
            Minimum temperature (K).
        T_max : float | None, optional
            Maximum temperature (K).
        nT : int, optional
            Number of temperature grid points, default 64.
        err_tol : float, optional
            Maximum allowed relative interpolation error, default 0.01.
        rate_min : float, optional
            Lower rate clamp, default ``1e-30``.
        rate_max : float, optional
            Upper rate clamp, default ``1e100``.
        fast_log : bool, optional
            Use fast logarithm approximation, default ``False``.
        include_all : bool, optional
            Include temperature-bounded reactions, default ``False``.
        verbose : bool, optional
            Print tabulation progress, default ``False``.
        """
        if isinstance(fname, str):
            fname = Path(fname)

        if fname.suffix != ".txt":
            fname = fname.with_suffix(".txt")

        write_data_table(
            reactions=self.reactions,
            logger=self.logger,
            fname=fname,
            label=label or self.label,
            T_min=T_min,
            T_max=T_max,
            nT=nT,
            err_tol=err_tol,
            rate_min=rate_min,
            rate_max=rate_max,
            fast_log=fast_log,
            format="txt",
            include_all=include_all,
            verbose=verbose,
        )

    def write_table(
        self,
        fname: str | Path,
        label: str | None = None,
        T_min=None,
        T_max=None,
        nT=64,
        err_tol=0.01,
        rate_min=1e-30,
        rate_max=1e100,
        fast_log=False,
        format="auto",
        include_all=False,
        verbose=False,
    ) -> None:
        """Write a pre-tabulated rate coefficient table in the specified format.

        This is the general-purpose version of ``to_hdf5`` / ``to_txt``.

        Parameters
        ----------
        fname : str | Path
            Output file path.
        label : str | None, optional
            Dataset label.  Defaults to ``self.label``.
        T_min : float | None, optional
            Minimum temperature (K).
        T_max : float | None, optional
            Maximum temperature (K).
        nT : int, optional
            Number of temperature grid points, default 64.
        err_tol : float, optional
            Maximum allowed relative interpolation error, default 0.01.
        rate_min : float, optional
            Lower rate clamp, default ``1e-30``.
        rate_max : float, optional
            Upper rate clamp, default ``1e100``.
        fast_log : bool, optional
            Use fast logarithm approximation, default ``False``.
        format : str, optional
            Output format: ``"hdf5"``, ``"txt"``, or ``"auto"`` (inferred
            from the file extension), default ``"auto"``.
        include_all : bool, optional
            Include temperature-bounded reactions, default ``False``.
        verbose : bool, optional
            Print tabulation progress, default ``False``.
        """
        write_data_table(
            reactions=self.reactions,
            logger=self.logger,
            fname=fname,
            label=label or self.label,
            T_min=T_min,
            T_max=T_max,
            nT=nT,
            err_tol=err_tol,
            rate_min=rate_min,
            rate_max=rate_max,
            fast_log=fast_log,
            format=format,
            include_all=include_all,
            verbose=verbose,
        )
