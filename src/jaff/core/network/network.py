# ABOUTME: Network class: loads, validates and queries a chemical reaction network
# ABOUTME: Builds Species/Reactions catalogues and symbolic ODE/flux expressions

"""Top-level ``Network`` class for loading and querying JAFF reaction networks.

A ``Network`` is the main entry point for users of the JAFF library.  It
reads a reaction network file in any supported format (KROME, PRIZMO, UDFA,
KIDA, UCLChem), builds typed
``Species`` and ``Reactions`` catalogues, validates the network, and exposes
symbolic ODE and flux expressions for downstream code generation.

Auxiliary ``.jfunc`` files may accompany the network file.  They can supply
custom rate coefficients, chemical heating/cooling rates, and radiation
moment contributions using the ``@var`` / ``@function`` syntax parsed by
``AuxiliaryFunctionParser``.

Number densities are represented as ``nden[i]`` — an entry of the SymPy
``IndexedBase("nden", shape=(n_species,))`` returned by ``net.symbols.ndens``,
where ``i`` is the integer index of the species within the network.
"""

from __future__ import annotations

import logging
import sys
from dataclasses import dataclass
from functools import cached_property, lru_cache
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, Mapping

import numpy as np
from sympy import (
    Basic,
    Expr,
    parse_expr,
    symbols,
)

from ...common import is_jaff_file, print_motd, resolve_dependencies
from ...config.base import as_layer
from ...drivers import HDF5
from ...errors import ParserError
from ...io import JaffLogger, jaff_progress
from ...io._io import to_jaff_file, write_data_table
from ...physics import (
    Dust,
    Radiation,
    Thermodynamics,
)
from ...physics.dust.config import DustConfig
from ...physics.photo_reactions.config import RadiationConfig
from ...physics.thermodynamics.config import ThermodynamicsConfig
from ..chemistry import element_table
from ..elements import Elements
from ..parsers import NetworkParser
from ..parsers.network._formats import ParsedRecord
from ..reaction import (
    RateSegment,
    Reaction,
    Reactions,
)
from ..reaction.config import ReactionsConfig
from ..species import Specie, Species
from ._functions import NetworkFunctions
from ._spec import NetworkSpec
from ._symbols import NetworkSymbols
from .config import NetworkConfig, NetworkConfigSource, load_network_config

if TYPE_CHECKING:
    from .._typing import ElementProps
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

    The keywords are the keys of the ``[network]`` table, so the layer is that
    table's keys with the arguments given (not ``None``); an omitted one never
    overrides the ``jaff.toml`` or ``config``.  A section given as a config
    model becomes the dict of its set keys (:func:`~jaff.config.base.as_layer`).
    """
    return {key: as_layer(value) for key, value in kwargs.items() if value is not None}


@dataclass(frozen=True)
class ConservationReport:
    """Structured result of :meth:`Network.check_conservation`.

    Attributes
    ----------
    mass : tuple[str, ...]
        Reactions (``verbatim`` form, network order) whose reactant and
        product masses differ by more than one electron mass.
    charge : tuple[str, ...]
        Reactions (``verbatim`` form, network order) whose net charge changes.
    """

    mass: tuple[str, ...]
    charge: tuple[str, ...]

    @property
    def ok(self) -> bool:
        """``True`` when every reaction conserves mass and charge."""
        return not self.mass and not self.charge


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

    ``Network`` is the main entry point of JAFF.  It reads a reaction network
    file, auto-detects its format (KIDA, UDFA, PRIZMO, KROME, UCLCHEM), validates mass and charge
    conservation, and assembles the species and reaction catalogues together
    with the stoichiometry matrices.  It also handles the optional radiation
    field, photochemistry cross-sections, dust, and the auxiliary ``.jfunc``
    functions file.

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
        Unique elements derived from all species; used for composition
        matrices.
    reactant_matrix : np.ndarray
        Integer stoichiometry matrix, shape ``(n_reactions, n_species)``.
        Entry ``[i, j]`` counts how many times species *j* appears as a
        reactant in reaction *i*.
    product_matrix : np.ndarray
        Integer stoichiometry matrix for products, same shape.
    thermodynamics : Thermodynamics
        Equation of state and thermal rates, built lazily: ``eos`` (internal
        energy for the configured EOS, an ideal gas by default),
        ``dEdt_chemical`` (accumulated over all reactions), ``dEdt_extra``
        (from the ``heatingcoolingrate`` auxiliary function, if present) and
        ``dEdt_tot``, all in erg cm⁻³ s⁻¹, plus the matching temperature rates
        ``dTdt_chemical``, ``dTdt_extra`` and ``dTdt_tot`` in K s⁻¹.
    symbols : NetworkSymbols
        Canonical symbols, density expressions (``ndens``, ``ntot``, ``rho``,
        ``n_hnuc``), introspection sets and symbol standardization for this
        network.
    radiation : Radiation or None
        Radiation field object; ``None`` when no ``[network.radiation]`` table
        (in the config or the ``radiation`` argument) is given.
    dust : Dust or None
        Dust module; ``None`` when no ``[network.dust]`` table (in the config or
        the ``dust`` argument) is given.
    mass_dict : dict[str, ElementProps]
        The element table (:func:`~jaff.core.chemistry.element_table`) the
        species and elements are built from; shared by every network.
    spec : NetworkSpec
        Facts about the source files (not settings): the resolved ``fname``,
        the ``.jfunc`` file used (``jfunc_path``).
    functions : NetworkFunctions
        The ``.jfunc`` file's functions, standardized on use.
    config : NetworkConfig
        The validated config: the ``jaff.toml`` next to the network, the
        ``config`` layers, then the keyword arguments.
    logger : logging.Logger
        The JAFF logger this network reports to.

    Examples
    --------
    Load a predefined network by name and inspect it:

    >>> from jaff import Network
    >>> net = Network("COthin")
    >>> net.label
    'react_COthin'
    >>> net.species.count, net.reactions.count
    (37, 238)
    >>> net.reactant_matrix.shape
    (238, 37)

    ``config`` and the keyword arguments are merged over the ``jaff.toml``:

    >>> net = Network(
    ...     "h_photoionization",
    ...     config={"label": "hphoto"},
    ...     radiation={"bands": [13.6, 20.0, 50.0]},
    ... )
    >>> net.config.label, net.radiation.nbands
    ('hphoto', 2)
    """

    def __init__(
        self,
        fname: str | Path,
        config: NetworkConfigSource | None = None,
        strict: bool | None = None,
        label: str | None = None,
        funcfile: bool | str | Path | None = None,
        expand_nuclei: bool | None = None,
        reactions: ReactionsConfig | Mapping[str, Any] | None = None,
        radiation: RadiationConfig | Mapping[str, Any] | None = None,
        dust: DustConfig | Mapping[str, Any] | None = None,
        thermodynamics: ThermodynamicsConfig | Mapping[str, Any] | None = None,
    ):
        """Load a reaction network from *fname*.

        Parameters
        ----------
        fname : str | Path
            Path to the network file, or the name of a predefined network (a
            sub-directory of :data:`~jaff.config.NETWORKS_DIR` containing a
            single ``.jet`` file).  A predefined network name wins over a
            same-named path on disk.  Supported formats: any text format
            auto-detected by ``NetworkParser`` (KIDA, UDFA, PRIZMO, KROME,
            UCLCHEM).  ``.jaff`` files written by :meth:`to_jaff` cannot be
            loaded at the moment (``ParserError``).
        config : str | Path | dict | NetworkConfig | None, optional
            Settings merged over the network's own defaults, the ``jaff.toml``
            next to the network file (if any): a TOML file (its ``[network]``
            table), a dict of the ``[network]`` table's keys
            (``{"strict": True, "reactions": {...}}``) or a
            :class:`~jaff.core.network.config.NetworkConfig` (only the keys
            explicitly set on it count); see
            :func:`~jaff.core.network.config.load_network_config`.  The keyword
            arguments below are the same ``[network]`` keys, named and nested
            as in the config file; given ones (not ``None``) are merged over
            *config*.  Tables merge key by key, except an ``eos`` or a
            shielding table, which replaces the lower one whole.
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
        expand_nuclei : bool | None, optional
            When ``True`` (default), an ``n_<element>_nuc`` symbol (e.g.
            ``n_H_nuc``, ``n_He_nuc``) in rate expressions is expanded to the
            element-nucleus sum — ``Σ atom-count · nden[i]`` over all species
            bearing that element.  Set to ``False`` to keep it as a free symbol
            ``n<element>_nuc`` (e.g. ``nh_nuc``).  (Plain ``n_X`` always resolves
            to species ``X``; only ``n_X_nuc`` is affected by this flag.)
        reactions : ReactionsConfig | dict | None, optional
            The ``[network.reactions]`` table: e.g. ``duplicate_policy``,
            ``t_cutoff``, ``photo`` (``pi_database``, ``use_proxy_photoreaction``,
            ...) and per-reaction ``reaction`` settings, with the keys of the
            config file.
        radiation : RadiationConfig | dict | None, optional
            The ``[network.radiation]`` table; enables the network's
            :class:`Radiation` and photochemistry.  Without it (here or in the
            config), no radiation field is built and photochemistry is disabled.
        dust : DustConfig | dict | None, optional
            The ``[network.dust]`` table; enables the network's :class:`Dust`.
            Without it (here or in the config), the dust module is disabled.
        thermodynamics : ThermodynamicsConfig | dict | None, optional
            The ``[network.thermodynamics]`` table, e.g. ``{"eos": {"type":
            "ideal", "gamma": 1.4}}``.  An ``eos`` given here replaces the
            config's whole; without one, an ideal gas with the default adiabatic
            index is used.

        Raises
        ------
        FileNotFoundError
            If *fname* is neither an existing file nor a predefined network
            name, or a predefined network directory has no ``.jet`` file.
        ParserError
            If *fname* is a ``.jaff`` file, a predefined network directory has
            more than one ``.jet`` file, the config is invalid (including a per-reaction config entry that
            names no reaction of the network, or photo settings on a non-photo
            reaction), a duplicate rate is found under ``duplicate_policy =
            "error"``, or the network or ``.jfunc`` file is
            malformed.
        TypeError
            If *config* is not a path, dict or ``NetworkConfig``.
        SystemExit
            If *strict* is enabled and a reaction does not conserve mass or
            charge, or :meth:`check_sink_sources`, :meth:`check_recombinations`
            or :meth:`check_isomers` finds a problem.
        """
        self.logger: logging.Logger = JaffLogger().get_logger()
        if is_jaff_file(Path(fname)):
            raise ParserError(
                f"{fname}: loading .jaff files is not supported; "
                "load the original network file"
            )

        overrides = _keyword_layer(
            strict=strict,
            label=label,
            funcfile=funcfile,
            expand_nuclei=expand_nuclei,
            reactions=reactions,
            radiation=radiation,
            dust=dust,
            thermodynamics=thermodynamics,
        )
        self.config: NetworkConfig = load_network_config(
            Path(fname),
            config,
            overrides=overrides,
        )
        self.spec: NetworkSpec = NetworkSpec(fname, self.config.funcfile)

        print_motd()

        self.species: Species = Species()
        self.reactions: Reactions = Reactions(self.config.reactions)
        self.symbols: NetworkSymbols = NetworkSymbols(self)
        self.functions: NetworkFunctions = NetworkFunctions(self, self.spec.jfunc_path)
        self.reactant_matrix: np.ndarray | None = None
        self.product_matrix: np.ndarray | None = None
        self.radiation: Radiation | None = (
            Radiation(self, self.config.radiation)
            if self.config.radiation is not None
            else None
        )
        self.dust: Dust | None = (
            Dust(self, self.config.dust) if self.config.dust is not None else None
        )

        self.logger.info(f"Loading network from {self.spec.fname}")
        self.logger.info(f"Network label: [yellow]{self.label}[/]")

        parsers = self.__load_network()

        self.reactions.check_config()
        self.__post_process(parsers)
        self.reactions.link()
        if self.radiation is not None:
            self.radiation.compute()
        self.thermodynamics: Thermodynamics = Thermodynamics(
            self, self.config.thermodynamics
        )
        self.symbols.log_summary()

        self.check_conservation(self.config.strict)
        self.check_sink_sources(self.config.strict)
        self.check_recombinations(self.config.strict)
        self.check_isomers(self.config.strict)

        self.__generate_reaction_matrices()

        self.logger.info("[green]Network loaded successfully![/]")

    @cached_property
    def filename(self) -> Path:
        """Absolute path to the source network file (from :attr:`spec`)."""
        return self.spec.fname

    @property
    def mass_dict(self) -> dict[str, ElementProps]:
        """The element table (:func:`~jaff.core.chemistry.element_table`)."""
        return element_table()

    @property
    def label(self) -> str:
        """Human-readable label: the ``label`` setting, else the file stem."""
        return self.config.label or self.spec.fname.stem

    def __load_network(self) -> list[Parser]:
        """Parse the network file, then build the species and the reactions.

        The species catalogue and :attr:`elements` are complete before the
        first reaction is built.

        Returns
        -------
        list[Parser]
            The parsers that produced reactions, for :meth:`__post_process`.
        """
        with NetworkParser(self.spec.fname, self.logger) as netp:
            rows, global_vars = netp.get_parsed()
            parsers = netp.get_parsers()

        special_species = self.__build_species(rows)
        self.elements: Elements = Elements(self.species._list)
        self.__build_reactions(rows, global_vars, special_species)

        self.logger.info(f"Loaded {self.reactions.count} reactions")
        self.logger.info(
            f"Loaded {len(self.reactions.with_type('photo'))} photo-chemistry reactions"
        )

        return parsers

    def __build_species(self, rows: list[ParsedRecord]) -> dict[str, Specie]:
        """Add every species named in *rows* to :attr:`species`, in file order.

        Pseudo-species (names starting with ``_``, e.g. ``_PHOTON``) are not
        added to the catalogue.

        Parameters
        ----------
        rows : list[ParsedRecord]
            The parsed reaction rows.

        Returns
        -------
        dict[str, Specie]
            The pseudo-species, keyed by name.
        """
        seen: set[str] = set()
        special_species: dict[str, Specie] = {}

        for row in rows:
            for name in row.reactants + row.products:
                if name in seen:
                    continue
                seen.add(name)
                if name.startswith("_"):
                    special_species[name] = Specie(name, -1)
                else:
                    self.species.add(Specie(name, self.species.count))

        return special_species

    def __build_reactions(
        self,
        rows: list[ParsedRecord],
        global_vars: dict[str, Basic],
        special_species: dict[str, Specie],
    ) -> None:
        """Build a reaction from each row, merging rows of the same reaction.

        Each row's rate is resolved (global variables, ``.jfunc`` calls) and
        standardized before its reaction is built; ``rc_<N>`` references are
        left for :meth:`Reactions.link`.

        Parameters
        ----------
        rows : list[ParsedRecord]
            The parsed reaction rows.
        global_vars : dict[str, Basic]
            The network file's global variables.
        special_species : dict[str, Specie]
            The pseudo-species, keyed by name (see :meth:`__build_species`).
        """
        functions = self.functions

        resolved_vars: dict[str, Expr] = {
            var: resolve_dependencies(expr, {}, functions)
            for var, expr in global_vars.items()
        }
        subs_dict: dict[Basic, Basic] = {
            symbols(var.lower()): expr for var, expr in resolved_vars.items()
        }

        for reaction in jaff_progress.track(
            rows,
            description=f"Creating {self.label} network",
        ):
            reactants: list[str] = reaction.reactants
            products: list[str] = reaction.products
            tmin: float | None = reaction.tmin
            tmax: float | None = reaction.tmax
            rate: str = reaction.rate
            si = reaction.source_index
            aux_chem_rate = f"chemrate{si}"

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

            rate_expr = self.__parse_rate(aux_chem_rate, rate, resolved_vars)
            rate_expr = resolve_dependencies(rate_expr, local_subs_dict, functions)
            rate_expr = self.symbols.standardize(rate_expr)

            rtype = reaction.type
            if (srxn, rtype) in self.reactions:
                self.reactions.add_rate(
                    self.reactions[srxn, rtype], RateSegment(rate_expr, tmin, tmax)
                )
                continue

            rea = Reaction(
                self,
                reactants=rr,
                products=pp,
                rate_coefficient=rate_expr,
                tmin=tmin,
                tmax=tmax,
                original_string=reaction.line,
                index=si,
                type=rtype,
                format=reaction.format,
                metadata=reaction.metadata,
            )
            self.reactions.add(rea)

    def __post_process(self, parsers: list[Parser]) -> None:
        """Freeze the catalogues and run each parser's ``post_process`` hook.

        Species and reactions are final from here on: hooks may change
        reactions but not add any (``add()`` raises once frozen).  Rates are
        already standardised, so a hook writing a shorthand symbol (``n_X``,
        ...) must standardise it itself (``network.symbols.standardize``).
        Runs before :meth:`Reactions.link`, so ``rc_<N>`` in what a
        hook writes is resolved like every other rate.

        Parameters
        ----------
        parsers : list[Parser]
            Parsers that produced reactions.
        """
        self.reactions._freeze()
        self.species._freeze()

        for parser in parsers:
            parser.post_process(self)

    def __parse_rate(
        self,
        aux_chem_rate: str,
        rate: str,
        global_vars: dict[str, Expr],
    ) -> Expr:
        """Convert a raw rate string to a SymPy expression.

        Checks, in priority order:

        1. Whether a ``.jfunc`` function named *aux_chem_rate* exists (custom
           rate, read with :meth:`NetworkFunctions.value`).
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
        global_vars : dict[str, Expr]
            Resolved global variable map from the network file.

        Returns
        -------
        Expr
            The rate expression.
        """
        if aux_chem_rate in self.functions:
            rate_expr = self.functions.value(aux_chem_rate)
        elif rate in global_vars:
            rate_expr = symbols(rate)
        else:
            rate_expr = _parse_rate_expr(rate)

        if not isinstance(rate_expr, Expr):
            raise ParserError(f"Rate expression is not an Expr: {rate_expr}")

        return rate_expr

    def write_shielding_tables(self, directory: str | Path) -> list[Path]:
        """Write the data tables of the reactions' shielding models into *directory*.

        Models such as Leiden shielding build a per-reaction table the
        generated code interpolates over (``photo.shielding_tables``).

        Parameters
        ----------
        directory : str | Path
            Output directory; must exist.

        Returns
        -------
        list[Path]
            The files written, in reaction order.
        """
        written: list[Path] = []
        for reaction in self.reactions:
            tables = reaction.photo.shielding_tables if reaction.photo is not None else {}
            for name, table in tables.items():
                path = Path(directory) / name
                HDF5().from_dict(str(path), table, mode="w")
                written.append(path)

        return written

    def to_jaff(self, filename: str | Path):
        """Serialise this network to a gzip-compressed ``.jaff`` JSON file.

        It uses the ``jaff.network_json`` format marker and
        ``schema_version = 1``.  Loading ``.jaff`` files is not supported at
        the moment.

        Parameters
        ----------
        filename : str | Path
            Output file path.  Unless it already ends with ``.jaff`` or
            ``.jaff.gz``, its suffix is replaced by ``.jaff`` (``net.txt``
            becomes ``net.jaff``).

        Raises
        ------
        ValueError
            If a reaction rate contains an undefined SymPy function, which
            cannot be serialised.
        TypeError
            If a rate is neither a SymPy expression, a ``str`` nor ``None``.
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
            unique reactions.  Other values log only the counts.

        Notes
        -----
        Reactions are matched by their serialized form, so the same reaction
        with a different mechanism or rate still counts as common.
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
            When ``1`` (default), print the common and unique species sets.
            Other values log only the counts.
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

    def check_conservation(self, strict: bool) -> ConservationReport:
        """Report (and warn, or abort) reactions that do not conserve mass or charge.

        Uses :meth:`~jaff.core.reaction.Reaction.check_mass` and
        :meth:`~jaff.core.reaction.Reaction.check_charge` on every reaction and
        logs a warning naming each violation.  Run automatically on
        construction with ``strict`` set to the network's ``strict`` setting.

        Parameters
        ----------
        strict : bool
            If ``True``, call ``sys.exit(1)`` when a violation is found.

        Returns
        -------
        ConservationReport
            The reactions that violate mass or charge conservation.

        Examples
        --------
        >>> net = Network("COthin")
        >>> net.check_conservation(strict=False).ok
        True
        """
        mass = [r for r in self.reactions if not r.check_mass()]
        charge = [r for r in self.reactions if not r.check_charge()]

        for r in mass:
            self.logger.warning(f"Mass not conserved in reaction {r.index}: {r.verbatim}")
        for r in charge:
            self.logger.warning(
                f"Charge not conserved in reaction {r.index}: {r.verbatim}"
            )

        report = ConservationReport(
            mass=tuple(r.verbatim for r in mass),
            charge=tuple(r.verbatim for r in charge),
        )

        if not report.ok and strict:
            self.logger.error("Conservation errors found")
            sys.exit(1)

        return report

    def check_sink_sources(self, strict: bool) -> SinkSourceReport:
        """Report (and warn, or abort) species never produced or never consumed.

        A *sink* species appears as a reactant in at least one reaction but
        is never produced.  A *source* species is produced but never consumed.
        Special species (names starting with ``_``, such as ``_DUMMY``) are not
        network species, so they never appear in the report.  Each sink
        and source is logged, followed by a warning.  Run automatically on
        construction with ``strict`` set to the network's ``strict`` setting.

        Parameters
        ----------
        strict : bool
            If ``True``, call ``sys.exit()`` (exit status 0) when sinks or
            sources are found.

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

        if not report.ok and strict:
            self.logger.error("Exiting since strict is enabled")
            sys.exit()

        return report

    def check_recombinations(self, strict: bool) -> RecombinationReport:
        """Report (and warn, or abort) cations lacking an electron recombination.

        Every positively charged species should be a reactant, alongside
        ``e-``, in at least one reaction; a warning is logged for each one that
        is not.  Run automatically on construction with ``strict`` set to the
        network's ``strict`` setting.

        Parameters
        ----------
        strict : bool
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

        if not report.ok and strict:
            self.logger.error("Recombination errors found")
            sys.exit(1)

        return report

    def check_isomers(self, strict: bool) -> IsomerReport:
        """Report (and warn, or abort) species sharing the same atomic composition.

        Isomers are detected by comparing ``Specie.exploded`` tuples.  For
        example, HCO+ and HOC+ both have ``["+", "C", "H", "O"]`` and would
        be reported here.  A warning is logged for each group.  Run
        automatically on construction with ``strict`` set to the network's
        ``strict`` setting.

        Parameters
        ----------
        strict : bool
            If ``True``, call ``sys.exit(1)`` when isomers are found.

        Returns
        -------
        IsomerReport
            Groups of two or more species names sharing a composition.

        Examples
        --------
        >>> net = Network("COthin")
        >>> net.check_isomers(strict=False).groups
        (('HOC+', 'HCO+'),)
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

        if not report.ok and strict:
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
        where the product runs over all reactant species.  A reactant that
        appears several times enters with its stoichiometric power:

        ```
        F_i = k_i * prod_j nden[R_j] ** alpha_j
        ```

        with ``k_i`` the rate coefficient, ``R_j`` the core reactants and
        ``alpha_j`` their stoichiometric coefficients.

        Returns
        -------
        list[Expr]
            One SymPy expression per reaction, in reaction-index order.

        Examples
        --------
        >>> net = Network("COthin")
        >>> str(net.reactions[0])
        'H + e- -> H+ + e- + e-'
        >>> net.sfluxes()[0].has(net.symbols.ndens[0] * net.symbols.ndens[1])
        True
        """
        return list(self.reactions.rates())

    def sodes(self) -> list[Basic]:
        """Return symbolic ODE right-hand sides for all species.

        Each entry is the net rate of change of a species number density
        (cm⁻³ s⁻¹), formed by summing flux contributions over all reactions
        in which the species participates as a reactant or product:

        ```
        dn_s/dt = sum_{i: s in products(i)} F_i - sum_{i: s in reactants(i)} F_i
        ```

        where ``F_i`` are the fluxes from :meth:`sfluxes`; see
        :meth:`Reactions.odes`.  A species appearing
        several times on one side contributes the flux once per occurrence.

        Returns
        -------
        list[Basic]
            One SymPy expression per species, in species-index order.
        """
        return self.reactions.odes(self.species)

    def sradodes(self, order: int = 0) -> list[Expr]:
        """Return symbolic radiation moment ODE right-hand sides.

        Each radiation band has two moments: the energy or photon density
        ``rd_i`` (depending on ``radiation.config.mode``) and the flux ``f_i``.  Each
        density equation removes what photo reactions absorb and adds what
        reactions with a ``spectralYield<N>`` auxiliary function emit, plus the
        network ``globalEmissivity``, and removes what the network
        ``globalAbsorptionCoefficient`` absorbs.  The flux equations carry the
        photo-reaction absorption and, unless
        ``network.radiation.absorption_attenuates_flux`` is off, the
        ``globalAbsorptionCoefficient`` absorption too.  With dust enabled the
        dust reductions are applied to both.

        Parameters
        ----------
        order : int, optional
            Layout of the returned list, one of ``0`` (default), ``1``, ``2``
            or ``3``:

            - ``0``: ``[rd_0, f_0, rd_1, f_1, ...]``
            - ``1``: ``[f_0, rd_0, f_1, rd_1, ...]``
            - ``2``: ``[rd_0, rd_1, ..., f_0, f_1, ...]``
            - ``3``: ``[f_0, f_1, ..., rd_0, rd_1, ...]``

        Returns
        -------
        list[Expr]
            Two SymPy expressions (density and flux) per radiation band, so
            ``2 * radiation.nbands`` in total, laid out as *order* says.

        Raises
        ------
        RuntimeError
            If the network has no radiation field.
        ValueError
            If *order* is not one of ``0``, ``1``, ``2``, ``3``.

        Examples
        --------
        >>> net = Network("h_photoionization", radiation={"bands": [13.6, 20.0]})
        >>> len(net.sradodes(order=2))
        2
        """
        if self.radiation is None:
            raise RuntimeError(
                "No radiation bands found. Radiation odes cannot be generated"
            )

        return self.radiation.odes(order)

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

        The table is tabulated on a temperature grid and written in the layout
        read by the quokka hydro code.

        Parameters
        ----------
        fname : str | Path
            Output file path.  Unless it ends with ``.hdf5`` or ``.hdf``, its
            suffix is replaced by ``.hdf5``.
        label : str | None, optional
            Network label.  Unused in HDF5 output (see :meth:`to_txt`).
        T_min : float | None, optional
            Minimum temperature of the grid (K).  Defaults to the smallest
            ``reaction.tmin`` of the network.
        T_max : float | None, optional
            Maximum temperature of the grid (K).  Defaults to the largest
            ``reaction.tmax`` of the network.
        nT : int, optional
            Initial number of temperature grid points, default 64.  The
            adaptive refinement may increase it significantly.
        err_tol : float | None, optional
            Maximum permitted relative interpolation error, default 0.01
            (1 %).  ``None`` skips the adaptive refinement and returns exactly
            ``nT`` points.
        rate_min : float, optional
            Small positive floor added to the relative-error denominator, so
            near-zero rates do not divide by zero, default ``1e-30``.  Rates
            are not clamped to it.
        rate_max : float, optional
            Upper clamp applied to each rate before storing, preventing
            overflow in the log-rate arithmetic, default ``1e100``.
        fast_log : bool, optional
            Make the temperature grid uniform in ``fast_log2`` space rather
            than ``log10`` space, default ``False``.
        include_all : bool, optional
            Include every reaction, even those whose tabulated rate contains
            ``NaN`` or is constant, default ``False``.
        verbose : bool, optional
            Log the grid size and worst-case error after each refinement step,
            default ``False``.

        Notes
        -----
        Only rates that depend solely on ``tgas`` can be tabulated (``av`` is
        set to ``0`` and ``crate`` to ``1`` first); other rates come out as
        ``NaN`` and are dropped unless *include_all* is set.  Starting from
        ``nT`` points, the grid is doubled by inserting midpoints until the
        maximum relative interpolation error falls below *err_tol*.

        The file holds one ``reaction_coeff`` group with the attributes
        ``input_names``, ``input_units``, ``xlo``, ``xhi`` and ``spacing``
        and the datasets ``output_names``, ``output_units`` and ``data``.
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
        a whitespace-separated text table compatible with quokka.

        Parameters
        ----------
        fname : str | Path
            Output file path.  Unless it ends with ``.txt``, its suffix is
            replaced by ``.txt``.
        label : str | None, optional
            Network label written into the file header.  Defaults to
            the network's ``label``.
        T_min : float | None, optional
            Minimum temperature of the grid (K).  Defaults to the smallest
            ``reaction.tmin`` of the network.
        T_max : float | None, optional
            Maximum temperature of the grid (K).  Defaults to the largest
            ``reaction.tmax`` of the network.
        nT : int, optional
            Initial number of temperature grid points, default 64.  The
            adaptive refinement may increase it significantly.
        err_tol : float | None, optional
            Maximum permitted relative interpolation error, default 0.01
            (1 %).  ``None`` skips the adaptive refinement and returns exactly
            ``nT`` points.
        rate_min : float, optional
            Small positive floor added to the relative-error denominator, so
            near-zero rates do not divide by zero, default ``1e-30``.  Rates
            are not clamped to it.
        rate_max : float, optional
            Upper clamp applied to each rate before storing, preventing
            overflow in the log-rate arithmetic, default ``1e100``.
        fast_log : bool, optional
            Make the temperature grid uniform in ``fast_log2`` space rather
            than ``log10`` space, default ``False``.
        include_all : bool, optional
            Include every reaction, even those whose tabulated rate contains
            ``NaN`` or is constant, default ``False``.
        verbose : bool, optional
            Log the grid size and worst-case error after each refinement step,
            default ``False``.

        Notes
        -----
        Only rates that depend solely on ``tgas`` can be tabulated (``av`` is
        set to ``0`` and ``crate`` to ``1`` first); other rates come out as
        ``NaN`` and are dropped unless *include_all* is set.  Starting from
        ``nT`` points, the grid is doubled by inserting midpoints until the
        maximum relative interpolation error falls below *err_tol*.

        After a ``#``-commented header naming the network and listing the
        tabulated reactions, the file has the quokka 1-D table layout:

        ```
        1              # table dimensionality
        <nReact>       # number of outputs per entry
        2 or 3         # axis spacing (2 = log10, 3 = fast_log)
        <nTemp>        # number of temperature points
        <T_min> <T_max>
        <coeff row 0>
        ...
        ```
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
        Unlike them it does not change the suffix of *fname*.

        Parameters
        ----------
        fname : str | Path
            Output file path.
        label : str | None, optional
            Network label written into a text file header.  Defaults to
            the network's ``label``; unused in HDF5 output.
        T_min : float | None, optional
            Minimum temperature of the grid (K).  Defaults to the smallest
            ``reaction.tmin`` of the network.
        T_max : float | None, optional
            Maximum temperature of the grid (K).  Defaults to the largest
            ``reaction.tmax`` of the network.
        nT : int, optional
            Initial number of temperature grid points, default 64.  The
            adaptive refinement may increase it significantly.
        err_tol : float | None, optional
            Maximum permitted relative interpolation error, default 0.01
            (1 %).  ``None`` skips the adaptive refinement and returns exactly
            ``nT`` points.
        rate_min : float, optional
            Small positive floor added to the relative-error denominator, so
            near-zero rates do not divide by zero, default ``1e-30``.  Rates
            are not clamped to it.
        rate_max : float, optional
            Upper clamp applied to each rate before storing, preventing
            overflow in the log-rate arithmetic, default ``1e100``.
        fast_log : bool, optional
            Make the temperature grid uniform in ``fast_log2`` space rather
            than ``log10`` space, default ``False``.
        format : str, optional
            Output format: ``"hdf5"``, ``"txt"``, or ``"auto"`` (inferred
            from the file extension: ``.txt`` gives text, ``.hdf`` / ``.hdf5``
            give HDF5), default ``"auto"``.
        include_all : bool, optional
            Include every reaction, even those whose tabulated rate contains
            ``NaN`` or is constant, default ``False``.
        verbose : bool, optional
            Log the grid size and worst-case error after each refinement step,
            default ``False``.

        Raises
        ------
        ValueError
            If *format* is not one of the supported strings, or if
            ``format="auto"`` and the file extension is not recognised.

        Notes
        -----
        Only rates that depend solely on ``tgas`` can be tabulated (``av`` is
        set to ``0`` and ``crate`` to ``1`` first); other rates come out as
        ``NaN`` and are dropped unless *include_all* is set.  Starting from
        ``nT`` points, the grid is doubled by inserting midpoints until the
        maximum relative interpolation error falls below *err_tol*.
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
