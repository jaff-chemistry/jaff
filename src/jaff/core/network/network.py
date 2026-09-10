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
import re
import sys
from dataclasses import dataclass
from functools import cached_property, lru_cache, reduce
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from sympy import (
    Basic,
    Expr,
    Float,
    Function,
    Idx,
    MatrixSymbol,
    parse_expr,
    symbols,
)
from sympy.core.function import AppliedUndef, UndefinedFunction

from ...common import is_jaff_file, load_mass_dict, motd, resolve_dependencies
from ...errors import ParserError
from ...io import JaffLogger, jaff_progress
from ...io._io import JaffProps, from_jaff_file, to_jaff_file, write_data_table
from ...physics import (
    Dust,
    Photochemistry,
    Radiation,
    constants,
    get_eos,
    get_sfluxes,
    get_sodes,
    get_sradodes,
)
from ..elements import Elements
from ..parsers import NetworkParser
from ..reaction import RateSegment, RateSegments, Reaction, Reactions
from ..species import Specie, Species
from ._spec import NetworkSpec

if TYPE_CHECKING:
    from .._typing import ElementProps
    from ..parsers.auxiliary_func._typing import AuxiliaryFunctionsDict


@lru_cache(maxsize=200000)
def _parse_rate_expr(rate: str) -> Expr:
    """Parse a rate string into a (non-evaluated) SymPy expression, memoized.

    Large networks contain many reactions with identical rate strings (≈40% of
    KIDA-2024 rates repeat), and the parsed expression depends only on the
    string — species substitution happens later in ``_standardize_symbols`` —
    so results are cached across reactions and networks.  SymPy expressions are
    immutable, making the shared objects safe to reuse.
    """
    return parse_expr(rate, evaluate=False)


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
    dEdt_chem : Basic
        SymPy expression for the total chemical heating/cooling rate
        (erg cm⁻³ s⁻¹), accumulated over all reactions.
    dEdt_other : Basic
        Additional heating/cooling rate from the ``heatingcoolingrate``
        auxiliary function, if present.
    dRad_dt_extra : Basic
        Extra radiation moment source terms from ``@function`` definitions.
    radiation : Radiation | None
        Radiation field object; ``None`` when no radiation bands are specified.
    mass_dict : dict[str, ElementProps]
        Element mass dictionary used during species construction.
    spec : NetworkSpec
        The normalized construction parameters this network was built from.
    """

    #: Valid temperature cutoff behaviours for rate expressions.
    _valid_tcutoffs: list[str] = ["clip", "extrapolate"]

    _simple_map: dict[str, str] = {
        "nh0": "H",
        "nh2": "H2",
        "ne": "e-",
        "nhp": "H+",
    }

    _n_suffixes: dict[str, str] = {"p": "+", "m": "-", "0": ""}

    def __init__(
        self,
        fname: str | Path,
        config: str | Path | None = None,
        errors: bool = False,
        label: str | None = None,
        funcfile: bool | str | Path = True,
        duplicate_policy: str | None = None,  # preserve-first, preserve-last, error
        replace_nH: bool = True,
        rad_bands: list[str | int | float | Basic] = [],
        rad_powerlaw_index: int | float = 0,
        rad_energy_density: bool = False,
        dust: bool = False,
        background_field: str = "draine",
        c: float | str = constants.c.cgs.value,  # Speed of light in cgs unit
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
        errors : bool, optional
            If ``True``, treat conservation violations and duplicate reactions
            as fatal errors (process exits).  Default ``False`` (warnings
            only).
        config : str | Path | None, optional
            Path to a ``jaff.toml`` network config file.  When ``None``
            (default), JAFF looks for ``jaff.toml`` in the network file's
            directory.  Only the ``[network.rates]`` global and
            ``[network.reactions."<serialized>"].T_cutoff`` per-reaction
            temperature-cutoff keys are read from it.
        label : str | None, optional
            Human-readable name for this network.  Defaults to the file stem.
        funcfile : bool | str | Path, optional
            Path to a ``.jfunc`` auxiliary function file.  When ``True``
            (default), JAFF scans the network's directory for
            ``<network>.jfunc``.  Pass ``False`` to skip auxiliary-function
            loading entirely.
        duplicate_policy : str | None, optional
            How to resolve two rate coefficients sharing a reaction, mechanism,
            and temperature range: ``"preserve-first"`` (keep the first, drop
            later duplicates), ``"preserve-last"`` (keep the last), or
            ``"error"`` (raise a ``ParserError``).  When ``None`` (default),
            JAFF reads the network ``jaff.toml`` ``[network].duplicate_policy``
            key, falling back to ``"preserve-first"``.
        replace_nH : bool, optional
            When ``True`` (default), the shorthand symbol ``nh`` (and ``n_H``,
            ``n_He``) in rate expressions is expanded to a sum of
            ``nden[i]`` terms over all H-bearing (He-bearing) species.  Set
            to ``False`` to keep ``nh`` / ``nhe`` as free symbols.
        rad_bands : list, optional
            Radiation band boundaries used to construct the ``Radiation``
            object.  An empty list (default) disables radiation transport.
        rad_powerlaw_index : int | float, optional
            Power-law spectral index for the radiation field, default ``0``.
        rad_energy_density : bool, optional
            If ``True``, radiation moments are energy densities rather than
            number densities, default ``False``.
        c : float | str, optional
            Speed of light in CGS units (cm s⁻¹) or a symbol.  Defaults to
            ``constants.c.cgs.value``.
        _from_cli : bool, optional
            Internal flag: suppresses the MOTD banner when ``True``.

        Raises
        ------
        FileNotFoundError
            If *fname* is neither an existing file nor a predefined network
            name, or a predefined network directory has no ``.jet`` file.
        ParserError
            If a predefined network directory has more than one ``.jet`` file.
        """
        self.logger: logging.Logger = JaffLogger().get_logger()

        self.spec: NetworkSpec = NetworkSpec(
            fname,
            config,
            errors,
            label,
            funcfile,
            duplicate_policy,
            replace_nH,
            rad_bands,
            rad_powerlaw_index,
            rad_energy_density,
            c,
            _from_cli,
            _metadata,
        )

        jaff_props: JaffProps = {}  # type: ignore
        loaded_from_jaff_file = is_jaff_file(self.spec.fname)
        if loaded_from_jaff_file:
            jaff_props = from_jaff_file(self.spec.fname, errors)
            self.spec.fname = jaff_props.get("file_name", self.spec.fname)
        # Finalize the label: a ``.jaff`` payload label wins, else the supplied
        # label, else the file stem.
        self.spec.label = jaff_props.get("label", label or self.spec.fname.stem)
        if not _from_cli:
            print(motd())

        self.mass_dict: dict[str, ElementProps] = {}
        self.species: Species = Species()
        self.reactions: Reactions = Reactions()
        self.reactant_matrix: np.ndarray | None = None
        self.product_matrix: np.ndarray | None = None
        self.dEdt_chem: Basic = Float(0.0)
        self.dEdt_other: Basic = Float(0.0)
        self.dRad_dt_extra: Basic = Float(0.0)
        self._dust_enabled = dust
        self.radiation: Radiation | None = (
            Radiation(
                self,
                rad_bands,
                rad_powerlaw_index,
                rad_energy_density,
                c,
                background_field,
            )
            if len(rad_bands) > 0
            else None
        )
        self.__photochemistry: None | Photochemistry = None
        self.dust: Dust | None = Dust(self) if dust else None
        self.__element_sums: dict[str, Expr | None] = {}
        self.__tgas_clamp_cache: dict[tuple[float | None, float | None], Expr] = {}

        self.logger.info(f"Loading network from {self.spec.fname}")
        self.logger.info(f"Network label: [yellow]{self.label}[/]")

        self.mass_dict: dict[str, ElementProps] = load_mass_dict()
        Species.configure(self.mass_dict)

        if not loaded_from_jaff_file:
            self.__load_network()
        else:
            self.__load_network_from_jaff_file(jaff_props)

        self.__normalize_network_extras(replace_nH, loaded_from_jaff_file)

        self.check_sink_sources(errors)
        self.check_recombinations(errors)
        self.check_isomers(errors)

        self.__generate_reaction_matrices()

        self.elements: Elements = Elements(self.species._list)

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

    def __load_network(self):
        """Parse the network file and build species, reactions, and auxiliary quantities."""
        specie_names = set()
        special_species: dict[str, Specie] = {}
        free_symbols = set()
        undef_funcs = set()
        interp_funcs = set()

        n_photo = 0
        default_tcutoff: str = "clip"
        config = self.spec.config
        reactions_config: dict = config.get("reactions", {})
        reaction_props: dict = self.spec._metadata.get("reaction_props", {})
        jaff_global_tcutoff = config.get("rates", {}).get("T_cutoff")
        jaffgen_global_tcutoff = self.spec._metadata.get("rate_props", {}).get("T_cutoff")

        with NetworkParser(self.spec.fname, self.logger) as netp:
            reactions_list, global_vars = netp.get_parsed()

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
            aux_chem_rate = f"chemrate{i}"
            aux_delta_rad = f"deltarad{i}"
            aux_delta_e = f"deltae{i}"

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
            local_tcutoff: str = (
                reaction_props.get(srxn, {}).get("T_cutoff")
                or reactions_config.get(srxn, {}).get("T_cutoff")
                or jaffgen_global_tcutoff
                or jaff_global_tcutoff
                or default_tcutoff
            ).lower()
            if local_tcutoff not in self._valid_tcutoffs:
                raise ParserError(
                    f"Invalid temperature cutoff for reaction: {srxn}\n"
                    f"valid cutoffs are: {','.join(self._valid_tcutoffs)}"
                )

            rate_expr, n_photo = self.__parse_rate(
                aux_chem_rate, rate, aux_funcs, global_vars, n_photo
            )
            rate_expr = resolve_dependencies(rate_expr, local_subs_dict, aux_funcs)

            rtype = reaction.get("type", "unknown")
            if (srxn, rtype) in self.reactions:
                rea = self.reactions[srxn, rtype]
                if (tmin, tmax) in rea.rate_segments._by_prop:
                    # Same reaction, mechanism, and temperature range with a
                    # second rate coefficient — resolved per duplicate_policy.
                    existing = rea.rate_segments._by_prop[(tmin, tmax)]
                    if self.spec.duplicate_policy == "error":
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
                        if self.spec.duplicate_policy == "preserve-first"
                        else (rate_expr, existing.rate)
                    )
                    self.logger.warning(
                        f"Duplicate rate for [cyan]{rea.verbatim}[/] "
                        f"[type={rtype}] over T=({tmin}, {tmax}) resolved by "
                        f"policy '{self.spec.duplicate_policy}': kept {kept}, "
                        f"dropped {dropped}"
                    )
                    if self.spec.duplicate_policy == "preserve-first":
                        continue
                    # "preserve-last": fall through; add() replaces the segment.

                rea.rate_segments.add(RateSegment(rate_expr, tmin, tmax))
                continue

            # deltarad{i}: radiation energy emission per photon energy (eV)
            # per reaction added to the moment-0 equations at codegen time
            deltaRad: Basic = Float(0.0)
            if aux_delta_rad in aux_funcs:
                deltaRad = aux_funcs[aux_delta_rad]["def"]

            # deltae{i}: chemical energy change per reaction, accumulates into dEdt_chem
            deltaE: Basic = Float(0.0)
            if aux_delta_e in aux_funcs:
                deltaE = aux_funcs[aux_delta_e]["def"]

            for expr in [rate_expr, deltaE, deltaRad]:
                free_symbols |= self.free_symbols(expr)
                self.__detect_undefined_functions(expr, undef_funcs, interp_funcs)

            rea = Reaction(
                reactants=rr,
                products=pp,
                rate=rate_expr,
                tmin=tmin,
                tmax=tmax,
                dE=deltaE,
                dRad=deltaRad,
                original_string=reaction["string"],
                index=i,
                type=rtype,
                t_cutoff=local_tcutoff,
            )
            if "reaction_props" in self.spec._metadata:
                self.__parse_reaction_metadata(rea)

            self.reactions.add(rea)

            if rea.type == "photo":
                if self.__photochemistry is None:
                    self.__photochemistry = Photochemistry()

                rea.xsecs_dict = self.__photochemistry.get_xsec(rea)

            if rea.type == "photo" and self.radiation is not None:
                if aux_chem_rate not in aux_funcs:
                    self.radiation.set_reaction_rate_coefficient(rea)
                elif aux_chem_rate in aux_funcs and aux_delta_rad:
                    rea.custom_rad_rate = True
                    self.radiation.set_custom_rate(rea)
                else:
                    raise ParserError(
                        "If radiation is enabled and a custom rate is supplied\n"
                        "for a photo reaction, the auxilary deltaRad function is\n"
                        "necessary to weigh the first moment radiation equations\n"
                        f"Please add a custom deltaRad function for reaction {i}"
                    )

        if "heatingcoolingrate" in aux_funcs:
            self.dEdt_other = aux_funcs["heatingcoolingrate"]["def"]
            self.dEdt_other = self._standardize_symbols(
                self.dEdt_other, self.spec.replace_nH
            )
            free_symbols |= self.free_symbols(self.dEdt_other)
            self.__detect_undefined_functions(self.dEdt_other, undef_funcs, interp_funcs)

        self.logger.info(
            f"Variables found: {', '.join(sorted(f'[cyan]{s}[/]' for s in free_symbols))}"
        )
        self.logger.info(f"Loaded {self.reactions.count} reactions")
        self.logger.info(f"Loaded {n_photo} photo-chemistry reactions")

        if interp_funcs:
            self.logger.info(
                f"Found the following interpolation functions: {', '.join([f'[cyan]{func}[/]' for func in interp_funcs])}"
            )
        if undef_funcs:
            self.logger.warning(
                f"Found undefined functions {', '.join([f'[red]{func}[/]' for func in undef_funcs])}"
            )

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
            rea = Reaction(
                reactants=reaction["reactants"],
                products=reaction["products"],
                rate=reaction["rate"],
                dE=reaction["dE"],
                dRad=reaction["dRad"],
                tmin=reaction["tmin"],
                tmax=reaction["tmax"],
                t_cutoff=reaction.get("t_cutoff", "clip"),
                original_string=reaction["original_string"],
                index=i,
                type=reaction.get("reaction_type", "unknown"),
            )
            rea.custom_rad_rate = reaction["custom_rad_rate"]
            segments = reaction.get("rate_segments")
            if segments:
                rea.rate_segments = RateSegments(
                    [RateSegment(s["rate"], s["tmin"], s["tmax"]) for s in segments],
                    rea.t_cutoff,
                )
            self.reactions.add(rea)

            if rea.type == "photo":
                if self.__photochemistry is None:
                    self.__photochemistry = Photochemistry()
                rea.xsecs_dict = self.__photochemistry.get_xsec(rea) or reaction.get(
                    "xsecs_dict"
                )

            if rea.type == "photo" and self.radiation is not None:
                if rea.custom_rad_rate:
                    self.radiation.set_custom_rate(rea)
                    continue

                self.radiation.set_reaction_rate_coefficient(rea)

    def __normalize_network_extras(
        self, replace_nH: bool, loaded_from_jaff: bool = False
    ):
        """Standardize convenience symbols in all rate and auxiliary expressions.

        Replaces shorthand symbols (``nh``, ``ne``, ``ntot``, ``n_X``, …) with
        ``nden[i]`` references in every reaction rate, the chemical heating/cooling
        sum :attr:`dEdt_chem`, and the extra radiation source term
        :attr:`dRad_dt_extra`.

        Parameters
        ----------
        replace_nH : bool
            When ``True``, expand hydrogen-density shorthands to sums over
            H-bearing species.
        loaded_from_jaff : bool, optional
            When ``True``, the network was restored from a ``.jaff`` file whose
            stored ``rate`` is already the final (piecewise-collapsed,
            standardized) expression.  Re-evaluating the rate segments would
            wrap the already-collapsed rate in a second piecewise, so the stored
            rate is used as-is instead.
        """
        nden = self.ndens
        for r in self.reactions:
            if loaded_from_jaff:
                r.rate = self._standardize_symbols(r.rate, replace_nH)
            elif r.type == "photo" and self.radiation is not None:
                r.rate = self._standardize_symbols(r.rate, replace_nH)
                r.rate_segments[0].rate = r.rate
            else:
                for seg in r.rate_segments:
                    seg.rate = self._standardize_symbols(seg.rate, replace_nH)
                r.rate = r.rate_segments.sort().evaluate_equivalent_rate()

            r.tmin, r.tmax = r.rate_segments[0].tmin, r.rate_segments[-1].tmax
            dE_dt = r.dE * r.rate
            dRad_dt = r.dRad * r.rate
            for s in r.reactants.core:
                dE_dt *= nden[self.species[s.name].index]
                dRad_dt *= nden[self.species[s.name].index]
            self.dEdt_chem += dE_dt
            self.dRad_dt_extra += dRad_dt
        self.dEdt_chem = self._standardize_symbols(self.dEdt_chem, replace_nH)
        self.dRad_dt_extra = self._standardize_symbols(self.dRad_dt_extra, replace_nH)

    @staticmethod
    def __parse_rate(
        aux_chem_rate: str,
        rate: str,
        aux_funcs: dict[str, AuxiliaryFunctionsDict],
        global_vars: dict[str, Expr],
        n_photo: int,
    ) -> tuple[Expr, int]:
        """Convert a raw rate string to a SymPy expression.

        Checks, in priority order:

        1. Whether an auxiliary function named *aux_chem_rate* exists (custom rate).
        2. Whether *rate* is a global variable name.
        3. Whether *rate* describes a photo-reaction (contains ``"photo"``).
        4. Falls back to ``sympy.parse_expr``.

        Parameters
        ----------
        aux_chem_rate : str
            Key for the optional custom-rate auxiliary function (e.g. ``"chemrate0"``).
        rate : str
            Raw rate string from the network file.
        aux_funcs : dict[str, AuxiliaryFunctionsDict]
            Parsed auxiliary functions dictionary.
        global_vars : dict[str, Basic]
            Resolved global variable map from the network file.
        n_photo : int
            Running counter of photo-reactions seen so far.

        Returns
        -------
        tuple[Expr, int]
            ``(rate_expr, n_photo)`` where *n_photo* is
            incremented by 1 for photo-reactions.
        """
        if aux_chem_rate in aux_funcs:
            rate_expr = aux_funcs[aux_chem_rate]["def"]
        elif rate in global_vars:
            rate_expr = symbols(rate)
        elif "photo" in rate.lower():
            f: UndefinedFunction = Function("photorates")  # type: ignore
            n_photo += 1

            match = re.match(r"(?i:photo)\((.*?)\)", rate)
            if match:
                args_str = match.group(1)
                photo_args: list[str | float] = [
                    arg.strip() for arg in args_str.split(",") if arg.strip()
                ]
                while len(photo_args) < 2:
                    photo_args.append(1.0e99)

                rate_expr = f(n_photo, photo_args[0], photo_args[1])
            else:
                photo_args: list[str | float] = [
                    arg.strip() for arg in rate.split(",") if arg.strip()
                ]
                while len(photo_args) < 3:
                    photo_args.append(1.0e99)

                rate_expr = f(n_photo, photo_args[1], photo_args[2])
        else:
            rate_expr = _parse_rate_expr(rate)

        if not isinstance(rate_expr, Expr):
            raise ParserError(f"Rate expression is not an Expr: {rate_expr}")

        return rate_expr, n_photo

    def __parse_reaction_metadata(self, reaction: Reaction) -> None:
        if reaction.serialized not in self.spec._metadata["reaction_props"]:
            return

        rprops = self.spec._metadata["reaction_props"][reaction.serialized]
        if "shielding" in rprops:
            if reaction.type != "photo":
                raise ParserError(f"{reaction} is not a photo reaction")

            shielding_props = rprops["shielding"]
            if "type" not in shielding_props:
                shielding_props["type"] = "leiden"

            reaction._metadata["shielding"] = {
                k: (v.lower() if isinstance(v, str) else v)
                for k, v in shielding_props.items()
            }
            reaction._metadata["jaffgen"] = {
                "jaffgen_object": self.spec._metadata["jaffgen_object"]
            }

    @staticmethod
    def __detect_undefined_functions(
        expr: Expr | Basic, undef_funcs: set, interp_funcs: set
    ) -> None:
        """Scan *expr* for undefined function calls and categorise them.

        Functions whose names contain ``"interp"`` are added to *interp_funcs*;
        all others are added to *undef_funcs*.

        Parameters
        ----------
        expr : Expr | Basic
            SymPy expression to scan.
        undef_funcs : set
            Accumulator for unrecognised undefined function names.
        interp_funcs : set
            Accumulator for interpolation function names.
        """
        for f in expr.atoms(AppliedUndef):
            if "interp" in f.func.__name__:
                interp_funcs |= {f.func.__name__}
                continue
            undef_funcs |= {f.func.__name__}

    def to_jaff(self, filename: str | Path):
        """Serialise this network to a binary ``.jaff`` file.

        Parameters
        ----------
        filename : str | Path
            Output file path.  The ``.jaff`` extension is conventional but
            not enforced.
        """
        to_jaff_file(filename, self)

    @staticmethod
    def free_symbols(expr: Basic) -> set[Basic]:
        """Return the free symbols of *expr*, excluding ``nden`` matrix entries.

        ``nden[i]`` references are excluded because they are internal index
        variables, not user-visible physical symbols.

        Parameters
        ----------
        expr : Basic
            A SymPy expression.

        Returns
        -------
        set[Basic]
            Free symbols that do not involve ``"nden"``.
        """
        return {fs for fs in expr.free_symbols if "nden" not in str(fs)}

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

    @cached_property
    def ndens(self) -> MatrixSymbol:
        """Symbolic ``nden`` column vector of species number densities.

        A SymPy :class:`~sympy.matrices.expressions.MatrixSymbol` of shape
        ``(species.count, 1)``.  Entry ``nden[i]`` is the number density of the
        species with index ``i``.  Cached so every consumer shares one symbol.

        Returns
        -------
        sympy.MatrixSymbol
            The ``nden`` matrix symbol.
        """
        return MatrixSymbol("nden", self.species.count, 1)

    @cached_property
    def ntot(self) -> Expr:
        """Total number density ``Σ_i nden[i]`` over all species.

        Returns
        -------
        sympy.Expr
            Symbolic sum of every entry of :attr:`ndens`.
        """
        return sum(self.ndens[Idx(i)] for i in range(self.species.count))

    @cached_property
    def rho(self) -> Expr:
        """Mass density ``Σ_i m_i · nden[i]`` over all species.

        Each species contributes its mass ``m_i`` times its number density.
        Species with an unset mass (``mass is None``) contribute ``0``.

        Returns
        -------
        sympy.Expr
            Symbolic mass density.
        """
        return reduce(
            lambda x, y: x + y,
            [(s.mass or 0.0) * self.ndens[Idx(s.index)] for s in self.species],
        )

    def eos(self, gamma: float = 1.6666666666667) -> Expr:
        """Symbolic ideal-gas specific internal energy of the network.

        Thin wrapper around :func:`jaff.physics.get_eos`, which builds the
        expression from this network's :attr:`ntot` and :attr:`rho`.  Used by
        the code generator to form the temperature column of the Jacobian via
        the chain rule ``∂ẋ/∂e = (∂ẋ/∂T) / (∂e/∂T)``.

        Parameters
        ----------
        gamma : float, optional
            Adiabatic index.  Default ``5/3 ≈ 1.6̄`` (monoatomic ideal gas).

        Returns
        -------
        sympy.Expr
            Symbolic specific internal energy in CGS units (erg/g).
        """
        return get_eos(self, gamma)

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

    def _standardize_symbols(self, expr: Basic, replace_nH: bool) -> Expr:
        """Replace convenience symbols (nh, ne, ntot, n_X, …) with nden[i] references.

        When replace_nH is False, H/He element sums become ``nh``/``nhe`` symbols
        instead of being expanded over all species.
        """
        if expr == Float(0.0):
            return Float(0.0)

        nden = self.ndens
        reps = {}

        def get_element_sum(element):
            if element not in self.__element_sums:
                terms = []
                for i, spec in enumerate(self.species):
                    count = spec.exploded.count(element)
                    if count > 0:
                        terms.append(count * nden[Idx(i)])
                self.__element_sums[element] = sum(terms) if terms else None
            return self.__element_sums[element]

        simple_map = self._simple_map
        n_suffixes = self._n_suffixes

        for fs in expr.free_symbols:
            name = str(fs)
            low_name = name.lower()
            repl = None

            if low_name == "ntot":
                repl = self.ntot

            elif low_name == "nh":
                repl = get_element_sum("H") if replace_nH else symbols("nh")

            elif low_name in simple_map:
                spec_name = simple_map[low_name]
                repl = nden[Idx(self.species[spec_name].index)]

            elif low_name == "chi_pe":
                if self.radiation is None:
                    raise ParserError(
                        "In order to replace the 'chi_pe' symbol, radiation must be enabled"
                    )

                if self.dust is None:
                    raise ParserError(
                        "In order to replace the 'chi_pe' symbol, dust must be enabled"
                    )

                repl = self.dust.pe.chi

            elif low_name.startswith("n_"):
                core = name[2:]

                if core in ["H", "He"]:
                    if replace_nH:
                        repl = get_element_sum(core)
                    else:
                        repl = symbols(f"n{core.lower()}")

                else:
                    if core == "e":
                        core = "e-"
                    elif core[-1] in n_suffixes:
                        core = core[:-1] + n_suffixes[core[-1]]

                    if core in self.species:
                        repl = nden[Idx(self.species[core].index)]

            if repl is not None:
                reps[fs] = repl

        return expr.xreplace(reps)

    def sfluxes(self) -> list[Expr]:
        """Return symbolic flux expressions for all reactions.

        The flux of reaction *i* is ``rate_i * nden[r1] * nden[r2] * ...``,
        where the product runs over all reactant species.

        Returns
        -------
        list[Expr]
            One SymPy expression per reaction, in reaction-index order.
        """
        return get_sfluxes(self.reactions, self.species)

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
        return get_sodes(self.reactions, self.species)

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
        return get_sradodes(self.radiation, self.species, order)

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
            fname.with_suffix(".hdf5")

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
            fname.with_suffix(".txt")

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
