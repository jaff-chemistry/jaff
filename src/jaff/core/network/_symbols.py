# ABOUTME: NetworkSymbols — canonical symbols and symbolic quantities of a Network
# ABOUTME: (fixed physics symbols, densities, introspection, standardization)
"""Canonical symbols and symbolic quantities of a :class:`~jaff.Network`.

Reached as ``net.symbols``.  Fixed physics symbols are class attributes, so code
without a network can use ``NetworkSymbols.tgas``.  This module's own runtime imports
are limited to SymPy, the standard library and :mod:`jaff.errors`, but importing it
still initialises the ``jaff.core.network`` package; modules imported during that
package's initialisation (``io``, ``core.reaction``, ``physics``) must therefore import
:class:`NetworkSymbols` inside the function that uses it.
"""

from __future__ import annotations

from functools import cached_property, reduce
from typing import TYPE_CHECKING, ClassVar, Sequence

from sympy import Basic, Expr, Float, Function, IndexedBase, S, Symbol
from sympy.core.function import AppliedUndef, UndefinedFunction

from ...errors import ParserError

if TYPE_CHECKING:
    from ..species import Specie
    from .network import Network


class NetworkSymbols:
    """Canonical symbols and symbolic quantities of a network.

    Attributes
    ----------
    tgas, tdust : Symbol
        Gas and dust temperature [K].
    av : Symbol
        Visual extinction [mag].
    crate : Symbol
        Cosmic-ray ionisation rate [s⁻¹].
    chi : Symbol
        UV field scaling (Draine units).
    chi_pe : Symbol
        Photoelectric-band field placeholder, replaced by :meth:`standardize`.
    zd : Symbol
        Grain charge (``Zd``).
    vdisp : Symbol
        Velocity dispersion [cm s⁻¹].
    photorates : UndefinedFunction
        Placeholder function for photo-reaction rates.

    Notes
    -----
    The fixed symbols carry no SymPy assumptions on purpose: rate strings are
    parsed with ``parse_expr``, which creates assumption-free symbols, and SymPy
    only treats two symbols as equal when name *and* assumptions match.  Adding
    e.g. ``positive=True`` here would silently break substitution and ``diff``.
    """

    tgas: ClassVar[Symbol] = Symbol("tgas")
    tdust: ClassVar[Symbol] = Symbol("tdust")
    av: ClassVar[Symbol] = Symbol("av")
    crate: ClassVar[Symbol] = Symbol("crate")
    chi: ClassVar[Symbol] = Symbol("chi")
    chi_pe: ClassVar[Symbol] = Symbol("chi_pe")
    zd: ClassVar[Symbol] = Symbol("Zd")
    vdisp: ClassVar[Symbol] = Symbol("vdisp")
    photorates: ClassVar[UndefinedFunction] = Function("photorates")  # type: ignore

    def __init__(self, net: Network) -> None:
        """Bind to *net*; derived quantities are computed lazily from it.

        Parameters
        ----------
        net : Network
            Network whose species, reactions and settings the symbols describe.
        """
        self._net: Network = net
        self._element_sums: dict[str, Expr | None] = {}

    @staticmethod
    def ncol(name: str) -> Symbol:
        """Column-density symbol ``ncol_<name>`` [cm⁻²].

        Parameters
        ----------
        name : str
            Species name as used in shielding configuration (e.g. ``"H2"``).

        Returns
        -------
        Symbol
            ``Symbol(f"ncol_{name}")``.
        """
        return Symbol(f"ncol_{name}")

    @staticmethod
    def free_symbols(expr: Basic) -> set[Basic]:
        """Free symbols of *expr*, excluding ``nden`` entries.

        ``nden[i]`` references are internal index variables, not user-visible
        physical symbols.

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

    @cached_property
    def ndens(self) -> IndexedBase:
        """Symbolic ``nden`` indexed base for species number densities.

        A SymPy :class:`~sympy.tensor.indexed.IndexedBase` that provides
        scalar-indexed access. Entry ``nden[i]`` is the number density of the
        species with index ``i``.  Cached so every consumer shares one symbol.

        Returns
        -------
        sympy.IndexedBase
            The ``nden`` indexed base symbol.
        """
        return IndexedBase("nden", shape=(self._net.species.count,))

    @cached_property
    def ntot(self) -> Expr:
        """Total number density ``Σ_i nden[i]`` over all species.

        Returns
        -------
        sympy.Expr
            Symbolic sum of every entry of :attr:`ndens`.
        """
        return sum(self.ndens[i] for i in range(self._net.species.count))

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
            [(s.mass or 0.0) * self.ndens[s.index] for s in self._net.species],
        )

    @cached_property
    def n_hnuc(self) -> Expr:
        """Total hydrogen-nuclei number density ``Σ_i n_H(i) · nden[i]``.

        Each species contributes its hydrogen-atom count (``H2`` counts twice,
        ``H+`` once, ...) times its number density, so the sum is the total H
        nuclei density rather than a molecular count.  Equivalent to the
        ``n_H_nuc`` grammar token; used directly by the dust radiation-moment
        source terms (see :mod:`jaff.physics._equations`), cached so every
        consumer shares one expression.

        Returns
        -------
        sympy.Expr
            Symbolic total hydrogen-nuclei number density.  ``Float(0.0)`` when
            the network contains no H-bearing species.
        """
        total = self.element_sum("H")

        return total if total is not None else Float(0.0)

    def element_sum(self, element: str) -> Expr | None:
        """Nucleus density of *element*: ``Σ_i count_i · nden[i]`` (memoised).

        Parameters
        ----------
        element : str
            Canonical element symbol (e.g. ``"H"``, ``"He"``).

        Returns
        -------
        Expr | None
            The sum, or ``None`` when no species bears *element*.
        """
        if element not in self._element_sums:
            terms = [
                count * self.ndens[i]
                for i, spec in enumerate(self._net.species)
                if (count := spec.exploded.count(element)) > 0
            ]
            self._element_sums[element] = sum(terms) if terms else None

        return self._element_sums[element]

    def weighted_rate(self, weights: Sequence[Basic | float], rtol: float = 0.0) -> Expr:
        """Rate of ``Σ_i w_i n_i``: ``Σ_r (Σ_i w_i ν_ri) F_r``.

        ``ν = product_matrix − reactant_matrix`` (core species) and ``F_r`` are the
        reaction fluxes (:meth:`Network.sfluxes`).  Flux expressions are substituted,
        so the result can be differentiated and shares CSE with the species rows.
        Reactions whose coefficient is exactly zero are skipped.

        Parameters
        ----------
        weights : Sequence[Basic | float]
            One weight per species, in species-index order.
        rtol : float, optional
            Numeric weights only: a reaction coefficient ``c_r = Σ_i w_i ν_ri`` is
            treated as zero when ``|c_r| <= rtol · Σ_i |w_i ν_ri|``, i.e. when it is a
            floating-point rounding residue of terms that cancel exactly in
            principle.  Default ``0.0`` keeps every non-zero coefficient.

        Returns
        -------
        Expr
            The weighted rate [weight units · cm⁻³ s⁻¹].

        Raises
        ------
        ValueError
            If the number of weights differs from the number of species.
        """
        weights = list(weights)
        if len(weights) != self._net.species.count:
            raise ValueError(
                f"weighted_rate needs one weight per species "
                f"({self._net.species.count}), got {len(weights)}"
            )

        net = self._net
        nu = net.product_matrix - net.reactant_matrix
        total = S.Zero
        for row, flux in zip(nu, net.sfluxes()):
            terms = [w * int(n) for w, n in zip(weights, row) if n]
            coeff = sum(terms)
            if rtol and abs(coeff) <= rtol * sum(abs(t) for t in terms):
                continue

            if coeff != 0:
                total += coeff * flux

        return total

    @cached_property
    def dntot_dt(self) -> Expr:
        """Chemical rate of total number density, ``Σ_r Δn_r F_r`` [cm⁻³ s⁻¹]."""
        return self.weighted_rate([1] * self._net.species.count)

    @cached_property
    def drho_dt(self) -> Expr:
        """Chemical rate of mass density, ``Σ_r Δm_r F_r`` [g cm⁻³ s⁻¹].

        Never forced to zero: every reaction contributes its actual mass change,
        so the per-cell rate is kept even for reactions that pass
        :meth:`Reaction.check_mass`.  Only floating-point rounding residues
        (``|Δm_r| <= 1e-12 · Σ_i |m_i ν_ri|``, far below one electron mass) are
        dropped, so the expression -- and the Jacobian sparsity -- does not depend
        on the last bits of the species masses.  Unset masses count as ``0``.
        """
        masses = [s.mass or 0.0 for s in self._net.species]
        return self.weighted_rate(masses, rtol=1e-12)

    # The introspection caches below are filled on first access (Network.__init__,
    # after loading).  Mutating rates or thermodynamics afterwards leaves them stale.

    def _expressions(self) -> list[Expr]:
        """Standardized network expressions: rates and energy/radiation sources.

        The aggregated ``dEdt_chemical`` and ``dRad_dt_extra`` are used instead of
        each reaction's raw ``dE``/``dRad`` because only the aggregates have their
        convenience symbols (``n_X``, ...) resolved to ``nden`` entries.
        """
        net = self._net
        return [
            *(r.rate for r in net.reactions),
            net.thermodynamics.dEdt_chemical.volumetric,
            net.thermodynamics.dEdt_extra.volumetric,
            net.dRad_dt_extra,
        ]

    @cached_property
    def variables(self) -> frozenset[Basic]:
        """Free symbols across all network expressions, excluding ``nden`` entries."""
        return frozenset().union(*(self.free_symbols(e) for e in self._expressions()))

    @cached_property
    def _applied_functions(self) -> frozenset[str]:
        """Names of all undefined (applied) functions across network expressions."""
        return frozenset(
            call.func.__name__
            for e in self._expressions()
            for call in e.atoms(AppliedUndef)
        )

    @cached_property
    def interp_functions(self) -> frozenset[str]:
        """Names of interpolation functions (containing ``"interp"``)."""
        return frozenset(name for name in self._applied_functions if "interp" in name)

    @cached_property
    def undefined_functions(self) -> frozenset[str]:
        """Names of undefined, non-interpolation functions."""
        return frozenset(name for name in self._applied_functions if "interp" not in name)

    @cached_property
    def _element_lookup(self) -> dict[str, str]:
        """Lower-cased element token -> canonical element symbol.

        Valid only once ``net.mass_dict`` is loaded; first touched by
        :meth:`standardize`, which runs after loading.
        """
        return {s.lower(): s for s in self._net.mass_dict}

    @cached_property
    def _charge_reverse(self) -> dict[str, Specie]:
        """j/k-normalized species identifier -> ``Specie`` (see ``charge_reverse_map``).

        Assumes a collision-free network; raises ValueError on case-distinct
        colliding species (e.g. CO / Co).
        """
        return self._net.species.charge_reverse_map()

    def standardize(self, expr: Basic) -> Expr:
        """Replace convenience symbols with ``nden``-based expressions.

        ``ntot`` → :attr:`ntot`; ``n_<species>`` → that species' density;
        ``n_e`` → electron density; ``n_<element>_nuc`` → :meth:`element_sum`
        (or the free symbol ``n<element>_nuc`` when the network was built with
        ``expand_nuclei=False``); ``rc_<int>`` → the rate of the reaction whose
        file-side number is ``<int>`` (looked up via ``reactions.by_source_index``,
        not the catalogue position, so it stays dedup-safe and consistent with
        ``chemRateN``); ``chi_pe`` → ``net.dust.pe.chi`` (needs radiation and
        dust).  Name matching is case-insensitive.

        Parameters
        ----------
        expr : Basic
            Expression to standardize.

        Returns
        -------
        Expr
            *expr* with every recognised convenience symbol replaced.

        Raises
        ------
        ParserError
            On an unknown species/element, a charged nucleus alias, a missing
            ``rc_`` target, or ``chi_pe`` without radiation or dust.
        ValueError
            If two species collide case-insensitively (e.g. ``CO`` / ``Co``).
        """
        if expr == Float(0.0):
            return Float(0.0)

        reps = {}
        for fs in expr.free_symbols:
            repl = self._replacement(str(fs))
            if repl is not None:
                reps[fs] = repl

        return expr.xreplace(reps)

    def _replacement(self, name: str) -> Basic | None:
        """Replacement for the convenience symbol *name*, or ``None`` to keep it.

        Raises
        ------
        ParserError
            For ``chi_pe`` without radiation/dust, or an ``rc_`` target that is
            not in the network (see also :meth:`_density_replacement`).
        """
        net = self._net
        low_name = name.lower()

        if low_name == "ntot":
            return self.ntot

        if low_name == self.chi_pe.name:
            if net.radiation is None:
                raise ParserError(
                    "In order to replace the 'chi_pe' symbol, radiation must be enabled"
                )
            if net.dust is None:
                raise ParserError(
                    "In order to replace the 'chi_pe' symbol, dust must be enabled"
                )
            return net.dust.pe.chi

        if low_name.startswith("n_"):
            return self._density_replacement(name)

        if low_name.startswith("rc_"):
            try:
                num = int(name[3:])
            except ValueError:
                net.logger.error(
                    f"The 'rc_' keyword in {net.spec.funcfile} must be followed by "
                    f"an integer\ndenoting the reaction number. Found {name}"
                )
                return None

            rxn = net.reactions.by_source_index(num)
            if rxn is None:
                raise ParserError(
                    f"'{name}' references reaction {num}, which is not in the network"
                )
            return rxn.rate

        return None

    def _density_replacement(self, name: str) -> Basic | None:
        """Replacement for an ``n_*`` density symbol (species, electron, nucleus).

        Raises
        ------
        ParserError
            For a charged nucleus alias, an unknown element, an element no species
            bears, or a species name not in the network.
        """
        net = self._net
        core = name[2:]
        core_low = core.lower()

        if core_low.endswith("_nuc"):
            base = core_low[:-4]
            if base.endswith(("j", "k")):
                raise ParserError(
                    f"'{name}' is invalid: a nucleus sum is per-element, "
                    f"so a charged nucleus alias is meaningless"
                )
            element = self._element_lookup.get(base)
            if element is None:
                raise ParserError(
                    f"'{name}' requests a nucleus sum for unknown element '{base}'"
                )
            if not net.spec.expand_nuclei:
                return Symbol(f"n{base}_nuc")

            total = self.element_sum(element)
            if total is None:
                raise ParserError(
                    f"'{name}': no species in the network bears element '{element}'"
                )
            return total

        if core == "e":
            if "e-" in net.species:
                return self.ndens[net.species["e-"].index]
            return None

        sp = self._charge_reverse.get(core_low)
        if sp is None:
            raise ParserError(
                f"Density symbol '{name}' does not match any species in this network"
            )
        return self.ndens[sp.index]

    def log_summary(self) -> None:
        """Log the network's free variables, interpolation and undefined functions."""
        logger = self._net.logger
        logger.info(
            "Variables found: "
            f"{', '.join(sorted(f'[cyan]{s}[/]' for s in self.variables))}"
        )

        if self.interp_functions:
            names = sorted(self.interp_functions)
            logger.info(
                "Found the following interpolation functions: "
                f"{', '.join(f'[cyan]{name}[/]' for name in names)}"
            )

        if self.undefined_functions:
            names = sorted(self.undefined_functions)
            logger.warning(
                "Found undefined functions "
                f"{', '.join(f'[red]{name}[/]' for name in names)}"
            )
