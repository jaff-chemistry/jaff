"""
Symbolic ODE and flux generators for astrochemical reaction networks.

This module builds SymPy symbolic expressions for:

- **Chemical fluxes** -- the reaction rate multiplied by the number densities
  of each reactant (``get_sfluxes``).
- **Chemical ODEs** -- the net rate of change of each species number density,
  obtained by summing fluxes that produce or destroy each species
  (``get_sodes``).
- **Radiation-moment ODEs** -- the zeroth-moment (energy/photon density) and
  first-moment (energy/photon flux) equations for each frequency band, taking
  into account photoionisation/photodissociation sinks and any user-supplied
  radiation source/sink terms (``get_sradodes``).
- **Equation of state** -- the ideal-gas specific internal energy used to
  couple the gas temperature into the Jacobian (``get_eos``).

The symbolic expressions are later code-generated (via SymPy's code printers)
into efficient numerical kernels.
"""

from __future__ import annotations

from functools import cache
from typing import TYPE_CHECKING

from sympy import Basic, Expr, Float, Idx, MatrixSymbol, symbols

from ..io._logger import jaff_progress
from .constants import k_B

if TYPE_CHECKING:
    from .. import Network, Reactions, Species
    from .photo_reactions._radiation import Radiation


def get_sfluxes(reactions: "Reactions", species: Species) -> list[Expr]:
    """
    Build the symbolic reaction flux for every reaction in the network.

    For a reaction with rate coefficient *k* and reactants A, B the flux is::

        flux_i = k_i * nden[idx_A] * nden[idx_B]

    The number densities are represented as entries of the SymPy
    ``MatrixSymbol`` ``nden`` (shape ``(species.core.count, 1)`` — only core
    species enter the integrated state), so the returned
    expressions reference ``nden[j]`` symbolically and can be differentiated or
    printed by any SymPy backend.

    Parameters
    ----------
    reactions : Reactions
        Collection of all reactions in the network.  Each element must expose
        ``.rate`` (a SymPy ``Expr``) and ``.reactants`` (an iterable of
        reactant objects with a string representation matching a key in
        *species*).
    species : Species
        Collection of all species.  Used to look up the numeric index of each
        reactant via ``species[str(reactant)].index``.

    Returns
    -------
    list of sympy.Expr
        List of length ``reactions.count``.  ``fluxes[i]`` is the symbolic
        flux expression for the *i*-th reaction, or ``Float(0.0)`` if the
        reaction has no rate (should not occur in practice).

    Notes
    -----
    The flux is purely a *loss* term from the reactants' perspective; signs
    are applied in :func:`get_sodes`.
    """
    fluxes: list[Expr] = [Float(0.0) for _ in range(reactions.count)]
    nden_matrix = MatrixSymbol("nden", species.core.count, 1)

    for i, reaction in enumerate(reactions):
        flux = reaction.rate
        for reactant in reaction.reactants.core:
            flux *= nden_matrix[species[str(reactant)].index]

        fluxes[i] = flux

    return fluxes


def get_sodes(reactions: "Reactions", species: Species) -> list[Basic]:
    """
    Assemble the symbolic ODE right-hand sides for all species.

    For each species *s* the ODE is::

        d nden[s] / dt = sum_{i: s in products(i)} flux_i
                       - sum_{i: s in reactants(i)} flux_i

    Parameters
    ----------
    reactions : Reactions
        Collection of all reactions in the network.
    species : Species
        Collection of all species, used to resolve array indices.

    Returns
    -------
    list of sympy.Basic
        List of length ``species.core.count`` (special pseudo-species are not
        integrated).  ``sodes[j]`` is the symbolic time-derivative of the
        *j*-th core species number density.

    Notes
    -----
    The index used for each participant is determined by ``fidx``:

    - If ``fidx`` is a string starting with ``"idx_"`` the *runtime* index
      attribute of the participant object is used (dynamic lookup, e.g. for
      named network slots).
    - Otherwise ``fidx`` is cast to ``int`` and used directly (static index,
      e.g. a literal position in a fixed-order array).

    This dual-path allows the same code to handle both named-species networks
    and fixed-layout networks produced by certain code-generation backends.
    """
    fluxes = get_sfluxes(reactions, species)
    sodes: list[Basic] = [Float(0.0) for _ in range(species.core.count)]

    for i, reaction in enumerate(reactions):
        for rr in reaction.reactants.core:
            # Choose the output-array slot: either the species' runtime index
            # (when fidx is a "idx_*" tag) or a literal integer position.
            idx = (
                rr.index
                if isinstance(rr.fidx, str) and rr.fidx.startswith("idx_")
                else int(rr.fidx)
            )
            sodes[idx] -= fluxes[i]

        # Add flux to every product (creation term)
        for pp in reaction.products.core:
            idx = (
                pp.index
                if isinstance(pp.fidx, str) and pp.fidx.startswith("idx_")
                else int(pp.fidx)
            )
            sodes[idx] += fluxes[i]

    return sodes


def get_sradodes(
    radiation: "Radiation" | None, species: Species, order: int = 0
) -> list[Expr]:
    """
    Build symbolic radiation-moment ODE right-hand sides for all frequency bands.

    The radiation field is described by two moments per band:

    - **Energy/photon density** ``den[i]`` (``radeden`` in erg/cm³, or
      ``photden`` in cm⁻³, depending on ``radiation.energy_density``).
    - **Energy/photon flux** ``rflux[i]``.

    For each band *i* the function computes:

    - ``grate[i]``: the *source/sink* term for the density moment, including
      photoionisation/photodissociation losses and any user-defined
      ``dRad`` source terms.
    - ``gflux[i]``: the spatial-gradient term for the flux moment (obtained
      by substituting ``den[i] → rflux[i]`` in ``grate``).

    The two moments for each band are interleaved in the output list according
    to the *order* convention (see :meth:`Radiation.ordered_index
    <jaff.physics._radiation.Radiation.ordered_index>`).

    Parameters
    ----------
    radiation : Radiation or None
        Radiation field descriptor containing band definitions and per-band
        per-reaction rate coefficients.  Must not be ``None``.
    species : Species
        Collection of all chemical species (used for number-density indexing).
    order : {0, 1, 2, 3}, optional
        Layout convention for the output array:

        - ``0`` (default): ``[den_0, flux_0, den_1, flux_1, ...]``
          (energy-density first, interleaved).
        - ``1``: ``[flux_0, den_0, flux_1, den_1, ...]``
          (flux first, interleaved).
        - ``2``: ``[den_0, den_1, ..., flux_0, flux_1, ...]``
          (all densities then all fluxes, energy-density block first).
        - ``3``: ``[flux_0, flux_1, ..., den_0, den_1, ...]``
          (all fluxes then all densities, flux block first).

    Returns
    -------
    list of sympy.Expr
        List of length ``2 * radiation.nbands``.  Even and odd slots (or
        front/back halves for order 2/3) contain density and flux terms
        according to the chosen *order*.

    Raises
    ------
    RuntimeError
        If *radiation* is ``None`` (no bands have been configured).
    ValueError
        If *order* is not one of ``{0, 1, 2, 3}``.

    Notes
    -----
    Each reaction's ``dRad`` contribution is its band-integrated ``delta_rad``
    (the radiation energy added per reaction event, erg) times the reaction
    flux ``k * prod(nden)``, then divided by the band-average photon energy
    ``group.eavg`` (erg) in **both** modes.  In energy-density mode ``props["k"]``
    already carries an extra ``eavg`` (``radeden`` is an energy density), so the
    division recovers the true flux and the result is erg cm⁻³ s⁻¹ added to
    ``radeden``; in photon-density mode the division converts the per-event
    energy to a photon count, giving cm⁻³ s⁻¹ added to ``photden``.

    The substitution ``den[i] → rflux[i]`` (via ``xreplace``) yields the
    flux-divergence term needed in the first-moment (flux) equation of the
    two-moment radiation transport system.
    """
    if radiation is None:
        raise RuntimeError("No radiation bands found. Radiation odes cannot be generated")

    if order not in [0, 1, 2, 3]:
        raise ValueError("Invalid order: Supported orders are 0, 1, 2, 3")

    rad_groups = radiation.groups
    nden = MatrixSymbol("nden", species.core.count, 1)

    # Choose the symbolic name for the radiation density variable based on
    # whether the field is tracked as energy density (erg/cm³) or photon
    # number density (cm⁻³).
    den = MatrixSymbol(
        "radeden" if radiation.energy_density else "photden",
        radiation.nbands,
        1,
    )
    rflux = MatrixSymbol("rflux", radiation.nbands, 1)
    # Mapping used to obtain the flux-moment equation from the density-moment
    # equation: replace each density symbol den[i] with the flux rflux[i].
    flux_map = {den[Idx(i)]: rflux[Idx(i)] for i in range(radiation.nbands)}
    grate, gflux = (
        [Float(0.0)] * radiation.nbands,
        [Float(0.0)] * radiation.nbands,
    )

    for group in jaff_progress.track(
        rad_groups, description="Generating radiation equations"
    ):
        group_rate: Basic = Float(0.0)
        group_dRad_dt_extra = Float(0.0)
        for reaction, props in group.props.items():
            rrate = props["k"]
            # Multiply by all reactant number densities (mass-action kinetics)
            # so rrate becomes the full reaction flux (k * prod(nden)).
            for reactant in reaction.reactants.core:
                rrate *= nden[Idx(species[str(reactant)].index)]

            # Photochemical reactions *remove* radiation, hence the minus sign.
            group_rate -= rrate
            group_dRad_dt_extra += props["delta_rad"] * rrate

        # The flux-moment equation is obtained by substituting den → rflux in
        # the density-moment equation (two-moment closure).
        flux = group_rate.xreplace(flux_map)

        # Add the user-supplied dRad term to the density equation.  Divide by
        # the band-average photon energy in BOTH modes:
        group_rate += group_dRad_dt_extra / (group.eavg or 1)

        grate[group.index] = group_rate
        gflux[group.index] = flux

    # Allocate output array: 2 slots per band (one density, one flux).
    radodes: list[Expr] = [Float(0.0) for _ in range(2 * radiation.nbands)]

    # Place each (rate, flux) pair at the positions dictated by the chosen
    # ordering convention.
    for i, (rate, flux) in enumerate(zip(grate, gflux)):
        ei, fi = radiation.ordered_index(i, order)
        radodes[ei] = rate
        radodes[fi] = flux

    return radodes


@cache
def get_eos(net: "Network", gamma: float = 1.6666666666667) -> Expr:
    """Return the symbolic ideal-gas specific internal energy.

    Uses the ideal-gas equation of state::

        e = n_tot · k_B · T_gas / (ρ · (γ − 1))   [erg / g]

    where ``n_tot`` is the total number density (:attr:`Network.ntot`), ``ρ``
    is the mass density (:attr:`Network.rho`), ``k_B`` is the Boltzmann
    constant in CGS (erg K⁻¹) and ``tgas`` is a SymPy symbol for the gas
    temperature in Kelvin.

    This expression drives the temperature column of the Jacobian via the
    chain rule ``∂ẋ/∂e = (∂ẋ/∂T) / (∂e/∂T)``.  The result is cached (via
    :func:`functools.cache`) since it depends only on *net* and *gamma*.

    Parameters
    ----------
    net : Network
        Network supplying the symbolic ``n_tot`` and ``ρ`` sums.
    gamma : float, optional
        Adiabatic index.  Default ``5/3 ≈ 1.6̄`` (monoatomic ideal gas).

    Returns
    -------
    sympy.Expr
        Symbolic specific internal energy in CGS units (erg/g).
    """
    tgas = symbols("tgas")

    return net.ntot * k_B.cgs.value * tgas / net.rho * 1.0 / (gamma - 1.0)
