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

The equation of state lives in :mod:`jaff.physics.thermodynamics.eos`.

The symbolic expressions are later code-generated (via SymPy's code printers)
into efficient numerical kernels.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from sympy import Basic, Expr, Float, Idx, IndexedBase

from ..io._logger import jaff_progress

if TYPE_CHECKING:
    from .. import Network, Reactions, Species
    from ..physics import RadiationGroup


def get_sfluxes(
    reactions: "Reactions",
    species: Species,
    nden: IndexedBase,
) -> list[Expr]:
    """
    Build the symbolic reaction flux for every reaction in the network.

    For a reaction with rate coefficient *k* and reactants A, B the flux is::

        flux_i = k_i * nden[idx_A] * nden[idx_B]

    The number densities are represented as indexed-base symbols ``nden`` that
    support scalar indexing (``nden[i]`` for species *i*).

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
    nden : IndexedBase
        The network's density base (``net.symbols.ndens``).

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

    for i, reaction in enumerate(reactions):
        flux = reaction.rate
        for reactant in reaction.reactants.core:
            flux *= nden[species[str(reactant)].index]

        fluxes[i] = flux

    return fluxes


def get_sodes(
    reactions: "Reactions",
    species: Species,
    nden: IndexedBase,
) -> list[Basic]:
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
    nden : IndexedBase
        The network's density base (``net.symbols.ndens``).

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
    fluxes = get_sfluxes(reactions, species, nden)
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


def get_sradodes(net: "Network", order: int = 0) -> list[Expr]:
    """
    Build symbolic radiation-moment ODE right-hand sides for all frequency bands.

    The radiation field is described by two moments per band:

    - **Energy/photon density** ``den[i]`` (``radeden`` in erg/cm³, or
      ``photden`` in cm⁻³, depending on ``radiation.mode``).
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
    net : Network
        The network whose :attr:`~Network.radiation` field supplies the band
        definitions and per-band per-reaction rate coefficients, and whose
        species provide number-density indexing.  ``net.radiation`` must not
        be ``None``.
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
        If ``net.radiation`` is ``None`` (no bands have been configured).
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
    if net.radiation is None:
        raise RuntimeError("No radiation bands found. Radiation odes cannot be generated")

    if order not in [0, 1, 2, 3]:
        raise ValueError("Invalid order: Supported orders are 0, 1, 2, 3")

    rad_groups = net.radiation.groups
    nden = net.symbols.ndens

    rflux = IndexedBase("rflux", shape=(net.radiation.nbands,))
    # Mapping used to obtain the flux-moment equation from the density-moment
    # equation: replace each density symbol den[i] with the flux rflux[i].
    flux_map = {g.sym: rflux[i] for i, g in enumerate(net.radiation.groups)}
    grate: list[Expr | float] = [Float(0.0) for _ in range(net.radiation.nbands)]
    gflux: list[Expr | float] = [Float(0.0) for _ in range(net.radiation.nbands)]

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
                rrate *= nden[Idx(net.species[str(reactant)].index)]

            # Photochemical reactions *remove* radiation, hence the minus sign.
            group_rate -= rrate
            group_dRad_dt_extra += props["delta_rad"] * rrate

        # The flux-moment equation is obtained by substituting den → rflux in
        # the density-moment equation (two-moment closure).
        flux = group_rate.xreplace(flux_map)

        # Add the user-supplied dRad term to the density equation.  Divide by
        # the band-average photon energy in BOTH modes:
        group_rate += group_dRad_dt_extra / (group.eavg or 1)

        if net.dust is not None:
            group_rate, flux = handle_dust_reduction(net, group, group_rate, flux, rflux)

        grate[group.index] = group_rate
        gflux[group.index] = flux

    # Allocate output array: 2 slots per band (one density, one flux).
    radodes: list[Expr] = [Float(0.0) for _ in range(2 * net.radiation.nbands)]

    # Place each (rate, flux) pair at the positions dictated by the chosen
    # ordering convention.
    for i, (rate, flux) in enumerate(zip(grate, gflux)):
        ei, fi = net.radiation.ordered_index(i, order)
        radodes[ei] = rate
        radodes[fi] = flux

    return radodes


def handle_dust_reduction(
    net: Network, group: RadiationGroup, grate: Expr, gflux: Expr, rflux: IndexedBase
) -> tuple[Expr, Expr]:
    """Subtract dust absorption/transport reductions from a band's ODE terms.

    For a single radiation band, this removes the radiation lost to dust from
    the energy-density source term ``grate`` and the flux source term
    ``gflux``.  Each reduction term is
    ``Zd * c * <moment> * n_hnuc * avg_cross_section_per_hnuc(kind, band)``,
    where ``<moment>`` is the band symbol ``group.sym`` for the energy density
    and ``rflux[group.index]`` for the flux.  The energy-density and flux
    reductions use the dust's ``u_reduction`` and ``f_reduction`` kinds
    respectively; a kind that is ``None`` or ``"none"`` is skipped.

    Parameters
    ----------
    net : Network
        The network supplying the radiation field (``net.radiation``) and the
        dust module (``net.dust``); both must be enabled.
    group : RadiationGroup
        The radiation band whose symbol, index, and energy bounds
        (``group.lower``, ``group.upper``) select the averaged cross-section.
    grate : Expr
        The band's energy-density source term to reduce.
    gflux : Expr
        The band's flux source term to reduce.
    rflux : IndexedBase
        Flux moment symbol, indexed by band to form the flux reduction term.

    Returns
    -------
    tuple of sympy.Expr
        The updated ``(grate, gflux)`` pair with dust reductions subtracted.
    """
    assert net.radiation is not None
    assert net.dust is not None

    u_reduction = net.dust.u_reduction
    f_reduction = net.dust.f_reduction
    if u_reduction not in (None, "none"):
        grate -= (
            net.symbols.zd
            * net.radiation.c
            * group.sym
            * net.symbols.n_hnuc
            * net.dust.tabular.avg_cross_section_per_hnuc(
                u_reduction, (group.lower, group.upper)
            )
        )
    if f_reduction not in (None, "none"):
        gflux -= (
            net.symbols.zd
            * net.radiation.c
            * rflux[group.index]
            * net.symbols.n_hnuc
            * net.dust.tabular.avg_cross_section_per_hnuc(
                f_reduction, (group.lower, group.upper)
            )
        )

    return grate, gflux
