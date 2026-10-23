# ABOUTME: Symbolic chemical fluxes, species ODEs and radiation-moment ODEs
# ABOUTME: (get_sfluxes, get_sodes, get_sradodes) built from a network
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

from sympy import Basic, Expr, Float, IndexedBase

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

    For a reaction with rate coefficient *k* and reactants A, B the flux is

    ```
    flux_i = k_i * nden[idx_A] * nden[idx_B]
    ```

    The number densities are represented as indexed-base symbols ``nden`` that
    support scalar indexing (``nden[i]`` for species *i*).  Each flux is the
    reaction's :attr:`Reaction.rate`
    (``rate_coefficient`` times the density of each core reactant;
    pseudo-species such as ``_PHOTON`` carry no density).
    :meth:`Network.sfluxes` calls this
    with the network's own arguments.

    Parameters
    ----------
    reactions : Reactions
        Collection of all reactions in the network; the flux of each is its
        ``.rate``.
    species : Species
        Collection of all species.  Not used: each reaction resolves its
        reactant indices from its own network.
    nden : IndexedBase
        The network's density base (``net.symbols.ndens``).  Not used, for
        the same reason.

    Returns
    -------
    list of sympy.Expr
        List of length ``reactions.count``.  ``fluxes[i]`` is the symbolic
        flux expression for the *i*-th reaction.

    Notes
    -----
    The flux is purely a *loss* term from the reactants' perspective; signs
    are applied in :func:`get_sodes`.

    Examples
    --------
    >>> from jaff import Network
    >>> from jaff.physics import get_sfluxes
    >>> net = Network("h_photoionization")
    >>> get_sfluxes(net.reactions, net.species, net.symbols.ndens)[1]
    1.65941781598291e-10*nden[1]*nden[2]/tgas**0.7
    """
    fluxes: list[Expr] = [Float(0.0) for _ in range(reactions.count)]

    for i, reaction in enumerate(reactions):
        fluxes[i] = reaction.rate

    return fluxes


def get_sodes(
    reactions: "Reactions",
    species: Species,
    nden: IndexedBase,
) -> list[Basic]:
    """
    Assemble the symbolic ODE right-hand sides for all species.

    For each species *s* the ODE is

    ```
    d nden[s] / dt = sum_{i: s in products(i)} flux_i
                   - sum_{i: s in reactants(i)} flux_i
    ```

    with ``flux_i`` from :func:`get_sfluxes`.  A species appearing several
    times on one side of a reaction contributes the flux once per occurrence.

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
    Each participant's slot is its species ``index``: a species' ``fidx`` is
    always an ``"idx_<name>"`` string (see
    :meth:`~jaff.core.species.Specie.get_fidx`).  The code also has a branch
    that casts a ``fidx`` not starting with ``"idx_"`` to ``int`` and uses it
    as the slot; no species currently takes it.

    Examples
    --------
    >>> from jaff import Network
    >>> from jaff.physics import get_sodes
    >>> net = Network("h_photoionization")  # H + γ -> H+ + e-, H+ + e- -> H
    >>> [s.name for s in net.species.core]
    ['H', 'H+', 'e-']
    >>> sodes = get_sodes(net.reactions, net.species, net.symbols.ndens)
    >>> sodes[1] == -sodes[0]
    True
    """
    fluxes = get_sfluxes(reactions, species, nden)
    sodes: list[Basic] = [Float(0.0) for _ in range(species.core.count)]

    for i, reaction in enumerate(reactions):
        for rr in reaction.reactants.core:
            # Output-array slot: the species' index (fidx is always an "idx_*"
            # tag; the int(fidx) branch is not reached by any species).
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

    - ``grate[i]``: the *source/sink* term for the density moment: the
      photoionisation/photodissociation losses (each photo reaction's band
      loss coefficient ``rgrouped[i].k`` times its reactant densities), the
      absorption by the medium ``-c * <α>_i * den[i]`` (from a
      ``globalAbsorptionCoefficient``), the radiation emitted by reactions
      with a spectral yield (``reaction.rad.group_yields`` times the reaction
      rate) and the band's ``globalEmissivity``.
    - ``gflux[i]``: the sink term for the flux moment, obtained by
      substituting ``den[i] → rflux[i]`` in the losses.  Emission is
      isotropic and feeds the density moment only; the absorption term is
      included only when ``radiation.absorption_attenuates_flux`` is set.

    With dust enabled (``net.dust``), the dust attenuation of each moment is
    subtracted too (see :class:`~jaff.physics.Dust`).  Photo reactions whose
    rate comes from a ``.jfunc`` ``chemRate<N>`` are not coupled to the bands
    and are skipped.

    The two moments for each band are placed in the output list according to
    the *order* convention (see :meth:`Radiation.ordered_index`).

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
        If ``net.radiation`` is ``None`` (no bands have been configured), or a
        photo reaction without a custom rate has no per-band data.
    ValueError
        If *order* is not one of ``{0, 1, 2, 3}``.

    Notes
    -----
    The substitution ``den[i] → rflux[i]`` (via ``xreplace``) turns the
    density-moment losses into the attenuation term of the first-moment (flux)
    equation of the two-moment radiation transport system.  Transport terms
    (the divergence of the flux, the radiation pressure) are not included.

    Examples
    --------
    >>> from jaff import Network
    >>> from jaff.physics import get_sradodes
    >>> net = Network(
    ...     "h_photoionization",
    ...     radiation_props={"bands": [1, 13.6, 100, "inf"]},
    ... )
    >>> get_sradodes(net)  # [den_0, flux_0, den_1, flux_1, den_2, flux_2]
    [0, 0, -5.85154082319921e-8*nden[0]*photden[1], -5.85154082319921e-8*nden[0]*rflux[1], -1.1420465158197e-10*nden[0]*photden[2], -1.1420465158197e-10*nden[0]*rflux[2]]
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

    # Per-band totals over all photo reactions: the radiation removed.
    group_rates: list[Expr] = [Float(0.0) for _ in rad_groups]

    for reaction in jaff_progress.track(
        net.reactions.with_type("photo"), description="Generating radiation equations"
    ):
        # A custom-rate reaction is not coupled to the bands.
        if reaction.photo is not None and reaction.photo.custom_rad_rate:
            continue

        # With radiation, Reaction.__init__ fills rgrouped for every other photo
        # reaction or raises, so missing data is a broken invariant.
        if reaction.photo is None or not reaction.photo.rgrouped:
            raise RuntimeError(
                f"Photo reaction {reaction} has no per-band radiation data"
            )

        # rgrouped is built by iterating rad_groups in order, one entry per band.
        for group, props in zip(rad_groups, reaction.photo.rgrouped, strict=True):
            rrate = props.k
            # Multiply by all reactant number densities (mass-action kinetics)
            # so rrate becomes the full reaction flux (k * prod(nden)).
            for reactant in reaction.reactants.core:
                rrate *= nden[net.species[str(reactant)].index]

            # Photochemical reactions *remove* radiation, hence the minus sign.
            group_rates[group.index] -= rrate

    # Per-group radiation emitted by any reaction with a spectral yield.
    group_emissions: list[Expr] = [Float(0.0) for _ in rad_groups]
    for reaction in net.reactions.radiation_emitters():
        for group, group_yield in zip(rad_groups, reaction.rad.group_yields, strict=True):
            group_emissions[group.index] += group_yield * reaction.rate

    for group in rad_groups:
        group_rate = group_rates[group.index]

        # Absorption by the medium: -c <α>_i den_i.  Added before the den → rflux
        # substitution it also attenuates the flux; afterwards only the density.
        absorption: Expr = Float(0.0)
        if group.absorption_coefficient is not None:
            absorption = -net.radiation.c * group.absorption_coefficient * group.sym

        if net.radiation.absorption_attenuates_flux:
            group_rate += absorption
            absorption = Float(0.0)

        # The flux-moment equation is obtained by substituting den → rflux in
        # the density-moment equation (two-moment closure).
        flux = group_rate.xreplace(flux_map)
        group_rate += absorption

        # Emission is isotropic, so it feeds the density equation only.
        group_rate += group_emissions[group.index]
        if group.emissivity is not None:
            group_rate += group.emissivity

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
    respectively; a kind that is ``"none"`` is skipped.

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
    if u_reduction != "none":
        grate -= (
            net.symbols.zd
            * net.radiation.c
            * group.sym
            * net.symbols.n_hnuc
            * net.dust.tabular.avg_cross_section_per_hnuc(
                u_reduction, (group.lower, group.upper)
            )
        )
    if f_reduction != "none":
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
