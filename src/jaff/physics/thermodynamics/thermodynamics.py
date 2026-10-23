# ABOUTME: Thermodynamics: per-network dE/dt and dT/dt expressions plus the lazily
# ABOUTME: built equation of state (eos) they are derived from

from __future__ import annotations

from functools import cached_property
from typing import TYPE_CHECKING

from sympy import Expr, Float, diff

from .config import ThermodynamicsConfig
from .eos import EosFactory
from .internal_energy import DEDt, InternalEnergy

if TYPE_CHECKING:
    from ...core import Network


class Thermodynamics:
    """Thermal state equations of a network: heating rates, EOS and ``dT/dt``.

    Bound to one :class:`~jaff.core.network.Network`; reachable as
    ``net.thermodynamics``. All attributes are built lazily on first access and
    cached.

    Attributes
    ----------
    config : ThermodynamicsConfig
        The ``[network.thermodynamics]`` settings it was built with (read-only).
    eos : InternalEnergy
        Internal energy for the configured equation of state (``config.eos``;
        an ideal gas with ``gamma = 1.6666666666667`` by default).
    dEdt_chemical : DEDt
        Chemical heating/cooling rate ``Σ_r dE_r F_r`` [erg cm⁻³ s⁻¹].
    dEdt_extra : DEDt
        Non-reactive heating/cooling rate from the ``heatingcoolingrate``
        auxiliary function (``0`` when absent) [erg cm⁻³ s⁻¹].
    dEdt_tot : DEDt
        ``dEdt_chemical + dEdt_extra`` [erg cm⁻³ s⁻¹].
    dTdt_chemical : sympy.Expr
        ``(Ė_chemical − Σ_i ∂E/∂n_i · ṅ_i) / (∂E/∂T)`` [K s⁻¹].
    dTdt_extra : sympy.Expr
        ``Ė_extra / (∂E/∂T)`` [K s⁻¹].
    dTdt_tot : sympy.Expr
        ``dTdt_chemical + dTdt_extra``, built as one fraction [K s⁻¹].

    Notes
    -----
    ``E = E(T, n)`` is the volumetric internal energy :attr:`eos`. The
    temperature rates follow from ``Ė = ∂E/∂T · Ṫ + Σ_i ∂E/∂n_i · ṅ_i``:

    ```
    dTdt_chemical = (Ė_chemical − Σ_i ∂E/∂n_i · ṅ_i) / (∂E/∂T)
    dTdt_extra    =  Ė_extra / (∂E/∂T)
    dTdt_tot      = (Ė_tot − Σ_i ∂E/∂n_i · ṅ_i) / (∂E/∂T)
    ```

    Only reactions change the particle numbers, so the composition term
    ``Σ_i ∂E/∂n_i · ṅ_i`` belongs to the chemical part.

    :class:`~jaff.physics.thermodynamics.internal_energy.DEDt` wraps a
    volumetric rate and offers the same forms as
    :class:`~jaff.physics.thermodynamics.internal_energy.InternalEnergy`
    (``volumetric``, ``specific``, ``per_particle``, ``molar``, or
    ``normaliser(form)``); its normalised forms are time derivatives of the
    normalised energies (quotient rule with ``ρ̇`` and ``ṅ_tot``). Energies or
    rates of exactly the same class and network support ``+`` and ``-``.

    The code generator uses :attr:`eos` to form the energy column of the
    Jacobian by the chain rule ``∂ẋ/∂e = (∂ẋ/∂T) / (∂e/∂T)``, with ``e`` the
    EOS in the chosen normalisation.

    Examples
    --------
    >>> from jaff import Network
    >>> net = Network("COthin")
    >>> th = net.thermodynamics
    >>> net.config.thermodynamics.eos
    IdealEosConfig(type='ideal', gamma=1.6666666666667)
    >>> th.eos.per_particle.simplify()  # 3/2 k_B T [erg]
    2.0709734999999e-16*tgas
    >>> th.dEdt_tot.volumetric  # COthin has no reaction energies or heatingcoolingrate
    0

    Choose another EOS with ``thermodynamics={"eos": ...}`` (or a
    ``[network.thermodynamics.eos]`` table in the config):

    >>> net = Network("COthin", thermodynamics={"eos": {"type": "ideal", "gamma": 1.4}})
    >>> net.thermodynamics.eos.per_particle.simplify()  # k_B T / (gamma - 1)
    3.4516225e-16*tgas
    """

    def __init__(
        self,
        net: Network,
        config: ThermodynamicsConfig,
    ):
        """Bind the thermodynamics to a network and its settings.

        Parameters
        ----------
        net : Network
            Network whose reactions and symbols are used.
        config : ThermodynamicsConfig
            The ``[network.thermodynamics]`` table; :class:`~jaff.core.Network`
            passes ``net.config.thermodynamics``.
        """
        self.net: Network = net
        self.config: ThermodynamicsConfig = config

    @cached_property
    def dEdt_chemical(self) -> DEDt:
        """Chemical heating/cooling rate ``Σ_r dE_r F_r`` [erg cm⁻³ s⁻¹]."""
        return self._get_dEdt_chemical()

    @cached_property
    def dEdt_extra(self) -> DEDt:
        """Non-reactive heating/cooling rate (``heatingcoolingrate``) [erg cm⁻³ s⁻¹]."""
        return self._get_dEdt_extra()

    @cached_property
    def dEdt_tot(self) -> DEDt:
        """Total heating/cooling rate, ``dEdt_chemical + dEdt_extra``."""
        return self.dEdt_chemical + self.dEdt_extra

    @cached_property
    def eos(self) -> InternalEnergy:
        """Internal energy of the network for its configured EOS (built once).

        Returns
        -------
        InternalEnergy
            Built by :class:`EosFactory` from ``config.eos``.
        """
        return EosFactory(self.net, self.config.eos).generate()

    @cached_property
    def _dE_dT(self) -> Expr:
        """``∂E/∂T`` of the EOS volumetric internal energy [erg cm⁻³ K⁻¹]."""
        return diff(self.eos.volumetric, self.net.symbols.tgas)

    @cached_property
    def _composition_rate(self) -> Expr:
        """``Σ_i ∂E/∂n_i · ṅ_i``: energy needed to keep T fixed as n changes.

        Only reactions change particle numbers, so this term belongs to the
        chemical part of ``Ṫ``.
        """
        sym = self.net.symbols
        energy = self.eos.volumetric
        de_dn = [diff(energy, sym.ndens[i]) for i in range(self.net.species.count)]

        return sym.weighted_rate(de_dn)

    @cached_property
    def dTdt_chemical(self) -> Expr:
        """Temperature rate from the chemistry [K s⁻¹].

        ``(Ė_chemical − Σ_i ∂E/∂n_i · ṅ_i) / (∂E/∂T)``: the heat released by
        reactions plus the effect of the reactions changing the particle number
        sharing the thermal energy ``E = E(T, n)``.
        """
        return (self.dEdt_chemical.volumetric - self._composition_rate) / self._dE_dT

    @cached_property
    def dTdt_extra(self) -> Expr:
        """Temperature rate from non-reactive heating/cooling [K s⁻¹].

        ``Ė_extra / (∂E/∂T)``: these mechanisms depend on the local densities but
        do not change them, so there is no composition term.
        """
        return self.dEdt_extra.volumetric / self._dE_dT

    @cached_property
    def dTdt_tot(self) -> Expr:
        """Total temperature rate ``(Ė_tot − Σ_i ∂E/∂n_i · ṅ_i) / (∂E/∂T)`` [K s⁻¹].

        Equal to ``dTdt_chemical + dTdt_extra``, built as a single fraction.
        """
        return (self.dEdt_tot.volumetric - self._composition_rate) / self._dE_dT

    def _get_dEdt_chemical(self) -> DEDt:
        """Build the chemical heating rate from the reactions.

        Returns
        -------
        DEDt
            ``Σ_r dE_r · k_r · Π n_reactants`` as a volumetric rate.
        """
        dEdt = Float(0.0)
        for r in self.net.reactions:
            dEdt += r.dE * r.rate

        return DEDt(dEdt, self.net)

    def _get_dEdt_extra(self) -> DEDt:
        """Build the non-reactive heating rate from ``heatingcoolingrate``.

        Returns
        -------
        DEDt
            The standardised auxiliary-function definition, or zero if absent.
        """
        dEdt = Float(0.0)
        if "heatingcoolingrate" in self.net.functions:
            dEdt = self.net.reactions.resolve_rate_references(
                self.net.functions.value("heatingcoolingrate")
            )

        return DEDt(dEdt, self.net)
