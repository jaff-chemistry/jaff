# ABOUTME: Thermodynamics: per-network dE/dt and dT/dt expressions plus the lazily
# ABOUTME: built equation of state (eos) they are derived from

from __future__ import annotations

from functools import cached_property
from typing import TYPE_CHECKING

from sympy import Expr, Float, diff

from .eos import EosFactory
from .internal_energy import DEDt, InternalEnergy

if TYPE_CHECKING:
    from ...core import Network


class Thermodynamics:
    """Thermal state equations of a network: heating rates, EOS and ``dT/dt``.

    Bound to one :class:`~jaff.core.network.Network`; reachable as
    ``net.thermodynamics``. All attributes are built lazily and cached.

    Parameters
    ----------
    net : Network
        Network whose reactions, symbols and EOS properties are used.
    dEdt_extra : sympy.Expr, optional
        Stored volumetric non-reactive heating/cooling rate (e.g. restored from a
        ``.jaff`` file). When ``None`` it is built from the ``heatingcoolingrate``
        auxiliary function.

    Attributes
    ----------
    eos : InternalEnergy
        Internal energy for the configured equation of state.
    dEdt_chemical, dEdt_extra, dEdt_tot : DEDt
        Chemical, non-reactive and total heating/cooling rates.
    dTdt_chemical, dTdt_extra, dTdt_tot : sympy.Expr
        Corresponding temperature rates [K s⁻¹].
    """

    def __init__(self, net: Network, dEdt_extra: Expr | None = None):
        self.net: Network = net
        # Stored volumetric heating/cooling rate (e.g. restored from a .jaff file);
        # when None, dEdt_extra is built from the heatingcoolingrate aux function.
        self._stored_dEdt_extra: Expr | None = dEdt_extra

    @cached_property
    def dEdt_chemical(self) -> DEDt:
        """Chemical heating/cooling rate ``Σ_r dE_r F_r`` [erg cm⁻³ s⁻¹]."""
        return self._get_dEdt_chemical()

    @cached_property
    def dEdt_extra(self) -> DEDt:
        """Non-reactive heating/cooling rate (``heatingcoolingrate``) [erg cm⁻³ s⁻¹]."""
        if self._stored_dEdt_extra is not None:
            return DEDt(self._stored_dEdt_extra, self.net)

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
            Built by :class:`EosFactory` from ``net.eos_props``.
        """
        return EosFactory(self.net, self.net.eos_props).generate()

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
            _dEdt = r.dE * r.rate
            for s in r.reactants.core:
                _dEdt *= self.net.symbols.ndens[self.net.species[s.name].index]

            dEdt += _dEdt

        return DEDt(self.net.symbols.standardize(dEdt), self.net)

    def _get_dEdt_extra(self) -> DEDt:
        """Build the non-reactive heating rate from ``heatingcoolingrate``.

        Returns
        -------
        DEDt
            The standardised auxiliary-function definition, or zero if absent.
        """
        dEdt = Float(0.0)
        if "heatingcoolingrate" in self.net.spec.aux_funcs:
            dEdt = self.net.symbols.standardize(
                self.net.spec.aux_funcs["heatingcoolingrate"]["def"]
            )

        return DEDt(dEdt, self.net)
