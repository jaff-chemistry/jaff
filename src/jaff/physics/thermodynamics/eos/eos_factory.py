# ABOUTME: EosFactory: builds the symbolic InternalEnergy of a network from its
# ABOUTME: EOS config (ideal, multi_gamma and Fermi-degenerate types)

from __future__ import annotations

from typing import TYPE_CHECKING, Dict, cast

from sympy import Float

from ...constants import k_B
from ..internal_energy import InternalEnergy
from .config import EosConfig, IdealEosConfig, MultiGammaEosConfig

if TYPE_CHECKING:
    from ....core import Network


class EosFactory:
    """Build the symbolic :class:`InternalEnergy` selected by an EOS config.

    Each EOS type maps to a builder method through :attr:`_BUILDERS`; the
    valid types must match :data:`~jaff.physics.thermodynamics.eos.config.EOS_TYPES`.
    """

    _BUILDERS: Dict[str, str] = {
        "ideal": "ideal",
        "multi_gamma": "multi_gamma",
        "fermi_degenerate": "fermi_degenerate",
        "relativistic_fermi_degenerate": "relativistic_fermi_degenerate",
    }

    def __init__(self, net: Network, config: EosConfig) -> None:
        """Bind the factory to a network and an EOS configuration.

        Parameters
        ----------
        net : Network
            Network supplying the symbolic densities (``ntot``, ``ndens``,
            ``rho``) and species list.
        config : EosConfig
            The ``[network.thermodynamics.eos]`` config; its ``type`` selects the
            builder.
        """
        self.config: EosConfig = config
        self._net: Network = net

    def generate(self) -> InternalEnergy:
        """Build the EOS selected by ``config.type``.

        Returns
        -------
        InternalEnergy
            Symbolic internal energy of the bound network.
        """
        return getattr(self, self._BUILDERS[self.config.type])()

    def ideal(self) -> InternalEnergy:
        """Volumetric ideal-gas internal energy with a single adiabatic index.

        ``E = n_tot · k_B · T_gas / (γ − 1)`` [erg cm⁻³].

        Returns
        -------
        InternalEnergy
            InternalEnergy wrapping the volumetric internal energy [erg cm⁻³].
        """
        ntot = self._net.symbols.ntot
        tgas = self._net.symbols.tgas
        gamma = cast(IdealEosConfig, self.config).gamma
        e = ntot * k_B.cgs.value * tgas / (gamma - 1.0)

        return InternalEnergy(e, self._net)

    def multi_gamma(self) -> InternalEnergy:
        """Internal energy summed over species with per-species adiabatic indices.

        Each species uses ``config.gamma_map[name]``, falling back to
        ``config.default_gamma``.

        Returns
        -------
        InternalEnergy
            InternalEnergy wrapping the volumetric internal energy [erg cm⁻³].
        """
        e = Float(0.0)
        config = cast(MultiGammaEosConfig, self.config)
        gamma_map = config.gamma_map
        default_gamma = config.default_gamma
        for sp in self._net.species:
            gamma = gamma_map.get(sp.name, default_gamma)
            e += (
                self._net.symbols.ndens[sp.index]
                * k_B.cgs.value
                * self._net.symbols.tgas
                / (gamma - 1.0)
            )

        return InternalEnergy(e, self._net)

    def fermi_degenerate(self) -> InternalEnergy:
        """Non-relativistic degenerate Fermi gas EOS (not implemented).

        Raises
        ------
        NotImplementedError
            Always.
        """
        raise NotImplementedError()

    def relativistic_fermi_degenerate(self) -> InternalEnergy:
        """Relativistic degenerate Fermi gas EOS (not implemented).

        Raises
        ------
        NotImplementedError
            Always.
        """
        raise NotImplementedError()
