# ABOUTME: EosFactory: builds the symbolic InternalEnergy of a network from its
# ABOUTME: EosProps (ideal, multi_gamma and Fermi-degenerate types)

from __future__ import annotations

from typing import TYPE_CHECKING, Dict

from sympy import Float

from ...constants import k_B
from ..internal_energy import InternalEnergy
from .eos_props import EosProps

if TYPE_CHECKING:
    from ....core import Network


class EosFactory:
    """Build the symbolic :class:`InternalEnergy` selected by an :class:`EosProps`.

    Each EOS type maps to a builder method through :attr:`_BUILDERS`; the
    valid types must match :attr:`EosProps._REQUIRED`.
    """

    _BUILDERS: Dict[str, str] = {
        "ideal": "ideal",
        "multi_gamma": "multi_gamma",
        "fermi_degenerate": "fermi_degenerate",
        "relativistic_fermi_degenerate": "relativistic_fermi_degenerate",
    }

    def __init__(self, net: Network, props: EosProps) -> None:
        """Bind the factory to a network and an EOS configuration.

        Parameters
        ----------
        net : Network
            Network supplying the symbolic densities (``ntot``, ``ndens``,
            ``rho``) and species list.
        props : EosProps
            Validated EOS configuration selecting the builder.
        """
        self.props: EosProps = props
        self._net: Network = net

    def generate(self) -> InternalEnergy:
        """Build the EOS selected by ``props.type``.

        Returns
        -------
        InternalEnergy
            Symbolic internal energy of the bound network.

        Raises
        ------
        ValueError
            If ``props.type`` has no entry in :attr:`_BUILDERS`.
        """
        if self.props.type not in self._BUILDERS:
            raise ValueError(
                f"Invalid eos: '{self.props.type}'. "
                f"Valid eos types are: {', '.join(self._BUILDERS)}"
            )

        return getattr(self, self._BUILDERS[self.props.type])()

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
        e = ntot * k_B.cgs.value * tgas / (self.props.gamma - 1.0)  # type: ignore

        return InternalEnergy(e, self._net)

    def multi_gamma(self) -> InternalEnergy:
        """Internal energy summed over species with per-species adiabatic indices.

        Each species uses ``props.gamma_map[name]``, falling back to
        ``props.default_gamma``.

        Returns
        -------
        InternalEnergy
            InternalEnergy wrapping the volumetric internal energy [erg cm⁻³].
        """
        e = Float(0.0)
        gamma_map = self.props.gamma_map  # type: ignore
        default_gamma = self.props.default_gamma  # type: ignore
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
