# ABOUTME: EosFactory: builds the symbolic InternalEnergy of a network from its
# ABOUTME: EOS config (ideal, multi_gamma and Fermi-degenerate types)

from __future__ import annotations

from typing import TYPE_CHECKING, cast

from sympy import Float

from ...constants import k_B
from ..internal_energy import InternalEnergy
from .config import EosConfig, IdealEosConfig, MultiGammaEosConfig

if TYPE_CHECKING:
    from ....core import Network


class EosFactory:
    """Build the symbolic :class:`InternalEnergy` selected by an EOS config.

    Each EOS type is built by the method of the same name; there is one for
    every type in :data:`~jaff.physics.thermodynamics.eos.config.EOS_TYPES`.
    A network builds its EOS through this factory on first access to
    ``net.thermodynamics.eos``.

    Notes
    -----
    The EOS is configured by the ``[network.thermodynamics.eos]`` table, or
    ``thermodynamics={"eos": {...}}`` passed to :class:`~jaff.core.network.Network`; its
    ``type`` key (case-insensitive, ``"ideal"`` when omitted) selects the model
    and its other keys. Omitted keys take the model's default. All energies are
    volumetric [erg cm⁻³]:

    | ``type`` | Keys | Volumetric energy |
    | --- | --- | --- |
    | ``ideal`` | ``gamma`` (float > 1, default 1.6666666666667) | ``E = n_tot k_B T_gas / (gamma − 1)`` |
    | ``multi_gamma`` | ``default_gamma`` (float > 1), ``gamma_map`` (dict[str, float > 1]) | ``E = Σ_i n_i k_B T_gas / (gamma_i − 1)``, ``gamma_i`` from ``gamma_map`` by species name, else ``default_gamma`` |
    | ``fermi_degenerate`` | -- | not implemented |
    | ``relativistic_fermi_degenerate`` | -- | not implemented |

    The config is validated when the network loads: an unknown ``type``, a
    missing or unexpected key, or an adiabatic index ``<= 1`` is an invalid
    config (:class:`~jaff.errors.ParserError` from ``Network``). The
    not-implemented types raise ``NotImplementedError`` only when the EOS is
    built.

    Examples
    --------
    >>> from jaff import Network
    >>> from jaff.physics import EosFactory
    >>> eos = {
    ...     "type": "multi_gamma",
    ...     "default_gamma": 5 / 3,
    ...     "gamma_map": {"H2": 1.4},
    ... }
    >>> net = Network("COthin", thermodynamics={"eos": eos})
    >>> eos = EosFactory(net, net.config.thermodynamics.eos).generate()
    >>> nden = net.symbols.ndens
    >>> eos.volumetric.coeff(nden[net.species["H2"].index])  # k_B T / (1.4 - 1)
    3.4516225e-16*tgas
    >>> eos.volumetric.coeff(nden[net.species["CO"].index])  # k_B T / (5/3 - 1)
    2.0709735e-16*tgas
    """

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
        return getattr(self, self.config.type)()

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
