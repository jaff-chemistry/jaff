# ABOUTME: SpecieThermodynamics: one species' thermodynamic properties under the network's
# ABOUTME: EOS (its adiabatic index), attached as specie.thermodynamics by Thermodynamics

from __future__ import annotations

from typing import TYPE_CHECKING

from sympy import Expr, Float

from .eos.config import EosConfig, IdealEosConfig, MultiGammaEosConfig

if TYPE_CHECKING:
    from ...core.species import Specie

__all__ = ["SpecieThermodynamics"]


class SpecieThermodynamics:
    """Thermodynamic properties of one species under the network's EOS.

    :class:`~jaff.physics.thermodynamics.Thermodynamics` attaches one to every
    species as ``specie.thermodynamics`` when the EOS describes species
    individually (:meth:`supports`: ``ideal`` and ``multi_gamma``), else
    ``None``. Host codes with their own EOS read the adiabatic indices from here
    (templates with ``REPEAT specie_gammas``), and the ``multi_gamma`` EOS builds
    its energy from them, so the two cannot disagree.

    Attributes
    ----------
    specie : Specie
        The species described.
    config : EosConfig
        The ``[network.thermodynamics.eos]`` table the properties come from.
    """

    def __init__(self, specie: Specie, config: EosConfig) -> None:
        """Bind the properties to a species and the EOS settings.

        Parameters
        ----------
        specie : Specie
            The species described.
        config : EosConfig
            An EOS config this class supports (see :meth:`supports`).

        Raises
        ------
        ValueError
            If the EOS type has no per-species properties.
        """
        if not self.supports(config):
            raise ValueError(f"EOS type {config.type!r} has no per-species properties")
        self.specie: Specie = specie
        self.config: EosConfig = config

    @staticmethod
    def supports(config: EosConfig) -> bool:
        """Whether an EOS config describes species individually (ideal, multi_gamma).

        Parameters
        ----------
        config : EosConfig
            The ``[network.thermodynamics.eos]`` table.

        Returns
        -------
        bool
        """
        return isinstance(config, IdealEosConfig | MultiGammaEosConfig)

    @property
    def gamma(self) -> Expr:
        """Adiabatic index of the species, as an expression (constant for now).

        ``ideal`` gives its ``gamma``; ``multi_gamma`` takes it from
        ``gamma_map`` by species name, else ``default_gamma``. An expression so
        that a temperature-dependent index (e.g. H2's rotational levels) can be
        one of ``tgas`` later.
        """
        config = self.config
        if isinstance(config, MultiGammaEosConfig):
            return Float(config.gamma_map.get(self.specie.name, config.default_gamma))

        assert isinstance(config, IdealEosConfig)
        return Float(config.gamma)

    @property
    def is_constant(self) -> bool:
        """Whether :attr:`gamma` is a number (no dependence on ``tgas`` or others)."""
        return bool(self.gamma.is_number)
