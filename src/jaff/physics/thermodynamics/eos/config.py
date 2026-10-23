"""Config of the equation of state (``[network.thermodynamics.eos]``).

One model per EOS ``type``; :data:`EosConfig` picks the model from ``type``
(case-insensitive, ``"ideal"`` when omitted).
"""

from __future__ import annotations

from typing import Annotated, Any, ClassVar, Dict, Literal, Union

from pydantic import BeforeValidator, Discriminator, Field, Tag

from ....config.base import ConfigModel, lower


#: Adiabatic index; must exceed one.
AdiabaticIndex = Annotated[float, Field(gt=1.0)]


class IdealEosConfig(ConfigModel):
    """Ideal-gas EOS (``type = "ideal"``).

    Attributes
    ----------
    gamma : float
        Adiabatic index, > 1.
    """

    _atomic: ClassVar[bool] = True

    type: Annotated[Literal["ideal"], BeforeValidator(lower)] = "ideal"
    gamma: AdiabaticIndex = 1.6666666666667


class MultiGammaEosConfig(ConfigModel):
    """Ideal gas with a per-species adiabatic index (``type = "multi_gamma"``).

    Attributes
    ----------
    default_gamma : float
        Adiabatic index of species absent from ``gamma_map``, > 1.
    gamma_map : dict[str, float]
        Adiabatic index per species name, each > 1.
    """

    _atomic: ClassVar[bool] = True

    type: Annotated[Literal["multi_gamma"], BeforeValidator(lower)] = "multi_gamma"
    default_gamma: AdiabaticIndex
    gamma_map: Dict[str, AdiabaticIndex]


class FermiDegenerateEosConfig(ConfigModel):
    """Non-relativistic degenerate Fermi gas (``type = "fermi_degenerate"``)."""

    _atomic: ClassVar[bool] = True

    type: Annotated[Literal["fermi_degenerate"], BeforeValidator(lower)] = (
        "fermi_degenerate"
    )


class RelativisticFermiDegenerateEosConfig(ConfigModel):
    """Relativistic degenerate Fermi gas (``type = "relativistic_fermi_degenerate"``)."""

    _atomic: ClassVar[bool] = True

    type: Annotated[Literal["relativistic_fermi_degenerate"], BeforeValidator(lower)] = (
        "relativistic_fermi_degenerate"
    )


def _eos_type(value: Any) -> Any:
    """Return the EOS ``type`` of *value*, lower-cased; ``"ideal"`` when omitted."""
    if isinstance(value, dict):
        eos_type = value.get("type", "ideal")

        return eos_type.lower() if isinstance(eos_type, str) else eos_type

    return getattr(value, "type", None)


EosConfig = Annotated[
    Union[
        Annotated[IdealEosConfig, Tag("ideal")],
        Annotated[MultiGammaEosConfig, Tag("multi_gamma")],
        Annotated[FermiDegenerateEosConfig, Tag("fermi_degenerate")],
        Annotated[
            RelativisticFermiDegenerateEosConfig, Tag("relativistic_fermi_degenerate")
        ],
    ],
    Discriminator(_eos_type),
]
