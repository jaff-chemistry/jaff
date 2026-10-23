# ABOUTME: Config models of the equation of state ([network.thermodynamics.eos]),
# ABOUTME: one per EOS type, selected by type
"""Config of the equation of state (``[network.thermodynamics.eos]``).

One model per EOS ``type``; :data:`EosConfig` picks the model from ``type``
(case-insensitive, ``"ideal"`` when omitted).  The two Fermi-degenerate types
are accepted but not implemented yet: building their EOS raises
``NotImplementedError``.
"""

from __future__ import annotations

from typing import Annotated, Any, ClassVar, Dict, Literal, Union

from pydantic import BeforeValidator, Discriminator, Field, Tag

from ....config.base import ConfigModel, lower

__all__ = [
    "AdiabaticIndex",
    "EOS_TYPES",
    "EosConfig",
    "FermiDegenerateEosConfig",
    "IdealEosConfig",
    "MultiGammaEosConfig",
    "RelativisticFermiDegenerateEosConfig",
]


#: Adiabatic index; must exceed one.
AdiabaticIndex = Annotated[float, Field(gt=1.0)]


class IdealEosConfig(ConfigModel):
    """Ideal-gas EOS (``type = "ideal"``).

    Attributes
    ----------
    type : str
        ``"ideal"`` (case-insensitive); the EOS used when ``type`` is omitted.
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
    type : str
        ``"multi_gamma"`` (case-insensitive); required.
    default_gamma : float
        Adiabatic index of species absent from ``gamma_map``, > 1.
    gamma_map : dict[str, float]
        Adiabatic index per species name, each > 1.
    """

    _atomic: ClassVar[bool] = True

    @classmethod
    def _extend_json_schema(cls, schema: Dict[str, Any]) -> None:
        """Mark ``type`` required: only the ideal EOS may omit it."""
        schema["required"] = [*schema.get("required", []), "type"]

    type: Annotated[Literal["multi_gamma"], BeforeValidator(lower)] = "multi_gamma"
    default_gamma: AdiabaticIndex
    gamma_map: Dict[str, AdiabaticIndex]


class FermiDegenerateEosConfig(ConfigModel):
    """Non-relativistic degenerate Fermi gas (``type = "fermi_degenerate"``).

    Not implemented yet: the config is accepted, but building the EOS raises
    ``NotImplementedError``.

    Attributes
    ----------
    type : str
        ``"fermi_degenerate"`` (case-insensitive); required.
    """

    _atomic: ClassVar[bool] = True

    @classmethod
    def _extend_json_schema(cls, schema: Dict[str, Any]) -> None:
        """Mark ``type`` required: only the ideal EOS may omit it."""
        schema["required"] = [*schema.get("required", []), "type"]

    type: Annotated[Literal["fermi_degenerate"], BeforeValidator(lower)] = (
        "fermi_degenerate"
    )


class RelativisticFermiDegenerateEosConfig(ConfigModel):
    """Relativistic degenerate Fermi gas (``type = "relativistic_fermi_degenerate"``).

    Not implemented yet: the config is accepted, but building the EOS raises
    ``NotImplementedError``.

    Attributes
    ----------
    type : str
        ``"relativistic_fermi_degenerate"`` (case-insensitive); required.
    """

    _atomic: ClassVar[bool] = True

    @classmethod
    def _extend_json_schema(cls, schema: Dict[str, Any]) -> None:
        """Mark ``type`` required: only the ideal EOS may omit it."""
        schema["required"] = [*schema.get("required", []), "type"]

    type: Annotated[Literal["relativistic_fermi_degenerate"], BeforeValidator(lower)] = (
        "relativistic_fermi_degenerate"
    )


def _eos_type(value: Any) -> Any:
    """Return the EOS ``type`` of *value*, lower-cased; ``"ideal"`` when omitted."""
    if isinstance(value, dict):
        eos_type = value.get("type", "ideal")

        return eos_type.lower() if isinstance(eos_type, str) else eos_type

    return getattr(value, "type", None)


#: EOS ``type`` values, one per model.
EOS_TYPES: tuple[str, ...] = (
    "ideal",
    "multi_gamma",
    "fermi_degenerate",
    "relativistic_fermi_degenerate",
)


def _check_eos_type(value: Any) -> Any:
    """Reject an unknown EOS ``type`` naming the choices, before the union is tried."""
    eos_type = _eos_type(value)
    if isinstance(value, dict) and eos_type not in EOS_TYPES:
        raise ValueError(
            f"unknown EOS type {eos_type!r}; choose one of {', '.join(EOS_TYPES)}"
        )

    return value


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
    BeforeValidator(_check_eos_type),
]
