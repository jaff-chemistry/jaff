"""Config of the dust module (``[network.dust]``)."""

from __future__ import annotations

from typing import Annotated, Literal, Self

from pydantic import BeforeValidator, Field, field_validator, model_validator

from ...config.base import ConfigModel, lower

__all__ = [
    "DUST_RV",
    "DustConfig",
    "DustReduction",
    "PhotoelectricEmissionConfig",
]


DustReduction = Annotated[
    Literal["extinction", "absorption", "scattering", "transport", "none"],
    BeforeValidator(lower),
]


class PhotoelectricEmissionConfig(ConfigModel):
    """Photoelectric emission of the dust (``[network.dust.pe]``).

    Attributes
    ----------
    threshold_low : float
        Lower edge of the photoelectric band in eV (grain work function).
    threshold_high : float
        Upper edge of the photoelectric band in eV, above ``threshold_low``.
    """

    threshold_low: float = 6.0
    threshold_high: float = 13.6

    @model_validator(mode="after")
    def _check_band(self) -> Self:
        """The band is not empty: ``threshold_low`` below ``threshold_high``."""
        if self.threshold_low >= self.threshold_high:
            raise ValueError(
                f"threshold_low ({self.threshold_low}) must be below "
                f"threshold_high ({self.threshold_high})"
            )

        return self


#: R_V values with a tabulated dust model.
DUST_RV: tuple[float, ...] = (3.1, 4.0, 5.5)


class DustConfig(ConfigModel):
    """Dust (``[network.dust]``); the table enables the dust module.

    Attributes
    ----------
    rv : float
        Total-to-selective extinction ratio R_V selecting the tabulated dust
        model: ``3.1``, ``4.0`` or ``5.5``.
    u_reduction : str
        Dust cross-section attenuating the radiation density: ``"extinction"``,
        ``"absorption"``, ``"scattering"``, ``"transport"`` or ``"none"``.
    f_reduction : str
        Dust cross-section attenuating the radiation flux (same choices).
    pe : PhotoelectricEmissionConfig
        Photoelectric emission settings.
    """

    rv: float = Field(3.1, json_schema_extra={"enum": list(DUST_RV)})
    u_reduction: DustReduction = "absorption"
    f_reduction: DustReduction = "transport"
    pe: PhotoelectricEmissionConfig = PhotoelectricEmissionConfig()

    @field_validator("rv", mode="after")
    @classmethod
    def _check_rv(cls, rv: float) -> float:
        """R_V must have a tabulated dust model."""
        if rv not in DUST_RV:
            raise ValueError(f"no dust model for R_V = {rv}; choose one of {DUST_RV}")

        return rv
