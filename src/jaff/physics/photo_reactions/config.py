# ABOUTME: Config model of the radiation field ([network.radiation]): bands, spectral
# ABOUTME: index, mode, speed of light, background field, emissivity / absorption plugins
"""Config of the radiation field (``[network.radiation]``)."""

from __future__ import annotations

from typing import Annotated, List, Literal, Self

import astropy.constants as _ac
from pydantic import BeforeValidator, ValidationInfo, field_validator, model_validator

from ...config.base import ConfigModel, lower
from ...config.plugin_table import check_no_copies
from .absorption_coefficient.config import AbsorptionCoefficientConfig
from .emissivity.config import EmissivityConfig

__all__ = [
    "BackgroundField",
    "RadiationConfig",
    "RadiationMode",
]


RadiationMode = Annotated[Literal["nph", "u"], BeforeValidator(lower)]
BackgroundField = Annotated[
    Literal[
        "bb_4000",
        "bb_10000",
        "bb_20000",
        "draine",
        "habing",
        "mathis",
        "solar",
        "tw_hydra",
    ],
    BeforeValidator(lower),
]


class RadiationConfig(ConfigModel):
    """Radiation field (``[network.radiation]``); the table enables radiation.

    Only the shape of the settings is checked here; whether the band-averaged
    spectrum converges is checked by the radiation field itself.

    Attributes
    ----------
    bands : list of (float or "inf")
        Increasing photon-energy band edges in eV, at least two; ``"inf"`` is
        allowed as the last edge only.
    profile_index : float or list of float
        Spectral index *α* of the photon-number spectrum ``n(E) ∝ E^(α-2)``:
        one for every band, or a list with one per band.
    mode : str
        Tracked quantity: ``"nph"`` (photon number density) or ``"u"``
        (energy density).
    c : float or str
        Speed of light in cm/s, or a symbol name (e.g. ``"c_hat"``) for a
        reduced speed of light.
    background_field : str
        Reference radiation field used to scale ``chi_pe``.
    absorption_attenuates_flux : bool
        Whether the ``globalAbsorptionCoefficient`` also attenuates the flux
        equations, not only the density equations.
    emissivity : list[EmissivityConfig]
        Emission sources, each ``{ type = "<plugin>", ... }``: a volume emissivity
        spectrum in ``ph_E`` (erg cm⁻³ s⁻¹ eV⁻¹).  Their spectra are summed with the
        ``.jfunc`` ``globalEmissivity`` and integrated over each band; an entry
        repeating another exactly (same type and settings) is rejected.  Empty by
        default.
    absorption_coefficient : list[AbsorptionCoefficientConfig]
        Absorbers, each ``{ type = "<plugin>", ... }``: an absorption coefficient
        spectrum in ``ph_E`` (cm⁻¹).  Their spectra are summed with the ``.jfunc``
        ``globalAbsorptionCoefficient`` and averaged over each band; an entry repeating
        another exactly is rejected.  Empty by default.
    """

    bands: List[float | Literal["inf"]]
    profile_index: float | List[float] = 0.0
    mode: RadiationMode = "nph"
    c: float | str = float(_ac.c.cgs.value)  # ty: ignore[unresolved-attribute]
    background_field: BackgroundField = "draine"
    absorption_attenuates_flux: bool = True
    emissivity: list[EmissivityConfig] = []
    absorption_coefficient: list[AbsorptionCoefficientConfig] = []

    @field_validator("emissivity", "absorption_coefficient", mode="after")
    @classmethod
    def _no_copies(cls, entries: list, info: ValidationInfo) -> list:
        """Reject an entry repeating another exactly (it would be counted twice)."""
        check_no_copies(entries, str(info.field_name))

        return entries

    @field_validator("bands", mode="after")
    @classmethod
    def _check_bands(cls, bands: List[float | str]) -> List[float | str]:
        """At least two increasing edges, with ``"inf"`` only as the last one."""
        if len(bands) < 2:
            raise ValueError(f"need at least two edges (one band), got {bands}")

        if "inf" in bands[:-1]:
            raise ValueError(f"'inf' is allowed as the last edge only, got {bands}")

        edges = [float("inf") if edge == "inf" else float(edge) for edge in bands]
        if any(low >= high for low, high in zip(edges, edges[1:])):
            raise ValueError(f"band edges must increase, got {bands}")

        return bands

    @model_validator(mode="after")
    def _check_profile_index_per_band(self) -> Self:
        """A list ``profile_index`` has one entry per band."""
        nbands = len(self.bands) - 1
        if isinstance(self.profile_index, list) and len(self.profile_index) != nbands:
            raise ValueError(
                f"profile_index has {len(self.profile_index)} entries but there are "
                f"{nbands} bands; give one per band or a single value"
            )

        return self
