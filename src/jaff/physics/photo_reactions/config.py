# ABOUTME: Config model of the radiation field ([network.radiation]): bands,
# ABOUTME: spectral index, mode, speed of light, background field
"""Config of the radiation field (``[network.radiation]``)."""

from __future__ import annotations

from typing import Annotated, List, Literal, Self

import astropy.constants as _ac
from pydantic import BeforeValidator, field_validator, model_validator

from ...config.base import ConfigModel, lower

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
    """

    bands: List[float | Literal["inf"]]
    profile_index: float | List[float] = 0.0
    mode: RadiationMode = "nph"
    c: float | str = float(_ac.c.cgs.value)  # ty: ignore[unresolved-attribute]
    background_field: BackgroundField = "draine"
    absorption_attenuates_flux: bool = True

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
