# ABOUTME: Tabular: the tabulated Milky Way dust cross sections for the configured
# ABOUTME: R_V, read from dust.hdf5, and their band averages
from __future__ import annotations

from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np
from sympy import Basic, Expr

from ...common import arr_integrate
from ...config import DATA_DIR
from ...drivers import HDF5
from ...errors import ParserError
from ...types import HDF5Dict

if TYPE_CHECKING:
    from .dust import Dust


class Tabular:
    """Tabulated per-H-nucleus dust cross sections for radiation attenuation.

    Loads the ``mw_rv{rv}`` table from ``dust.hdf5`` (selected by the parent
    :class:`Dust`'s R_V) and exposes the per-hydrogen-nucleus extinction,
    absorption, scattering and transport cross sections together with the
    supporting grids (wavelength, photon energy, albedo, scattering moments).
    :meth:`avg_cross_section_per_hnuc` provides band-averaged cross sections
    weighted by the network's radiation photon-density profile.
    """

    def __init__(self, dust: Dust):
        """Load the tabulated cross sections for the parent dust model.

        Parameters
        ----------
        dust : Dust
            The parent dust model; supplies the R_V used to select the
            ``mw_rv{rv}`` table and the network radiation used for band
            averaging.
        """
        self.dust: Dust = dust
        self._data: HDF5Dict = self._load_cross_sections()

    @cached_property
    def extinction_cross_section_per_hnuc(self) -> np.ndarray:
        """Extinction cross section per H nucleus (cm^2/H)."""
        # units in cm^2/H
        return self._data["C_ext_per_H"]["_data"]

    @cached_property
    def scattering_cross_section_per_hnuc(self) -> np.ndarray:
        """Scattering cross section per H nucleus (cm^2/H)."""
        # units in cm^2/H
        return self.albedo * self.extinction_cross_section_per_hnuc

    @cached_property
    def absorption_cross_section_per_hnuc(self) -> np.ndarray:
        """Absorption cross section per H nucleus (cm^2/H)."""
        # units in cm^2/H
        return (1 - self.albedo) * self.extinction_cross_section_per_hnuc

    @cached_property
    def transport_cross_section_per_hnuc(self) -> np.ndarray:
        """Transport (momentum-transfer) cross section per H nucleus (cm^2/H)."""
        # units in cm^2/H
        return (
            1 - self.albedo * self.scattering_cos
        ) * self.extinction_cross_section_per_hnuc

    @cached_property
    def wavelength(self) -> np.ndarray:
        """Wavelength grid for the tabulated cross sections (micrometers)."""
        # units in micrometer
        return self._data["wavelength"]["_data"]

    @cached_property
    def photon_energy(self) -> np.ndarray:
        """Photon-energy grid (eV), sorted ascending and aligned with the
        cross-section arrays."""
        # units in eV, sorted ascending (aligned with all cross-section arrays)
        return self._data["photon_energy"]["_data"]

    @cached_property
    def albedo(self) -> np.ndarray:
        """Single-scattering albedo (dimensionless)."""
        return self._data["albedo"]["_data"]

    @cached_property
    def scattering_cos(self) -> np.ndarray:
        """Mean scattering-angle cosine <cos theta> (dimensionless)."""
        return self._data["cos"]["_data"]

    def avg_cross_section_per_hnuc(
        self, kind: str, bounds: tuple[float | Basic, float | Basic]
    ) -> Expr | float:
        """Radiation-weighted average cross section per H nucleus over a band.

        Computes the photon-density-weighted mean of the selected cross section
        over the energy band, using the network radiation's photon-density
        profile as the weight.

        Parameters
        ----------
        kind : str
            Cross-section kind to average: ``"extinction"``, ``"absorption"``,
            ``"scattering"`` or ``"transport"``.
        bounds : tuple of (float or sympy.Basic, float or sympy.Basic)
            Lower and upper band edges, in eV.  May be symbolic.

        Returns
        -------
        sympy.Expr or float
            The band-averaged cross section per H nucleus (cm^2/H); symbolic
            when the bounds or radiation profile are symbolic.

        Raises
        ------
        ValueError
            If ``kind`` is not one of the recognised cross-section kinds.
        ParserError
            If the network has no radiation field configured.
        """
        # bounds must be in eV
        types = {
            "extinction": self.extinction_cross_section_per_hnuc,
            "absorption": self.absorption_cross_section_per_hnuc,
            "scattering": self.scattering_cross_section_per_hnuc,
            "transport": self.transport_cross_section_per_hnuc,
        }

        if kind not in types:
            raise ValueError(f"Invalid type of cross-section: {kind}")

        if self.dust.net.radiation is None:
            raise ParserError(
                "Radiation must be enabled to caculate average dust cross-section"
            )

        lower, upper = bounds
        ph_profile = self.dust.net.radiation.get_photden_profile(self.photon_energy)

        return arr_integrate(
            types[kind] * ph_profile,
            self.photon_energy,
            (lower, upper),
        ) / arr_integrate(
            ph_profile,
            self.photon_energy,
            (lower, upper),
        )

    def _load_cross_sections(self) -> HDF5Dict:
        # rv is validated by DustConfig._check_rv; :.1f matches the dust.hdf5 group
        # names (mw_rv3.1 / mw_rv4.0 / mw_rv5.5).
        data: HDF5Dict = HDF5().to_dict(
            f"{DATA_DIR / 'dust' / 'dust.hdf5'}::mw_rv{self.dust.config.rv:.1f}"
        )
        return data
