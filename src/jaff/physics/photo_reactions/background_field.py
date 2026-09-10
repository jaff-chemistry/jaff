"""Reference background radiation field loaded from the bundled HDF5 store.

A :class:`BackgroundField` is the tabulated spectrum of a named reference
radiation field (``"draine"``, an ISRF variant, a blackbody, ...) used as the
normalisation for scaled-field quantities such as the photoelectric-band
``chi_pe`` (see :class:`~jaff.physics.dust.PhotoelectricEmission`).

The fields live in ``data/background_radiation/radiation.hdf5``, one group per
field, each with co-sorted ``wavelength`` (ascending, nm) and ``intensity``
(photon flux, photons s⁻¹ cm⁻² nm⁻¹) datasets.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ...config import DATA_DIR
from ...drivers import HDF5
from ...drivers.pooch import download_background_radiation

if TYPE_CHECKING:
    import numpy as np


class BackgroundField:
    """Tabulated reference radiation field.

    Parameters
    ----------
    type : str, optional
        Field name -- the HDF5 group to load (default ``"draine"``).

    Attributes
    ----------
    type : str
        The field name that was loaded.
    wavelength : numpy.ndarray
        Wavelength grid in nm, ascending.
    intensity : numpy.ndarray
        Photon-flux spectrum aligned with :attr:`wavelength`, in
        photons s⁻¹ cm⁻² nm⁻¹.
    """

    def __init__(self, type: str = "draine"):
        """Load the named background field from the bundled HDF5 store.

        Parameters
        ----------
        type : str, optional
            Field name / HDF5 group to load, by default ``"draine"``.
        """
        self.type: str = type

        download_background_radiation()

        bgrad = HDF5().to_dict(
            f"{DATA_DIR / 'background_radiation' / 'radiation.hdf5'}::{type}"
        )
        self.wavelength: np.ndarray = bgrad["wavelength"]["_data"]
        self.intensity: np.ndarray = bgrad["intensity"]["_data"]
