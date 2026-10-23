# ABOUTME: PhotoelectricEmission: the radiation field in the dust photoelectric
# ABOUTME: band scaled to the reference background field (chi_pe)
"""Photoelectric emission: the scaled radiation field in the photoelectric band.

Photoelectric heating by dust grains is driven by far-UV photons in the
photoelectric band ``[E_low, E_high]``, which defaults to ``[6, 13.6] eV``
(from the grain work function up to the hydrogen ionisation edge) but is
configurable via ``[network.dust.pe]`` ``threshold_low`` / ``threshold_high``.
Rate expressions parametrise this by ``chi_pe`` -- the local radiation field
scaled to a reference (Draine/ISRF) field, integrated over that band.

:class:`PhotoelectricEmission` computes ``chi_pe`` as the ratio of two energy
densities in the photoelectric band:

- **numerator**: the energy density carried by the network's own radiation
  bands (:attr:`Radiation.groups`), summing the fraction of each band that
  overlaps ``[E_low, E_high]``;
- **denominator**: the energy density of the reference background field
  (:attr:`Radiation.background_field`) in the same band.

Both are in erg/cm³, so their ratio is dimensionless; it is then multiplied
by ``Radiation.c / c``, the configured speed of light (``[network.radiation]
c``, possibly a reduced-speed symbol such as ``c_hat``) over the true one,
which is 1 for the default ``c``.  The resulting symbol is
substituted into rate expressions when a network is built with
``[network.dust]`` and radiation enabled.
"""

from __future__ import annotations

from functools import cached_property
from typing import TYPE_CHECKING

from astropy import units as u
from sympy import Expr

from ...common import arr_integrate, smart_integrate
from ..constants import c, h

if TYPE_CHECKING:
    from ...core.network import Network
    from .config import PhotoelectricEmissionConfig


class PhotoelectricEmission:
    """Scaled photoelectric-band radiation field (``chi_pe``) for a network.

    Parameters
    ----------
    network : Network
        The parent network; its :attr:`~jaff.core.network.Network.radiation`
        supplies both the band densities and the reference background field.
    config : PhotoelectricEmissionConfig
        The ``[network.dust.pe]`` config: the photoelectric band edges
        ``threshold_low`` and ``threshold_high``, in eV.

    Attributes
    ----------
    net : Network
        Back-reference to the parent network.
    E_low : astropy.units.Quantity
        Lower edge of the photoelectric band, ``threshold_low`` (grain work
        function, default 6 eV).
    E_high : astropy.units.Quantity
        Upper edge of the photoelectric band, ``threshold_high`` (H ionisation
        edge, default 13.6 eV).
    """

    def __init__(self, network: Network, config: PhotoelectricEmissionConfig):
        """Initialise the photoelectric-emission model.

        Parameters
        ----------
        network : Network
            The parent network.
        config : PhotoelectricEmissionConfig
            The ``[network.dust.pe]`` config: the photoelectric band edges
            ``threshold_low`` and ``threshold_high``, in eV.
        """
        # photoelectric emission activation energy in eV
        self.E_low: u.Quantity = config.threshold_low * u.eV
        # photoelectric emission cutoff energy in eV
        self.E_high: u.Quantity = config.threshold_high * u.eV
        self.net: Network = network

    @cached_property
    def chi(self) -> None | Expr | float:
        """Scaled radiation field in the photoelectric band (``chi_pe``).

        Computes the dimensionless ratio

        ```text
        chi_pe = u_bands(E_low-E_high) / u_background(E_low-E_high) * c_rad / c
        ```

        over the photoelectric band ``[E_low, E_high]`` (``threshold_low`` /
        ``threshold_high``, default 6-13.6 eV), where ``u_bands`` is the energy
        density summed over the network's radiation bands (each band weighted
        by the fraction of its energy that falls in the photoelectric band) and
        ``u_background`` is the energy density of the reference background
        field in the same band.  ``c_rad`` is the radiation's speed of light
        ``Radiation.c`` (``[network.radiation] c``, possibly a symbol such as
        ``c_hat`` for a reduced speed of light) and ``c`` the true one, so the
        factor is 1 for the default ``c``.

        The background field is tabulated as a photon-flux spectrum
        ``I(lambda)`` (photons s⁻¹ cm⁻² nm⁻¹).  Its energy density follows from
        ``u = integral I(lambda) * E_photon / c dlambda = h integral I/lambda
        dlambda`` (the speed of light cancels), evaluated over the band; the
        ``nm -> cm`` factor makes ``h/lambda`` a photon energy in CGS.

        Returns
        -------
        None or sympy.Expr or float
            The ``chi_pe`` expression (symbolic in the band density variables),
            or ``None`` when the network has no radiation field configured.

        Notes
        -----
        Numerator and denominator are both energy densities in erg/cm³, so the
        result is dimensionless.  In photon-density mode each band's symbol is a
        photon number density, converted to an energy density by multiplying by
        the band-average photon energy ``eavg`` (erg).
        """
        rad = self.net.radiation
        if rad is None:
            return None

        num_tot = 0.0
        # Background radiation intensity integrated within the photoelectric band
        # and normalized to get the energy density
        den: Expr | float = arr_integrate(
            rad.background_field.intensity / rad.background_field.wavelength,
            rad.background_field.wavelength,  # wavelength is in nm
            (
                self.E_high.to(u.nm, equivalencies=u.spectral()).value,
                self.E_low.to(u.nm, equivalencies=u.spectral()).value,
            ),
        ) * (h.cgs.value / u.nm.to(u.cm))

        for grp in rad.groups:
            lower = max(grp.lower, self.E_low.value)
            upper = (
                min(grp.upper, self.E_high.value)
                if isinstance(grp.upper, (int, float))
                else max(self.E_high.value, grp.lower)
            )
            if upper < lower:
                upper = lower

            # Each band uses its own spectral profile (per-band profile index).
            energy_frac = smart_integrate(
                grp.energy_profile,
                rad.ph_E,
                (lower, upper),
            ) / smart_integrate(
                grp.energy_profile,
                rad.ph_E,
                (grp.lower, grp.upper),
            )

            num = grp.sym * energy_frac  # type: ignore

            # multiply with average energy in group if number densities are enabled
            if rad.mode == "nph":
                num *= grp.eavg or 0.0

            num_tot += num

        return (num_tot / den) * (rad.c / c.cgs.value)  # chi
