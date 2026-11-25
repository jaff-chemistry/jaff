# ABOUTME: Regression tests for per-band average photon energy (eavg) in Radiation.
# ABOUTME: Each band mean must divide by its OWN photon integral, not the total.

import math
from unittest.mock import patch

import astropy.units as u
import numpy as np
import pytest
import sympy as sp

from jaff.physics.photo_reactions.config import RadiationConfig
from jaff.physics.photo_reactions._radiation import Radiation
from tests.helpers.xsec_props import make_xsec_props, stub_photo

EV2ERG = u.eV.to(u.erg)


def _radiation(bands, profile_index, mode="u"):
    """Build a Radiation with the background field mocked out."""
    with patch("jaff.physics.photo_reactions._radiation.BackgroundField"):
        props = RadiationConfig(
            bands=list(bands),
            profile_index=profile_index,
            mode=mode,
            c=1.0,
            background_field="draine",
        )
        return Radiation(None, props)


def _means_eV(rad):
    return [float(g.eavg) / EV2ERG for g in rad.groups]


def test_flat_photon_spectrum_analytic_means():
    """profile_index=2 => n(E)∝E^0; band mean is the band midpoint."""
    rad = _radiation([1.0, 2.0, 3.0], profile_index=2)
    # Buggy code (÷ total photon integral) gave [0.75, 1.25]; correct is midpoints.
    assert _means_eV(rad) == pytest.approx([1.5, 2.5], rel=1e-9)


def test_flat_energy_spectrum_analytic_mean():
    """profile_index=1 => n(E)∝E^-1; <E> = (E_hi-E_lo)/ln(E_hi/E_lo)."""
    rad = _radiation([1.0, 2.0], profile_index=1)
    expected = (2.0 - 1.0) / math.log(2.0 / 1.0)
    assert _means_eV(rad)[0] == pytest.approx(expected, rel=1e-6)


def test_band_means_lie_within_band_edges():
    """Every band mean must lie within its own [lower, upper] (bug broke this)."""
    rad = _radiation([1.0, 2.0, 3.0, 5.0], profile_index=2)
    for g, mean in zip(rad.groups, _means_eV(rad)):
        assert g.lower <= mean <= g.upper


def test_single_vs_multiband_consistency():
    """A band's mean is independent of how many other bands share the spectrum."""
    solo = _radiation([1.0, 2.0], profile_index=2)
    multi = _radiation([1.0, 2.0, 3.0, 4.0], profile_index=2)
    # The [1,2] band is the first group in both; its mean must match.
    assert _means_eV(solo)[0] == pytest.approx(_means_eV(multi)[0], rel=1e-12)


def test_eavg_independent_of_mode():
    """eavg is a property of the spectrum, identical in 'nph' and 'u' modes."""
    u_mode = _means_eV(_radiation([1.0, 2.0, 3.0], profile_index=2, mode="u"))
    nph_mode = _means_eV(_radiation([1.0, 2.0, 3.0], profile_index=2, mode="nph"))
    assert u_mode == pytest.approx(nph_mode, rel=1e-12)


class _FakeReaction:
    """Minimal reaction with a flat tabulated decay cross section."""

    def __init__(self):
        E = np.linspace(1.0, 3.0, 201)
        sigma = np.full_like(E, 1.0e-18)  # flat 1e-18 cm^2
        self.photo = stub_photo(xsec=make_xsec_props(photon_energy=E, photodecay=sigma))


def test_photon_and_energy_modes_equivalent():
    """For the same physical field, u-mode and nph-mode rates must agree.

    The physical relation between modes is ``u_i = n_i * <E>_i`` with the TRUE
    band mean.  We fix that mean analytically (band midpoints for a flat photon
    spectrum) — independent of the model's ``eavg`` — so the identity only holds
    when the model's ``eavg`` equals the true mean.  A wrong ``eavg`` (the bug)
    would leave a residual factor and break the equality.
    """
    bands = [1.0, 2.0, 3.0]
    true_mean_eV = {0: 1.5, 1: 2.5}  # flat photon spectrum band midpoints

    rad_n = _radiation(bands, profile_index=2, mode="nph")
    rad_u = _radiation(bands, profile_index=2, mode="u")

    k_n, _, _ = _FakeReaction().photo.band_rates(rad_n)
    k_u, _, _ = _FakeReaction().photo.band_rates(rad_u)

    # Map radeden_i -> photden_i * <E>_i using the INDEPENDENT analytic mean
    # (converted to erg to match eavg's units).
    subs = {
        rad_u.den[sp.Idx(g.index)]: rad_n.den[sp.Idx(g.index)]
        * (true_mean_eV[g.index] * EV2ERG)
        for g in rad_u.groups
    }
    converted = k_u.xreplace(subs)
    assert sp.simplify(converted - k_n) == 0


def test_photden_tot_independent_of_band_count():
    """Splitting the same range into more bands must not scale photden_tot."""
    solo = _radiation([1.0, 4.0], profile_index=2)
    multi = _radiation([1.0, 2.0, 3.0, 4.0], profile_index=2)
    assert float(multi.photden_tot) == pytest.approx(float(solo.photden_tot))
    assert float(multi.photden_tot) == pytest.approx(3.0)


def test_photden_tot_is_sum_of_band_integrals_for_per_band_index():
    """Per-band indices: photden_tot = sum of each band's own integral."""
    with patch("jaff.physics.photo_reactions._radiation.BackgroundField"):
        props = RadiationConfig(bands=[1.0, 2.0, 4.0], profile_index=[2.0, 1.0], c=1.0)
        rad = Radiation(None, props)
    # ∫_1^2 E^0 dE + ∫_2^4 E^-1 dE = 1 + ln 2
    assert float(rad.photden_tot) == pytest.approx(1.0 + math.log(2.0))


class TestBandHelpers:
    """Radiation offers the mode-dependent band pieces; photo reactions assemble."""

    def test_rate_assembly_is_not_on_radiation(self):
        assert not hasattr(Radiation, "get_reaction_rate_props")

    def test_photon_density_in_photon_mode(self):
        rad = _radiation([1.0, 2.0, 3.0], profile_index=2, mode="nph")

        assert rad.photon_density(rad.groups[1]) == rad.den[1]

    def test_photon_density_in_energy_mode(self):
        rad = _radiation([1.0, 2.0, 3.0], profile_index=2, mode="u")
        group = rad.groups[1]

        assert rad.photon_density(group) == rad.den[1] / group.eavg

    def test_loss_cross_section_in_photon_mode_is_the_photon_average(self):
        rad = _radiation([1.0, 2.0, 3.0], profile_index=2, mode="nph")
        E = np.linspace(1.0, 3.0, 201)

        assert rad.loss_cross_section(E, E, rad.groups[0], None, 7.0) == 7.0

    def test_loss_cross_section_in_energy_mode_is_energy_weighted(self):
        rad = _radiation([1.0, 2.0, 3.0], profile_index=2, mode="u")
        E = np.linspace(1.0, 3.0, 201)
        group = rad.groups[0]
        expected = rad.get_group_average_cross_section(E, E, group, energy_weighted=True)

        assert rad.loss_cross_section(E, E, group, None, 7.0) == pytest.approx(expected)


class TestSymbols:
    def test_flux_symbol(self):
        rad = _radiation([1.0, 2.0, 3.0], profile_index=2, mode="nph")

        assert rad.flux.label.name == "rflux"
        assert rad.flux.shape == (2,)

    def test_energy_density_in_energy_mode(self):
        rad = _radiation([1.0, 2.0, 3.0], profile_index=2, mode="u")

        assert rad.energy_density(rad.groups[1]) == rad.den[1]

    def test_energy_density_in_photon_mode(self):
        rad = _radiation([1.0, 2.0, 3.0], profile_index=2, mode="nph")
        group = rad.groups[1]

        assert rad.energy_density(group) == rad.den[1] * group.eavg
