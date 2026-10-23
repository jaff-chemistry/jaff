# ABOUTME: Radiation.get_group/total_average_cross_section: photon-number-weighted
# ABOUTME: band and full-spectrum mean cross sections, for tabulated and symbolic σ(E)

from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest
import sympy as sp

from jaff.physics.photo_reactions.config import RadiationConfig
from jaff.physics.photo_reactions._radiation import Radiation
from tests.xsec_props import make_xsec_props, stub_photo

BANDS = [1.0, 2.0, 3.0]


def _radiation(profile_index=2, mode="nph"):
    """Radiation over BANDS with the background field mocked out."""
    with patch("jaff.physics.photo_reactions._radiation.BackgroundField"):
        props = RadiationConfig(
            bands=BANDS, profile_index=profile_index, c=1.0, mode=mode
        )
        return Radiation(None, props)


def _tabulated(fn):
    """σ(E) sampled densely over the band range, as (σ, E) arrays."""
    energy = np.linspace(BANDS[0], BANDS[-1], 20001)
    return fn(energy), energy


def test_flat_tabulated_xsec_averages_to_itself():
    rad = _radiation()
    sigma, energy = _tabulated(lambda e: np.full_like(e, 1e-18))

    groups = [rad.get_group_average_cross_section(sigma, energy, g) for g in rad.groups]
    assert [float(x) for x in groups] == pytest.approx([1e-18, 1e-18], rel=1e-6, abs=0)
    assert float(rad.get_total_average_cross_section(sigma, energy)) == pytest.approx(
        1e-18, rel=1e-6, abs=0
    )


def test_flat_symbolic_xsec_averages_to_itself():
    rad = _radiation()
    sigma = sp.Float(1e-18) + 0 * rad.ph_E

    groups = [rad.get_group_average_cross_section(sigma, rad.ph_E, g) for g in rad.groups]
    assert [float(x) for x in groups] == pytest.approx([1e-18, 1e-18], rel=1e-6, abs=0)
    total = rad.get_total_average_cross_section(sigma, rad.ph_E)
    assert float(total) == pytest.approx(1e-18, rel=1e-6, abs=0)


def test_linear_xsec_band_average_is_band_midpoint():
    """Flat photon spectrum (α=2) and σ ∝ E: <σ>_i is σ at the band midpoint."""
    rad = _radiation(profile_index=2)
    sigma, energy = _tabulated(lambda e: 1e-18 * e)

    tab = [rad.get_group_average_cross_section(sigma, energy, g) for g in rad.groups]
    sym = [
        rad.get_group_average_cross_section(1e-18 * rad.ph_E, rad.ph_E, g)
        for g in rad.groups
    ]
    assert [float(x) for x in tab] == pytest.approx([1.5e-18, 2.5e-18], rel=1e-6, abs=0)
    assert [float(x) for x in sym] == pytest.approx([1.5e-18, 2.5e-18], rel=1e-9, abs=0)


def test_tabulated_profile_built_once_per_reaction():
    """The decay total, band and absorption averages share one profile array."""
    rad = _radiation()
    sigma, energy = _tabulated(lambda e: np.full_like(e, 1e-18))
    xsec = make_xsec_props(
        absorption=True, photon_energy=energy, photo_absorption=sigma, photodecay=sigma
    )
    rxn = SimpleNamespace(photo=stub_photo(xsec=xsec))

    with patch.object(
        rad, "get_photden_profile", wraps=rad.get_photden_profile
    ) as profile:
        rxn.photo.band_rates(rad)

    assert profile.call_count == 1


# Flat photon spectrum (α=2) and σ = 1e-18·E over BANDS [1, 2, 3]:
# photon-weighted <σ>_i = σ at the band midpoint; energy-weighted
# <σ>_i = ∫E·σ dE / ∫E dE = 1e-18·(14/9) and 1e-18·(38/15).
PHOTON_AVG = [1.5e-18, 2.5e-18]
ENERGY_AVG = [1e-18 * 14 / 9, 1e-18 * 38 / 15]


def _linear_reaction():
    """Reaction with a tabulated decay cross section σ = 1e-18·E."""
    sigma, energy = _tabulated(lambda e: 1e-18 * e)
    xsec = make_xsec_props(photon_energy=energy, photodecay=sigma)

    return SimpleNamespace(photo=stub_photo(xsec=xsec))


def test_energy_weighted_group_average():
    rad = _radiation(mode="u")
    sigma, energy = _tabulated(lambda e: 1e-18 * e)

    tab = [
        rad.get_group_average_cross_section(sigma, energy, g, energy_weighted=True)
        for g in rad.groups
    ]
    sym = [
        rad.get_group_average_cross_section(
            1e-18 * rad.ph_E, rad.ph_E, g, energy_weighted=True
        )
        for g in rad.groups
    ]
    assert [float(x) for x in tab] == pytest.approx(ENERGY_AVG, rel=1e-6, abs=0)
    assert [float(x) for x in sym] == pytest.approx(ENERGY_AVG, rel=1e-9, abs=0)


def test_energy_mode_band_loss_is_energy_weighted():
    """In "u" mode the energy a band loses weights σ by photon energy."""
    rad = _radiation(mode="u")
    _, _, grouped = _linear_reaction().photo.band_rates(rad)

    losses = [float(g.k / rad.den[i]) for i, g in enumerate(grouped)]
    assert losses == pytest.approx(ENERGY_AVG, rel=1e-6, abs=0)


def test_energy_mode_reaction_rate_stays_photon_weighted():
    """The reaction rate counts photons: photon-weighted <σ>, n_i = u_i / <E>_i."""
    rad = _radiation(mode="u")
    k_tot, _, _ = _linear_reaction().photo.band_rates(rad)

    for i, grp in enumerate(rad.groups):
        coeff = float(k_tot.coeff(rad.den[i]))
        assert coeff == pytest.approx(PHOTON_AVG[i] / grp.eavg, rel=1e-6, abs=0)


def test_photon_mode_band_loss_is_photon_weighted():
    rad = _radiation(mode="nph")
    _, _, grouped = _linear_reaction().photo.band_rates(rad)

    losses = [float(g.k / rad.den[i]) for i, g in enumerate(grouped)]
    assert losses == pytest.approx(PHOTON_AVG, rel=1e-6, abs=0)
