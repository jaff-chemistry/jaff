# ABOUTME: Radiation.absorbed_average / emitted_average: a quantity in ph_E averaged over the
# ABOUTME: band photons absorbed (σ·n_ph, per band) or emitted (Y/E, over all bands)

import math

import numpy as np
import pytest
import sympy as sp
from scipy.integrate import quad

BANDS = {"bands": [13.6, 20.0, 100.0]}
H_PHOTO = "H -> H+ + E  []  PHOTO, 13.6\nH+ + E -> H  []  1e-12\n"


@pytest.fixture
def radiation(make_network):
    return make_network(H_PHOTO, name="h.jet", radiation=BANDS).radiation


class TestAbsorbed:
    def test_flat_symbolic_cross_section_weights_by_the_photon_spectrum(self, radiation):
        ph_E = radiation.ph_E
        values = radiation.absorbed_average(ph_E, sp.Float(1e-18), ph_E)

        for group, value in zip(radiation.groups, values):
            profile = sp.lambdify(ph_E, group.nph_profile)
            lo, hi = float(group.lower), float(group.upper)
            expected = (
                quad(lambda e: e * profile(e), lo, hi)[0] / quad(profile, lo, hi)[0]
            )
            assert float(value) == pytest.approx(expected, rel=1e-6)

    def test_tabulated_polynomial_keeps_other_symbols(self, radiation):
        ph_E, tgas = radiation.ph_E, sp.Symbol("tgas")
        energy = np.linspace(13.6, 100.0, 400)
        xsec = np.full_like(energy, 1e-18)

        plain = radiation.absorbed_average(ph_E, xsec, energy)
        mixed = radiation.absorbed_average(ph_E + 2 * tgas, xsec, energy)

        for p, m in zip(plain, mixed):
            assert float(sp.expand(m - 2 * tgas)) == pytest.approx(float(p), rel=1e-9)

    def test_tabulated_non_polynomial_mix_is_a_value_error(self, radiation):
        ph_E, tgas = radiation.ph_E, sp.Symbol("tgas")
        energy = np.linspace(13.6, 100.0, 50)

        with pytest.raises(ValueError, match="polynomial in ph_E"):
            radiation.absorbed_average(sp.exp(ph_E / tgas), np.ones_like(energy), energy)

    def test_free_of_ph_E_is_every_band_value(self, radiation):
        assert radiation.absorbed_average(
            sp.Float(3.0), sp.Float(1e-18), radiation.ph_E
        ) == [
            3.0,
            3.0,
        ]


class TestEmitted:
    def test_flat_yield_photons_fall_as_one_over_E(self, radiation):
        ph_E = radiation.ph_E
        flat = sp.Piecewise((1e-11, (ph_E >= 13.6) & (ph_E <= 30)), (0, True))
        log = math.log(30 / 13.6)

        value = radiation.emitted_average(ph_E, flat)

        assert float(value) == pytest.approx(16.4 / log, rel=1e-6)

    def test_no_photons_in_the_bands_is_a_value_error(self, radiation):
        ph_E = radiation.ph_E
        outside = sp.Piecewise((1e-11, ph_E >= 200), (0, True))

        with pytest.raises(ValueError, match="no photons inside the radiation bands"):
            radiation.emitted_average(ph_E, outside)
