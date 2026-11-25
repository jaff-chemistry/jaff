# ABOUTME: reaction.energy (ReactionEnergy): the energy change per event and its dE/dt term;
# ABOUTME: a deltaE in ph_E is averaged over the band photons absorbed or emitted

import math

import pytest
import sympy as sp
from scipy.integrate import quad

from jaff.errors import ParserError

H_PHOTO = "H -> H+ + E  []  PHOTO, 13.6\nH+ + E -> H  []  1e-12\n"
BANDS = {"bands": [13.6, 20.0, 100.0]}
EV = 1.602176634e-12
KB = 1.380649e-16

FLAT_XSEC = "@function photoCrossSection0()\n@return 1e-18\n"
ABOVE_THRESHOLD = "@function deltaE0(ph_E)\n@return (ph_E - 13.6) * 1.602176634e-12\n"


def _flat_yield(index, top):
    # Energy emitted per eV, flat on [13.6, top] eV: photons per eV fall as 1/E.
    return (
        f"@function spectralYield{index}()\n"
        f"@return Piecewise((1e-11, (ph_E >= 13.6) & (ph_E <= {top})), (0, True))\n"
    )


# The recombination (reaction 1) loses the photon's energy above 13.6 eV.
RECOMBINATION_LOSS = "@function deltaE1(ph_E)\n@return -(ph_E - 13.6) * 1.602176634e-12\n"


def _mean_loss(top):
    """<-(E - 13.6)> over photons N(E) ∝ 1/E on [13.6, top] (eV), in erg."""
    log = math.log(top / 13.6)

    return -((top - 13.6) - 13.6 * log) / log * EV


def _network(make_network, tmp_path, jfunc, **kwargs):
    (tmp_path / "h.jfunc").write_text(jfunc)

    return make_network(H_PHOTO, name="h.jet", **kwargs)


def _photo_energy(net):
    return net.reactions[0].energy


class TestPerEvent:
    def test_constant_deltaE(self, make_network, tmp_path):
        net = _network(make_network, tmp_path, "@function deltaE0()\n@return 6.4e-12\n")
        energy = _photo_energy(net)

        assert float(energy.per_event) == pytest.approx(6.4e-12)
        assert energy.photon == "none"
        assert (
            sp.simplify(energy.change() - 6.4e-12 * net.reactions[0].rate.expanded) == 0
        )

    def test_no_energy_change_is_zero(self, make_network, tmp_path):
        net = _network(make_network, tmp_path, "")

        assert float(net.reactions[1].energy.per_event) == 0.0
        assert net.reactions[1].energy.change() == 0

    def test_reaction_dE_is_gone(self, make_network, tmp_path):
        net = _network(make_network, tmp_path, "")

        assert not hasattr(net.reactions[0], "dE")
        assert not hasattr(net.reactions, "dE")


class TestAbsorbed:
    def test_flat_cross_section_averages_over_the_photon_spectrum(
        self, make_network, tmp_path
    ):
        net = _network(
            make_network, tmp_path, FLAT_XSEC + ABOVE_THRESHOLD, radiation=BANDS
        )
        energy = _photo_energy(net)
        ph_E = net.symbols.ph_E

        assert energy.photon == "absorbed"
        for band, value in zip(net.reactions[0].photo.bands, energy.band_values):
            profile = sp.lambdify(ph_E, band.group.nph_profile)
            lo, hi = float(band.group.lower), float(band.group.upper)
            mean_E = quad(lambda e: e * profile(e), lo, hi)[0] / quad(profile, lo, hi)[0]

            assert float(value) == pytest.approx((mean_E - 13.6) * EV, rel=1e-6)

    def test_tabulated_cross_section_average_lies_in_the_band(
        self, make_network, tmp_path
    ):
        net = _network(make_network, tmp_path, ABOVE_THRESHOLD, radiation=BANDS)
        bands = net.reactions[0].photo.bands

        for band, value in zip(bands, _photo_energy(net).band_values):
            lo, hi = float(band.group.lower), float(band.group.upper)
            assert (lo - 13.6) * EV <= float(value) <= (hi - 13.6) * EV

    def test_other_symbols_stay_symbolic(self, make_network, tmp_path):
        # Polynomial in ph_E with a tgas term: the tgas part passes through.
        plain = _network(make_network, tmp_path, ABOVE_THRESHOLD, radiation=BANDS)
        mixed = _network(
            make_network,
            tmp_path,
            "@function deltaE0(ph_E, tgas)\n"
            "@return (ph_E - 13.6) * 1.602176634e-12 + 1.380649e-16 * tgas\n",
            radiation=BANDS,
        )
        tgas = mixed.symbols.tgas

        for p, m in zip(
            _photo_energy(plain).band_values, _photo_energy(mixed).band_values
        ):
            assert float(sp.expand(m - KB * tgas)) == pytest.approx(float(p), rel=1e-9)

    def test_change_sums_band_rates_times_band_energies(self, make_network, tmp_path):
        net = _network(make_network, tmp_path, ABOVE_THRESHOLD, radiation=BANDS)
        reaction = net.reactions[0]
        energy = reaction.energy
        expected = (
            sum(b.rate * v for b, v in zip(reaction.photo.bands, energy.band_values))
            * reaction.rate.densities
        )

        for form in ("expanded", "reduced"):
            assert sp.simplify(energy.change(form) - expected) == 0

    def test_ph_E_reaches_neither_dEdt_nor_the_solver_inputs(
        self, make_network, tmp_path
    ):
        net = _network(make_network, tmp_path, ABOVE_THRESHOLD, radiation=BANDS)
        ph_E = net.symbols.ph_E

        for form in ("expanded", "reduced"):
            assert (
                ph_E not in net.thermodynamics.dEdt_chemical(form).volumetric.free_symbols
            )
        system = net.ode_system(thermal="dedt", radiation=True)
        assert ph_E not in set(system.parameters)

    def test_band_rates_sum_to_the_rate_coefficient(self, make_network, tmp_path):
        net = _network(make_network, tmp_path, ABOVE_THRESHOLD, radiation=BANDS)
        photo = net.reactions[0].photo

        assert sp.simplify(sum(b.rate for b in photo.bands) - photo.band_rate) == 0


class TestNoPhotonDistribution:
    def test_without_radiation(self, make_network, tmp_path):
        with pytest.raises(ParserError, match="ph_E"):
            _network(make_network, tmp_path, ABOVE_THRESHOLD)

    def test_on_a_reaction_without_photons(self, make_network, tmp_path):
        jfunc = "@function deltaE1(ph_E)\n@return ph_E * 1.602176634e-12\n"
        with pytest.raises(ParserError, match="ph_E"):
            _network(make_network, tmp_path, jfunc, radiation=BANDS)


class TestEmitted:
    def test_averages_over_the_emitted_photons(self, make_network, tmp_path):
        net = _network(
            make_network,
            tmp_path,
            _flat_yield(1, 30) + RECOMBINATION_LOSS,
            radiation=BANDS,
        )
        energy = net.reactions[1].energy

        assert energy.photon == "emitted"
        assert float(energy.mean) == pytest.approx(_mean_loss(30), rel=1e-6)

    def test_only_photons_inside_the_bands_count(self, make_network, tmp_path):
        # The yield runs to 150 eV; the bands stop at 100 eV.
        net = _network(
            make_network,
            tmp_path,
            _flat_yield(1, 150) + RECOMBINATION_LOSS,
            radiation=BANDS,
        )

        assert float(net.reactions[1].energy.mean) == pytest.approx(
            _mean_loss(100), rel=1e-6
        )

    def test_change_is_the_mean_times_the_rate(self, make_network, tmp_path):
        net = _network(
            make_network,
            tmp_path,
            _flat_yield(1, 30) + RECOMBINATION_LOSS,
            radiation=BANDS,
        )
        reaction = net.reactions[1]
        ph_E = net.symbols.ph_E

        for form in ("expanded", "reduced"):
            change = reaction.energy.change(form)
            assert change == reaction.energy.mean * reaction.rate.form(form)
            assert (
                ph_E not in net.thermodynamics.dEdt_chemical(form).volumetric.free_symbols
            )

    def test_no_photons_inside_the_bands(self, make_network, tmp_path):
        jfunc = (
            "@function spectralYield1()\n"
            "@return Piecewise((1e-11, ph_E >= 200), (0, True))\n" + RECOMBINATION_LOSS
        )
        with pytest.raises(ParserError, match="no photons inside the radiation bands"):
            _network(make_network, tmp_path, jfunc, radiation=BANDS)

    def test_a_photo_reaction_that_also_emits_is_ambiguous(self, make_network, tmp_path):
        with pytest.raises(ParserError, match="absorbs and emits"):
            _network(
                make_network,
                tmp_path,
                _flat_yield(0, 30) + ABOVE_THRESHOLD,
                radiation=BANDS,
            )

    def test_an_emitting_reaction_needs_radiation(self, make_network, tmp_path):
        with pytest.raises(ParserError, match="no radiation field"):
            _network(make_network, tmp_path, _flat_yield(1, 30) + RECOMBINATION_LOSS)
