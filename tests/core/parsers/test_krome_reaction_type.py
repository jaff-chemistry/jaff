# ABOUTME: KROME reaction-type inference from the rate: whole-word Av / user_Av match
# ABOUTME: and the case-insensitive, word-bounded RATE_REPS aliases

from __future__ import annotations

import pytest

from jaff.core.parsers.network._formats.krome.reaction import KromeReaction

STATE = {"idx": True, "nreact": 2, "nprod": 2, "tmin": False, "tmax": False, "rate": True}


def _parse(rate: str) -> dict:
    return KromeReaction().parse(f"1,H2,H+,H2+,H,{rate}", 1, STATE, "net.jet")


class TestAvMarksPhoto:
    @pytest.mark.parametrize(
        "rate",
        [
            "7.1d-7*exp(-.5*user_Av)",
            "7.1d-7*exp(-.5*user_av)",
            "7.1d-7*exp(-.5*USER_AV)",
            "1e-10*exp(-2*Av)",
            "1e-10*av",
        ],
    )
    def test_av_word_is_photo(self, rate):
        assert _parse(rate)["type"] == "photo"

    @pytest.mark.parametrize(
        "rate",
        [
            "sumsav*exp(-asav*invT)",  # av inside a word
            "1e-10*user_sigmavel",
            "1e-10*user_avg",
            "1e-10*avx",
        ],
    )
    def test_av_inside_a_word_is_not_photo(self, rate):
        assert _parse(rate)["type"] == "unknown"


class TestRateAliases:
    @pytest.mark.parametrize("name", ["user_Av", "user_av", "USER_AV"])
    def test_user_av_becomes_av_in_any_case(self, name):
        assert _parse(f"1e-10*exp(-{name})")["rate"] == "1e-10*exp(-av)"

    @pytest.mark.parametrize(
        "name", ["user_crflux", "user_CRflux", "user_crate", "User_Crate"]
    )
    def test_cosmic_ray_flux_becomes_crate_in_any_case(self, name):
        assert _parse(f"1e-17*{name}")["rate"] == "1e-17*crate"

    def test_alias_inside_a_longer_name_is_kept(self):
        assert _parse("1e-10*user_avg")["rate"] == "1e-10*user_avg"


class TestPhotonToken:
    """KROME's ``g`` (any case) is the photon: ``_PHOTON``, never a grain."""

    @staticmethod
    def _line(line: str) -> dict:
        state = {**STATE, "nreact": 3, "nprod": 4, "tmin": True, "tmax": True}
        return KromeReaction().parse(line, 1, state, "net.jet")

    @pytest.mark.parametrize("token", ["g", "G"])
    def test_photon_product_of_an_electron_capture_is_a_recombination(self, token):
        fields = self._line(f"1,H+,E,,H,{token},,,NONE,NONE,3.92d-13*invTe**0.6353d0")

        assert fields["products"] == ["H", "_PHOTON"]
        assert fields["type"] == "recombination"

    def test_photon_product_without_an_electron_is_a_radiative_association(self):
        fields = self._line("1,H,H+,,H2+,g,,,NONE,NONE,1d-18")

        assert fields["products"] == ["H2+", "_PHOTON"]
        assert fields["type"] == "radiative_association"

    @pytest.mark.parametrize("token", ["g", "G"])
    def test_photon_reactant_is_photo(self, token):
        fields = self._line(f"1,HE+,{token},,HE++,E,,,NONE,NONE,5d1*user_Trad")

        assert fields["reactants"] == ["He+", "_PHOTON"]
        assert fields["type"] == "photo"

    def test_no_grain_from_g(self):
        fields = self._line("1,H+,E,,H,g,,,NONE,NONE,1d-12")

        assert "_GRAIN" not in fields["reactants"] + fields["products"]
