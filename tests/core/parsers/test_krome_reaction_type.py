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
