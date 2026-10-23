# ABOUTME: Radiation built from its RadiationConfig: band edges, per-band profile indices,
# ABOUTME: divergence checks, and radiation enabled from the network config

import logging
from unittest.mock import patch

import pytest
import sympy as sp

from jaff.physics.photo_reactions._radiation import Radiation
from jaff.physics.photo_reactions.config import RadiationConfig

BACKGROUND = "jaff.physics.photo_reactions._radiation.BackgroundField"


def _radiation(**config):
    """Build a Radiation from RadiationConfig(**config), background field mocked."""
    with patch(BACKGROUND):
        return Radiation(None, RadiationConfig(**config))


class TestRadiationFromConfig:
    def test_inf_edge_becomes_infinity(self):
        rad = _radiation(bands=[6.0, 13.6, "inf"], profile_index=-1.0)

        assert rad.bands == [6.0, 13.6, sp.oo]
        assert rad.nbands == 2

    def test_scalar_profile_index_applies_to_every_band(self):
        rad = _radiation(bands=[1.0, 2.0, 4.0], profile_index=2.0)

        assert [grp.profile_idx for grp in rad.groups] == [2.0, 2.0]

    def test_list_profile_index_sets_each_band(self):
        rad = _radiation(bands=[1.0, 2.0, 4.0], profile_index=[2.0, 1.0])

        assert [grp.profile_idx for grp in rad.groups] == [2.0, 1.0]

    def test_string_c_becomes_symbol(self):
        assert _radiation(bands=[6.0, 13.6], c="c_hat").c == sp.Symbol("c_hat")

    def test_numeric_c_kept(self):
        assert _radiation(bands=[6.0, 13.6], c=3e9).c == pytest.approx(3e9)

    def test_mode_and_absorption_flag(self):
        rad = _radiation(bands=[6.0, 13.6], mode="u", absorption_attenuates_flux=False)

        assert rad.config.mode == "u"
        assert rad.config.absorption_attenuates_flux is False
        assert rad.den.name == "radeden"

    def test_keeps_its_config(self):
        config = RadiationConfig(bands=[6.0, 13.6])
        with patch(BACKGROUND):
            rad = Radiation(None, config)

        assert rad.config is config

    @pytest.mark.parametrize("name", ["mode", "absorption_attenuates_flux"])
    def test_settings_are_not_copied(self, name):
        assert not hasattr(_radiation(bands=[6.0, 13.6]), name)


class TestDivergence:
    def test_zero_lower_edge_checks_first_band_index(self):
        # α0 < 0 → <E> diverges at E→0; last band's index is irrelevant.
        with pytest.raises(RuntimeError, match="bands\\[0\\]"):
            _radiation(bands=[0.0, 1.0, 2.0], profile_index=[-1.0, 2.0], mode="u")

    def test_zero_lower_edge_ok_when_first_band_converges(self):
        _radiation(bands=[0.0, 1.0, 2.0], profile_index=[2.0, -1.0], mode="u")

    def test_inf_upper_edge_checks_last_band_index(self):
        # α_last > 0 → <E> diverges at E→∞; first band's index is irrelevant.
        with pytest.raises(RuntimeError, match="inf"):
            _radiation(bands=[1.0, 2.0, "inf"], profile_index=[-1.0, 2.0], mode="u")

    def test_inf_upper_edge_ok_when_last_band_converges(self):
        _radiation(bands=[1.0, 2.0, "inf"], profile_index=[2.0, -1.0], mode="u")

    def test_low_edge_warning_uses_first_band_index(self):
        logger = logging.getLogger()
        with (
            patch.object(logger, "warning") as warn,
            patch("jaff.physics.photo_reactions._radiation.JaffLogger") as jl,
        ):
            jl.return_value.get_logger.return_value = logger
            _radiation(bands=[0.5, 1.0, 2.0], profile_index=[2.0, 0.0])
            warn.assert_not_called()
            _radiation(bands=[0.5, 1.0, 2.0], profile_index=[0.0, 2.0])
            warn.assert_called_once()


class TestNetworkRadiation:
    NET = "H + H -> H2 [10,1000] 1e-10\nH2 -> H + H [10,1000] 1e-10\n"

    def test_jaff_toml_radiation_table_enables_radiation(self, make_network, tmp_path):
        (tmp_path / "jaff.toml").write_text(
            '[network.radiation]\nbands = [6.0, 13.6, "inf"]\nprofile_index = -1\n'
        )

        with patch(BACKGROUND):
            net = make_network(self.NET, name="n.dat", funcfile=False)

        assert net.radiation is not None
        assert net.radiation.nbands == 2

    def test_no_radiation_table_means_no_radiation(self, make_network):
        assert make_network(self.NET, name="n.dat", funcfile=False).radiation is None

    def test_radiation_keyword_accepts_a_dict(self, make_network):
        with patch(BACKGROUND):
            net = make_network(
                self.NET,
                name="n.dat",
                funcfile=False,
                radiation={"bands": [6.0, 13.6]},
            )

        assert net.radiation is not None
        assert net.radiation.bands == [6.0, 13.6]

    def test_radiation_keyword_beats_jaff_toml(self, make_network, tmp_path):
        (tmp_path / "jaff.toml").write_text("[network.radiation]\nbands = [6.0, 13.6]\n")

        with patch(BACKGROUND):
            net = make_network(
                self.NET,
                name="n.dat",
                funcfile=False,
                radiation=RadiationConfig(bands=[1.0, 2.0, 3.0]),
            )

        assert net.radiation is not None
        assert net.radiation.nbands == 2
