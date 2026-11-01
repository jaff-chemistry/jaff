# ABOUTME: Dust built from its DustConfig (R_V, reductions, photoelectric band) and
# ABOUTME: enabled from the network config ([network.dust] table or the dust keyword)

from unittest.mock import MagicMock, patch

import astropy.units as u
import pytest

from jaff.physics import Dust
from jaff.physics.dust.config import DustConfig

TABULAR = "jaff.physics.dust.dust.Tabular"
NET = "H + H -> H2 [10,1000] 1e-10\nH2 -> H + H [10,1000] 1e-10\n"


def _dust(**config):
    """Build a Dust from DustConfig(**config), its tabulated data mocked out."""
    with patch(TABULAR):
        return Dust(MagicMock(), DustConfig.model_validate(config))


class TestDustFromConfig:
    def test_defaults(self):
        dust = _dust()

        assert dust.config.rv == pytest.approx(3.1)
        assert dust.config.u_reduction == "absorption"
        assert dust.config.f_reduction == "transport"

    def test_settings_are_used(self):
        dust = _dust(rv=4.0, u_reduction="none", f_reduction="extinction")

        assert dust.config.rv == pytest.approx(4.0)
        assert dust.config.u_reduction == "none"
        assert dust.config.f_reduction == "extinction"

    def test_pe_gets_its_section(self):
        dust = _dust(pe={"threshold_low": 5.0})

        assert dust.pe.config is dust.config.pe

    @pytest.mark.parametrize("name", ["rv", "u_reduction", "f_reduction"])
    def test_settings_are_not_copied(self, name):
        assert not hasattr(_dust(), name)

    def test_photoelectric_band_from_pe_section(self):
        pe = _dust(pe={"threshold_low": 5.0}).pe

        assert pe.E_low.to_value(u.eV) == pytest.approx(5.0)
        assert pe.E_high.to_value(u.eV) == pytest.approx(13.6)


class TestNetworkDust:
    def test_dust_table_enables_dust(self, make_network, tmp_path):
        (tmp_path / "jaff.toml").write_text("[network.dust]\nrv = 5.5\n")

        with patch(TABULAR):
            net = make_network(NET, name="n.dat", funcfile=False)

        assert isinstance(net.dust, Dust)
        assert net.dust.config.rv == pytest.approx(5.5)

    def test_no_dust_table_means_no_dust(self, make_network):
        assert make_network(NET, name="n.dat", funcfile=False).dust is None

    def test_dust_keyword_accepts_a_dict(self, make_network):
        with patch(TABULAR):
            net = make_network(
                NET,
                name="n.dat",
                funcfile=False,
                dust={"pe": {"threshold_high": 12.0}},
            )

        assert net.dust is not None
        assert net.dust.pe.E_high.to_value(u.eV) == pytest.approx(12.0)

    def test_dust_keyword_beats_jaff_toml(self, make_network, tmp_path):
        (tmp_path / "jaff.toml").write_text(
            '[network.dust]\nrv = 5.5\nf_reduction = "none"\n'
        )

        with patch(TABULAR):
            net = make_network(NET, name="n.dat", funcfile=False, dust=DustConfig(rv=4.0))

        assert net.dust is not None
        assert net.dust.config.rv == pytest.approx(4.0)
        assert net.dust.config.f_reduction == "none"


class TestAttenuation:
    BANDS = {"bands": [6.0, 13.6, 100.0]}

    def _net(self, make_network, **dust):
        return make_network(
            NET, name="d.jet", funcfile=False, radiation=self.BANDS, dust=dust
        )

    def test_loss_terms(self, make_network):
        net = self._net(make_network)
        rad, group = net.radiation, net.radiation.groups[1]
        bounds = (group.lower, group.upper)
        avg = net.dust.tabular.avg_cross_section_per_hnuc
        base = net.symbols.zd * rad.c * net.symbols.n_hnuc

        density_loss, flux_loss = net.dust.attenuation(group)

        assert density_loss == base * rad.den[1] * avg("absorption", bounds)
        assert flux_loss == base * rad.flux[1] * avg("transport", bounds)

    def test_none_drops_the_term(self, make_network):
        net = self._net(make_network, u_reduction="none", f_reduction="none")

        losses = net.dust.attenuation(net.radiation.groups[0])

        assert [float(loss) for loss in losses] == pytest.approx([0.0, 0.0])

    def test_radiation_subtracts_it(self, make_network):
        net = self._net(make_network)
        group = net.radiation.groups[1]
        density_loss, flux_loss = net.dust.attenuation(group)
        ei, fi = net.radiation.ordered_index(group.index, 0)
        odes = net.radiation.odes()

        # No photo reactions or emission: dust is the only term.
        assert odes[ei] == -density_loss
        assert odes[fi] == -flux_loss
