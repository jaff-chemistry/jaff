# ABOUTME: override_custom_rate: under radiation, recompute a photo reaction's rate from
# ABOUTME: its cross sections instead of its .jfunc chemRate (global + per-reaction)

from types import SimpleNamespace
from unittest.mock import patch

import pytest

from jaff.cli import JaffGen
from jaff.errors import ParserError
from jaff.physics import RadiationProps
from jaff.physics.photo_reactions._radiation import Radiation

BANDS = [6.0, 13.6, 100.0]
PHOTO = "C -> C+ + E  []  PHOTO, 11.3\n"
CHEMRATE = "@function chemrate0()\n@return 1e-10\n"
KEY = "C._PHOTON__C+.e-"


def _photo(net):
    return next(r for r in net.reactions if r.type == "photo")


@pytest.fixture
def custom_network(make_network, tmp_path):
    """C photoionisation with a .jfunc chemRate; kwargs go to ``Network``."""
    (tmp_path / "p.jfunc").write_text(CHEMRATE)

    def _make(**kwargs):
        return make_network(PHOTO, name="p.jet", **kwargs)

    return _make


class TestRadiationProps:
    def test_default_is_false(self):
        assert RadiationProps(bands=BANDS).override_custom_rate is False

    def test_accepts_true(self):
        props = RadiationProps(bands=BANDS, override_custom_rate=True)
        assert props.override_custom_rate is True

    @pytest.mark.parametrize("value", ["true", 1, None])
    def test_non_bool_rejected(self, value):
        with pytest.raises(ParserError, match="override_custom_rate"):
            RadiationProps(bands=BANDS, override_custom_rate=value)

    def test_radiation_stores_flag(self):
        props = RadiationProps(bands=BANDS, override_custom_rate=True)
        with patch("jaff.physics.photo_reactions._radiation.BackgroundField"):
            rad = Radiation(None, props)
        assert rad.override_custom_rate is True


class TestGlobalFlag:
    def test_override_recomputes_rate_from_xsecs(self, custom_network):
        props = RadiationProps(bands=BANDS, override_custom_rate=True)
        rxn = _photo(custom_network(radiation_props=props))

        assert rxn.photo.custom_rad_rate is False
        assert rxn.rate_coefficient.has(rxn.net.radiation.den[0])
        assert len(rxn.photo.rgrouped) == len(BANDS) - 1

    def test_without_radiation_custom_rate_is_kept(self, custom_network, tmp_path):
        cfg = tmp_path / "jaff.toml"
        cfg.write_text(
            f'[network.reactions.reaction."{KEY}".photo]\noverride_custom_rate = true\n'
        )
        rxn = _photo(custom_network(config=str(cfg)))

        assert rxn.photo.custom_rad_rate is True
        assert float(rxn.rate_coefficient) == pytest.approx(1e-10, rel=1e-6, abs=0)


class TestPerReaction:
    def test_defaults_to_none(self, custom_network):
        assert _photo(custom_network()).photo.override_custom_rate is None

    def test_jaff_toml_override_beats_global(self, custom_network, tmp_path):
        cfg = tmp_path / "jaff.toml"
        cfg.write_text(
            f'[network.reactions.reaction."{KEY}".photo]\noverride_custom_rate = true\n'
        )
        rxn = _photo(
            custom_network(config=str(cfg), radiation_props=RadiationProps(bands=BANDS))
        )

        assert rxn.photo.override_custom_rate is True
        assert rxn.photo.custom_rad_rate is False

    def test_jaffgen_override_beats_jaff_toml(self, custom_network, tmp_path):
        cfg = tmp_path / "jaff.toml"
        cfg.write_text(
            f'[network.reactions.reaction."{KEY}".photo]\noverride_custom_rate = false\n'
        )
        meta = {"reaction_props": {KEY: {"override_custom_rate": True}}}
        rxn = _photo(custom_network(config=str(cfg), _metadata=meta))

        assert rxn.photo.override_custom_rate is True

    def test_non_photo_reaction_raises(self, make_network, tmp_path):
        cfg = tmp_path / "jaff.toml"
        cfg.write_text(
            '[network.reactions.reaction."H.H__H2".photo]\noverride_custom_rate = true\n'
        )
        with pytest.raises(ParserError, match="photo"):
            make_network("H + H -> H2  []  1e-10\n", name="t.jet", config=str(cfg))

    def test_non_bool_raises(self, custom_network, tmp_path):
        cfg = tmp_path / "jaff.toml"
        cfg.write_text(
            f'[network.reactions.reaction."{KEY}".photo]\noverride_custom_rate = "yes"\n'
        )
        with pytest.raises(ParserError, match="override_custom_rate"):
            custom_network(config=str(cfg))


def test_jaffgen_global_flag_reaches_radiation(tmp_path):
    cfg = tmp_path / "jaffgen.toml"
    cfg.write_text(
        "[network.radiation]\nbands = [6, 13.6, 100]\noverride_custom_rate = true\n"
    )
    gen = JaffGen(
        SimpleNamespace(
            network="h_photoionization",
            config=str(cfg),
            label=None,
            funcfile=None,
            duplicate_policy=None,
            expand_nuclei=None,
            errors=None,
            network_config=None,
            outdir=str(tmp_path / "out"),
            indir=None,
            files=None,
            template="microphysics",
        )
    )

    assert gen.net.radiation is not None
    assert gen.net.radiation.override_custom_rate is True
