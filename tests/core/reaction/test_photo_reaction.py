# ABOUTME: PhotoReaction: a photo reaction's component (its config.photo section plus
# ABOUTME: computed state and the photo logic); Reaction only builds it

import pytest

from jaff.core.reaction import PhotoReaction

H_PHOTO = "H -> H+ + E  []  PHOTO, 13.6\nH+ + E -> H  []  1e-12\n"
BANDS = {"bands": [13.6, 20.0, 100.0]}


def _net(make_network, **kwargs):
    return make_network(H_PHOTO, name="h.jet", funcfile=False, **kwargs)


def _photo_reaction(net):
    return next(r for r in net.reactions if r.type == "photo")


class TestComponent:
    def test_photo_reaction_gets_one(self, make_network):
        rxn = _photo_reaction(_net(make_network))

        assert isinstance(rxn.photo, PhotoReaction)
        assert rxn.photo.reaction is rxn

    def test_takes_the_photo_section(self, make_network):
        rxn = _photo_reaction(_net(make_network))

        assert rxn.photo.config is rxn.config.photo

    def test_non_photo_reaction_has_none(self, make_network):
        net = _net(make_network)

        assert next(r for r in net.reactions if r.type != "photo").photo is None

    def test_computes_band_rates_with_radiation(self, make_network):
        rxn = _photo_reaction(_net(make_network, radiation=BANDS))

        assert rxn.photo.bands is not None and len(rxn.photo.bands) == 2
        assert rxn.rate_coefficient.has(rxn.net.radiation.den)


class TestPhotoMethods:
    def test_band_xsecs_on_the_component(self, make_network):
        rxn = _photo_reaction(_net(make_network, radiation=BANDS))

        assert rxn.photo.band_xsecs.columns.tolist() == [
            "lower",
            "upper",
            "eavg",
            "xsec",
            "xsec_frac",
        ]
        assert len(rxn.photo.band_xsecs) == 2

    @pytest.mark.parametrize("name", ["band_xsecs", "plot_xsecs"])
    def test_not_on_reaction(self, make_network, name):
        assert not hasattr(_photo_reaction(_net(make_network)), name)

    def test_plot_xsecs(self, make_network):
        rxn = _photo_reaction(_net(make_network))

        fig, axes = rxn.photo.plot_xsecs(show=False, close=True)

        assert fig is not None


class TestExports:
    def test_photo_reaction_props_is_gone(self):
        import jaff.core
        import jaff.core.reaction

        assert not hasattr(jaff.core, "PhotoReactionProps")
        assert not hasattr(jaff.core.reaction, "PhotoReactionProps")
        assert jaff.core.PhotoReaction is PhotoReaction


class TestPhotochemistryOwner:
    def test_reactions_own_one(self, make_network):
        from jaff.physics import Photochemistry

        net = _net(make_network)

        assert isinstance(net.reactions.photochemistry, Photochemistry)
        assert net.reactions.photochemistry is net.reactions.photochemistry
        assert not hasattr(net, "_photochemistry")

    def test_built_from_the_reactions_config(self, make_network):
        from jaff.physics import Photochemistry

        net = _net(make_network)
        rxn = _photo_reaction(net)

        assert Photochemistry(net.reactions.config).get_xsec(rxn).database == "norad"


class TestXsecTotal:
    def test_set_with_radiation(self, make_network):
        rxn = _photo_reaction(_net(make_network, radiation=BANDS))

        assert float(rxn.photo.xsec_total) > 0.0

    def test_none_without_radiation(self, make_network):
        rxn = _photo_reaction(_net(make_network))

        assert rxn.photo.xsec_total is None

    def test_not_on_the_database_record(self, make_network):
        rxn = _photo_reaction(_net(make_network, radiation=BANDS))

        assert not hasattr(rxn.photo.xsec, "total")


class TestBands:
    def test_each_band_holds_its_group(self, make_network):
        net = _net(make_network, radiation=BANDS)
        bands = _photo_reaction(net).photo.bands

        assert [band.group for band in bands] == list(net.radiation.groups)
        assert all(
            a is b for a, b in zip([band.group for band in bands], net.radiation.groups)
        )

    def test_band_is_frozen(self, make_network):
        from dataclasses import FrozenInstanceError

        from jaff.physics import PhotoBand

        band = _photo_reaction(_net(make_network, radiation=BANDS)).photo.bands[0]

        assert isinstance(band, PhotoBand)
        with pytest.raises(FrozenInstanceError):
            band.k = 0

    def test_rgrouped_is_gone(self, make_network):
        rxn = _photo_reaction(_net(make_network, radiation=BANDS))

        assert not hasattr(rxn.photo, "rgrouped")

    def test_reaction_props_module_is_gone(self):
        import importlib.util

        assert importlib.util.find_spec("jaff.core.reaction.reaction_props") is None

    def test_substitute_replaces_bands(self, make_network):
        import sympy as sp

        rxn = _photo_reaction(_net(make_network, radiation=BANDS))
        before = rxn.photo.bands
        den = before[0].k.free_symbols

        rxn.photo.substitute({s: sp.Integer(1) for s in den})

        assert rxn.photo.bands[0] is not before[0]
        assert rxn.photo.bands[0].group is before[0].group


class TestBandRates:
    def test_match_the_reference_values(self):
        import sympy as sp

        from jaff import Network

        net = Network("h_photoionization", radiation={"bands": [1, 13.6, 100, "inf"]})
        photo = net.reactions.photo_reactions()[0].photo
        band_rate, xsec_total, bands = photo.band_rates(net.radiation)
        den = net.radiation.den

        assert float(band_rate.coeff(den[1])) == pytest.approx(5.85154082319921e-8)
        assert float(band_rate.coeff(den[2])) == pytest.approx(1.1420465158197e-10)
        assert float(xsec_total) == pytest.approx(1.2403886126668026e-19)
        assert [band.group for band in bands] == list(net.radiation.groups)
        assert photo.band_rate == band_rate
        assert isinstance(band_rate, sp.Expr)
