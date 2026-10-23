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

        assert rxn.photo.rgrouped is not None and len(rxn.photo.rgrouped) == 2
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
