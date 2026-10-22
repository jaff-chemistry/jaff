# ABOUTME: Photo-reaction line shielding: the configured shielding model multiplies
# ABOUTME: the band rate coefficients of the reaction it is attached to

from unittest.mock import patch

import pytest

from jaff.physics import RadiationProps

H2_PHOTO = "H2 -> H + H  []  PHOTO, 11.2\n"
H2_KEY = "H2._PHOTON__H.H"


@pytest.fixture
def shielded_network(make_network):
    """H2 photodissociation shielded by *shielding*; returns its band rates."""

    def _make(shielding):
        meta = {"reaction_props": {H2_KEY: {"shielding": shielding}}}
        with patch("jaff.physics.photo_reactions._radiation.BackgroundField"):
            net = make_network(
                H2_PHOTO,
                name="h2.jet",
                funcfile=False,
                radiation_props=RadiationProps(bands=[6.0, 13.6]),
                _metadata=meta,
            )

        rxn = net.reactions[H2_KEY, "photo"]
        assert rxn.photo is not None and rxn.photo.rgrouped is not None

        return [band.k for band in rxn.photo.rgrouped]

    return _make


def test_value_key_does_not_bypass_the_shielding_model(shielded_network):
    modelled = shielded_network({"type": "hg2015"})

    assert shielded_network({"type": "hg2015", "value": 0.5}) == modelled
