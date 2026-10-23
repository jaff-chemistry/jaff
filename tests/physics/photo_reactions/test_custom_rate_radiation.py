# ABOUTME: Photo reactions with a .jfunc chemRate under radiation: the custom rate is
# ABOUTME: kept but not coupled to the bands, and a warning names each such reaction

import logging

import pytest

from jaff.physics.photo_reactions.config import RadiationConfig

BANDS = [6.0, 13.6, 100.0]
PHOTO = "C -> C+ + E  []  PHOTO, 11.3\n"
CHEMRATE = "@function chemrate0()\n@return 1e-10\n"


def _photo(net):
    return next(r for r in net.reactions if r.type == "photo")


@pytest.fixture
def custom_network(make_network, tmp_path):
    """C photoionisation with a .jfunc chemRate; kwargs go to ``Network``."""
    (tmp_path / "p.jfunc").write_text(CHEMRATE)

    def _make(**kwargs):
        return make_network(PHOTO, name="p.jet", **kwargs)

    return _make


@pytest.fixture
def warnings(caplog):
    """Messages of the JAFF warnings logged while the test runs."""
    logging.disable(logging.NOTSET)  # conftest globally disables logging
    try:
        with caplog.at_level(logging.WARNING, logger="JAFF"):
            yield lambda: [r.message for r in caplog.records]
    finally:
        logging.disable(logging.CRITICAL)


def test_without_radiation_custom_rate_is_kept(custom_network):
    rxn = _photo(custom_network())

    assert rxn.photo.custom_rad_rate is True
    assert float(rxn.rate_coefficient) == pytest.approx(1e-10, rel=1e-6, abs=0)


def test_under_radiation_custom_rate_is_kept(custom_network):
    rxn = _photo(custom_network(radiation_props=RadiationConfig(bands=BANDS)))

    assert rxn.photo.custom_rad_rate is True
    assert float(rxn.rate_coefficient) == pytest.approx(1e-10, rel=1e-6, abs=0)


def test_under_radiation_has_no_band_data(custom_network):
    rxn = _photo(custom_network(radiation_props=RadiationConfig(bands=BANDS)))

    assert rxn.photo.rgrouped is None
    assert rxn.band_xsecs.empty


def test_under_radiation_removes_no_photons(custom_network):
    net = custom_network(radiation_props=RadiationConfig(bands=BANDS))

    assert all(ode.is_zero for ode in net.sradodes(order=2))


def test_under_radiation_warns_naming_the_reaction(custom_network, warnings):
    custom_network(radiation_props=RadiationConfig(bands=BANDS))

    msgs = [m for m in warnings() if "chemRate0" in m]
    assert len(msgs) == 1
    assert "C + _PHOTON -> C+ + e-" in msgs[0]
    assert "not coupled" in msgs[0]
    assert "override_custom_rate" in msgs[0]


def test_without_radiation_no_warning(custom_network, warnings):
    custom_network()

    assert not [m for m in warnings() if "chemRate0" in m]


def test_overridden_custom_rate_does_not_warn(custom_network, warnings):
    custom_network(
        radiation_props=RadiationConfig(bands=BANDS),
        config={"network": {"reactions": {"photo": {"override_custom_rate": True}}}},
    )

    assert not [m for m in warnings() if "chemRate0" in m]
