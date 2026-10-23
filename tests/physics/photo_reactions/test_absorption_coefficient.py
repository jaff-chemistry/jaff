# ABOUTME: globalAbsorptionCoefficient aux function: network-wide absorption coefficient
# ABOUTME: α(ph_E) [cm^-1], averaged per group and removed from the density equations

from types import SimpleNamespace

import pytest
from pydantic import ValidationError
import sympy as sp

from jaff import Network
from jaff.cli import JaffGen
from jaff.physics.photo_reactions.config import RadiationConfig

BANDS = [6.0, 13.6, 100.0]
NET = "H + H -> H2  []  1e-17\n"
FLAT = "@function globalAbsorptionCoefficient()\n@return 1e-3\n"
LINEAR = "@function globalAbsorptionCoefficient()\n@return 1e-3*ph_E\n"
# Flat photon spectrum (profile_index=2): n(E) ∝ 1 and u(E) ∝ E in every group.
FLAT_PHOTONS = {"profile_index": 2}


def _photon_mean(a: float, b: float) -> float:
    """Photon-weighted mean of α = 1e-3·E over [a, b] for n(E) ∝ 1."""
    return 1e-3 * (a + b) / 2


def _energy_mean(a: float, b: float) -> float:
    """Energy-weighted mean of α = 1e-3·E over [a, b] for u(E) ∝ E."""
    return 1e-3 * (2 / 3) * (b**3 - a**3) / (b**2 - a**2)


@pytest.fixture
def absorbing_network(make_network, tmp_path):
    """Network with *jfunc* as its aux file; kwargs go to ``Network``."""

    def _make(jfunc=FLAT, **kwargs):
        (tmp_path / "a.jfunc").write_text(jfunc)
        return make_network(NET, name="a.jet", **kwargs)

    return _make


def _means(rad):
    return [float(g.absorption_coefficient) for g in rad.groups]


def test_without_radiation_ignored(absorbing_network):
    assert absorbing_network().radiation is None


def test_without_function_is_none(absorbing_network):
    rad = absorbing_network(
        jfunc="", radiation_props=RadiationConfig(bands=BANDS)
    ).radiation

    assert rad.absorption_coefficient is None
    assert all(g.absorption_coefficient is None for g in rad.groups)


def test_spectral_expression_stored(absorbing_network):
    rad = absorbing_network(radiation_props=RadiationConfig(bands=BANDS)).radiation

    assert float(rad.absorption_coefficient) == pytest.approx(1e-3, rel=1e-6, abs=0)


@pytest.mark.parametrize("mode", ["nph", "u"])
def test_flat_coefficient_averages_to_itself(absorbing_network, mode):
    props = RadiationConfig(bands=BANDS, mode=mode, profile_index=0.5)
    rad = absorbing_network(radiation_props=props).radiation

    assert _means(rad) == pytest.approx([1e-3, 1e-3], rel=1e-6, abs=0)


def test_photon_mode_is_photon_weighted(absorbing_network):
    props = RadiationConfig(bands=BANDS, mode="nph", **FLAT_PHOTONS)
    rad = absorbing_network(jfunc=LINEAR, radiation_props=props).radiation

    expected = [_photon_mean(6.0, 13.6), _photon_mean(13.6, 100.0)]
    assert _means(rad) == pytest.approx(expected, rel=1e-6, abs=0)


def test_energy_mode_is_energy_weighted(absorbing_network):
    props = RadiationConfig(bands=BANDS, mode="u", **FLAT_PHOTONS)
    rad = absorbing_network(jfunc=LINEAR, radiation_props=props).radiation

    expected = [_energy_mean(6.0, 13.6), _energy_mean(13.6, 100.0)]
    assert _means(rad) == pytest.approx(expected, rel=1e-6, abs=0)


def test_density_symbols_are_standardized(absorbing_network):
    jfunc = "@function globalAbsorptionCoefficient(n_H)\n@return 1e-21*n_H\n"
    net = absorbing_network(jfunc=jfunc, radiation_props=RadiationConfig(bands=BANDS))

    nden = net.symbols.ndens[net.species["H"].index]
    for group in net.radiation.groups:
        assert group.absorption_coefficient.has(nden)
        assert not group.absorption_coefficient.has(sp.Symbol("n_H"))


def test_symbols_reported_in_network_variables(absorbing_network):
    jfunc = "@function globalAbsorptionCoefficient(tgas)\n@return 1e-6*tgas\n"
    net = absorbing_network(jfunc=jfunc, radiation_props=RadiationConfig(bands=BANDS))

    assert net.symbols.tgas in net.symbols.variables


def _absorption_losses(net):
    """Expected density and flux losses -c<α>_i·den_i and -c<α>_i·rflux_i."""
    rad = net.radiation
    rflux = sp.IndexedBase("rflux", shape=(rad.nbands,))
    density = [
        -rad.c * g.absorption_coefficient * rad.den[i] for i, g in enumerate(rad.groups)
    ]
    flux = [
        -rad.c * g.absorption_coefficient * rflux[i] for i, g in enumerate(rad.groups)
    ]

    return density, flux


def test_absorbs_from_density_and_flux_by_default(absorbing_network):
    net = absorbing_network(radiation_props=RadiationConfig(bands=BANDS))
    assert net.radiation.absorption_attenuates_flux is True
    odes = net.sradodes(order=2)
    nb = net.radiation.nbands
    density, flux = _absorption_losses(net)

    for i in range(nb):
        assert sp.simplify(odes[i] - density[i]) == 0
        assert sp.simplify(odes[nb + i] - flux[i]) == 0


def test_flag_off_leaves_flux_unchanged(absorbing_network):
    props = RadiationConfig(bands=BANDS, absorption_attenuates_flux=False)
    net = absorbing_network(radiation_props=props)
    odes = net.sradodes(order=2)
    nb = net.radiation.nbands
    density, _ = _absorption_losses(net)

    for i in range(nb):
        assert sp.simplify(odes[i] - density[i]) == 0
        assert odes[nb + i].is_zero


@pytest.mark.parametrize("value", ["false", 0, None])
def test_flag_must_be_bool(value):
    with pytest.raises(ValidationError, match="absorption_attenuates_flux"):
        RadiationConfig(bands=BANDS, absorption_attenuates_flux=value)


def test_jaffgen_flag_reaches_radiation(tmp_path):
    cfg = tmp_path / "jaffgen.toml"
    cfg.write_text(
        "[network.radiation]\nbands = [6, 13.6, 100]\n"
        "absorption_attenuates_flux = false\n"
    )
    gen = JaffGen(
        SimpleNamespace(
            network="h_photoionization",
            config=str(cfg),
            label=None,
            funcfile=None,
            duplicate_policy=None,
            expand_nuclei=None,
            strict=None,
            network_config=None,
            outdir=str(tmp_path / "out"),
            indir=None,
            files=None,
            template="microphysics",
        )
    )

    assert gen.net.radiation is not None
    assert gen.net.radiation.absorption_attenuates_flux is False


def test_jaff_round_trip_keeps_coefficient(absorbing_network, tmp_path):
    net = absorbing_network(radiation_props=RadiationConfig(bands=BANDS))
    net.to_jaff(tmp_path / "a.jaff")

    back = Network(
        str(tmp_path / "a.jaff"),
        funcfile=False,
        radiation_props=RadiationConfig(bands=BANDS),
    )

    assert back.radiation is not None
    assert back.radiation.absorption_coefficient is not None
    assert float(back.radiation.absorption_coefficient) == pytest.approx(
        1e-3, rel=1e-6, abs=0
    )
    assert _means(back.radiation) == pytest.approx(_means(net.radiation), rel=1e-6, abs=0)
