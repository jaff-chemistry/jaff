# ABOUTME: End-to-end test that the dust + radiation modules enable and build
# ABOUTME: on a minimal H photoionization network driven by a jaffgen.toml config

from types import SimpleNamespace

import pytest
import sympy as sp

from jaff.cli import JaffGen
from jaff.errors import ParserError
from jaff.physics import Dust, Radiation
from tests._paths import REPO_ROOT

REPO = REPO_ROOT
H_PHOTO = REPO / "networks" / "h_photoionization" / "h_photo.jet"

# Mirrors the commented [network.radiation] / [network.dust] blocks shipped in
# src/jaff/templates/generator/microphysics/jaffgen.toml.
CONFIG_TOML = """\
[network.radiation]
bands = [6, 11.2, 13.6]
profile_index = 0
mode = "nph"
c = "c_hat"
background_field = "draine"

[network.reactions.photo]
use_proxy_photoreaction = true

[network.dust]
rv = 3.1
u_reduction = "extinction"
f_reduction = "extinction"

[network.dust.pe]
threshold_low = 6
threshold_high = 13.6
"""


def _args(config, outdir):
    """CLI arg namespace driving jaffgen in-process (see test_codegen_e2e)."""
    return SimpleNamespace(
        network=str(H_PHOTO),
        config=str(config),
        label=None,
        funcfile=None,
        duplicate_policy=None,
        expand_nuclei=None,
        strict=None,
        network_config=None,
        outdir=str(outdir),
        indir=None,
        files=None,
        template="microphysics",
    )


@pytest.fixture
def gen(tmp_path):
    """Run the full jaffgen pipeline with dust + radiation enabled."""
    config = tmp_path / "jaffgen.toml"
    config.write_text(CONFIG_TOML)
    return JaffGen(_args(config, tmp_path / "out"))


def test_radiation_enabled(gen):
    """[network.radiation] with bands builds a Radiation model on the network."""
    assert isinstance(gen.net.radiation, Radiation)


def test_dust_enabled(gen):
    """Presence of [network.dust] builds a Dust model on the network."""
    assert isinstance(gen.net.dust, Dust)


def test_dust_has_photoelectric_emission(gen):
    """The enabled dust model carries its photoelectric-emission sub-model."""
    assert gen.net.dust.pe is not None


def test_photoelectric_chi_weights_bands_by_eavg(gen):
    """Bands fully inside the PE window: each band's chi weight is its eavg (nph)."""
    rad = gen.net.radiation
    chi = sp.expand(gen.net.dust.pe.chi)
    # Coefficients carry the c_hat symbol; it cancels in the ratio.
    coeffs = [chi.coeff(grp.sym) for grp in rad.groups]
    eavgs = [float(grp.eavg) for grp in rad.groups]
    assert float(coeffs[0] / coeffs[1]) == pytest.approx(eavgs[0] / eavgs[1])


def _gen_with_profile_index(tmp_path, profile_index):
    """Run jaffgen with CONFIG_TOML's profile_index replaced."""
    config = tmp_path / "jaffgen.toml"
    config.write_text(
        CONFIG_TOML.replace("profile_index = 0", f"profile_index = {profile_index}")
    )
    return JaffGen(_args(config, tmp_path / "out"))


def test_list_profile_index_sets_per_band_index(tmp_path):
    """A TOML list profile_index gives each radiation band its own index."""
    gen = _gen_with_profile_index(tmp_path, "[0, 1]")
    assert [grp.profile_idx for grp in gen.net.radiation.groups] == [0, 1]


def test_list_profile_index_length_mismatch_raises(tmp_path):
    """profile_index list must have one entry per band (2 bands here)."""
    with pytest.raises(ParserError, match="profile_index"):
        _gen_with_profile_index(tmp_path, "[0, 1, 2]")
