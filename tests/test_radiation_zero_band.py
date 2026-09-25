# ABOUTME: Regression tests for zero-cross-section bands in Radiation rate assembly.
# ABOUTME: A band where both cross sections vanish must contribute 0, never NaN.

import math
from unittest.mock import patch

import numpy as np
import sympy as sp

from jaff.physics import RadiationProps
from jaff.physics.photo_reactions._radiation import Radiation


class _Reaction:
    """Reaction with tabulated absorption + decay cross sections (pa enabled)."""

    def __init__(self, energies, absorption, decay):
        E = np.array(energies, dtype=float)
        self.xsecs_dict = {
            "_equations": {"pa": True},
            "photon_energy": E,
            "photo_absorption": np.array(absorption, dtype=float),
            "photodecay": np.array(decay, dtype=float),
        }
        self.dRad = sp.Float(0.0)
        self._metadata = {}
        self.rad_groups = []
        self.rad_xsecs = None
        self.rate = None


def _radiation(bands, profile_index, mode="nph"):
    with patch("jaff.physics.photo_reactions._radiation.BackgroundField"):
        props = RadiationProps(
            bands=list(bands),
            profile_index=profile_index,
            mode=mode,
            c=1.0,
            background_field="draine",
        )
        return Radiation(None, props)


def _numeric(rate, rad):
    """Evaluate a symbolic rate with every density symbol set to 1.0."""
    subs = {rad.den[sp.Idx(i)]: 1.0 for i in range(rad.nbands)}
    return float(rate.xreplace(subs))


def test_zero_band_before_nonzero_no_nan():
    """Band 0 (zero xsec) then band 1 (nonzero): rate finite, band 0 gives 0."""
    rad = _radiation([1.0, 2.0, 3.0], profile_index=2, mode="nph")
    rxn = _Reaction([1.0, 2.0, 3.0], absorption=[0.0, 0.0, 1.0], decay=[0.0, 0.0, 1.0])
    rad.set_reaction_rate_coefficient(rxn)

    assert not rxn.rate.has(sp.nan)
    assert math.isfinite(_numeric(rxn.rate, rad))
    # Band 0's stored coefficient is exactly zero; only band 1 contributes.
    assert rad.groups[0].props[rxn]["k"] == 0
    assert rad.groups[1].props[rxn]["k"] != 0


def test_energy_density_mode_zero_band_no_nan():
    """The /eavg energy-density path must also stay finite through a zero band."""
    rad = _radiation([1.0, 2.0, 3.0], profile_index=2, mode="u")
    rxn = _Reaction([1.0, 2.0, 3.0], absorption=[0.0, 0.0, 1.0], decay=[0.0, 0.0, 1.0])
    rad.set_reaction_rate_coefficient(rxn)

    assert not rxn.rate.has(sp.nan)
    assert math.isfinite(_numeric(rxn.rate, rad))


def test_all_zero_cross_section_table():
    """An all-zero cross-section table yields a zero rate and total, no NaN."""
    rad = _radiation([1.0, 2.0, 3.0], profile_index=2, mode="nph")
    rxn = _Reaction([1.0, 2.0, 3.0], absorption=[0.0, 0.0, 0.0], decay=[0.0, 0.0, 0.0])
    rad.set_reaction_rate_coefficient(rxn)

    assert not rxn.rate.has(sp.nan)
    assert _numeric(rxn.rate, rad) == 0.0
    assert rxn.rad_xsecs == 0.0


def test_xsec_frac_zero_total_is_zero():
    """When the total cross section is zero, every band fraction is defined as 0."""
    rad = _radiation([1.0, 2.0, 3.0], profile_index=2, mode="nph")
    rxn = _Reaction([1.0, 2.0, 3.0], absorption=[0.0, 0.0, 0.0], decay=[0.0, 0.0, 0.0])
    rad.set_reaction_rate_coefficient(rxn)

    for g in rad.groups:
        assert g.props[rxn]["xsec_frac"] == 0.0
