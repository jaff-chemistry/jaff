# ABOUTME: Tests for temperature-cutoff (clip/extrapolate) resolution
# ABOUTME: Covers jaff.toml + jaffgen metadata sources and Model-A precedence

from pathlib import Path

import pytest
from sympy import Piecewise, symbols

from jaff import Network
from jaff.errors import ParserError

KIDA = str(Path(__file__).parent / "fixtures" / "sample_kida.dat")

# Unique, tgas-dependent reaction in sample_kida.dat (bounded below at tmin,
# open above).  Under the default "clip" behaviour its rate is wrapped in a
# Piecewise that holds the boundary value below tmin; "extrapolate" removes the
# wrapper and returns the bare rate.
RXN = "H+.H2__H.H2+"


def _load(config=None, metadata=None):
    """Load sample_kida with an optional jaff.toml path and Network metadata."""
    return Network(
        KIDA,
        funcfile=False,
        _from_cli=True,
        config=config,
        _metadata=metadata or {},
    )


def _is_clipped(net, srxn=RXN):
    """True when the reaction's rate still carries the boundary clamp (clip)."""
    return net.reactions[srxn].rate_coefficient.has(Piecewise)


def _write_toml(tmp_path, body):
    path = tmp_path / "jaff.toml"
    path.write_text(body)
    return str(path)


class TestTemperatureCutoffResolution:
    """T_cutoff resolution from jaff.toml and jaffgen metadata."""

    def test_default_is_clip(self):
        # No config, no metadata -> built-in default clamps tgas.
        assert _is_clipped(_load()) is True
        assert symbols("tgas") in _load().reactions[RXN].rate_coefficient.free_symbols

    def test_jafftoml_global_extrapolate(self, tmp_path):
        cfg = _write_toml(tmp_path, '[network.reactions]\nt_cutoff = "extrapolate"\n')
        assert _is_clipped(_load(config=cfg)) is False

    def test_jafftoml_per_reaction_overrides_own_global(self, tmp_path):
        # Global clip, per-reaction extrapolate -> per-reaction wins.
        cfg = _write_toml(
            tmp_path,
            '[network.reactions]\nt_cutoff = "clip"\n'
            f'[network.reactions.reaction."{RXN}"]\nt_cutoff = "extrapolate"\n',
        )
        assert _is_clipped(_load(config=cfg)) is False

    def test_case_insensitive(self, tmp_path):
        cfg = _write_toml(tmp_path, '[network.reactions]\nt_cutoff = "EXTRAPOLATE"\n')
        assert _is_clipped(_load(config=cfg)) is False

    def test_invalid_value_raises(self, tmp_path):
        cfg = _write_toml(tmp_path, '[network.reactions]\nt_cutoff = "bogus"\n')
        with pytest.raises(ParserError):
            _load(config=cfg)


class TestPrecedenceModelA:
    """jaffgen-local > jaff-local > jaffgen-global > jaff-global > default."""

    def test_jaffgen_local_beats_jaff_local(self, tmp_path):
        # jaff.toml per-reaction clip, jaffgen per-reaction extrapolate -> extrapolate.
        cfg = _write_toml(
            tmp_path, f'[network.reactions.reaction."{RXN}"]\nt_cutoff = "clip"\n'
        )
        meta = {"reaction_props": {RXN: {"T_cutoff": "extrapolate"}}}
        assert _is_clipped(_load(config=cfg, metadata=meta)) is False

    def test_jaff_local_beats_jaffgen_global(self, tmp_path):
        # Model A: a per-reaction setting is never overridden by any global.
        cfg = _write_toml(
            tmp_path, f'[network.reactions.reaction."{RXN}"]\nt_cutoff = "extrapolate"\n'
        )
        meta = {"rate_props": {"T_cutoff": "clip"}}
        assert _is_clipped(_load(config=cfg, metadata=meta)) is False

    def test_jaffgen_global_beats_jaff_global(self, tmp_path):
        cfg = _write_toml(tmp_path, '[network.reactions]\nt_cutoff = "clip"\n')
        meta = {"rate_props": {"T_cutoff": "extrapolate"}}
        assert _is_clipped(_load(config=cfg, metadata=meta)) is False

    def test_reaction_props_tcutoff_without_jaffgen_object(self):
        # T_cutoff-only reaction_props (no shielding) must not require jaffgen_object.
        meta = {"reaction_props": {RXN: {"T_cutoff": "extrapolate"}}}
        assert _is_clipped(_load(metadata=meta)) is False
