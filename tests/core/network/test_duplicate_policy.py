# ABOUTME: Tests for duplicate_policy resolution and duplicate rate-segment handling
# ABOUTME: duplicate_policy config precedence + end-to-end dedup behaviour (CLI wiring lives elsewhere)


import pytest

from jaff import Network
from jaff.errors import ParserError
from tests._paths import FIXTURES_DIR

FIXTURES = FIXTURES_DIR
KIDA = str(FIXTURES / "sample_kida.dat")
# Two "H+ + e- -> H" lines over the same (10, 10000) range, alpha 1e-13 then
# 9.99e-13.  The surviving coefficient reveals which policy was applied.
DUP = str(FIXTURES / "duplicate_temp_range.dat")
DUP_RXN = "H+.e-__H"


def _load(fname=KIDA, config=None, duplicate_policy=None):
    return Network(
        fname,
        funcfile=False,
        config=config,
        reactions=None
        if duplicate_policy is None
        else {"duplicate_policy": duplicate_policy},
    )


def _write_toml(tmp_path, body):
    path = tmp_path / "jaff.toml"
    path.write_text(body)
    return str(path)


class TestPolicyResolution:
    """duplicate_policy resolution in the network config."""

    def test_default_is_preserve_first(self):
        assert _load().config.reactions.duplicate_policy == "preserve-first"

    @pytest.mark.parametrize("policy", ["preserve-first", "preserve-last", "error"])
    def test_explicit_arg(self, policy):
        assert _load(duplicate_policy=policy).config.reactions.duplicate_policy == policy

    def test_invalid_arg_raises(self):
        with pytest.raises(ParserError, match="duplicate_policy"):
            _load(duplicate_policy="bogus")

    def test_jafftoml_source(self, tmp_path):
        cfg = _write_toml(
            tmp_path, '[network.reactions]\nduplicate_policy = "preserve-last"\n'
        )
        assert _load(config=cfg).config.reactions.duplicate_policy == "preserve-last"

    def test_arg_beats_jafftoml(self, tmp_path):
        cfg = _write_toml(
            tmp_path, '[network.reactions]\nduplicate_policy = "preserve-last"\n'
        )
        assert (
            _load(config=cfg, duplicate_policy="error").config.reactions.duplicate_policy
            == "error"
        )

    def test_invalid_jafftoml_value_raises(self, tmp_path):
        cfg = _write_toml(tmp_path, '[network.reactions]\nduplicate_policy = "bogus"\n')
        with pytest.raises(ParserError, match="duplicate_policy"):
            _load(config=cfg)


class TestDuplicateResolution:
    """End-to-end dedup of a reaction duplicated over the same temperature range."""

    def test_preserve_first_keeps_first_coefficient(self):
        rea = _load(DUP, duplicate_policy="preserve-first").reactions[DUP_RXN]
        assert rea.rate_segments.count == 1
        assert "1.00000000000000e-13" in str(rea.rate_coefficient)

    def test_preserve_last_keeps_last_coefficient(self):
        # Regression: the duplicate carries a *different* rate, so removing the
        # stale segment must target the stored (old) segment, not the new one.
        rea = _load(DUP, duplicate_policy="preserve-last").reactions[DUP_RXN]
        assert rea.rate_segments.count == 1
        assert "9.99000000000000e-13" in str(rea.rate_coefficient)

    def test_error_raises_parsererror(self):
        with pytest.raises(ParserError, match="Duplicate reaction"):
            _load(DUP, duplicate_policy="error")


class TestSameRangeDuplicateInline:
    """Behaviour on an exact same-range duplicate written inline."""

    def test_exact_duplicate_collapses_to_one(self, make_network):
        net = make_network(
            [
                "H + H -> H2 [10,1000] 1e-10",
                "H + H -> H2 [10,1000] 1e-10",  # exact duplicate
            ]
        )
        assert net.reactions.count == 1

    def test_error_policy_rejects_inline_duplicate(self, make_network):
        with pytest.raises(ParserError, match="Duplicate reaction"):
            make_network(
                [
                    "H + H -> H2 [10,1000] 1e-10",
                    "H + H -> H2 [10,1000] 1e-10",
                ],
                reactions={"duplicate_policy": "error"},
            )

    def test_different_temperature_ranges_not_duplicates(self, make_network):
        # Same reaction over disjoint temperature ranges is kept as two segments,
        # not collapsed.
        net = make_network(
            [
                "H + H -> H2 [10,1000] 1e-10",
                "H + H -> H2 [2000,5000] 1e-10",
            ]
        )
        assert net.reactions.count == 1
        assert net.reactions[0].rate_segments.count == 2


class TestPerReactionPolicy:
    """A per-reaction duplicate_policy wins over the global one."""

    def test_per_reaction_policy_applies(self):
        config = {
            "reactions": {
                "duplicate_policy": "error",
                "reaction": {DUP_RXN: {"duplicate_policy": "preserve-last"}},
            }
        }

        rea = _load(DUP, config=config).reactions[DUP_RXN]

        assert "9.99000000000000e-13" in str(rea.rate_coefficient)
