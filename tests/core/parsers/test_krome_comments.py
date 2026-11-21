# ABOUTME: KROME end-of-line comments on reaction lines: stripped from the rate (as KROME
# ABOUTME: does), and a comment mentioning "grain(s)" marks the reaction as type "grain"

import pytest

from jaff.core.parsers.network._formats.krome.reaction import KromeReaction

STATE = {"idx": True, "nreact": 2, "nprod": 2, "tmin": False, "tmax": False, "rate": True}


def _parse(line: str) -> dict:
    return KromeReaction().parse(line, 1, STATE, "net.jet")


@pytest.mark.parametrize("marker", ["#", "!"])
def test_comment_is_stripped_from_the_rate(marker):
    fields = _parse(f"1,H2,H+,H2+,H,1.d-17*Tgas {marker} H2 + H -> 3H")

    assert fields["rate"] == "1.e-17*Tgas"


@pytest.mark.parametrize(
    "comment",
    ["# grain", "! grain", "#GRAIN", "# H2 formation on grains", "! on Grains, approx"],
)
def test_grain_comment_marks_a_grain_reaction(comment):
    fields = _parse(f"1,H,H,H2,,3d-18*sqrt(Tgas) {comment}")

    assert fields["type"] == "grain"
    assert fields["rate"] == "3e-18*sqrt(Tgas)"


@pytest.mark.parametrize("comment", ["# migraine", "# graining", "# no surface"])
def test_other_comments_keep_the_inferred_type(comment):
    assert _parse(f"1,H,H,H2,,3d-18 {comment}")["type"] == "unknown"


def test_grain_comment_overrides_the_inferred_type():
    # exp(-av) would make it photo.
    assert _parse("1,H,H,H2,,1d-10*exp(-user_Av) # grain")["type"] == "grain"


def test_without_a_comment_nothing_changes():
    fields = _parse("1,H,H,H2,,3d-18")

    assert (fields["rate"], fields["type"]) == ("3e-18", "unknown")


def test_the_source_line_keeps_its_comment():
    assert _parse("1,H,H,H2,,3d-18 # grain")["line"] == "1,H,H,H2,,3d-18 # grain"


def test_grain_reaction_in_a_network(make_network):
    net = make_network(["@format:idx,R,R,P,rate", "1,H,H,H2,3d-18*sqrt(Tgas) # grain"])
    (reaction,) = list(net.reactions)

    assert reaction.type == "grain"
    assert str(reaction.rate.coefficient) == "3.0e-18*sqrt(tgas)"
