# ABOUTME: Format detection: UDFA claims only "<index>:" lines, PRIZMO only "A -> B [..]" ones,
# ABOUTME: so KROME lines with ':' / '->' stay KROME; look-alike lines still raise errors

import pytest
import sympy as sp

from jaff.core.parsers.network._formats.prizmo.reaction import PrizmoReaction
from jaff.core.parsers.network._formats.udfa.reaction import UdfaReaction
from jaff.errors import ParserError

KROME_WITH_COLON = [
    "1,H+,E,H,1d-12*get_Hnuclei(n(:))",
    "1,H+,E,H,1d-12*sum(n(1:nmols))",
    "1,H+,E,H,if(Tgas>1d2): 1d-12",
    "24,HE+,g,,HE++,E,,,NONE,NONE,if(user_Trad.gt.1d3): 5.0d1*user_Trad**(1.63d0)",
]
KROME_WITH_ARROW = [
    "3,H2,H,H,H,1.d-17*Tgas !H2 + H -> 3H",
    "3,H2,H,H,H,1.d-17*Tgas # H2 + H -> 3 H [10, 100]",
]
UDFA = [
    "1:RR:H:e-:H-::::1:1e-16:0:0:10:10000",
    '1:AD:C-:C:C2:e-:::1:5.00e-10:0.00:0.0:10:41000:L:C:"10.1086/190665":"Prasad"',
    "  12 :PH:H2::H:H:::1:5.00e-11:0.00:0.0:10:3000",
]
PRIZMO = [
    "H + H -> H2 [10,1000] 1e-10",
    "H + H -> H2 [,] 1e-10",
    "H -> H+ + E        []           PHOTO, 13.60",
    "  CO -> C + O [] 1.0e-10 * exp(-3.53 * av)",
]


@pytest.mark.parametrize("line", KROME_WITH_COLON)
def test_udfa_does_not_claim_a_krome_line_with_a_colon(line):
    assert not UdfaReaction.global_re.match(line)


@pytest.mark.parametrize("line", UDFA)
def test_udfa_claims_udfa_lines(line):
    assert UdfaReaction.global_re.match(line)


@pytest.mark.parametrize("line", KROME_WITH_ARROW)
def test_prizmo_does_not_claim_a_krome_line_with_an_arrow(line):
    assert not PrizmoReaction.global_re.match(line)


@pytest.mark.parametrize("line", PRIZMO)
def test_prizmo_claims_prizmo_lines(line):
    assert PrizmoReaction.global_re.match(line)


@pytest.mark.parametrize("line", ["# H + H -> H2 [10,1000] 1e-10", "! 1:RR:H:e-:H-::::1"])
def test_comment_lines_are_claimed_by_neither(line):
    assert not PrizmoReaction.global_re.match(line)
    assert not UdfaReaction.global_re.match(line)


def test_a_krome_line_with_a_colon_is_parsed_as_krome(make_network):
    net = make_network(["@format:idx,R,R,P,rate", "1,H+,E,H,1d-12*get_Hnuclei(n(:))"])
    (reaction,) = list(net.reactions)

    assert reaction.source.format == "krome"
    assert not reaction.rate.coefficient.has(sp.Function("get_hnuclei"))


@pytest.mark.parametrize(
    ("line", "message"),
    [
        ("H + H + H -> H2 + H 1e-10", "Invalid PRIZMO reaction"),  # no [tmin, tmax]
        ("RR:H:e-:H-::::1:1e-16:0:0:10:10000", "Invalid UDFA reaction"),  # no index
    ],
)
def test_a_line_that_only_looks_like_a_format_is_still_an_error(
    make_network, line, message
):
    # Not skipped: the last-resort handlers report it in that format's terms.
    with pytest.raises(ParserError, match=message):
        make_network([line])
