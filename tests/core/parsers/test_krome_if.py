# ABOUTME: KROME conditional rates `if(cond): rate`: the rate while the Fortran condition
# ABOUTME: holds, else 0 (Piecewise), via the krome_if placeholder resolved in post_process

import pytest
import sympy as sp

from jaff.errors import ParserError

HEADER = "@format:idx,R,R,P,P,tmin,tmax,rate"
TRAD, TGAS = sp.Symbol("user_trad"), sp.Symbol("tgas")


def _rate(make_network, rate, tmin="NONE", tmax="NONE"):
    net = make_network([HEADER, f"1,HE+,g,HE++,E,{tmin},{tmax},{rate}"])
    (reaction,) = list(net.reactions)

    return reaction.rate.coefficient


def test_rate_applies_while_the_condition_holds(make_network):
    rate = _rate(make_network, "if(user_Trad.gt.1d3): 5.0d1*user_Trad**(1.63d0)")

    assert float(rate.subs(TRAD, 500.0)) == 0.0
    assert float(rate.subs(TRAD, 2000.0)) == pytest.approx(50.0 * 2000.0**1.63)
    assert "krome_if" not in str(rate)


@pytest.mark.parametrize(
    ("condition", "true_at", "false_at"),
    [
        ("user_Trad.ge.1d3", 1000.0, 999.0),
        ("user_Trad.lt.1d3", 999.0, 1000.0),
        ("user_Trad.le.1d3", 1000.0, 1001.0),
        ("user_Trad>1d3", 1001.0, 1000.0),
        ("USER_TRAD .GT. 1D3", 1001.0, 1000.0),
    ],
)
def test_fortran_and_python_comparisons(make_network, condition, true_at, false_at):
    rate = _rate(make_network, f"if({condition}): 2d0")

    assert float(rate.subs(TRAD, true_at)) == 2.0
    assert float(rate.subs(TRAD, false_at)) == 0.0


def test_and_binds_tighter_than_or(make_network):
    # a .or. b .and. c  ==  a .or. (b .and. c)
    rate = _rate(
        make_network,
        "if(user_Trad.gt.1d4 .or. user_Trad.gt.1d2 .and. Tgas.lt.1d2): 3d0",
    )
    at = {TRAD: 2e4, TGAS: 500.0}  # a true
    assert float(rate.subs(at)) == 3.0
    at = {TRAD: 500.0, TGAS: 50.0}  # b and c true
    assert float(rate.subs(at)) == 3.0
    at = {TRAD: 500.0, TGAS: 500.0}  # only b true
    assert float(rate.subs(at)) == 0.0


def test_condition_with_temperature_limits(make_network):
    rate = _rate(make_network, "if(user_Trad.gt.1d3): 2d0", tmin="1d1", tmax="1d4")

    assert float(rate.subs({TRAD: 2000.0, TGAS: 100.0})) == 2.0
    assert float(rate.subs({TRAD: 500.0, TGAS: 100.0})) == 0.0


def test_a_colon_inside_the_rate_is_not_a_condition(make_network):
    net = make_network(["@format:idx,R,R,P,rate", "1,H+,E,H,1d-12*get_Hnuclei(n(:))"])
    (reaction,) = list(net.reactions)

    assert not reaction.rate.coefficient.has(sp.Piecewise)


@pytest.mark.parametrize("operator", [".eq.", ".ne.", "==", "/="])
def test_unsupported_comparison_raises(make_network, operator):
    with pytest.raises(ParserError, match="unsupported"):
        _rate(make_network, f"if(user_Trad{operator}1d3): 2d0")


def test_not_raises(make_network):
    with pytest.raises(ParserError, match="unsupported"):
        _rate(make_network, "if(.not. user_Trad.gt.1d3): 2d0")
