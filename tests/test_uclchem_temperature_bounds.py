# ABOUTME: UCLCHEM temperature bounds: a T_min of 0 means "no lower bound" (UCLCHEM
# ABOUTME: only switches a rate off below T_min, and gas temperature is never < 0)

import sympy

# H+ + e- -> H with UCLCHEM's placeholder range [0.0, 10000.0] (makerates default).
RECOMBINATION = "H+,E-,NAN,H,NAN,NAN,NAN,2.59e-13,-0.70,0.0,0.0,10000.0,0.0,False"


def test_zero_tmin_is_unbounded(make_network):
    net = make_network([RECOMBINATION])

    assert net.reactions[0].tmin is None


def test_zero_tmin_adds_no_singular_clip_branch(make_network):
    """Clipping at T = 0 evaluated (tgas/300)**-0.7 at zero: a zoo (or nan) branch."""
    rate = make_network([RECOMBINATION]).reactions[0].rate

    assert not rate.has(sympy.zoo, sympy.nan)
