# ABOUTME: Characterization tests pinning current parser shorthand-replacement behavior
# ABOUTME: Regression guard now that BASE_GLOBALS is owned by the KROME parser

import sympy

TGAS = sympy.symbols("tgas")


def _rate(make_network, rate_str, reactants="H,C", product="CH", **kw):
    net = make_network(
        [
            "@format:idx,R,R,P,tmin,tmax,rate",
            f"1,{reactants},{product},10,1000,{rate_str}",
        ],
        **kw,
    )
    return net, net.reactions[0].rate


def _at(rate, t=100):
    """Evaluate a (possibly clip-Piecewise-wrapped) rate inside its valid range."""
    return sympy.nsimplify(rate.subs(TGAS, t)) if hasattr(rate, "subs") else rate


# --------------------------------------------------------------------------- #
# temperature shorthands (KROME), incl. order-dependent compounds              #
# Rates are evaluated at T=100 (inside [10,1000]) to bypass the clip Piecewise. #
# --------------------------------------------------------------------------- #
def test_t32_expands_to_tgas_over_300(make_network):
    _, rate = _rate(make_network, "t32")
    assert _at(rate) == sympy.Rational(100, 300)


def test_te_expands(make_network):
    _, rate = _rate(make_network, "te")
    assert abs(float(rate.subs(TGAS, 100)) - 100 * 8.617343e-5) < 1e-12


def test_invt32_resolves_through_t32(make_network):
    # invt32 = 1/t32 = 1/(tgas/300) = 300/tgas.  Guards resolution ORDER:
    # if t32 is not substituted first, invt32 keeps a stray t32 symbol.
    _, rate = _rate(make_network, "invt32")
    assert sympy.Symbol("t32") not in rate.free_symbols
    assert _at(rate) == 3  # 300/100


def test_invte_resolves_through_te(make_network):
    _, rate = _rate(make_network, "invte")
    assert sympy.Symbol("te") not in rate.free_symbols
    assert abs(float(rate.subs(TGAS, 100)) - 1 / (100 * 8.617343e-5)) < 1e-6


def test_invtgas_expands(make_network):
    _, rate = _rate(make_network, "invtgas")
    assert _at(rate) == sympy.Rational(1, 100)


def test_sqrtgas_expands(make_network):
    _, rate = _rate(make_network, "sqrtgas")
    assert abs(float(rate.subs(TGAS, 100)) - 10.0) < 1e-9


# --------------------------------------------------------------------------- #
# user_ shorthands -> bare symbols                                             #
# --------------------------------------------------------------------------- #
def test_user_tdust_maps_to_tdust(make_network):
    _, rate = _rate(make_network, "user_tdust")
    assert rate == sympy.Symbol("tdust")


def test_user_av_maps_to_av(make_network):
    _, rate = _rate(make_network, "user_av")
    assert rate == sympy.Symbol("av")


# --------------------------------------------------------------------------- #
# KROME density accessors -> JAFF density namespace                            #
# --------------------------------------------------------------------------- #
def test_n_idx_h_is_atomic_H_density(make_network):
    # n(idx_H) -> n_H -> species H number density.
    net, rate = _rate(make_network, "n(idx_H)")
    assert rate == net.symbols.ndens[net.species["H"].index]


def test_n_idx_h2_is_H2_density(make_network):
    net, rate = _rate(make_network, "n(idx_H2)", reactants="H2,C", product="CH2")
    assert rate == net.symbols.ndens[net.species["H2"].index]


def test_n_global_idx_h2_is_H2_density(make_network):
    net, rate = _rate(
        make_network, "n_global(idx_H2)", reactants="H2,C", product="CH2"
    )
    assert rate == net.symbols.ndens[net.species["H2"].index]


def test_get_hnuclei_is_H_nucleus_sum(make_network):
    # get_hnuclei(n) -> n_H_nuc -> total H-nuclei sum (H once, H2 twice).
    net, rate = _rate(
        make_network, "get_hnuclei(n)", reactants="H,H", product="H2"
    )
    nden = net.symbols.ndens
    expected = nden[sympy.Idx(net.species["H"].index)] + 2 * nden[
        sympy.Idx(net.species["H2"].index)
    ]
    assert sympy.simplify(rate - expected) == 0
