# ABOUTME: KROME BASE_GLOBALS emit the n_X grammar, not the nh/nh0 layer

from jaff.core.parsers.network._formats.krome.parser import KromeParser


def test_base_globals_use_n_x_grammar():
    bg = KromeParser.BASE_GLOBALS
    assert bg["get_hnuclei(n)"] == "n_H_nuc"
    assert bg["n(idx_h)"] == "n_H"
    assert bg["n(idx_h2)"] == "n_H2"
    assert bg["n_global(idx_h2)"] == "n_H2"
    # temperature/user shorthands unchanged
    assert bg["t32"] == "tgas/3e2"
    assert bg["user_av"] == "av"
