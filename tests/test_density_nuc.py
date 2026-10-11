# ABOUTME: n_X_nuc = element nucleus sum; guards raise. Magic n_H untouched (Stage 1).

import sympy
import pytest
from jaff import Network
from jaff.errors import ParserError


def _net(tmp_path, body, **kw):
    p = tmp_path / "n.dat"
    p.write_text(body)
    return Network(str(p), funcfile=False, **kw)


def test_n_H_nuc_is_element_sum(tmp_path):
    net = _net(tmp_path, "@format:idx,R,R,P,rate\n1,H,H,H2,1\n2,H2,H,H,n_H_nuc\n")
    nden = net.symbols.ndens
    expected = nden[sympy.Idx(net.species["H"].index)] + 2 * nden[
        sympy.Idx(net.species["H2"].index)
    ]
    assert sympy.simplify(net.reactions[1].rate - expected) == 0


def test_n_H_nuc_symbol_when_not_expanded(tmp_path):
    net = _net(
        tmp_path, "@format:idx,R,R,P,rate\n1,H,H,H2,1\n2,H2,H,H,n_H_nuc\n",
        expand_nuclei=False,
    )
    assert net.reactions[1].rate == sympy.symbols("nh_nuc")


def test_charged_nucleus_sum_raises(tmp_path):
    with pytest.raises(ParserError):
        _net(tmp_path, "@format:idx,R,R,P,rate\n1,C,C+,C,n_Cj_nuc\n")


def test_unknown_element_nucleus_sum_raises(tmp_path):
    with pytest.raises(ParserError):
        _net(tmp_path, "@format:idx,R,R,P,rate\n1,H,C,CH,n_Xx_nuc\n")
