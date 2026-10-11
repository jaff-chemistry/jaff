# ABOUTME: n_X species-path: n_H is a species; 0-marker and missing species raise

import sympy
import pytest
from jaff import Network
from jaff.errors import ParserError


def _net(tmp_path, body):
    p = tmp_path / "n.dat"
    p.write_text(body)
    return Network(str(p), funcfile=False)


def test_n_H_is_neutral_species_not_sum(tmp_path):
    net = _net(tmp_path, "@format:idx,R,R,P,rate\n1,H,H,H2,1\n2,H2,H,H,n_H\n")
    assert net.reactions[1].rate == net.symbols.ndens[net.species["H"].index]


def test_zero_marker_no_longer_special_raises(tmp_path):
    with pytest.raises(ParserError):
        _net(tmp_path, "@format:idx,R,R,P,rate\n1,H,C,CH,n_H0\n")


def test_missing_species_raises(tmp_path):
    with pytest.raises(ParserError):
        _net(tmp_path, "@format:idx,R,R,P,rate\n1,H,C,CH,n_Ne\n")
