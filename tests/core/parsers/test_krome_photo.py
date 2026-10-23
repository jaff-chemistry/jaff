# ABOUTME: KROME `auto` rates: the parse-time krome_photo placeholder and the
# ABOUTME: post_process hook numbering them network-wide into photorates()

import pytest

from jaff.core.network import NetworkSymbols

NETWORK = "@format:idx,R,R,P,rate\n1,H+,E,H,1e-12\n2,C,C+,E,auto\n"


def test_auto_becomes_numbered_photorates(make_network):
    net = make_network(NETWORK, name="react_test")

    rxn = next(r for r in net.reactions if r.type == "photo")
    assert rxn.rate_coefficient.func == NetworkSymbols.photorates
    assert [float(a) for a in rxn.rate_coefficient.args] == pytest.approx([1, 1e99, 1e99])
    assert all("krome_photo" not in str(s.rate) for s in rxn.rate_segments)
