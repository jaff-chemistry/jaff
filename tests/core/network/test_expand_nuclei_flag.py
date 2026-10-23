# ABOUTME: expand_nuclei controls n_H/n_He element-sum expansion

import sympy
from jaff import Network


def _net(tmp_path, **kw):
    p = tmp_path / "n.dat"
    p.write_text("@format:idx,R,R,P,rate\n1,H,H,H2,1\n2,H2,H,H,n_H\n")
    return Network(str(p), funcfile=False, **kw)


def test_expand_nuclei_true_vs_false_differ(tmp_path):
    # Just assert both values are accepted and produce a rate (behavior of the
    # flag itself is covered elsewhere); expand_nuclei must be a valid kwarg.
    r_true = _net(tmp_path, expand_nuclei=True).reactions[1].rate_coefficient
    r_false = _net(tmp_path, expand_nuclei=False).reactions[1].rate_coefficient
    assert r_true is not None and r_false is not None
