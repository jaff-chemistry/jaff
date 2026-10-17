# ABOUTME: Tests for get_table adaptive refinement termination (nT_max cap)
# ABOUTME: A true step discontinuity never meets err_tol and must stop at nT_max

import logging

from jaff import Network
from jaff.io._io import get_table


def _step_network(tmp_path):
    """One reaction whose rate jumps 1 -> 1000 at exactly 100 K (no gap)."""
    dat = tmp_path / "step.dat"
    dat.write_text(
        "@format:idx,R,R,P,tmin,tmax,rate\n1,H,H,H2,10,100,1\n2,H,H,H2,100,1000,1000\n"
    )
    return Network(str(dat), funcfile=False)


def test_refinement_stops_at_nT_max(tmp_path, caplog):
    net = _step_network(tmp_path)
    logger = logging.getLogger("test_get_table")
    with caplog.at_level(logging.WARNING, logger="test_get_table"):
        temp, coef = get_table(
            net.reactions, logger, T_min=10.0, T_max=1000.0, nT=9, err_tol=0.01, nT_max=65
        )
    assert len(temp) == 65
    assert coef.shape == (1, 65)
    assert "nT_max = 65" in caplog.text


def test_refinement_converges_without_warning(tmp_path, caplog):
    net = _step_network(tmp_path)
    logger = logging.getLogger("test_get_table")
    with caplog.at_level(logging.WARNING, logger="test_get_table"):
        temp, _ = get_table(
            net.reactions, logger, T_min=10.0, T_max=90.0, nT=9, err_tol=0.01, nT_max=65
        )
    assert len(temp) <= 65
    assert "nT_max" not in caplog.text
