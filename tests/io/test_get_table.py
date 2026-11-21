# ABOUTME: io.get_table temperature bounds: inferred from the reactions' tmin/tmax, a clear
# ABOUTME: error when none can be inferred (no NaN grid), and T_min < T_max enforced

import numpy as np
import pytest

from jaff.io._io import get_table

RANGED = ["H + H -> H2 [10,1000] 1e-10", "H2 -> H + H [20,2000] 1e-15"]
UNRANGED = ["H + H -> H2  []  1e-10", "H2 -> H + H  []  1e-15"]


def _table(reactions, **kwargs):
    return get_table(list(reactions), None, nT=8, err_tol=None, **kwargs)


def test_bounds_are_inferred_from_the_reactions(make_network):
    temp, _ = _table(make_network(RANGED).reactions)

    assert temp[0] == pytest.approx(10.0)
    assert temp[-1] == pytest.approx(2000.0)


def test_without_any_temperature_limit_it_asks_for_bounds(make_network):
    with pytest.raises(ValueError, match="set T_min and T_max"):
        _table(make_network(UNRANGED).reactions)


def test_given_bounds_make_a_finite_table(make_network):
    temp, coeff = _table(make_network(UNRANGED).reactions, T_min=10.0, T_max=100.0)

    assert np.isfinite(temp).all()
    assert np.isfinite(coeff).all()


def test_an_empty_reaction_list_asks_for_bounds():
    with pytest.raises(ValueError, match="set T_min and T_max"):
        _table([])


@pytest.mark.parametrize(("t_min", "t_max"), [(100.0, 100.0), (1000.0, 10.0)])
def test_t_min_must_be_below_t_max(make_network, t_min, t_max):
    with pytest.raises(ValueError, match="T_min"):
        _table(make_network(RANGED).reactions, T_min=t_min, T_max=t_max)
