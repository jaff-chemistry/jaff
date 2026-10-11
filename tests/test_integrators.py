# ABOUTME: Tests for the integrator helpers in jaff.common._integrators.
# ABOUTME: smart_integrate dispatches arrays to arr_integrate and SymPy to sym_integrate.
import numpy as np
import pytest
import sympy as sp

from jaff.common import _integrators, arr_integrate, smart_integrate, sym_integrate

E = sp.Symbol("E")


def test_array_input_matches_arr_integrate():
    x = np.linspace(1.0, 3.0, 201)
    y = x**2
    assert smart_integrate(y, x, (1.5, 2.5)) == arr_integrate(y, x, (1.5, 2.5))


def test_expression_input_matches_sym_integrate():
    expr = sp.Piecewise((E**2, (E >= 1) & (E <= 3)), (0, True))
    got = smart_integrate(expr, E, (0.0, sp.oo))
    assert got == sym_integrate(expr, E, (0.0, sp.oo))
    assert got == pytest.approx(26.0 / 3.0, rel=1e-10)


def test_array_with_symbolic_x_raises():
    with pytest.raises(TypeError, match="ndarray"):
        smart_integrate(np.ones(3), E, (0.0, 1.0))


def test_expression_with_array_x_raises():
    with pytest.raises(TypeError, match="Symbol"):
        smart_integrate(E**2, np.linspace(0, 1, 3), (0.0, 1.0))


# --------------------------------------------------------------------------
# get_bounds fast path vs. the original solveset-based algorithm
# --------------------------------------------------------------------------

X = sp.Symbol("X")


def _band_profile():
    return sp.Piecewise((E**-2, E < 8), (E**-1, E < 10), (E**0.5, E < 13.6), (1, True))


def _seven_band_profile():
    edges = [6, 8, 10, 11.2, 13.6, 20, 50, 100]
    pieces = [(E ** (-2 + 0.3 * i), E < e) for i, e in enumerate(edges)]
    return sp.Piecewise(*pieces, (E**2, True))


def _verner_fe():
    from jaff.drivers import JaffDb

    with JaffDb() as j:
        row = j.table("verner_cross_sections").rows(
            conditions="reaction = 'Fe._PHOTON__Fe+.e-'"
        )[0]
    return sp.sympify(row["xsecs"])


def _verner_like():
    fit = 1e-18 * E**-2
    return sp.Piecewise((fit, (E >= 7.9) & (E <= 66)), (0, True))


FAST_PATH_CASES = {
    "verner_like_x_profile": lambda: _verner_like() * _band_profile(),
    "profile_alone": _band_profile,
    "zero_piece_middle": lambda: sp.Piecewise(
        (E, E < 2), (0, E < 5), (E**2, E <= 9), (0, True)
    ),
    "or_condition": lambda: sp.Piecewise(
        (E**2, (E < 1) | ((E > 3) & (E < 4))), (0, True)
    ),
    "reversed_relational": lambda: sp.Piecewise(
        (E, sp.Lt(5, E, evaluate=False)), (0, True)
    ),
    "reversed_le": lambda: sp.Piecewise(
        (E, sp.Le(2, E, evaluate=False) & sp.Ge(7.5, E, evaluate=False)), (0, True)
    ),
    "not_condition": lambda: sp.Piecewise((E, sp.Not(E < 3) & (E < 8)), (0, True)),
    "eq_ne": lambda: sp.Piecewise(
        (0, sp.Eq(E, 3)), (E, sp.Ne(E, 4) & (E < 6)), (1, True)
    ),
    "half_infinite": lambda: sp.Piecewise((E, E >= 1), (0, True)),
    "false_piece": lambda: sp.Piecewise((E, False), (E**2, E < 2), (0, True)),
    "verner_fe": _verner_fe,
    "verner_fe_x_7band": lambda: _verner_fe() * _seven_band_profile(),
}


@pytest.mark.parametrize("name", sorted(FAST_PATH_CASES))
def test_get_bounds_matches_solve_algorithm(name):
    expr = FAST_PATH_CASES[name]()
    assert _integrators.get_bounds(expr, E) == _integrators._get_bounds_solve(expr, E)


@pytest.mark.parametrize("name", sorted(FAST_PATH_CASES))
def test_fast_path_applies(name):
    folded = sp.piecewise_fold(FAST_PATH_CASES[name]())
    assert _integrators._get_bounds_fast(folded, E) is not None


def _random_condition(rng, depth=0):
    kind = rng.integers(0, 9 if depth < 2 else 6)
    t = [1, 2, 2.5, 3, 4][rng.integers(0, 5)]
    if kind <= 3:
        rel = [sp.Lt, sp.Le, sp.Gt, sp.Ge][kind]
        return rel(E, t) if rng.integers(0, 2) else rel(t, E, evaluate=False)
    if kind == 4:
        return sp.Eq(E, t)
    if kind == 5:
        return sp.Ne(E, t)
    args = [_random_condition(rng, depth + 1) for _ in range(2)]

    return [sp.And, sp.Or, lambda a, b: sp.Not(a)][kind - 6](*args)


@pytest.mark.parametrize("seed", range(40))
def test_get_bounds_random_conditions_match_solve_algorithm(seed):
    rng = np.random.default_rng(seed)
    pieces = []
    for i in range(int(rng.integers(1, 4))):
        val = 0 if rng.integers(0, 3) == 0 else E ** (i + 1)
        pieces.append((val, _random_condition(rng)))
    expr = sp.Piecewise(*pieces, (int(rng.integers(0, 2)), True))
    assert _integrators.get_bounds(expr, E) == _integrators._get_bounds_solve(expr, E)


FALLBACK_CASES = {
    "nonlinear": lambda: sp.Piecewise((E, E**2 < 4), (0, True)),
    "other_symbol_value_only": lambda: sp.Piecewise((X * E, E < 3), (0, True)),
    "shifted": lambda: sp.Piecewise((E, E - 3 < 0), (0, True)),
    "not_piecewise": lambda: E**2 + 1,
    "abs_kink": lambda: sp.Abs(E - 2),
}


@pytest.mark.parametrize("name", sorted(FALLBACK_CASES))
def test_get_bounds_fallback_matches_solve_algorithm(name):
    expr = FALLBACK_CASES[name]()
    assert _integrators.get_bounds(expr, E) == _integrators._get_bounds_solve(expr, E)


def test_get_bounds_fallback_rejects_nonlinear_condition():
    folded = sp.piecewise_fold(FALLBACK_CASES["nonlinear"]())
    assert _integrators._get_bounds_fast(folded, E) is None


def test_get_bounds_other_symbol_condition_falls_back():
    expr = sp.Piecewise((E, X < 3), (0, True))
    folded = sp.piecewise_fold(expr)
    assert _integrators._get_bounds_fast(folded, E) is None

    def outcome(fn):
        try:
            return fn(expr, E)
        except Exception as e:  # noqa: BLE001 - compare whatever the old path does
            return type(e)

    assert outcome(_integrators.get_bounds) == outcome(_integrators._get_bounds_solve)


def test_get_bounds_fast_path_skips_solver(monkeypatch):
    # Behavioural rather than wall-clock: the fast path must handle this
    # expression without falling back to SymPy's (slow) inequality solver.
    expr = _verner_fe() * _seven_band_profile()
    expected = _integrators._get_bounds_solve(expr, E)

    def _solver_called(*args, **kwargs):
        raise AssertionError("get_bounds fell back to _get_bounds_solve")

    monkeypatch.setattr(_integrators, "_get_bounds_solve", _solver_called)
    assert _integrators.get_bounds(expr, E) == expected


# --------------------------------------------------------------------------
# sym_integrate preparation cache
# --------------------------------------------------------------------------
def test_sym_integrate_prepares_once_per_expression():
    _integrators._prepare.cache_clear()
    expr = _verner_like() * _band_profile()
    first = sym_integrate(expr, E, (6.0, 20.0))
    second = sym_integrate(expr, E, (10.0, sp.oo))
    info = _integrators._prepare.cache_info()
    assert info.misses == 1
    assert info.hits == 1
    assert first > 0 and second > 0


def test_sym_integrate_cached_result_matches_uncached():
    expr = _verner_fe() * _seven_band_profile()
    bounds = [(6.0, 8.0), (10.0, 13.6), (13.6, 100.0), (0.0, sp.oo)]
    cached = [sym_integrate(expr, E, b) for b in bounds]
    _integrators._prepare.cache_clear()
    fresh = [sym_integrate(expr, E, b) for b in bounds]
    _integrators._prepare.cache_clear()
    assert cached == fresh
