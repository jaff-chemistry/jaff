# ABOUTME: Integrates the generated chemistry with scipy's Radau and the analytic Jacobian
# ABOUTME: solve(): densities by species name, parameters checked by name, optional T
"""Stiff integration of the generated chemistry (``chemistry.py``)."""

from dataclasses import dataclass

import numpy as np
from chemistry import NSPEC, PARAMETERS, SPECIES, jacobian, rhs
from scipy.integrate import solve_ivp

DEFAULT_RTOL = 1e-12
"""Relative tolerance of the integration."""

DEFAULT_ATOL = 1e-20
"""Absolute tolerance of the integration [cm^-3, and K for the temperature]."""


@dataclass
class Result:
    """Outcome of :func:`solve`.

    Attributes
    ----------
    t : np.ndarray
        Output times [s].
    n : dict[str, np.ndarray]
        Number density of each species at ``t`` [cm^-3], by name.
    temperature : np.ndarray
        Gas temperature at ``t`` [K] (constant unless evolved).
    solution : object
        The ``scipy.integrate.solve_ivp`` result.
    """

    t: np.ndarray
    n: dict
    temperature: np.ndarray
    solution: object


def _check_params(params):
    """Raise ``ValueError`` naming every missing and unknown parameter."""
    missing = sorted(set(PARAMETERS) - set(params))
    unknown = sorted(set(params) - set(PARAMETERS))
    if missing or unknown:
        raise ValueError(
            f"parameters: missing {missing}, unknown {unknown}; "
            f"expected exactly {list(PARAMETERS)}"
        )


def _rates_or_nan(size, *args):
    """``rhs(*args)``, or NaNs when a math function fails at a trial state.

    Radau probes states where a density or the temperature may be slightly
    negative; ``math.log`` and friends raise there. A non-finite rhs makes Radau
    reject the trial and shrink the step, so the failure stays recoverable.
    """
    try:
        return rhs(*args)[:size]
    except (ValueError, OverflowError, ZeroDivisionError):
        return np.full(size, np.nan)


def _initial_densities(n0):
    """Density vector from ``{species name: density}``; unlisted species start at 0."""
    unknown = sorted(set(n0) - set(SPECIES))
    if unknown:
        raise ValueError(f"unknown species {unknown}; known: {list(SPECIES)}")
    nden = np.zeros(NSPEC)
    for name, value in n0.items():
        nden[SPECIES[name]] = value

    return nden


def solve(
    n0,
    temperature,
    t_end,
    params,
    evolve_temperature=False,
    t_eval=None,
    rtol=DEFAULT_RTOL,
    atol=DEFAULT_ATOL,
):
    """Integrate the chemistry from t = 0 to *t_end* with Radau.

    Parameters
    ----------
    n0 : dict[str, float]
        Initial number densities [cm^-3] by species name; others start at 0.
    temperature : float
        Initial gas temperature [K]; held fixed unless *evolve_temperature*.
    t_end : float
        End time [s].
    params : dict[str, float]
        A value for every name in ``chemistry.PARAMETERS``, and nothing else.
    evolve_temperature : bool, optional
        Evolve the gas temperature with dT/dt from the EOS. Default ``False``.
    t_eval : array-like, optional
        Output times [s]; by default the solver's own steps.
    rtol, atol : float, optional
        Tolerances (default ``DEFAULT_RTOL``, ``DEFAULT_ATOL``).

    Returns
    -------
    Result
        Times, densities by name, temperature and the raw scipy result.

    Raises
    ------
    ValueError
        If a parameter is missing or unknown, or a species name is unknown.
    RuntimeError
        If the integration fails.
    """
    _check_params(params)
    nden0 = _initial_densities(n0)

    if evolve_temperature:
        y0 = np.append(nden0, temperature)

        def f(t, y):
            return _rates_or_nan(NSPEC + 1, y[:NSPEC], y[NSPEC], params)

        def jac(t, y):
            return jacobian(y[:NSPEC], y[NSPEC], params)

    else:
        y0 = nden0

        def f(t, y):
            return _rates_or_nan(NSPEC, y, temperature, params)

        def jac(t, y):
            return jacobian(y, temperature, params)[:NSPEC, :NSPEC]

    sol = solve_ivp(
        f, (0.0, t_end), y0, method="Radau", jac=jac, t_eval=t_eval, rtol=rtol, atol=atol
    )
    if not sol.success:
        raise RuntimeError(f"integration failed: {sol.message}")

    if evolve_temperature:
        tgas = sol.y[NSPEC]
    else:
        tgas = np.full(sol.t.shape, float(temperature))
    n = {name: sol.y[i] for name, i in SPECIES.items()}

    return Result(t=sol.t, n=n, temperature=tgas, solution=sol)
