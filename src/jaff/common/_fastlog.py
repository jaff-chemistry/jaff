#!/usr/bin/env python3
# ABOUTME: Fast log2 approximation (from quokka) and its numerical inverse, used to build
# ABOUTME: uniform-in-fast_log2 temperature grids for rate tables.

"""
Fast log2 approximation and its numerical inverse.

The algorithm is adapted from the quokka radiation-hydrodynamics code
(https://github.com/quokka-astro/quokka/) and is designed to produce a
*uniform-in-fast_log2* temperature grid for rate-coefficient tables so that
the runtime table-lookup in quokka can use the same cheap approximation
instead of a transcendental ``log2`` call.

The approximation
-----------------
For any positive float *x* decomposed by :func:`numpy.frexp` as::

    x = mantissa * 2^exponent,   mantissa ∈ [0.5, 1.0)

the fast log2 is defined as::

    fast_log2(x) = 2 * (mantissa - 1) + exponent

This is a piecewise-linear approximation to the true ``log2`` that is exact
at every power of 2 and has a maximum absolute error of about 0.086 bits.
"""

import numpy as np
from scipy.optimize import root_scalar


def fast_log2(x):
    """
    Compute a fast piecewise-linear approximation to log2(*x*).

    Uses the :func:`numpy.frexp` mantissa/exponent decomposition to evaluate::

        fast_log2(x) = 2 * (mantissa - 1) + exponent

    This expression is exact at every power of two and continuous everywhere,
    but deviates from the true log2 by at most ~0.086 in the worst case
    (at the midpoints between consecutive powers of two).

    Parameters
    ----------
    x : array_like
        Input value(s).  Must be positive; non-positive values map to ``NaN``.

    Returns
    -------
    float or numpy.ndarray
        Approximated log2 value(s).  Returns a scalar when *x* is scalar,
        an array otherwise.

    Notes
    -----
    The approximation is designed so that a uniform grid in ``fast_log2``
    space can be reconstructed at runtime using the same cheap operation,
    avoiding a transcendental ``log2`` call in inner loops.
    """
    x_ = np.asarray(x)
    scalar_input = x_.ndim == 0
    x_ = np.atleast_1d(x_)

    # frexp gives mantissa ∈ [0.5, 1.0): x = mantissa * 2^exponent
    mantissa, exponent = np.frexp(x_)

    # approximation: log2(x) ≈ 2(mantissa - 1) + exponent
    res = 2 * (mantissa - 1) + exponent

    res[x_ <= 0.0] = np.nan

    if scalar_input:
        return res[0]
    else:
        return res


def inverse_fast_log2(y):
    """
    Compute the numerical inverse of :func:`fast_log2` to near machine precision.

    For each element of *y* the function finds *x* such that
    ``fast_log2(x) == y`` using Brent's root-finding method, starting from
    the initial bracket ``[0.5 * 2^y, 2.0 * 2^y]`` (which safely brackets
    the root because :func:`fast_log2` is strictly monotone).  If Brent's
    method leaves a residual above the target tolerance, up to five Newton
    polishing steps are applied using the derivative::

        d(fast_log2)/dx ≈ 1 / (x * ln(2))

    Parameters
    ----------
    y : array_like
        Target fast_log2 value(s) to invert.

    Returns
    -------
    float or numpy.ndarray
        The *x* value(s) satisfying ``fast_log2(x) ≈ y``.  Returns a
        scalar when *y* is scalar, an array otherwise.

    Notes
    -----
    The tolerance goal is ``20 * eps * |y|`` where ``eps`` is double
    machine epsilon (~2.2e-16).  A warning is printed to stdout if the
    residual cannot be reduced below this threshold after all polishing
    steps -- this should not occur in normal use.
    """
    y_arr = np.atleast_1d(y)
    res = np.zeros(y_arr.shape)

    eps = np.finfo(np.float64).eps

    for i, y_ in enumerate(y_arr):
        # initial guess: inverse of exact log2
        x_guess = 2.0**y_

        # [0.5, 2]*x_guess brackets the root safely because fast_log2 is monotonic
        x_bracket = np.array([0.5, 2]) * x_guess

        root = root_scalar(
            lambda x: fast_log2(x) - y_,
            method="brentq",
            bracket=x_bracket,
            xtol=eps * x_guess,
            rtol=4 * eps,
            maxiter=1000,
        )

        resid = abs(fast_log2(root.root) - y_)
        tol = eps * abs(y_) * 20
        if resid > tol:
            # Polish with up to five Newton steps when Brent did not converge
            for _ in range(5):
                dx = -(fast_log2(root.root) - y_) * root.root * np.log(2)
                root_new = root.root + dx
                resid_new = abs(fast_log2(root_new) - y_)
                if resid_new < resid:
                    root.root = root_new
                    resid = resid_new
                else:
                    break

            if resid > tol:
                print(
                    "Warning: could not reduce residual below {:e} for "
                    "inverse_fast_log2({:e})".format(resid, y_)
                )

        res[i] = root.root

    if np.asarray(y).ndim == 0:
        return res[0]
    else:
        return res
