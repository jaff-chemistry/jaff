"""Tidy-DataFrame builders for the seaborn objects plotting layer.

The seaborn objects interface (``seaborn.objects``) consumes long/tidy
DataFrames.  These helpers assemble them and convert to the caller's requested
units.

.. important::
   Every frame returned here has a unique ``RangeIndex`` (``reset_index`` /
   ``ignore_index``).  The objects interface fails on duplicate index labels
   (``ValueError: Must have equal len keys and value``) when unscaling
   coordinates, which happens if per-process frames are concatenated without
   resetting the index.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from sympy import Symbol, lambdify

from ..common._integrators import get_bounds
from . import _units


def long_frame(
    series: list[tuple[np.ndarray, np.ndarray, str]],
    x_col: str,
    y_col: str,
    group_col: str,
) -> pd.DataFrame:
    """Assemble a tidy long frame from labelled ``(x, y, label)`` triples.

    Parameters
    ----------
    series : list[tuple[numpy.ndarray, numpy.ndarray, str]]
        One ``(x, y, label)`` per curve.  The ``x``/``y`` arrays may differ in
        length and range across curves.
    x_col, y_col, group_col : str
        Column names for the x, y, and grouping/label columns.

    Returns
    -------
    pandas.DataFrame
        Long frame with a unique index (duplicate indices break seaborn
        objects).
    """
    frames = [
        pd.DataFrame(
            {
                x_col: np.asarray(x, dtype=float),
                y_col: np.asarray(y, dtype=float),
                group_col: label,
            }
        )
        for x, y, label in series
    ]
    if not frames:
        return pd.DataFrame({x_col: [], y_col: [], group_col: []})
    # ignore_index=True is required: duplicate indices break seaborn objects.
    return pd.concat(frames, ignore_index=True)


def line_frame(
    x: np.ndarray, y: np.ndarray, x_col: str = "x", y_col: str = "y"
) -> pd.DataFrame:
    """Build a two-column frame for a single line.

    Parameters
    ----------
    x, y : numpy.ndarray
        Coordinates.
    x_col, y_col : str, optional
        Column names.

    Returns
    -------
    pandas.DataFrame
        Frame ``{x_col, y_col}`` with a fresh ``RangeIndex``.
    """
    return pd.DataFrame({x_col: np.asarray(x), y_col: np.asarray(y)})


def band_frame(
    band_xsecs: pd.DataFrame,
    energy_unit: str,
    xsec_unit: str,
) -> pd.DataFrame:
    """Convert a reaction's band-averaged cross sections to plotting units.

    Parameters
    ----------
    band_xsecs : pandas.DataFrame
        As returned by :attr:`jaff.core.reaction.Reaction.band_xsecs`
        (``lower``/``upper``/``eavg`` in eV, ``xsec`` in cm²).
    energy_unit : str
        Target unit for the band-edge / mid-point columns.
    xsec_unit : str
        Target unit for the cross-section column.

    Returns
    -------
    pandas.DataFrame
        A copy with converted units and added ``mid`` (geometric mid-point,
        for log axes) and ``width`` columns, restricted to rows with a finite
        cross section.  Empty if no band has a finite cross section.

    Notes
    -----
    Bands with a non-finite edge (an open top band with ``upper = inf``) keep
    that edge; the plotter is responsible for clipping the drawn width to the
    axis range.
    """
    df = band_xsecs.copy()
    # Drop bands without a tabulated cross section (custom-rate reactions).
    df = df[np.isfinite(df["xsec"])].reset_index(drop=True)
    if df.empty:
        return df

    df["lower"] = _units.convert_energy(df["lower"].to_numpy(), "eV", energy_unit)
    df["upper"] = _units.convert_energy(df["upper"].to_numpy(), "eV", energy_unit)
    df["eavg"] = _units.convert_energy(df["eavg"].to_numpy(), "eV", energy_unit)
    df["xsec"] = _units.convert_xsec(df["xsec"].to_numpy(), "cm2", xsec_unit)

    # Geometric mid-point works on both linear and log energy axes; falls back
    # to the arithmetic mean if an edge is non-positive.
    lo, hi = df["lower"].to_numpy(), df["upper"].to_numpy()
    with np.errstate(invalid="ignore"):
        geo = np.sqrt(lo * hi)
    df["mid"] = np.where((lo > 0) & np.isfinite(hi), geo, (lo + hi) / 2.0)
    df["width"] = hi - lo
    return df


def _sampled_xsecs(xsecs: dict, n: int = 2000) -> dict:
    """Tabulate a symbolic ``photodecay_expr`` (Verner) for plotting.

    Returns *xsecs* unchanged when it already carries arrays; otherwise a copy
    with ``photon_energy`` log-spaced over the expression's breakpoint range
    and ``photodecay`` evaluated on it.
    """
    expr = xsecs.get("photodecay_expr")
    if expr is None or xsecs.get("photon_energy") is not None:
        return xsecs
    E = Symbol("E")
    pts = get_bounds(expr, E)
    grid = np.geomspace(pts[0], pts[-1], n)
    values = np.asarray(lambdify(E, expr, "numpy")(grid), dtype=float)
    return {**xsecs, "photon_energy": grid, "photodecay": values}
