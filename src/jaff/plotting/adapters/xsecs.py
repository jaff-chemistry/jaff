# ABOUTME: Photo cross-section adapter: reaction -> series dicts per process
# ABOUTME: Unit conversion, Verner sampling, optional band-average bars

"""Photo cross sections as plotting series.

Duck-typed: a reaction is anything with ``photo.xsec`` (and optionally
``photo.band_xsecs`` / ``get_latex``).
"""

from __future__ import annotations

from dataclasses import replace
from typing import TYPE_CHECKING, Any

import numpy as np
from sympy import lambdify

from ...common._integrators import get_bounds
from ...io import JaffLogger
from .. import _units
from ._common import label_of as _label_of

if TYPE_CHECKING:
    from ...physics import XsecProps

#: Valid photo cross-section process keys.
XSEC_PROCESSES: tuple[str, ...] = ("photo_absorption", "photodecay")


#: Energy units that are really wavelengths (axis direction is reversed).
_WAVELENGTH_UNITS: frozenset[str] = frozenset({"nm", "um"})


def normalize_processes(processes: str | list[str] | None) -> list[str]:
    """Normalise a process selection to a validated list of keys.

    Raises
    ------
    KeyError
        If any key is not in :data:`XSEC_PROCESSES`.
    """
    if processes is None or processes == "all":
        return list(XSEC_PROCESSES)
    procs = [processes] if isinstance(processes, str) else list(processes)
    invalid = [p for p in procs if p not in XSEC_PROCESSES]
    if invalid:
        raise KeyError(
            f"Invalid cross-section(s) {invalid}. Supported: {', '.join(XSEC_PROCESSES)}"
        )
    return procs


def _sampled_xsecs(xsecs: XsecProps, n: int = 2000) -> XsecProps:
    """Tabulate a symbolic Verner ``photodecay_expr`` on a log energy grid."""
    expr = xsecs.photodecay_expr
    if expr is None or xsecs.photon_energy is not None:
        return xsecs
    from ...core.network import NetworkSymbols

    E = NetworkSymbols.ph_E
    pts = get_bounds(expr, E)
    if len(pts) < 2:
        return xsecs
    grid = np.geomspace(pts[0], pts[-1], n)
    values = np.asarray(lambdify(E, expr, "numpy")(grid), dtype=float)
    return replace(xsecs, photon_energy=grid, photodecay=values)


def _band_series(reaction: Any, energy_unit: str, xsec_unit: str, label: str) -> list:
    """Band-average bars from ``reaction.photo.band_xsecs`` (empty when no bands)."""
    df = getattr(getattr(reaction, "photo", None), "band_xsecs", None)
    if df is None or df.empty:
        return []
    y = df["xsec"].to_numpy(dtype=float)
    keep = np.isfinite(y)
    if not keep.any():
        return []

    wavelength = energy_unit in _WAVELENGTH_UNITS

    def energy(col: str) -> np.ndarray:
        values = df[col].to_numpy(dtype=float)[keep]
        # An infinite photon energy is zero wavelength: an open *low* edge.
        out = np.full(values.shape, -np.inf if wavelength else np.inf)
        finite = np.isfinite(values)
        if finite.any():
            out[finite] = _units.convert_energy(values[finite], "eV", energy_unit)
        return out

    edge_a, edge_b = energy("lower"), energy("upper")

    return [
        {
            "x": energy("eavg"),
            "y": np.asarray(_units.convert_xsec(y[keep], "cm2", xsec_unit), dtype=float),
            "x_lo": np.minimum(edge_a, edge_b),
            "x_hi": np.maximum(edge_a, edge_b),
            "kind": "bars",
            "label": label,
        }
    ]


def xsec_series(
    reaction: Any,
    *,
    processes: str | list[str] | None = "all",
    energy_unit: str = "eV",
    xsec_unit: str = "Mb",
    prefix: bool = False,
    bands: bool = False,
) -> list[dict]:
    """Series dicts for one reaction's photo cross sections.

    Parameters
    ----------
    reaction
        Reaction-like object exposing ``photo.xsec``.
    processes : str, list[str] or None
        ``"all"`` (default), one key, or a list of keys.
    energy_unit, xsec_unit : str
        Axis units (see :mod:`jaff.plotting._units`).
    prefix : bool
        Prefix labels with the reaction equation (several reactions share
        the axes).
    bands : bool
        Add band-average bars from ``photo.band_xsecs``.

    Returns
    -------
    list[dict]
        One line series per process with data, plus optional bars; empty
        (with an info log) when the reaction has no cross sections.
    """
    procs = normalize_processes(processes)
    logger = JaffLogger().get_logger()
    photo = getattr(reaction, "photo", None)
    xsecs: XsecProps | None = photo.xsec if photo is not None else None

    if xsecs is None:
        logger.info(f"No cross sections available for: {reaction}")
        return []

    if "photodecay" in procs:
        xsecs = _sampled_xsecs(xsecs)

    energy = xsecs.photon_energy
    if energy is None:
        logger.info(f"No photon energy grid available for: {reaction}")
        return []
    x = np.asarray(_units.convert_energy(np.asarray(energy), "eV", energy_unit))
    head = f"{_label_of(reaction)} — " if prefix else ""
    out = []
    for p in procs:
        data = getattr(xsecs, p)
        if data is None:
            continue

        y = np.asarray(_units.convert_xsec(np.asarray(data), "cm2", xsec_unit))
        out.append({"x": x, "y": y, "label": head + _units.PROCESS_LABELS[p]})

    if not out:
        logger.info(f"No data for requested cross-section(s) {procs} in: {reaction}")
        return []
    if bands:
        out += _band_series(reaction, energy_unit, xsec_unit, f"{head}Band average")
    return out
