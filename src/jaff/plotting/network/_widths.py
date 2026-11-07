# ABOUTME: Rate-scaled arrow widths for net.visualize: log-normalised between a
# ABOUTME: min and max width, floored 12 decades below the largest; legend key

from __future__ import annotations

import math

import numpy as np

#: Rates more than this many decades below the largest get the minimum width.
DECADES = 12


def _range(values: np.ndarray) -> tuple[float, float] | None:
    """log10 of the smallest (floored) and largest positive finite value."""
    good = values[np.isfinite(values) & (values > 0)]
    if good.size == 0:
        return None

    top = float(np.log10(good.max()))

    return max(float(np.log10(good.min())), top - DECADES), top


def _scale(logs: np.ndarray, low: float, top: float, lo: float, hi: float) -> np.ndarray:
    """Widths for log10 values *logs* on the range [*low*, *top*]."""
    if top == low:
        return np.full(logs.shape, hi)

    t = np.clip((logs - low) / (top - low), 0.0, 1.0)

    return lo + t * (hi - lo)


def edge_widths(values: np.ndarray, lo: float, hi: float) -> np.ndarray:
    """Line widths for rate *values*, log-normalised between *lo* and *hi*.

    The largest value gets *hi* and the smallest *lo*, unless it is more than
    :data:`DECADES` decades below the largest, which is where *lo* starts then.
    Zero, negative and non-finite values get *lo*; equal values all get *hi*.
    """
    values = np.asarray(values, dtype=float)
    widths = np.full(values.shape, lo)
    span = _range(values)
    if span is None:
        return widths

    good = np.isfinite(values) & (values > 0)
    widths[good] = _scale(np.log10(values[good]), *span, lo, hi)

    return widths


def width_key(values: np.ndarray, lo: float, hi: float) -> list[tuple[float, float]]:
    """Up to three ``(value, width)`` legend entries, powers of ten, widest first."""
    span = _range(np.asarray(values, dtype=float))
    if span is None:
        return []

    low, top = span
    first, last = math.floor(top), math.ceil(low)
    if first < last:  # no power of ten inside the range: show the largest value
        return [(10.0**top, hi)]

    decades = sorted({first, round((first + last) / 2), last}, reverse=True)
    widths = _scale(np.array(decades, dtype=float), low, top, lo, hi)

    return [(10.0**d, float(w)) for d, w in zip(decades, widths)]
