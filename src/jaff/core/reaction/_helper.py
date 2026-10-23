from typing import Any


def to_float_or_none(value: Any) -> float | None:
    """Coerce a band quantity to ``float``, or ``None`` when not representable.

    Band edges and averages may be plain numbers, SymPy numeric objects, or
    ``sympy.oo`` (open upper band, which casts to ``float('inf')``).  A value
    of ``None`` or a still-symbolic expression (e.g. the band-average cross
    section of a ``photoCrossSection`` that depends on ``tgas``) maps to
    ``None`` so it becomes ``NaN`` in a :class:`pandas.DataFrame`.

    Parameters
    ----------
    value : Any
        The quantity to convert.

    Returns
    -------
    float or None
        ``float(value)``, or ``None`` when *value* is ``None`` or cannot be
        cast.
    """
    if value is None:
        return None

    try:
        return float(value)
    except (TypeError, ValueError):
        return None
