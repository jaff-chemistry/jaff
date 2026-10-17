from __future__ import annotations

from typing import TYPE_CHECKING

from sympy import Basic, Expr, Piecewise, sympify

from ....types import Catalogue
from . import RateSegment

if TYPE_CHECKING:
    from .. import Reaction


class RateSegments(Catalogue[RateSegment]):
    """Ordered catalogue of a reaction's rate pieces over temperature ranges.

    A reaction whose rate is defined piecewise across several disjoint
    temperature ranges holds one :class:`RateSegment` per range, keyed by its
    ``(tmin, tmax)`` bounds.  The pieces are collapsed into a single SymPy
    expression by :meth:`evaluate_equivalent_rate`.

    Parameters
    ----------
    segments : list[RateSegment]
        Initial rate segments.
    mode : str
        Out-of-range behaviour.  ``"clip"`` holds the boundary rate below the
        first / above the last range; any other value lets the outermost
        segments extend unbounded.
    """

    def __init__(self, segments: list[RateSegment], mode: str):
        _by_temp: dict[tuple[float | None, float | None], RateSegment] = {}

        if segments is not None:
            _by_temp = {(rs.tmin, rs.tmax): rs for rs in segments}

        super().__init__(segments, _by_temp)
        self.mode: str = mode

    def add(self, segment: RateSegment):
        """Append a rate segment, keyed by its ``(tmin, tmax)`` range.

        If a segment already exists for the same ``(tmin, tmax)`` range it is
        replaced by *segment*.

        Parameters
        ----------
        segment : RateSegment

        Raises
        ------
        ValueError
            If *segment* is not a :class:`RateSegment`.
        """
        if not isinstance(segment, RateSegment):
            raise ValueError(f"'{segment}' must be an instance of 'RateSegment'")

        tup = (segment.tmin, segment.tmax)
        if tup in self._by_prop:
            self._list.remove(self._by_prop[tup])

        self._by_prop[tup] = segment
        self._list.append(segment)
        self.count = len(self._list)

    def __repr__(self):
        return "<RateSegment Object>"

    def evaluate_equivalent_rate(self, r: Reaction) -> Expr:
        """Collapse the segments into a single SymPy ``Piecewise`` rate.

        Segments must be sorted by ascending temperature (call :meth:`sort`
        first).  Each range contributes its own rate; the gap between two
        adjacent ranges is bridged by interpolating between the two bounding
        rates: geometrically (linear in ``log k``) when both are positive,
        linearly otherwise (e.g. negative or symbolic boundary rates).  Out-of-range behaviour follows :attr:`mode` (``"clip"`` holds
        the boundary value, otherwise the end segments extend unbounded).

        A lone segment with a temperature-independent rate, or one that is fully
        unbounded (no ``tmin`` *and* no ``tmax``, i.e. valid at every
        temperature), is returned as-is without wrapping in a ``Piecewise``.

        Under ``"clip"`` each **defined** outer bound is held independently: a
        range open on one side (e.g. ``tmin`` set, ``tmax`` ``None``) clamps its
        closed side and extrapolates the open side.  Interior bounds of a
        multi-range rate must always be defined.

        Returns
        -------
        Expr
            The equivalent rate expression in ``tgas``.

        Raises
        ------
        ValueError
            If adjacent segments have undefined bounds or overlapping ranges.
        """
        # Local import: a module-level one would cycle via jaff.core.network.
        from ...network import NetworkSymbols

        tgas = NetworkSymbols.tgas
        ls = self._list
        first = ls[0]
        last = ls[-1]

        # No piecewise needed only for a lone segment: a single, temperature-
        # independent rate, or one fully unbounded (valid at every temperature).
        # With multiple segments the later pieces may differ, so never shortcut.
        if len(ls) == 1 and not first.rate.has(tgas):
            return first.rate
        if len(ls) == 1 and first.tmin is None and first.tmax is None:
            return first.rate

        def _body(seg: RateSegment) -> tuple[Expr, Basic | bool]:
            cond: Basic | bool = True if seg.tmax is None else tgas < seg.tmax
            return (seg.rate, cond)

        segs: list[tuple[Expr, Basic | bool]] = []

        # Lower edge: hold the boundary rate below tmin (clip, when a lower
        # bound exists); an open lower end just extrapolates.
        if self.mode == "clip" and first.tmin is not None:
            segs.append((first.rate.xreplace({tgas: first.tmin}), tgas < first.tmin))

        segs.append(_body(first))

        # Interpolation gaps + subsequent ranges (interior bounds must exist).
        for i, seg in enumerate(ls[1:]):
            prev = ls[i]
            if prev.tmax is None or seg.tmin is None:
                raise ValueError(
                    f"Reaction {r} should have a well defined temperature since it's a multi temperature range reaction"
                )
            if prev.tmax > seg.tmin:
                raise ValueError(f"Temperature ranges shouldn't overlap for reaction {r}")

            a = prev.tmax  # left boundary
            b = seg.tmin  # right boundary
            if a != b:
                left = sympify(prev.rate.xreplace({tgas: a}))
                right = sympify(seg.rate.xreplace({tgas: b}))
                w = (tgas - a) / (b - a)
                if left.is_positive and right.is_positive:
                    # Geometric: rates spanning decades stay smooth in log
                    # space instead of dipping toward zero near the larger end.
                    # Written as one power so no huge constant is factored out
                    # (exp(a + b*T) gets split into exp(a)*exp(b*T), and
                    # exp(a) overflows a double for steep gaps).
                    interp = left * (right / left) ** w
                else:
                    interp = left * (1 - w) + right * w
                segs.append((interp, tgas < seg.tmin))

            # Always emit this range's own body (touching ranges have no gap).
            segs.append(_body(seg))

        # Upper edge: hold the boundary rate above tmax (clip, when an upper
        # bound exists); otherwise the last segment extends unbounded.
        if self.mode == "clip" and last.tmax is not None:
            segs.append((last.rate.xreplace({tgas: last.tmax}), True))
        else:
            segs[-1] = (segs[-1][0], True)

        return Piecewise(*segs)

    def sort(self) -> "RateSegments":
        """Sort the segments in place by ascending lower temperature bound.

        A segment with no lower bound (``tmin`` is ``None``) is open toward
        low temperatures, so it sorts first.  Only the lowest range may be open
        below and only the highest may be open above.

        Returns
        -------
        RateSegments
            ``self``, to allow chaining (e.g. ``segs.sort().evaluate_...()``).

        Raises
        ------
        ValueError
            If a range other than the lowest has no ``tmin``, or a range other
            than the highest has no ``tmax``.
        """
        self._list = sorted(self._list, key=lambda s: (s.tmin is not None, s.tmin))

        ranges = [(s.tmin, s.tmax) for s in self._list]
        for i, (tmin, tmax) in enumerate(ranges):
            if tmin is None and i > 0:
                raise ValueError(
                    f"Only the lowest temperature range may be open below; got {ranges}"
                )

            if tmax is None and i < len(ranges) - 1:
                raise ValueError(
                    f"Only the highest temperature range may be open above; got {ranges}"
                )

        return self
