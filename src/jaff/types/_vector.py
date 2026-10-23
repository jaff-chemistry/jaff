# ABOUTME: Vector: list subclass with NumPy-style element-wise arithmetic and
# ABOUTME: comparison operators (scalars broadcast; equal-length lists pair up).
"""
Element-wise arithmetic vector type.

This module provides :class:`Vector`, a :class:`list` subclass that overrides
all standard Python arithmetic and comparison operators to work *element-wise*
(NumPy-style) rather than with Python's default list semantics (where ``+``
concatenates and ``*`` repeats).

Scalar operands broadcast across every element; same-length list operands
produce a new :class:`Vector` via a pairwise operation.

Examples
--------
>>> v = Vector([1, 2, 3])
>>> v + 10
[11, 12, 13]
>>> v * Vector([2, 3, 4])
[2, 6, 12]
>>> -v
[-1, -2, -3]
"""

import operator
from typing import Callable, Generic, TypeVar

T = TypeVar("T")
U = TypeVar("U")


class Vector(list, Generic[T]):
    """
    List subclass with element-wise arithmetic and comparison operators.

    All binary arithmetic operators (``+``, ``-``, ``*``, ``/``, ``//``,
    ``%``, ``**``) and comparison operators (``==``, ``!=``, ``<``, ``<=``,
    ``>``, ``>=``) are overridden to apply element-wise, consistent with
    NumPy array semantics.

    Scalar operands broadcast to every element.  List/Vector operands must
    have the same length; a :exc:`ValueError` is raised otherwise.

    In-place operators (``+=``, ``-=``, etc.) update the vector contents in
    place and return ``self``.

    Parameters
    ----------
    *args
        Forwarded to :class:`list`.

    Examples
    --------
    >>> v = Vector([1, 2, 3])
    >>> v + 1
    [2, 3, 4]
    >>> v * Vector([10, 20, 30])
    [10, 40, 90]
    >>> v == 2
    [False, True, False]
    """

    # ------------------------------------------------------------------
    # Core dispatch helpers
    # ------------------------------------------------------------------

    def _apply_op(
        self, other: "Vector[T] | T", op: Callable, reverse: bool = False
    ) -> "Vector[T]":
        """
        Apply a binary operator element-wise and return a new :class:`Vector`.

        Parameters
        ----------
        other : Vector[T] or T
            Right-hand operand.  If a list (or Vector), must be the same
            length as ``self``; scalars broadcast.
        op : Callable
            Binary operator function (e.g. ``operator.add``).
        reverse : bool, optional
            If ``True``, swap the operand order so that the operation becomes
            ``op(other_elem, self_elem)``; used to implement reflected
            operators such as ``__radd__``.

        Returns
        -------
        Vector[T]
            New vector with the results.

        Raises
        ------
        ValueError
            If *other* is a list and its length differs from ``len(self)``.
        """
        if isinstance(other, list):
            if len(self) != len(other):
                raise ValueError(
                    "Vectors must be of the same length for element-wise operations."
                )

            if reverse:
                return Vector(op(b, a) for a, b in zip(self, other))
            return Vector(op(a, b) for a, b in zip(self, other))

        else:
            # Scalar broadcast: apply the operator between each element and
            # the scalar operand.
            if reverse:
                return Vector(op(other, a) for a in self)
            return Vector(op(a, other) for a in self)

    def _apply_iop(self, other: "Vector[T] | T", op: Callable) -> "Vector[T]":
        """
        Apply a binary operator in-place and return ``self``.

        Computes the result via :meth:`_apply_op` and overwrites the vector's
        contents using slice assignment.

        Parameters
        ----------
        other : Vector[T] or T
            Right-hand operand (same rules as :meth:`_apply_op`).
        op : Callable
            Binary operator function.

        Returns
        -------
        Vector[T]
            ``self`` after in-place update.
        """
        result = self._apply_op(other, op)
        # Slice assignment updates the list contents without replacing the object.
        self[:] = result
        return self

    # ------------------------------------------------------------------
    # Arithmetic operators
    # ------------------------------------------------------------------

    def __add__(self, other):
        """Element-wise addition."""
        return self._apply_op(other, operator.add)

    def __radd__(self, other):
        """Reflected element-wise addition."""
        return self._apply_op(other, operator.add, reverse=True)

    def __iadd__(self, other):
        """In-place element-wise addition."""
        return self._apply_iop(other, operator.add)

    def __sub__(self, other):
        """Element-wise subtraction."""
        return self._apply_op(other, operator.sub)

    def __rsub__(self, other):
        """Reflected element-wise subtraction."""
        return self._apply_op(other, operator.sub, reverse=True)

    def __isub__(self, other):
        """In-place element-wise subtraction."""
        return self._apply_iop(other, operator.sub)

    def __mul__(self, other):
        """Element-wise multiplication."""
        return self._apply_op(other, operator.mul)

    def __rmul__(self, other):
        """Reflected element-wise multiplication."""
        return self._apply_op(other, operator.mul, reverse=True)

    def __imul__(self, other):
        """In-place element-wise multiplication."""
        return self._apply_iop(other, operator.mul)

    def __truediv__(self, other):
        """Element-wise true division."""
        return self._apply_op(other, operator.truediv)

    def __rtruediv__(self, other):
        """Reflected element-wise true division."""
        return self._apply_op(other, operator.truediv, reverse=True)

    def __itruediv__(self, other):
        """In-place element-wise true division."""
        return self._apply_iop(other, operator.truediv)

    def __floordiv__(self, other):
        """Element-wise floor division."""
        return self._apply_op(other, operator.floordiv)

    def __rfloordiv__(self, other):
        """Reflected element-wise floor division."""
        return self._apply_op(other, operator.floordiv, reverse=True)

    def __ifloordiv__(self, other):
        """In-place element-wise floor division."""
        return self._apply_iop(other, operator.floordiv)

    def __mod__(self, other):
        """Element-wise modulo."""
        return self._apply_op(other, operator.mod)

    def __rmod__(self, other):
        """Reflected element-wise modulo."""
        return self._apply_op(other, operator.mod, reverse=True)

    def __imod__(self, other):
        """In-place element-wise modulo."""
        return self._apply_iop(other, operator.mod)

    def __pow__(self, other):
        """Element-wise exponentiation."""
        return self._apply_op(other, operator.pow)

    def __rpow__(self, other):
        """Reflected element-wise exponentiation."""
        return self._apply_op(other, operator.pow, reverse=True)

    def __ipow__(self, other):
        """In-place element-wise exponentiation."""
        return self._apply_iop(other, operator.pow)

    # ------------------------------------------------------------------
    # Comparison operators
    # ------------------------------------------------------------------

    def __eq__(self, other):
        """Element-wise equality test; returns a Vector of booleans."""
        return self._apply_op(other, operator.eq)

    def __ne__(self, other):
        """Element-wise inequality test; returns a Vector of booleans."""
        return self._apply_op(other, operator.ne)

    def __lt__(self, other):
        """Element-wise less-than test; returns a Vector of booleans."""
        return self._apply_op(other, operator.lt)

    def __le__(self, other):
        """Element-wise less-than-or-equal test; returns a Vector of booleans."""
        return self._apply_op(other, operator.le)

    def __gt__(self, other):
        """Element-wise greater-than test; returns a Vector of booleans."""
        return self._apply_op(other, operator.gt)

    def __ge__(self, other):
        """Element-wise greater-than-or-equal test; returns a Vector of booleans."""
        return self._apply_op(other, operator.ge)

    # ------------------------------------------------------------------
    # Unary operators
    # ------------------------------------------------------------------

    def __neg__(self) -> "Vector[T]":
        """Element-wise negation."""
        return Vector(-a for a in self)

    def __pos__(self) -> "Vector[T]":
        """Element-wise unary plus (returns a copy)."""
        return Vector(+a for a in self)

    def __abs__(self) -> "Vector[T]":
        """Element-wise absolute value."""
        return Vector(abs(a) for a in self)

    def __invert__(self) -> "Vector[T]":
        """Element-wise bitwise inversion."""
        return Vector(~a for a in self)

    # ------------------------------------------------------------------
    # Conversion utilities
    # ------------------------------------------------------------------

    def as_string(self) -> "Vector[str]":
        """
        Convert every element to :class:`str` and return a new Vector.

        Returns
        -------
        Vector[str]
            New vector with each element cast to ``str``.
        """
        return Vector(str(a) for a in self)

    def as_type(self, typ: Callable[[T], U]) -> "Vector[U]":
        """
        Cast every element to a new type using a callable.

        Parameters
        ----------
        typ : Callable[[T], U]
            A callable (type constructor or conversion function) applied to
            each element.

        Returns
        -------
        Vector[U]
            New vector with all elements converted to type *U*.

        Examples
        --------
        >>> Vector(["1", "2", "3"]).as_type(int)
        [1, 2, 3]
        """
        return Vector(typ(x) for x in self)
