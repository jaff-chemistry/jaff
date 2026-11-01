# ABOUTME: Versioned JSON (de)serializer for SymPy expressions: compact list-tag
# ABOUTME: encoding used by .jaff files, plus a verbose dict encoding for debugging.
"""
Versioned JSON serializer / deserializer for SymPy expressions.

Public API
----------
``SCHEMA_VERSION``
    Integer version tag embedded in every serialized payload (currently ``2``).
    Consumers must check this value and reject payloads with an unrecognised
    version.
:func:`dumps` / :func:`loads`
    JSON string round-trip for a single SymPy expression with a full metadata
    envelope (``format``, ``schema_version``, ``sympy_version``).  The
    round-trip holds for the default compact encoding only: :func:`loads`
    (via :func:`from_jsonable`) decodes compact nodes, so a payload written
    with ``dumps(..., compact=False)`` currently raises
    :class:`~jaff.errors.SympyJsonError` on loading.
:func:`to_jsonable` / :func:`from_jsonable`
    Lower-level helpers that convert to/from a JSON-compatible Python object
    (list/number/dict) without wrapping it in the metadata envelope.  Used by
    the JAFF network serializer when expressions are embedded in a larger JSON
    document.

Encoding formats
----------------
Two encoder variants are provided:

* :class:`_Encoder` (verbose) -- each node is a ``{"type": ..., ...}`` dict.
  Useful for debugging; only the private :class:`_Decoder` reads it back.
* :class:`_EncoderCompact` (default) -- each node is a compact list whose
  first element is a short tag string.  Tag mapping:

  ======  ================
  Tag     SymPy type
  ======  ================
  ``T``   ``BooleanTrue``
  ``F``   ``BooleanFalse``
  ``S``   ``Symbol``
  ``I``   ``Integer``
  ``Q``   ``Rational``
  float   ``Float`` (raw number)
  ``Flt`` ``Float`` (with precision; decoded only, never written)
  ``Str`` internal ``Str``
  ``MS``  ``MatrixSymbol``
  ``ME``  ``MatrixElement``
  ``ECP`` ``ExprCondPair``
  ``LT``  ``StrictLessThan``
  ``GT``  ``StrictGreaterThan``
  ``PW``  ``Piecewise``
  ``Pow`` ``Pow``
  ``Add`` ``Add``
  ``Mul`` ``Mul``
  ``exp`` ``exp``
  ``log`` ``log``
  ``Max`` ``Max``
  ``Min`` ``Min``
  ``IB``  ``IndexedBase``
  ``Idx`` ``Indexed``
  ======  ================

Cross-version compatibility
---------------------------
Optional SymPy internals (``sympy.core.symbol.Str``,
``sympy.matrices.expressions.matexpr.MatrixElement``,
``sympy.functions.elementary.piecewise.ExprCondPair``) are imported with
``try/except`` at module load time.  If a build of SymPy does not expose
them the corresponding encoder/decoder branches are disabled gracefully.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, Dict, List, Mapping, Tuple

import sympy

from ..errors import SympyJsonError

try:
    from sympy.core.symbol import Str as _SympyStr
except Exception:  # pragma: no cover
    _SympyStr = None

try:
    from sympy.matrices.expressions.matexpr import MatrixElement as _MatrixElement
except Exception:  # pragma: no cover
    _MatrixElement = None

try:
    from sympy.functions.elementary.piecewise import ExprCondPair as _ExprCondPair
except Exception:  # pragma: no cover
    _ExprCondPair = None


SCHEMA_VERSION = 2
"""int: Current schema version for serialized SymPy expressions."""


def dumps(
    expr: sympy.Basic,
    *,
    indent: int = 2,
    sort_keys: bool = True,
    compact: bool = True,
    include_assumptions: bool = True,
) -> str:
    """
    Serialize a SymPy expression to a JSON string with a metadata envelope.

    The resulting string is a JSON object with keys ``format``,
    ``schema_version``, ``sympy_version``, and ``expr``.

    Parameters
    ----------
    expr : sympy.Basic
        The expression to serialize.
    indent : int, optional
        JSON indentation width (default ``2``).
    sort_keys : bool, optional
        Whether to sort JSON object keys (default ``True``).
    compact : bool, optional
        Use the compact list-based encoding (default ``True``).  Pass
        ``False`` to use the verbose dict-based encoding (useful for
        debugging); :func:`loads` cannot currently read that form back.
    include_assumptions : bool, optional
        Whether to embed symbol assumption flags (e.g. ``positive=True``) in
        the output (default ``True``).

    Returns
    -------
    str
        A JSON string representing the expression.
    """
    payload = {
        "format": "jaff.sympy_json",
        "schema_version": SCHEMA_VERSION,
        "sympy_version": sympy.__version__,
        "expr": to_jsonable(
            expr, compact=compact, include_assumptions=include_assumptions
        ),
    }
    return json.dumps(payload, indent=indent, sort_keys=sort_keys)


def loads(s: str) -> sympy.Basic:
    """
    Deserialize a SymPy expression from a JSON string produced by :func:`dumps`.

    Parameters
    ----------
    s : str
        A JSON string with a ``jaff.sympy_json`` envelope.

    Returns
    -------
    sympy.Basic
        The reconstructed SymPy expression.

    Raises
    ------
    SympyJsonError
        If the payload is not a ``jaff.sympy_json`` document or uses an
        unsupported ``schema_version``.
    """
    payload = json.loads(s)
    if not isinstance(payload, dict) or payload.get("format") != "jaff.sympy_json":
        raise SympyJsonError("Not a jaff.sympy_json payload")
    version = payload.get("schema_version")
    if version != SCHEMA_VERSION:
        raise SympyJsonError(f"Unsupported schema_version={version!r}")
    return from_jsonable(payload.get("expr"))


def to_jsonable(
    expr: sympy.Basic, *, compact: bool = True, include_assumptions: bool = True
) -> Any:
    """
    Convert a SymPy expression to a JSON-compatible Python object.

    Does not add the ``jaff.sympy_json`` envelope; use :func:`dumps` for a
    self-contained serialized form.

    Parameters
    ----------
    expr : sympy.Basic
        The expression to encode.
    compact : bool, optional
        Use the compact list-based encoding (default ``True``).
    include_assumptions : bool, optional
        Embed symbol assumption flags in the output (default ``True``).

    Returns
    -------
    list or float or dict
        A JSON-serializable object representing *expr*.

    Raises
    ------
    TypeError
        If *expr* is not a :class:`sympy.Basic` instance.
    SympyJsonError
        If the expression tree contains an unsupported SymPy node type.
    """
    if not isinstance(expr, sympy.Basic):
        raise TypeError(f"Expected sympy.Basic, got {type(expr)!r}")
    if compact:
        return _EncoderCompact(include_assumptions=include_assumptions).encode(expr)
    return _Encoder(include_assumptions=include_assumptions).encode(expr)


def from_jsonable(obj: Any) -> sympy.Basic:
    """
    Reconstruct a SymPy expression from a JSON-compatible Python object.

    Accepts output produced by :func:`to_jsonable` (compact form: list or
    number).

    Parameters
    ----------
    obj : list or int or float
        The JSON-compatible node to decode.

    Returns
    -------
    sympy.Basic
        The reconstructed SymPy expression.

    Raises
    ------
    SympyJsonError
        If *obj* is not a recognised compact node (must be a list or number).
    """
    if not isinstance(obj, (list, int, float)):
        raise SympyJsonError(f"Expected list/number node, got {type(obj)!r}")
    return _DecoderCompact().decode(obj)


@dataclass(frozen=True)
class _SymbolKey:
    """Hashable cache key for a :class:`sympy.Symbol` with its assumptions."""

    name: str
    assumptions: Tuple[Tuple[str, bool], ...]


@dataclass(frozen=True)
class _MatrixSymbolKey:
    """Hashable cache key for a :class:`sympy.MatrixSymbol`."""

    name: str
    rows: Any
    cols: Any


class _Encoder:
    """
    Verbose dict-based SymPy expression encoder.

    Produces ``{"type": "<TypeName>", ...}`` dicts.  Primarily useful for
    debugging; the default serialization path uses :class:`_EncoderCompact`.

    Parameters
    ----------
    include_assumptions : bool
        Whether to include symbol assumptions in the output.
    """

    def __init__(self, *, include_assumptions: bool) -> None:
        """Initialise the verbose dict-based encoder.

        Parameters
        ----------
        include_assumptions : bool
            Whether to include symbol assumptions in the encoded output.
        """
        self._include_assumptions = include_assumptions

    def encode(self, expr: sympy.Basic) -> Dict[str, Any]:
        """
        Encode *expr* to a JSON-compatible dict.

        Parameters
        ----------
        expr : sympy.Basic
            The expression node to encode.

        Returns
        -------
        dict
            A JSON-serializable dictionary representation of *expr*.

        Raises
        ------
        SympyJsonError
            If *expr* contains an unsupported SymPy node type.
        """
        if expr is sympy.true:
            return {"type": "BooleanTrue"}
        if expr is sympy.false:
            return {"type": "BooleanFalse"}

        if isinstance(expr, sympy.Symbol):
            payload = {"type": "Symbol", "name": expr.name}
            if self._include_assumptions:
                payload["assumptions"] = _encode_assumptions(expr)
            return payload

        if isinstance(expr, sympy.Integer):
            return {"type": "Integer", "value": int(expr)}

        if isinstance(expr, sympy.Rational) and not isinstance(expr, sympy.Integer):
            return {"type": "Rational", "p": int(expr.p), "q": int(expr.q)}

        if isinstance(expr, sympy.Float):
            return {
                "type": "Float",
                "value": _encode_float_17(expr),
                "prec": int(expr._prec),
            }

        if _SympyStr is not None and isinstance(expr, _SympyStr):
            return {"type": "Str", "value": str(expr)}

        if isinstance(expr, sympy.MatrixSymbol):
            rows, cols = expr.shape
            return {
                "type": "MatrixSymbol",
                "name": expr.name,
                "rows": self.encode(sympy.Integer(rows))
                if isinstance(rows, int)
                else self.encode(rows),
                "cols": self.encode(sympy.Integer(cols))
                if isinstance(cols, int)
                else self.encode(cols),
            }

        if _MatrixElement is not None and isinstance(expr, _MatrixElement):
            return {
                "type": "MatrixElement",
                "base": self.encode(expr.parent),
                "i": self.encode(expr.i),
                "j": self.encode(expr.j),
            }

        if _ExprCondPair is not None and isinstance(expr, _ExprCondPair):
            return {
                "type": "ExprCondPair",
                "expr": self.encode(expr.expr),
                "cond": self.encode(expr.cond),
            }

        if isinstance(expr, sympy.StrictLessThan):
            return {
                "type": "StrictLessThan",
                "lhs": self.encode(expr.lhs),
                "rhs": self.encode(expr.rhs),
            }

        if isinstance(expr, sympy.StrictGreaterThan):
            return {
                "type": "StrictGreaterThan",
                "lhs": self.encode(expr.lhs),
                "rhs": self.encode(expr.rhs),
            }

        if isinstance(expr, sympy.Piecewise):
            pairs = []
            for pair in expr.args:
                if _ExprCondPair is None or not isinstance(pair, _ExprCondPair):
                    raise SympyJsonError("Unexpected Piecewise arg type")
                pairs.append(self.encode(pair))
            return {"type": "Piecewise", "pairs": pairs}

        if isinstance(expr, sympy.Pow):
            base, exp = expr.args
            return {"type": "Pow", "base": self.encode(base), "exp": self.encode(exp)}

        if isinstance(expr, sympy.Add):
            args = [self.encode(a) for a in expr.args]
            return {"type": "Add", "args": args}

        if isinstance(expr, sympy.Mul):
            args = [self.encode(a) for a in expr.args]
            return {"type": "Mul", "args": args}

        if isinstance(expr, sympy.IndexedBase):
            shape = expr.shape if expr.shape else None
            return {
                "type": "IndexedBase",
                "data": {"name": expr.name, "shape": shape},
            }

        if isinstance(expr, sympy.Indexed):
            base = expr.base
            indices = [self.encode(i) for i in expr.indices]
            return {"type": "Indexed", "base": self.encode(base), "indices": indices}

        func = expr.func
        if func is sympy.exp:
            return {"type": "exp", "args": [self.encode(expr.args[0])]}
        if func is sympy.log:
            return {"type": "log", "args": [self.encode(a) for a in expr.args]}
        if func is sympy.Max:
            return {"type": "Max", "args": [self.encode(a) for a in expr.args]}
        if func is sympy.Min:
            return {"type": "Min", "args": [self.encode(a) for a in expr.args]}

        raise SympyJsonError(f"Unsupported SymPy node: {type(expr).__name__}")


class _EncoderCompact:
    """
    Compact list-based SymPy expression encoder.

    Produces short-tag lists (e.g. ``["S", "tgas"]`` for a Symbol) that are
    significantly smaller than the verbose dict form and are the default output
    of :func:`to_jsonable`.

    Parameters
    ----------
    include_assumptions : bool
        Whether to include symbol assumptions in the output.
    """

    def __init__(self, *, include_assumptions: bool) -> None:
        """Initialise the compact list-based encoder.

        Parameters
        ----------
        include_assumptions : bool
            Whether to include symbol assumptions in the encoded output.
        """
        self._include_assumptions = include_assumptions

    def encode(self, expr: sympy.Basic) -> List[Any]:
        """
        Encode *expr* to a compact JSON-compatible list.

        Parameters
        ----------
        expr : sympy.Basic
            The expression node to encode.

        Returns
        -------
        list or float
            A compact JSON-serializable representation of *expr*.

        Raises
        ------
        SympyJsonError
            If *expr* contains an unsupported SymPy node type.
        """
        if expr is sympy.true:
            return ["T"]
        if expr is sympy.false:
            return ["F"]

        if isinstance(expr, sympy.Symbol):
            assumptions = _encode_assumptions(expr) if self._include_assumptions else {}
            if assumptions:
                return ["S", expr.name, assumptions]
            return ["S", expr.name]

        if isinstance(expr, sympy.Integer):
            return ["I", int(expr)]

        if isinstance(expr, sympy.Rational) and not isinstance(expr, sympy.Integer):
            return ["Q", int(expr.p), int(expr.q)]

        if isinstance(expr, sympy.Float):
            return _encode_float_17(expr)

        if _SympyStr is not None and isinstance(expr, _SympyStr):
            return ["Str", str(expr)]

        if isinstance(expr, sympy.MatrixSymbol):
            rows, cols = expr.shape
            return [
                "MS",
                expr.name,
                self.encode(sympy.Integer(rows))
                if isinstance(rows, int)
                else self.encode(rows),
                self.encode(sympy.Integer(cols))
                if isinstance(cols, int)
                else self.encode(cols),
            ]

        if _MatrixElement is not None and isinstance(expr, _MatrixElement):
            return [
                "ME",
                self.encode(expr.parent),
                self.encode(expr.i),
                self.encode(expr.j),
            ]

        if _ExprCondPair is not None and isinstance(expr, _ExprCondPair):
            return ["ECP", self.encode(expr.expr), self.encode(expr.cond)]

        if isinstance(expr, sympy.StrictLessThan):
            return ["LT", self.encode(expr.lhs), self.encode(expr.rhs)]

        if isinstance(expr, sympy.StrictGreaterThan):
            return ["GT", self.encode(expr.lhs), self.encode(expr.rhs)]

        if isinstance(expr, sympy.Piecewise):
            pairs = []
            for pair in expr.args:
                if _ExprCondPair is None or not isinstance(pair, _ExprCondPair):
                    raise SympyJsonError("Unexpected Piecewise arg type")
                pairs.append(self.encode(pair))
            return ["PW", pairs]

        if isinstance(expr, sympy.Pow):
            base, exp = expr.args
            return ["Pow", self.encode(base), self.encode(exp)]

        if isinstance(expr, sympy.Add):
            args = [self.encode(a) for a in expr.args]
            return ["Add", args]

        if isinstance(expr, sympy.Mul):
            args = [self.encode(a) for a in expr.args]
            return ["Mul", args]

        if isinstance(expr, sympy.IndexedBase):
            name = expr.name
            # Convert shape to list of ints (may contain sympy Integers)
            shape = None
            if expr.shape:
                shape = [
                    int(s) if isinstance(s, sympy.Integer) else s for s in expr.shape
                ]
            return ["IB", name, shape]

        if isinstance(expr, sympy.Indexed):
            base = expr.base
            indices = [self.encode(i) for i in expr.indices]
            return ["Idx", self.encode(base), indices]

        func = expr.func
        if func is sympy.exp:
            return ["exp", self.encode(expr.args[0])]
        if func is sympy.log:
            return ["log", [self.encode(a) for a in expr.args]]
        if func is sympy.Max:
            return ["Max", [self.encode(a) for a in expr.args]]
        if func is sympy.Min:
            return ["Min", [self.encode(a) for a in expr.args]]

        raise SympyJsonError(f"Unsupported SymPy node: {type(expr).__name__}")


class _Decoder:
    """
    Verbose dict-based SymPy expression decoder.

    Decodes the output of :class:`_Encoder`.  Caches
    :class:`sympy.Symbol` and :class:`sympy.MatrixSymbol` objects so that
    identical symbols with the same assumptions share the same Python object.
    """

    def __init__(self) -> None:
        """Initialise the verbose dict-based decoder with empty symbol caches."""
        self._symbol_cache: Dict[_SymbolKey, sympy.Symbol] = {}
        self._matrix_symbol_cache: Dict[_MatrixSymbolKey, sympy.MatrixSymbol] = {}

    def decode(self, obj: Any) -> sympy.Basic:
        """
        Decode a verbose dict node to a SymPy expression.

        Parameters
        ----------
        obj : dict
            A JSON-compatible dict node produced by :class:`_Encoder`.

        Returns
        -------
        sympy.Basic
            The reconstructed SymPy expression.

        Raises
        ------
        SympyJsonError
            If *obj* is not a dict, has a missing/invalid ``type`` field,
            or uses an unsupported node type.
        """
        if not isinstance(obj, dict):
            raise SympyJsonError(f"Expected dict node, got {type(obj)!r}")
        t = obj.get("type")
        if not isinstance(t, str):
            raise SympyJsonError("Missing/invalid node type")

        if t == "BooleanTrue":
            return sympy.true
        if t == "BooleanFalse":
            return sympy.false

        if t == "Symbol":
            name = obj.get("name")
            if not isinstance(name, str):
                raise SympyJsonError("Symbol.name must be a string")
            assumptions = obj.get("assumptions") or {}
            if not isinstance(assumptions, dict):
                raise SympyJsonError("Symbol.assumptions must be a dict")
            cleaned = _decode_assumptions(assumptions)
            key = _SymbolKey(name=name, assumptions=tuple(sorted(cleaned.items())))
            sym = self._symbol_cache.get(key)
            if sym is None:
                sym = sympy.Symbol(name, **cleaned)
                self._symbol_cache[key] = sym
            return sym

        if t == "Integer":
            value = obj.get("value")
            if not isinstance(value, int):
                raise SympyJsonError("Integer.value must be an int")
            return sympy.Integer(value)

        if t == "Rational":
            p = obj.get("p")
            q = obj.get("q")
            if not isinstance(p, int) or not isinstance(q, int):
                raise SympyJsonError("Rational.p and Rational.q must be ints")
            return sympy.Rational(p, q)

        if t == "Float":
            prec = obj.get("prec")
            value = obj.get("value")
            if not isinstance(prec, int):
                raise SympyJsonError("Float.prec must be int")
            if not isinstance(value, (str, int, float)):
                raise SympyJsonError("Float.value must be str or number")
            return sympy.Float(value, prec)

        if t == "Str":
            value = obj.get("value")
            if not isinstance(value, str):
                raise SympyJsonError("Str.value must be a string")
            if _SympyStr is None:
                raise SympyJsonError("Str node unsupported in this SymPy build")
            return _SympyStr(value)

        if t == "MatrixSymbol":
            name = obj.get("name")
            if not isinstance(name, str):
                raise SympyJsonError("MatrixSymbol.name must be a string")
            rows = self.decode(obj.get("rows"))
            cols = self.decode(obj.get("cols"))
            key = _MatrixSymbolKey(name=name, rows=rows, cols=cols)
            msym = self._matrix_symbol_cache.get(key)
            if msym is None:
                msym = sympy.MatrixSymbol(name, rows, cols)
                self._matrix_symbol_cache[key] = msym
            return msym

        if t == "MatrixElement":
            base = self.decode(obj.get("base"))
            i = self.decode(obj.get("i"))
            j = self.decode(obj.get("j"))
            if _MatrixElement is None:
                raise SympyJsonError("MatrixElement node unsupported in this SymPy build")
            return _MatrixElement(base, i, j)

        if t == "ExprCondPair":
            expr = self.decode(obj.get("expr"))
            cond = self.decode(obj.get("cond"))
            if _ExprCondPair is None:
                raise SympyJsonError("ExprCondPair node unsupported in this SymPy build")
            return _ExprCondPair(expr, cond)

        if t == "StrictLessThan":
            lhs = self.decode(obj.get("lhs"))
            rhs = self.decode(obj.get("rhs"))
            return sympy.StrictLessThan(lhs, rhs)

        if t == "StrictGreaterThan":
            lhs = self.decode(obj.get("lhs"))
            rhs = self.decode(obj.get("rhs"))
            return sympy.StrictGreaterThan(lhs, rhs)

        if t == "Piecewise":
            pairs_obj = obj.get("pairs")
            if not isinstance(pairs_obj, list):
                raise SympyJsonError("Piecewise.pairs must be a list")
            pairs = []
            for p in pairs_obj:
                pair = self.decode(p)
                if _ExprCondPair is None or not isinstance(pair, _ExprCondPair):
                    raise SympyJsonError(
                        "Piecewise.pairs must contain ExprCondPair nodes"
                    )
                pairs.append((pair.expr, pair.cond))
            return sympy.Piecewise(*pairs, evaluate=False)

        if t == "Pow":
            base = self.decode(obj.get("base"))
            exp = self.decode(obj.get("exp"))
            return sympy.Pow(base, exp, evaluate=False)

        if t == "Add":
            args = _decode_args_list(obj.get("args"))
            return sympy.Add(*[self.decode(a) for a in args], evaluate=False)

        if t == "Mul":
            args = _decode_args_list(obj.get("args"))
            return sympy.Mul(*[self.decode(a) for a in args], evaluate=False)

        if t == "exp":
            args = _decode_args_list(obj.get("args"))
            if len(args) != 1:
                raise SympyJsonError("exp expects 1 arg")
            return sympy.exp(self.decode(args[0]))

        if t == "log":
            args = _decode_args_list(obj.get("args"))
            if len(args) not in (1, 2):
                raise SympyJsonError("log expects 1 or 2 args")
            return sympy.log(*[self.decode(a) for a in args])

        if t == "Max":
            args = _decode_args_list(obj.get("args"))
            return sympy.Max(*[self.decode(a) for a in args], evaluate=False)

        if t == "Min":
            args = _decode_args_list(obj.get("args"))
            return sympy.Min(*[self.decode(a) for a in args], evaluate=False)

        if t == "IndexedBase":
            data = obj.get("data")
            if not isinstance(data, dict):
                raise SympyJsonError("IndexedBase data must be a dict")
            name = data.get("name")
            if not isinstance(name, str):
                raise SympyJsonError("IndexedBase.name must be a string")
            shape = data.get("shape")
            # shape can be None or a tuple/list
            if shape is not None and isinstance(shape, list):
                shape = tuple(shape)
            return sympy.IndexedBase(name, shape=shape)

        if t == "Indexed":
            base = self.decode(obj.get("base"))
            indices_obj = obj.get("indices")
            if not isinstance(indices_obj, list):
                raise SympyJsonError("Indexed.indices must be a list")
            indices = [self.decode(i) for i in indices_obj]
            return sympy.Indexed(base, *indices)

        raise SympyJsonError(f"Unsupported node type: {t!r}")


class _DecoderCompact:
    """
    Compact list-based SymPy expression decoder.

    Decodes the output of :class:`_EncoderCompact` and is the default decoder
    used by :func:`from_jsonable`.  Symbol objects are cached so that
    identical symbols share the same Python object within a single decode call.
    """

    def __init__(self) -> None:
        """Initialise the compact list-based decoder with empty symbol caches."""
        self._symbol_cache: Dict[_SymbolKey, sympy.Symbol] = {}
        self._matrix_symbol_cache: Dict[_MatrixSymbolKey, sympy.MatrixSymbol] = {}

    def decode(self, obj: Any) -> sympy.Basic:
        """
        Decode a compact list (or number) node to a SymPy expression.

        Raw ``int`` / ``float`` values are decoded as
        ``sympy.Float(obj, 53)`` (double precision).

        Parameters
        ----------
        obj : list, int, or float
            A compact JSON node produced by :class:`_EncoderCompact`.

        Returns
        -------
        sympy.Basic
            The reconstructed SymPy expression.

        Raises
        ------
        SympyJsonError
            If *obj* is not a recognised compact node or contains invalid
            payload for a given tag.
        """
        if isinstance(obj, (int, float)):
            return sympy.Float(obj, 53)
        if not isinstance(obj, list) or not obj:
            raise SympyJsonError(f"Expected list node, got {type(obj)!r}")
        t = obj[0]
        if not isinstance(t, str):
            raise SympyJsonError("Missing/invalid node type")

        if t == "T":
            return sympy.true
        if t == "F":
            return sympy.false

        if t == "S":
            if len(obj) < 2 or not isinstance(obj[1], str):
                raise SympyJsonError("Symbol name missing/invalid")
            name = obj[1]
            assumptions = {}
            if len(obj) >= 3:
                if not isinstance(obj[2], dict):
                    raise SympyJsonError("Symbol assumptions must be a dict")
                assumptions = _decode_assumptions(obj[2])
            key = _SymbolKey(name=name, assumptions=tuple(sorted(assumptions.items())))
            sym = self._symbol_cache.get(key)
            if sym is None:
                sym = sympy.Symbol(name, **assumptions)
                self._symbol_cache[key] = sym
            return sym

        if t == "I":
            if len(obj) != 2 or not isinstance(obj[1], int):
                raise SympyJsonError("Integer value missing/invalid")
            return sympy.Integer(obj[1])

        if t == "Q":
            if (
                len(obj) != 3
                or not isinstance(obj[1], int)
                or not isinstance(obj[2], int)
            ):
                raise SympyJsonError("Rational values missing/invalid")
            return sympy.Rational(obj[1], obj[2])

        if t == "Flt":
            if len(obj) < 3:
                raise SympyJsonError("Float value/precision missing")
            value = obj[1]
            prec = obj[2]
            if not isinstance(value, str):
                raise SympyJsonError("Float.value must be str")
            if not isinstance(prec, int):
                raise SympyJsonError("Float.prec must be int")
            return sympy.Float(value, prec)

        if t == "Str":
            if len(obj) != 2 or not isinstance(obj[1], str):
                raise SympyJsonError("Str value missing/invalid")
            if _SympyStr is None:
                raise SympyJsonError("Str node unsupported in this SymPy build")
            return _SympyStr(obj[1])

        if t == "MS":
            if len(obj) != 4 or not isinstance(obj[1], str):
                raise SympyJsonError("MatrixSymbol name/shape missing/invalid")
            name = obj[1]
            rows = self.decode(obj[2])
            cols = self.decode(obj[3])
            key = _MatrixSymbolKey(name=name, rows=rows, cols=cols)
            msym = self._matrix_symbol_cache.get(key)
            if msym is None:
                msym = sympy.MatrixSymbol(name, rows, cols)
                self._matrix_symbol_cache[key] = msym
            return msym

        if t == "ME":
            if len(obj) != 4:
                raise SympyJsonError("MatrixElement payload missing/invalid")
            base = self.decode(obj[1])
            i = self.decode(obj[2])
            j = self.decode(obj[3])
            if _MatrixElement is None:
                raise SympyJsonError("MatrixElement node unsupported in this SymPy build")
            return _MatrixElement(base, i, j)

        if t == "ECP":
            if len(obj) != 3:
                raise SympyJsonError("ExprCondPair payload missing/invalid")
            expr = self.decode(obj[1])
            cond = self.decode(obj[2])
            if _ExprCondPair is None:
                raise SympyJsonError("ExprCondPair node unsupported in this SymPy build")
            return _ExprCondPair(expr, cond)

        if t == "LT":
            if len(obj) != 3:
                raise SympyJsonError("StrictLessThan payload missing/invalid")
            return sympy.StrictLessThan(self.decode(obj[1]), self.decode(obj[2]))

        if t == "GT":
            if len(obj) != 3:
                raise SympyJsonError("StrictGreaterThan payload missing/invalid")
            return sympy.StrictGreaterThan(self.decode(obj[1]), self.decode(obj[2]))

        if t == "PW":
            if len(obj) != 2 or not isinstance(obj[1], list):
                raise SympyJsonError("Piecewise pairs missing/invalid")
            pairs = []
            for p in obj[1]:
                pair = self.decode(p)
                if _ExprCondPair is None or not isinstance(pair, _ExprCondPair):
                    raise SympyJsonError(
                        "Piecewise pairs must contain ExprCondPair nodes"
                    )
                pairs.append((pair.expr, pair.cond))
            return sympy.Piecewise(*pairs, evaluate=False)

        if t == "Pow":
            if len(obj) != 3:
                raise SympyJsonError("Pow payload missing/invalid")
            return sympy.Pow(self.decode(obj[1]), self.decode(obj[2]), evaluate=False)

        if t == "Add":
            if len(obj) != 2 or not isinstance(obj[1], list):
                raise SympyJsonError("Add args missing/invalid")
            return sympy.Add(*[self.decode(a) for a in obj[1]], evaluate=False)

        if t == "Mul":
            if len(obj) != 2 or not isinstance(obj[1], list):
                raise SympyJsonError("Mul args missing/invalid")
            return sympy.Mul(*[self.decode(a) for a in obj[1]], evaluate=False)

        if t == "exp":
            if len(obj) != 2:
                raise SympyJsonError("exp expects 1 arg")
            return sympy.exp(self.decode(obj[1]))

        if t == "log":
            if len(obj) != 2 or not isinstance(obj[1], list):
                raise SympyJsonError("log args missing/invalid")
            if len(obj[1]) not in (1, 2):
                raise SympyJsonError("log expects 1 or 2 args")
            return sympy.log(*[self.decode(a) for a in obj[1]])

        if t == "Max":
            if len(obj) != 2 or not isinstance(obj[1], list):
                raise SympyJsonError("Max args missing/invalid")
            return sympy.Max(*[self.decode(a) for a in obj[1]], evaluate=False)

        if t == "Min":
            if len(obj) != 2 or not isinstance(obj[1], list):
                raise SympyJsonError("Min args missing/invalid")
            return sympy.Min(*[self.decode(a) for a in obj[1]], evaluate=False)

        if t == "IB":
            # IndexedBase: ["IB", name, shape_or_none]
            if len(obj) < 2 or not isinstance(obj[1], str):
                raise SympyJsonError("IndexedBase name missing/invalid")
            name = obj[1]
            shape = None
            if len(obj) >= 3:
                shape_obj = obj[2]
                if shape_obj is not None:
                    if isinstance(shape_obj, list):
                        shape = tuple(shape_obj)
                    elif isinstance(shape_obj, tuple):
                        shape = shape_obj
                    else:
                        shape = self.decode(shape_obj)
                        if isinstance(shape, (list, tuple)):
                            shape = tuple(shape) if isinstance(shape, list) else shape
            return sympy.IndexedBase(name, shape=shape)

        if t == "Idx":
            # Indexed: ["Idx", base, [indices...]]
            if len(obj) < 3:
                raise SympyJsonError("Indexed payload missing")
            base = self.decode(obj[1])
            indices_list = obj[2]
            if not isinstance(indices_list, list):
                raise SympyJsonError("Indexed indices must be a list")
            indices = [self.decode(i) for i in indices_list]
            return sympy.Indexed(base, *indices)

        raise SympyJsonError(f"Unsupported node type: {t!r}")


def _decode_args_list(value: Any) -> List[Dict[str, Any]]:
    """
    Validate and return the ``args`` list from a verbose dict node.

    Parameters
    ----------
    value : Any
        The ``args`` field extracted from a JSON node dict.

    Returns
    -------
    list of dict
        The validated list of argument dicts.

    Raises
    ------
    SympyJsonError
        If *value* is not a list or any element is not a dict.
    """
    if not isinstance(value, list):
        raise SympyJsonError("args must be a list")
    for item in value:
        if not isinstance(item, dict):
            raise SympyJsonError("args must contain dict nodes")
    return value


def _encode_float_17(value: sympy.Float) -> float:
    """
    Round-trip *value* through a 17-significant-digit string representation.

    17 decimal digits are sufficient to uniquely identify any IEEE 754 double,
    ensuring lossless round-trip through JSON's ``float`` type.

    Parameters
    ----------
    value : sympy.Float
        The SymPy float to encode.

    Returns
    -------
    float
        A Python float equal (to double precision) to *value*.
    """
    return float(str(sympy.Float(value, 17)))


def _encode_assumptions(sym: sympy.Symbol) -> Dict[str, bool]:
    """
    Extract bool-valued assumption flags from a :class:`sympy.Symbol`.

    Only ``str`` keys paired with ``bool`` values are included; internal
    SymPy assumptions with non-bool values are silently dropped.

    Parameters
    ----------
    sym : sympy.Symbol
        The symbol whose assumptions to extract.

    Returns
    -------
    dict[str, bool]
        Mapping of assumption name to value.
    """
    out: Dict[str, bool] = {}
    for k, v in (sym.assumptions0 or {}).items():
        if isinstance(k, str) and isinstance(v, bool):
            out[k] = v
    return out


def _decode_assumptions(assumptions: Mapping[str, Any]) -> Dict[str, bool]:
    """
    Filter a raw assumptions dict to keep only ``str``-keyed ``bool`` values.

    Parameters
    ----------
    assumptions : Mapping[str, Any]
        Raw assumptions mapping loaded from JSON.

    Returns
    -------
    dict[str, bool]
        Cleaned mapping safe to pass as keyword arguments to
        :class:`sympy.Symbol`.
    """
    out: Dict[str, bool] = {}
    for k, v in assumptions.items():
        if isinstance(k, str) and isinstance(v, bool):
            out[k] = v
    return out
