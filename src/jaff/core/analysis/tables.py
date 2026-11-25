# ABOUTME: The tables of net.analysis: a value and a label per parameter, species
# ABOUTME: density or unknown function, with set/get/pop and the value checks
"""Tables of numeric values for a network's unknowns."""

from __future__ import annotations

import builtins
import difflib
import math
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, replace
from typing import Any

import numpy as np
from sympy import Basic, Expr, Indexed, Symbol
from sympy.core.function import AppliedUndef

__all__ = ["Entry", "AnalysisTable", "Parameters", "Densities", "Functions"]

#: Keys an entry given as a mapping may have.
ENTRY_FIELDS = frozenset({"value", "label"})


@dataclass(frozen=True)
class Entry:
    """The value and legend label of one name; both ``None`` when unset.

    Attributes
    ----------
    value : Any
        The value, or ``None`` when unset.
    label : str or None
        Legend label; ``{}`` in it takes the formatted value. ``None`` keeps the
        entry out of legends.
    """

    value: Any = None
    label: str | None = None


class AnalysisTable:
    """A fixed set of names, each with an :class:`Entry`.

    The names come from the network and cannot change; values are checked by
    the subclass when they are set.

    Parameters
    ----------
    names : iterable of str
        The table's names, in order.
    """

    #: What one name is, for messages.
    kind = "name"

    def __init__(self, names: Iterable[str]):
        self._entries: dict[str, Entry] = {name: Entry() for name in names}

    def __contains__(self, name: object) -> bool:
        return name in self._entries

    def __iter__(self):
        return iter(self._entries)

    def __len__(self) -> int:
        return len(self._entries)

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.values()!r})"

    @property
    def names(self) -> tuple[str, ...]:
        """The table's names, in order."""
        return tuple(self._entries)

    def set(self, name: str, value: Any, label: str | bool | None = None) -> None:
        """Set the value of *name*, and its label if *label* is given.

        Parameters
        ----------
        name : str
            One of :attr:`names`.
        value : Any
            The value, checked by the table.
        label : str or bool, optional
            Legend label; ``{}`` in it takes the formatted value. ``True`` gives
            the table's default label, ``False`` removes the label and ``None``
            (the default) keeps the current one.

        Raises
        ------
        KeyError
            If *name* is not in the table.
        ValueError
            If *value* is not allowed.
        """
        key = self._key(name)
        entry = replace(self._entries[key], value=self._check(key, value))
        if label is True:
            entry = replace(entry, label=self.default_label(key))
        elif label is False:
            entry = replace(entry, label=None)
        elif label is not None:
            entry = replace(entry, label=label)
        self._entries[key] = entry

    def update(self, entries: Mapping[str, Any]) -> None:
        """Set several entries: a value, or a mapping with ``value`` and ``label``."""
        for name, given in entries.items():
            if not isinstance(given, Mapping):
                self.set(name, given)
                continue
            extra = set(given) - ENTRY_FIELDS
            if extra:
                raise ValueError(
                    f"{self.kind} {name!r}: unknown entry field(s) "
                    f"{', '.join(sorted(extra))}; expected value and label"
                )
            self.set(name, given.get("value"), given.get("label"))

    def get(self, name: str | None = None) -> Any:
        """The entry of *name*, or a copy of every entry without *name*."""
        if name is None:
            return dict(self._entries)

        return self._entries[self._key(name)]

    def pop(self, name: str | None = None) -> Entry | None:
        """Clear *name* and return its old entry; without *name* clear every entry."""
        if name is None:
            self._entries = {key: Entry() for key in self._entries}
            return None
        key = self._key(name)
        old = self._entries[key]
        self._entries[key] = Entry()

        return old

    def values(self) -> dict[str, Any]:
        """The value of every name that has one."""
        return {k: e.value for k, e in self._entries.items() if e.value is not None}

    def missing(self, names: Iterable[str]) -> builtins.set[str]:
        """The *names* without a value (names outside the table included)."""
        return {
            n for n in names if n not in self._entries or self._entries[n].value is None
        }

    def default_label(self, name: str) -> str:
        """The label ``label=True`` gives *name*."""
        return f"{name} = {{}}"

    def _check(self, name: str, value: Any) -> Any:
        """The value to store for *name*; raises ``ValueError`` if not allowed."""
        return value

    def _key(self, name: str) -> str:
        if name in self._entries:
            return name
        close = difflib.get_close_matches(name, self._entries, n=3)
        hint = f"; did you mean {', '.join(close)}?" if close else ""
        raise KeyError(f"{name!r} is not a {self.kind} of the network{hint}")


def _number(kind: str, name: str, value: Any) -> float:
    """*value* as a finite float."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{kind} {name!r} needs a number, got {value!r}") from None
    if not math.isfinite(number):
        raise ValueError(f"{kind} {name!r} needs a finite number, got {value!r}")

    return number


class Parameters(AnalysisTable):
    """Values of the symbols in the rate coefficients (``tgas``, ``av``, ...).

    Each value is a finite number.
    """

    kind = "parameter"

    def default_label(self, name: str) -> str:
        """``$name = {}$``."""
        return f"${name} = {{}}$"

    def _check(self, name: str, value: Any) -> float:
        return _number(self.kind, name, value)


class Densities(AnalysisTable):
    """Number density (cm⁻³) of each species, by species name.

    Each value is a finite number, zero or more.

    Parameters
    ----------
    names : iterable of str
        Species names, in species-index order.
    latex : mapping of str to str
        LaTeX of each species, for the default label.
    """

    kind = "species"

    def __init__(self, names: Iterable[str], latex: Mapping[str, str]):
        super().__init__(names)
        self._latex = dict(latex)

    def update(self, entries: Mapping[str, Any] | Sequence[float] | np.ndarray) -> None:
        """Set several densities: by name, or an array in species order."""
        if isinstance(entries, Mapping):
            super().update(entries)
            return
        values = np.asarray(entries, dtype=float)
        if values.shape != (len(self),):
            raise ValueError(
                f"A density array needs one value per species ({len(self)}), "
                f"got shape {values.shape}"
            )
        for name, value in zip(self.names, values):
            self.set(name, value)

    def default_label(self, name: str) -> str:
        """``$n(<species LaTeX>) = {}$``."""
        return f"$n({self._latex.get(name, name)}) = {{}}$"

    def _check(self, name: str, value: Any) -> float:
        number = _number("density of", name, value)
        if number < 0:
            raise ValueError(f"density of {name!r} cannot be negative, got {value!r}")

        return number


class Functions(AnalysisTable):
    """Value of each unknown function, one per function whatever its arguments.

    A value is a number, a sympy expression of the parameters, or a callable
    without arguments returning either.

    Parameters
    ----------
    names : iterable of str
        Unknown function names.
    parameters : Parameters
        The parameters an expression may use.
    """

    kind = "function"

    def __init__(self, names: Iterable[str], parameters: Parameters):
        super().__init__(names)
        self._parameters = parameters

    def _check(self, name: str, value: Any) -> Any:
        if callable(value) and not isinstance(value, Basic):
            return value
        if isinstance(value, Basic):
            return self.check_expression(name, value)

        return _number(self.kind, name, value)

    def check_expression(self, name: str, expr: Basic) -> Expr:
        """*expr* if it uses only parameters and no unknown function."""
        if not isinstance(expr, Expr):
            raise ValueError(f"function {name!r} needs an expression, got {expr!r}")
        unknown = sorted(f.func.__name__ for f in expr.atoms(AppliedUndef))
        if unknown:
            raise ValueError(
                f"function {name!r}: an expression may not use unknown functions "
                f"({', '.join(unknown)})"
            )
        indexed = {a.base.label for a in expr.atoms(Indexed)}
        symbols = {s.name for s in expr.free_symbols - indexed if isinstance(s, Symbol)}
        symbols |= {str(s) for s in indexed}
        foreign = sorted(symbols - set(self._parameters.names))
        if foreign:
            raise ValueError(
                f"function {name!r}: {', '.join(foreign)} is not a parameter of the network"
            )

        return expr


#: A value of :class:`Functions`.
FunctionValue = float | Expr | Callable[[], float | Expr]
