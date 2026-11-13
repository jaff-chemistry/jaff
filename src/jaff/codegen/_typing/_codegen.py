# ABOUTME: Structured values Codegen methods return: ExtrasDict (CSE temporaries etc.),
# ABOUTME: IndexedReturn (extras + expressions) and CseGroup (shared temporaries + parts)
"""Definitions of the structured values Codegen methods return."""

from dataclasses import dataclass
from typing import TypedDict

from ...types import IndexedList

ExtrasDict = TypedDict(
    "ExtrasDict",
    {
        "cse": IndexedList,
    },
)
"""Extra data returned alongside the main code-generation expressions.

The ``extras`` entry of an :class:`IndexedReturn`, a plain ``dict`` with the
keys below.

Attributes
----------
cse : IndexedList
    Common-subexpression temporaries found by :func:`sympy.cse`, as
    :class:`~jaff.types.IndexedValue` items ``(indices, expr_str)`` where
    ``indices`` is ``[i]`` for the ``i``-th temporary.  Empty when CSE is off.
"""


IndexedReturn = TypedDict(
    "IndexedReturn",
    {
        "extras": ExtrasDict,
        "expressions": IndexedList,
    },
)
"""Structured return type for indexed code-generation methods.

Returned, as a plain ``dict`` with the keys below, by the ``get_indexed_*``
methods of :class:`~jaff.codegen.Codegen` that support CSE
(:meth:`~jaff.codegen.Codegen.get_indexed_rates`,
:meth:`~jaff.codegen.Codegen.get_indexed_odes`,
:meth:`~jaff.codegen.Codegen.get_indexed_rhs`,
:meth:`~jaff.codegen.Codegen.get_indexed_radodes` and
:meth:`~jaff.codegen.Codegen.get_indexed_jacobian`).  The expressions refer to
the CSE temporaries by name, so both lists are needed to emit the code.

Attributes
----------
extras : ExtrasDict
    Supplementary data, including CSE temporaries under ``extras["cse"]``.
expressions : IndexedList
    Primary expressions as :class:`~jaff.types.IndexedValue` items
    ``(indices, expr_str)``, one per reaction, species, or Jacobian element
    (``indices`` is ``[i]``, or ``[i, j]`` for a Jacobian element).

Examples
--------
>>> from jaff import Network
>>> from jaff.codegen import Codegen
>>> out = Codegen(network=Network("COthin"), lang="python").get_indexed_rates()
>>> sorted(out)
['expressions', 'extras']
>>> print(out["extras"]["cse"][0])
[0] -> (tgas < 2.73)
>>> out["expressions"][5].indices, out["expressions"][5].value
([5], '3.00351426766844e-16*tgas**0.25')
"""


@dataclass(frozen=True)
class CseGroup:
    """Expressions of several kinds printed with one shared set of CSE temporaries.

    Returned by :meth:`~jaff.codegen.Codegen.cse_group`.  Emit :attr:`temporaries`
    first, then each member's expressions in the same scope: they read the
    temporaries by name.

    Attributes
    ----------
    temporaries : IndexedList
        The shared temporaries as :class:`~jaff.types.IndexedValue` items
        ``([i], expr_str)``, in evaluation order (each after the ones it reads).
    expressions : tuple of IndexedList
        One list per member, in the order the members were given, of
        ``(indices, expr_str)`` items as the member's ``get_indexed_*`` method
        would index them (``[i]``, or ``[i, j]`` for a Jacobian element).
    """

    temporaries: IndexedList
    expressions: tuple[IndexedList, ...]
