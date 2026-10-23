# ABOUTME: NetworkFunctions: a network's .jfunc functions, used by calling them with
# ABOUTME: arguments (call) or reading their parameters as network inputs (value)

"""The :class:`NetworkFunctions` table of a :class:`~jaff.Network` (``net.functions``)."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Iterator, Sequence, cast

from sympy import Basic, Expr

from ...errors import ParserError
from ..parsers import AuxiliaryFunctionParser

if TYPE_CHECKING:
    from ..parsers.auxiliary_func import AuxiliaryFunctionProps
    from .network import Network


class NetworkFunctions:
    """The functions of a network's ``.jfunc`` file, held as ``net.functions``.

    A function is only ever used by calling it, in one of two ways:

    - :meth:`call` with arguments, as a network-file rate or global variable
      does (``fselfH2(1e21*n_CO, 1e5)``): each parameter is a placeholder for
      its argument;
    - :meth:`value` without arguments, as the network does when it reads a
      function itself (``chemRate<N>``, ``deltaE<N>``, ``heatingCoolingRate``,
      ...): each parameter names the network quantity it stands for, so
      ``chemRate3(tgas, n_e)`` uses the gas temperature and the electron
      density.

    Bodies are standardized (``n_X``, ``ntot``, ... become ``nden``
    expressions; see :meth:`NetworkSymbols.standardize`) except for their own
    parameters, on first use, so a function the network never uses may name a
    species it does not have.  Calls between functions were already inlined by
    the ``.jfunc`` parser.  ``rc_<N>`` references are left for
    :meth:`NetworkSymbols.resolve_rate_references`.

    Names are case-insensitive.

    Examples
    --------
    >>> net = Network("GOW")
    >>> "heatingCoolingRate" in net.functions
    True
    """

    def __init__(self, net: Network, jfunc_path: Path | None):
        """Load the functions of *jfunc_path*.

        Parameters
        ----------
        net : Network
            The network whose symbols standardize the functions.
        jfunc_path : Path or None
            The ``.jfunc`` file (``net.spec.jfunc_path``); ``None`` for none.

        Raises
        ------
        ParserError
            If the ``.jfunc`` file is malformed.
        """
        self._net: Network = net
        self._functions: dict[str, AuxiliaryFunctionProps] = self._load(jfunc_path)
        self._bodies: dict[str, Expr] = {}
        self._values: dict[str, Expr] = {}

    @staticmethod
    def _load(jfunc_path: Path | None) -> dict[str, AuxiliaryFunctionProps]:
        """Parse *jfunc_path* into a dict keyed by lower-cased name."""
        if jfunc_path is None:
            return {}

        with AuxiliaryFunctionParser(jfunc_path) as afp:
            return afp.get_dict()

    def __contains__(self, name: object) -> bool:
        """Whether the ``.jfunc`` file defines a function called *name*."""
        return isinstance(name, str) and name.lower() in self._functions

    def __iter__(self) -> Iterator[str]:
        """The function names, lower-cased, in file order."""
        return iter(self._functions)

    def call(self, name: str, args: Sequence[Basic]) -> Expr:
        """The function's body with each parameter replaced by its argument.

        Parameters
        ----------
        name : str
            The function name.
        args : Sequence[Basic]
            One argument per parameter, in order.

        Returns
        -------
        Expr
            The body with the arguments substituted.

        Raises
        ------
        KeyError
            If no function is called *name*.
        ParserError
            If the number of arguments does not match the parameters.
        """
        key = self._key(name)
        params = self._functions[key].args
        if len(args) != len(params):
            raise ParserError(
                f"Function '{name}' called with {len(args)} argument(s) "
                f"but expects {len(params)}"
            )

        return cast(Expr, self._body(key).xreplace(dict(zip(params, args))))

    def value(self, name: str) -> Expr:
        """The function's body with its parameters read as network inputs.

        A parameter ``n_e`` becomes the electron density, ``tgas`` stays the gas
        temperature, and so on.  Cached.

        Parameters
        ----------
        name : str
            The function name.

        Returns
        -------
        Expr
            The standardized body.

        Raises
        ------
        KeyError
            If no function is called *name*.
        """
        key = self._key(name)
        if key not in self._values:
            standardize = self._net.symbols.standardize
            params = self._functions[key].args
            body = self._body(key)
            self._values[key] = cast(
                Expr, body.xreplace({p: standardize(p) for p in params})
            )

        return self._values[key]

    def _key(self, name: str) -> str:
        """The table key of *name*; ``KeyError`` when there is no such function."""
        key = name.lower()
        if key not in self._functions:
            raise KeyError(name)

        return key

    def _body(self, key: str) -> Expr:
        """The body standardized except for its own parameters; cached."""
        if key not in self._bodies:
            function = self._functions[key]
            self._bodies[key] = self._net.symbols.standardize(
                function.expr, keep=function.args
            )

        return self._bodies[key]
