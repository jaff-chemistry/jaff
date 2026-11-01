# ABOUTME: Template functions for $$JAFF CALL: the @jaffgen_function decorator and the
# ABOUTME: lazy loader of jaffgen_template_functions.py (run once, on first use)

from __future__ import annotations

import runpy
from contextvars import ContextVar
from pathlib import Path
from typing import TYPE_CHECKING, Callable, TypeVar, overload

if TYPE_CHECKING:
    from ...core.network.network import Network

FUNCTIONS_FILE = "jaffgen_template_functions.py"
"""Name of the file holding a template's ``$$JAFF CALL`` functions."""

TemplateFunction = Callable[["Network"], str]
"""A template function: takes the network, returns the text to insert."""

F = TypeVar("F", bound=Callable[..., str])

# (name, function) pairs registered while a functions file is being executed;
# None outside a load, so decorating elsewhere (tests, imports) registers nothing.
_collecting: ContextVar[list[tuple[str, Callable[..., str]]] | None] = ContextVar(
    "_collecting", default=None
)


@overload
def jaffgen_function(target: F) -> F: ...


@overload
def jaffgen_function(target: str) -> Callable[[F], F]: ...


def jaffgen_function(target: F | str) -> F | Callable[[F], F]:
    """Register a function for ``$$JAFF CALL`` in ``jaffgen_template_functions.py``.

    The function takes the network and returns the text that replaces the
    ``$$JAFF CALL`` line. It is called by its own name, or by the name given
    to the decorator. The function itself is returned unchanged.

    Parameters
    ----------
    target : callable or str
        The function to register, or the name to register the next one under.

    Returns
    -------
    callable
        The function, or a decorator registering one under *target*.

    Examples
    --------
    In ``jaffgen_template_functions.py``:

    ```python
    from jaff.codegen import jaffgen_function


    @jaffgen_function
    def species_count(net):
        return f"const int n = {len(net.species)};"


    @jaffgen_function("init")
    def _initial_values(net):
        return "\\n".join(f"y[{i}] = 0.0;" for i in range(len(net.species)))
    ```

    and in a template, ``$$JAFF CALL species_count`` or ``$$JAFF CALL init``.
    """
    if isinstance(target, str):
        name = target

        def register(function: F) -> F:
            return _register(name, function)

        return register

    return _register(_name(target), target)


def _name(function: Callable[..., str]) -> str:
    """The function's own name (its ``repr`` for a callable without one)."""
    return getattr(function, "__name__", repr(function))


def _register(name: str, function: F) -> F:
    """Add *function* under *name* to the load in progress, if any."""
    collected = _collecting.get()
    if collected is not None:
        collected.append((name, function))

    return function


class TemplateFunctions:
    """The ``$$JAFF CALL`` functions of one ``jaffgen_template_functions.py``.

    The file is executed the first time a function is looked up, and only
    then: templates without ``$$JAFF CALL`` never run it. One instance serves
    every template of a jaffgen run, so the file runs at most once per run.

    Parameters
    ----------
    path : Path
        The functions file.

    Attributes
    ----------
    path : Path
        The functions file.

    Notes
    -----
    The file is ordinary Python run by jaffgen: render only templates you trust.
    """

    def __init__(self, path: Path) -> None:
        """Remember the file; nothing is executed until a function is needed.

        Parameters
        ----------
        path : Path
            The functions file.
        """
        self.path = Path(path)
        self._functions: dict[str, Callable[..., str]] | None = None

    def names(self) -> list[str]:
        """The registered function names, sorted; loads the file if needed.

        Returns
        -------
        list of str
            Names usable in ``$$JAFF CALL``.

        Raises
        ------
        RuntimeError
            If the file fails to run or registers a name twice.
        """
        return sorted(self._load())

    def get(self, name: str) -> Callable[..., str] | None:
        """The function registered as *name*, or ``None``; loads the file if needed.

        Parameters
        ----------
        name : str
            The name used in ``$$JAFF CALL``.

        Returns
        -------
        callable or None
            The function, or ``None`` when no function has that name.

        Raises
        ------
        RuntimeError
            If the file fails to run or registers a name twice.
        """
        return self._load().get(name)

    def _load(self) -> dict[str, Callable[..., str]]:
        """Run the file once and keep what its decorators registered."""
        if self._functions is not None:
            return self._functions

        collected: list[tuple[str, Callable[..., str]]] = []
        token = _collecting.set(collected)
        try:
            runpy.run_path(str(self.path), run_name="jaffgen_template_functions")
        except Exception as e:
            raise RuntimeError(
                f"{self.path.name} failed to load: {type(e).__name__}: {e}"
            ) from e
        finally:
            _collecting.reset(token)

        functions: dict[str, Callable[..., str]] = {}
        for name, function in collected:
            if name in functions:
                raise RuntimeError(
                    f"{self.path.name} registers {name} more than once "
                    f"({_name(functions[name])} and {_name(function)})"
                )

            functions[name] = function
        self._functions = functions

        return functions
