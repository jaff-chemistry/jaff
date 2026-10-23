# ABOUTME: Shared helpers: file-extension groups, mass table loading, Fortran literal
# ABOUTME: conversion, symbolic dependency resolution and .jaff path detection.
"""
Miscellaneous helper utilities shared across JAFF subsystems.

This module provides:

* Extension-list constants for recognizing source/data file types.
* :func:`load_mass_dict` -- load atomic masses from the embedded SQLite DB.
* :func:`f90_convert` -- coerce Fortran-style numeric literals to Python/C.
* :func:`resolve_symbolic_dependencies` -- topological substitution of a
  symbol dependency graph.
* :func:`float_piecewise_branches` -- float the integer branch values of
  every ``Piecewise`` in an expression.
* :func:`resolve_dependencies` -- single-expression dependency resolution
  with special handling for built-in JAFF functions (``merge``, ``log10``).
* :func:`is_jaff_file` -- quick path-based check for ``.jaff`` files.
* :func:`load_module_from_path` and :func:`import_subpackages` -- dynamic
  import helpers (currently not used elsewhere in JAFF).
"""

from __future__ import annotations

import importlib.util
import re
from pathlib import Path
from typing import TYPE_CHECKING

from sympy import Basic, Expr, Float, Integer, Piecewise, srepr
from sympy.core.function import AppliedUndef

from ..errors import ParserError

if TYPE_CHECKING:
    from ..core._typing import ElementProps
    from ..core.parsers.auxiliary_func import AuxiliaryFunctionProps

# ---------------------------------------------------------------------------
# File-extension groups used by parsers and code generators
# ---------------------------------------------------------------------------

HDF_EXTENSIONS = [".hdf5", ".hdf", ".h5"]
CSV_EXTENSIONS = [".csv", ".txt", ".dat"]
C_EXTENSIONS = [".c", ".h"]
CPP_EXTENSIONS = [".cpp", ".hpp", ".cc", ".cxx", ".hh", ".hxx"]
FORTRAN_EXTENSIONS = [".f", ".f90", ".f95", ".f03", ".f08", ".for", ".f77"]
RUST_EXTENSIONS = [".rs"]

# Matches Fortran double-precision exponent notation, e.g. ``1.0d-3`` → ``1.0e-3``
F90_PATTERN = re.compile(r"([0-9_.])d([0-9_+-])")


def load_mass_dict() -> dict:
    """
    Load atomic mass data from the JAFF embedded SQLite database.

    Opens the bundled ``atomic_masses`` table via :class:`~jaff.drivers.sqlite.JaffDb`
    and returns a dict keyed by element symbol.

    Returns
    -------
    dict[str, ElementProps]
        Mapping from element symbol (e.g. ``"H"``) to an
        :class:`~jaff.core._typing.ElementProps` dict with keys
        ``name``, ``mass``, ``atomic_mass``, ``protons``, ``neutrons``,
        and ``electrons``.
    """
    from ..drivers.sqlite import JaffDb

    with JaffDb() as jdb:
        rows = jdb.table("atomic_masses").all_rows()

    mass_dict: dict[str, ElementProps] = {}
    for row in rows:
        mass_dict[row["element"]] = {
            "name": row["name"],
            "mass": row["mass"],
            "atomic_mass": row["atomic_mass"],
            "protons": row["protons"],
            "neutrons": row["neutrons"],
            "electrons": row["electrons"],
        }

    return mass_dict


def f90_convert(expr: str) -> str:
    """
    Convert a Fortran-style numeric expression string to Python/C syntax.

    Applies three transformations:

    1. Replaces ``dexp(`` with ``exp(`` (Fortran intrinsic → standard name).
    2. Replaces the Fortran double-precision exponent marker ``d`` with ``e``
       in numeric literals (e.g. ``1.0d-3`` → ``1.0e-3``), using
       :data:`F90_PATTERN`.
    3. Removes Fortran array subscript notation ``(:)`` which has no meaning
       in a scalar expression context.

    Parameters
    ----------
    expr : str
        The expression string to convert.

    Returns
    -------
    str
        The converted expression string.
    """
    expr = expr.replace("dexp(", "exp(").replace("(:)", "")

    return F90_PATTERN.sub(r"\1e\2", expr)


def resolve_symbolic_dependencies(
    dep_map: dict[str, Basic],
    external_refs: dict[str, Basic] | None = None,
    fname: Path | str = "",
) -> dict[str, Basic]:
    """
    Topologically resolve a symbol dependency graph into fully-substituted expressions.

    Given a mapping from symbol names to SymPy expressions that may reference
    other symbols in the same mapping, this function performs a depth-first
    traversal and returns a new mapping where every expression has all
    internal dependencies substituted away.  Symbols present in
    *external_refs* but not in *dep_map* are treated as resolved leaf values.

    Parameters
    ----------
    dep_map : dict[str, sympy.Basic]
        Mapping from symbol name to its defining expression.  Expressions may
        contain references to other keys in this dict or to keys in
        *external_refs*.
    external_refs : dict[str, sympy.Basic] or None, optional
        Additional pre-resolved symbol definitions.  When a name appears in
        both *external_refs* and *dep_map*, the local *dep_map* definition
        takes precedence.
    fname : pathlib.Path or str, optional
        Source file path shown in error messages (default ``""``).

    Returns
    -------
    dict[str, sympy.Basic]
        A new dict with the same keys as *dep_map* but with all internal
        symbol references expanded to their final values.

    Raises
    ------
    ParserError
        If a cyclic dependency is detected (e.g. ``A`` depends on ``B``
        which depends on ``A``).

    Notes
    -----
    Both free :class:`sympy.Symbol` objects and undefined applied functions
    (:class:`sympy.core.function.AppliedUndef`) are considered dependencies
    and resolved recursively.
    """
    resolved = {}
    visiting = set()
    external = external_refs or {}

    def dfs(sym: str):
        """Recursively resolve *sym*, memoizing results in *resolved*."""
        if sym in resolved:
            return resolved[sym]

        if sym in external and sym not in dep_map:
            return external[sym]

        if sym in visiting:
            raise ParserError(f"Cyclic dependency found for {sym}", fname=fname)

        visiting.add(sym)
        expr = dep_map[sym]

        replacements = {}

        # Resolve free symbol references (e.g. plain ``tgas`` in the expression)
        for s in expr.free_symbols:
            s_name = str(s)
            if s_name in dep_map or s_name in external:
                replacements[s] = dfs(s_name)

        # Resolve applied undefined function references (e.g. ``f(x)`` where
        # ``f`` is a user-defined auxiliary function)
        for f in expr.atoms(AppliedUndef):
            f_str = str(f)

            if f_str in dep_map or f_str in external:
                replacements[f] = dfs(f_str)
            else:
                f_name = f.func.__name__
                if f_name in dep_map or f_name in external:
                    replacements[f] = dfs(f_name)

        new_expr = expr.xreplace(replacements)

        visiting.remove(sym)
        resolved[sym] = new_expr

        return new_expr

    return {sym: dfs(sym) for sym in dep_map}


def float_piecewise_branches(expr: Basic) -> Basic:
    """Cast bare integer *value* branches of every ``Piecewise`` to ``Float``.

    A SymPy ``Piecewise`` prints to Fortran as ``merge(tsource, fsource,
    cond)``, and Fortran's ``merge`` requires both value branches to share the
    same type *and* kind.  When one branch is a real expression (``real*8``)
    and another is a bare integer literal -- e.g. the ``0`` fallback of a
    temperature-guarded rate, or the ``10``/``18`` regime constants in the H2
    cosmic-ray heating term -- gfortran rejects it with *"fsource argument of
    merge must be same type and kind as tsource"*.

    This walks *expr* and, for every ``Piecewise`` (including nested ones),
    replaces any branch **value** that is a bare :class:`sympy.Integer` with
    the equivalent :class:`sympy.Float`, so it prints as ``0.0d0`` / ``10.0d0``
    and matches the real-kind sibling branch.  Conditions are left untouched
    (integer comparisons are fine), as are integer literals used as exponents
    or array dimensions elsewhere in the expression -- only Piecewise branch
    values are affected.

    Parameters
    ----------
    expr : sympy.Basic
        Any expression that may contain one or more ``Piecewise`` nodes.

    Returns
    -------
    sympy.Basic
        The expression with integer Piecewise branch values floated.
    """

    def _fix(pw: Piecewise) -> Piecewise:
        new_args = [
            (Float(value) if isinstance(value, Integer) else value, cond)
            for value, cond in pw.args
        ]
        return Piecewise(*new_args, evaluate=False)

    return expr.replace(lambda node: isinstance(node, Piecewise), _fix)


def resolve_dependencies(
    expr: Basic,
    subs_dict: dict[Basic, Basic] | None = None,
    aux_funcs: dict[str, AuxiliaryFunctionProps] | None = None,
) -> Expr:
    """
    Resolve undefined SymPy function calls inside a single expression.

    Walks all :class:`~sympy.core.function.AppliedUndef` atoms in *expr* and
    replaces each with a concrete SymPy expression, drawing from three sources:

    * **Built-in shims** -- ``merge(a, b, cond)`` → ``Piecewise((a, cond), (b, True))``
      and ``log10(x)`` → ``log(x)/log(10)``.
    * **aux_funcs** -- user-defined auxiliary functions provided by the
      parser; arguments are recursively resolved before substitution.
    * **subs_dict** -- direct substitution table, matched case-insensitively
      on the call's string form.

    A ``subs_dict`` match is recorded first but is overwritten when the
    function is also a built-in or an *aux_funcs* entry, so in effect the
    built-ins win, then *aux_funcs*, then *subs_dict*.

    Parameters
    ----------
    expr : sympy.Basic
        The expression in which to resolve dependencies.
    subs_dict : dict[sympy.Basic, sympy.Basic] or None, optional
        Pre-populated substitution table.  Modified in-place as new
        substitutions are discovered.  Defaults to an empty dict.
    aux_funcs : dict[str, AuxiliaryFunctionProps] or None, optional
        Auxiliary function definitions keyed by lowercase function name.
        Each entry supplies ``expr`` (the body expression) and ``args``
        (the ordered parameter list).  Defaults to an empty dict.

    Returns
    -------
    sympy.Expr
        The expression with all recognized undefined function calls replaced.
        Unrecognized functions are left unchanged.
    """
    if subs_dict is None:
        subs_dict = {}
    if aux_funcs is None:
        aux_funcs = {}

    for f in sorted(expr.atoms(AppliedUndef), key=srepr):
        name = f.func.__name__.lower()

        # Check the direct substitution table first (case-insensitive key match)
        for k, v in subs_dict.items():
            if str(k).lower() == str(f).lower():
                subs_dict[f] = v
                break

        # Built-in: ``merge(a, b, cond)`` → Piecewise
        if name == "merge":
            subs_dict[f] = float_piecewise_branches(
                Piecewise(
                    (f.args[0], f.args[2]),
                    (f.args[1], True),
                )
            )

        # Built-in: ``log10(x)`` → sympy natural-log form
        elif name == "log10":
            import sympy

            subs_dict[f] = sympy.log(f.args[0]) / sympy.log(10)

        # User-defined auxiliary function: inline-expand with argument mapping
        elif name in aux_funcs:
            func_def = aux_funcs[name].expr
            func_args = aux_funcs[name].args

            arg_map = dict(
                zip(
                    func_args,
                    [resolve_dependencies(arg, subs_dict, aux_funcs) for arg in f.args],
                )
            )

            subs_dict[f] = resolve_dependencies(
                func_def.xreplace(arg_map), subs_dict, aux_funcs
            )

    return expr.xreplace(subs_dict)


def is_jaff_file(file: Path) -> bool:
    """
    Return ``True`` if *file* has a ``.jaff`` or ``.jaff.gz`` extension.

    The check is case-insensitive.

    Parameters
    ----------
    file : pathlib.Path
        The path to test.

    Returns
    -------
    bool
        ``True`` iff the path ends with ``.jaff`` or ``.jaff.gz``.
    """
    return file.suffix.lower() == ".jaff" or [fn.lower() for fn in file.suffixes] == [
        ".jaff",
        ".gz",
    ]


def load_module_from_path(path: str | Path, module_name: str):
    mpath: Path = path  # type: ignore
    if isinstance(path, str):
        mpath = Path(path).resolve()

    spec = importlib.util.spec_from_file_location(module_name, mpath)

    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load module from {mpath}")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    return module


def import_subpackages(package_name: str) -> None:
    """Import every non-private subpackage of *package_name*.

    Iterates the package directory and imports each non-private subpackage,
    which triggers any registration side-effects (e.g. the ``@register``
    decorator in ``_formats``).

    Uses :func:`Path.iterdir` (not :func:`pkgutil.walk_packages`) to avoid
    eagerly importing private subpackages while walking — a pitfall of
    :func:`pkgutil.walk_packages` which calls :func:`__import__` internally
    for every package it encounters.

    Parameters
    ----------
    package_name : str
        Fully-qualified package name
        (e.g. ``"jaff.core.parsers.network._formats"``).
    """
    from importlib import import_module

    pkg = import_module(package_name)
    pkg_path = Path(pkg.__path__[0])
    for entry in pkg_path.iterdir():
        if entry.is_dir() and not entry.name.startswith("_"):
            import_module(f"{package_name}.{entry.name}")
