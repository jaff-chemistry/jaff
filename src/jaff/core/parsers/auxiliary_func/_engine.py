# ABOUTME: AuxiliaryFunctionParser: parses .jfunc auxiliary function files
# ABOUTME: Handles @var globals and @function ... @return blocks

"""Parser for JAFF auxiliary function files (.jfunc).

A ``.jfunc`` file may define global symbolic constants and named functions
that augment a reaction network.  The format uses two directives:

``@var name = expression``
    Define a global symbolic constant.  The expression is parsed by SymPy.

``@function name(arg1, arg2, ...)``
    Begin a function block.  Lines inside the block are local variable
    assignments; the block ends with a ``@return <expression>`` statement.
    Comments (``# arg doc``) inside a block attach documentation to the
    named argument.

Continuation lines end with ``\\``.  Inline comments start with ``#``.

The parsed results are stored as :class:`AuxiliaryFunctionProps` objects keyed by
function name, holding their symbolic definitions, argument lists, and argument comments.  Global
variables are resolved into the function bodies so that callers receive
fully-substituted SymPy expressions.

Cross-function dependencies are resolved with a DFS to support one function
calling another, while detecting circular references.

Calls to the unit conversions of :mod:`jaff.physics.units` (``nm_to_eV(x)``,
``pc_to_cm(r)`` ...; any letter case) and the built-in CGS constants
(:data:`CONSTANTS`: ``kB``, ``c_light``, ``amu`` ...; exact case) are replaced by
their value once the file is parsed. A ``@var``, ``@function``, function
argument or local variable named like one of them is an error.
"""

from functools import cached_property
from pathlib import Path
from tokenize import TokenError
from typing import Any, Callable

import sympy as sp
from sympy.core.function import AppliedUndef

from ....common import float_piecewise_branches, resolve_symbolic_dependencies
from ....errors import ParserError
from ....physics import constants as _c
from ....physics import units
from .function_props import AuxiliaryFunctionProps

#: The unit conversions of :mod:`jaff.physics.units`, by lower-cased name: a
#: ``.jfunc`` file may call them in any case.
CONVERTERS: dict[str, Callable[[Any], Any]] = {
    name.lower(): getattr(units, name) for name in units.__all__
}

#: Built-in constants a ``.jfunc`` file may use by name (exact case), as CGS
#: values from :mod:`jaff.physics.constants`.
CONSTANTS: dict[str, float] = {
    "c_light": float(_c.c.cgs.value),
    "h_planck": float(_c.h.cgs.value),
    "hbar": float(_c.hbar.cgs.value),
    "kB": float(_c.k_B.cgs.value),
    "G_grav": float(_c.G.cgs.value),
    "e_charge": float(_c.e.esu.value),
    "m_e": float(_c.m_e.cgs.value),
    "m_p": float(_c.m_p.cgs.value),
    "m_n": float(_c.m_n.cgs.value),
    "amu": float(_c.u.cgs.value),
    "N_A": float(_c.N_A.value),
    "sigma_SB": float(_c.sigma_sb.cgs.value),
    "a_rad": float(_c.a_rad.cgs.value),
    "sigma_T": float(_c.sigma_T.cgs.value),
    "a_bohr": float(_c.a0.cgs.value),
    "alpha_fs": float(_c.alpha.value),
    "Ry": float(_c.Ry.to("erg").value),
    "pc": float(_c.pc.cgs.value),
    "au": float(_c.au.cgs.value),
    "M_sun": float(_c.M_sun.cgs.value),
    "L_sun": float(_c.L_sun.cgs.value),
    "R_sun": float(_c.R_sun.cgs.value),
}


class AuxiliaryFunctionParser:
    """Parser for JAFF ``.jfunc`` auxiliary function files.

    Parses the file on construction, resolving all global constants and
    cross-function dependencies.  Intended to be used as a context manager so
    that internal state is cleared on exit.

    Attributes
    ----------
    file : Path
        Absolute path to the parsed ``.jfunc`` file.
    globals : dict[str, sp.Basic]
        Global symbolic constants defined with ``@var``, keyed by name.
    func_dict : dict[str, AuxiliaryFunctionProps]
        Parsed functions, keyed by lower-cased function name.  Each entry
        is an :class:`AuxiliaryFunctionProps` with ``expr`` (SymPy
        expression), ``args`` (list of SymPy symbols), and ``arg_comments``
        (dict of argument doc strings).

    Examples
    --------
    >>> with AuxiliaryFunctionParser("network.jfunc") as afp:
    ...     funcs = afp.get_dict()
    """

    def __init__(self, file: Path | str):
        """Parse a ``.jfunc`` file and resolve all symbolic dependencies.

        Parameters
        ----------
        file : Path | str
            Path to the ``.jfunc`` file.

        Raises
        ------
        ValueError
            If *file* is neither a ``Path`` nor a ``str``.
        FileNotFoundError
            If the resolved path does not exist.
        ParserError
            On syntax errors in the file content.
        """
        if not isinstance(file, (Path, str)):
            raise ValueError(
                f"Invalid file type detected: {file}\n"
                f"file must be a pathlib.Path object or a str"
            )
        if isinstance(file, str):
            file = Path(file)

        file = file.resolve()
        if not file.exists():
            raise FileNotFoundError(file)

        self.file: Path = file
        self.og_line: str = ""  # raw line from file (before continuation merge)
        self.line: str = ""  # processed line ready for parsing
        self.cline: str = ""  # accumulator for continuation lines
        self.nline: int = 0  # current 1-based line number (for error messages)
        self.globals: dict[str, sp.Basic] = {}
        self.globals_parsed: bool = False  # True once global vars are resolved
        self.func_dict: dict[str, AuxiliaryFunctionProps] = {}
        self.scope: str = "global"  # "global" | "function"
        self.current_func: str = ""  # name of the function block being parsed

        self.__parse_file()
        self.__resolve_func_deps()
        self.__apply_builtins()

    def __repr__(self):
        """Return canonical string representation of this parser instance.

        Returns
        -------
        str
            String of form ``"AuxiliaryFunctionParserObject: <file>"``.
        """
        return f"AuxiliaryFunctionParserObject: {self.file}"

    def __str__(self):
        """Return human-readable string representation of this parser.

        Returns
        -------
        str
            String of form ``"AuxiliaryFunctionParserObject: <file>"``.
        """
        return f"AuxiliaryFunctionParserObject: {self.file}"

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Clear internal state on context manager exit."""
        self.globals.clear()
        self.__dict__.clear()
        self.og_line = self.line = self.cline = ""

    def __check_name(self, name: str, kind: str) -> None:
        """Raise if *name* is a unit conversion (any case) or a built-in constant.

        Raises
        ------
        ParserError
            If *name* matches a :mod:`jaff.physics.units` function or a key of
            :data:`CONSTANTS`.
        """
        name = name.strip()
        if name.lower() in CONVERTERS:
            taken = "a unit conversion of jaff.physics.units"
        elif name in CONSTANTS:
            taken = "a built-in constant"
        else:
            return

        raise ParserError(
            f"{kind} '{name}' has the name of {taken}; choose another name",
            self.og_line,
            self.nline,
            self.file,
        )

    def __apply_builtins(self) -> None:
        """Replace unit-conversion calls and built-in constants by their values.

        Runs after the file's own functions are inlined, so the only undefined
        calls left are converters and functions the network supplies (such as
        interpolators); the latter are left as they are. Conversions are
        replaced innermost first, then constants.

        Raises
        ------
        ParserError
            If a conversion is called with other than one argument.
        """

        def is_converter(e: sp.Basic) -> bool:
            return isinstance(e, AppliedUndef) and e.func.__name__.lower() in CONVERTERS

        def convert(e: Any) -> sp.Basic:
            name = e.func.__name__.lower()
            if len(e.args) != 1:
                raise ParserError(
                    f"Unit conversion '{name}' takes one argument, got {len(e.args)}",
                    fname=self.file,
                )

            return sp.sympify(CONVERTERS[name](e.args[0]))

        for func in self.func_dict.values():
            expr = func.expr.replace(is_converter, convert)
            values = {
                symbol: sp.Float(CONSTANTS[symbol.name])
                for symbol in expr.free_symbols
                if isinstance(symbol, sp.Symbol) and symbol.name in CONSTANTS
            }
            func.expr = expr.xreplace(values)

    def get_dict(self) -> dict[str, AuxiliaryFunctionProps]:
        """Return the parsed function dictionary.

        Returns
        -------
        dict[str, AuxiliaryFunctionProps]
            Maps lower-cased function names to their symbolic definitions,
            argument lists, and argument documentation strings.
        """
        return self.func_dict

    def __parse_file(self) -> None:
        """Iterate over the file, handling line continuations before dispatching."""
        with open(self.file, "r", encoding="utf-8") as f:
            for nline, line in enumerate(f, start=1):
                self.nline = nline
                self.og_line = line.strip()

                if len(self.og_line) == 0:
                    continue

                # A trailing backslash joins this line with the next.
                if self.og_line.endswith("\\"):
                    self.cline += f"{self.og_line[:-1]} "
                    continue

                self.cline += line
                self.line = self.cline.strip()
                self.cline = ""
                self.__parse_line()

        # A trailing backslash with no following line leaves an unterminated
        # continuation buffer; the pending directive would be silently dropped.
        if self.cline.strip():
            raise ParserError(
                "Unterminated line continuation at end of file",
                self.og_line,
                self.nline,
                self.file,
            )

        # A function block that never saw a ``@return`` leaves the scope open;
        # its ``expr`` would stay 0.0 with unresolved locals.
        if self.scope == "function":
            raise ParserError(
                f"Unterminated function block '{self.current_func}': "
                f"missing @return statement",
                self.og_line,
                self.nline,
                self.file,
            )

    def __parse_line(self) -> None:
        """Dispatch the current line to the appropriate handler.

        When inside a function block, delegates to :meth:`__parse_function`.
        Lines starting with ``@`` are dispatched to :meth:`__parse_token`.
        Any other line outside a function block (comments included) is
        ignored.
        """
        if self.scope == "function":
            self.__parse_function()
            return

        if self.line.startswith("@"):
            self.parse = True
            self.__parse_token()

            return

        if self.line.startswith("#"):
            return

    def __parse_token(self):
        """Extract the ``@token`` keyword and dispatch to its registered handler.

        Raises
        ------
        ParserError
            If the token is not in :attr:`tokens`.
        """
        line = self.line[1:].strip()
        split_line: list[str] = line.split(maxsplit=1)
        token, segment = split_line[0], split_line[1]

        if token not in self.tokens:
            raise ParserError(
                f"Invalid token detected: @{token}", self.og_line, self.nline, self.file
            )

        self.tokens[token](self.__strip_trailing_comment(segment.strip())[0])

    def __set_scope(self, scope: str):
        """Set the current parser scope.

        Parameters
        ----------
        scope : str
            New scope string — either ``"global"`` or ``"function"``.
        """
        self.scope = scope

    def __handle_var(self, segment: str) -> None:
        """Parse a ``@var name = expression`` directive and store in globals.

        Parameters
        ----------
        segment : str
            Text after the ``@var`` keyword (e.g. ``"alpha = 1.7e-9"``).

        Raises
        ------
        ParserError
            If the segment does not contain exactly one ``=`` separator, or if
            the expression cannot be parsed by SymPy.
        """
        tokens: list[str] = segment.split("=", maxsplit=1)
        if len(tokens) != 2:
            raise ParserError(
                "Invalid variable encountered", self.og_line, self.nline, self.file
            )
        var, expr = tokens[0].strip(), tokens[1].strip()
        self.__check_name(var, "Variable")

        try:
            self.globals[str(var)] = sp.parse_expr(expr)
        except (SyntaxError, TokenError, TypeError, sp.SympifyError) as e:
            raise ParserError(
                f"Invalid expression encountered: {e}",
                self.og_line,
                self.nline,
                self.file,
            )

    def __handle_function_decleration(self, segment: str):
        """Parse a ``@function name(args)`` declaration and open a new function scope.

        Resolves all pending global variables before entering the function scope
        so they are available inside the function body.

        Parameters
        ----------
        segment : str
            Text after the ``@function`` keyword
            (e.g. ``"alpha_H(T, nH)"``).

        Raises
        ------
        ParserError
            If the declaration is malformed (missing parentheses).  An empty
            argument list (``"f()"``) is accepted.
        """
        self.__set_scope("function")
        if not self.globals_parsed:
            self.globals_parsed = True
            self.globals = resolve_symbolic_dependencies(
                dep_map=self.globals, fname=self.file
            )

        parts = segment.split("(", maxsplit=1)
        if len(parts) != 2:
            raise ParserError(
                "Invalid function decleration", self.og_line, self.nline, self.file
            )
        name, rest = parts[0].lower(), parts[1]
        self.__check_name(parts[0].strip(), "Function")
        rest = rest.split(")", maxsplit=1)

        if len(rest) != 2:
            raise ParserError(
                "Invalid function decleration", self.og_line, self.nline, self.file
            )
        names = [arg.strip() for arg in rest[0].split(",") if len(arg.strip()) > 0]
        for arg in names:
            self.__check_name(arg, "Argument")

        args = [sp.parse_expr(arg) for arg in names]

        self.current_func = name
        self.func_dict[name] = AuxiliaryFunctionProps(expr=sp.Float(0.0), args=args)

    def __parse_function(self):
        """Dispatch a line inside a function body to the correct sub-handler.

        Handles three line types:
        - ``@return <expr>`` — delegates to :meth:`__handle_function_return`.
        - ``# arg doc``      — delegates to :meth:`__handle_function_comment`.
        - ``var = expr``     — stores a local variable assignment.

        Raises
        ------
        ParserError
            If the line is neither a comment, a return, nor a valid assignment.
        """
        line = self.line
        if line.startswith("@return"):
            self.__handle_function_return()
            return

        elif line.startswith("#"):
            self.__handle_function_comment()
            return

        line, _ = self.__strip_trailing_comment(line)
        splitline = line.split("=", maxsplit=1)

        if len(splitline) != 2:
            raise ParserError(
                "Invalid line detected", self.og_line, self.nline, self.file
            )

        var, expr = splitline[0].strip(), splitline[1].strip()
        self.__check_name(var, "Variable")
        try:
            self.func_dict[self.current_func]._locals[str(var)] = sp.parse_expr(expr)
        except (SyntaxError, TokenError, TypeError, sp.SympifyError) as e:
            raise ParserError(
                f"Invalid expression encountered: {e}",
                self.og_line,
                self.nline,
                self.file,
            )

    def __handle_function_comment(self) -> None:
        """Attach an argument documentation comment to the current function.

        A comment of the form ``# argname doc text`` inside a function block
        is stored as ``func_dict[name].arg_comments[argname] = doc text``
        when *argname* is a declared argument of the current function.
        """
        line = self.line[1:].strip()
        split_line = line.split(maxsplit=1)

        if len(split_line) < 2:
            return

        arg, comment = split_line[0], split_line[1]
        func = self.func_dict[self.current_func]
        if sp.Symbol(arg) in func.args:
            func.arg_comments[arg] = comment.strip()

    def __handle_function_return(self) -> None:
        """Parse the ``@return`` statement, resolve locals, and close the function scope.

        Resolves local variable dependencies, substitutes them into the return
        expression, and stores the fully-expanded SymPy expression in
        ``func_dict[name].expr`` and clears its ``_locals``.  Resets the scope
        to ``"global"``
        afterwards.

        Raises
        ------
        ParserError
            If the return expression cannot be parsed by SymPy.
        """
        line = self.line
        line, _ = self.__strip_trailing_comment(line)
        func = self.func_dict[self.current_func]
        func._locals = resolve_symbolic_dependencies(
            func._locals,
            external_refs=self.globals,
            fname=self.file,
        )

        try:
            funcdef = sp.parse_expr(
                line.removeprefix("@return").strip(),
                local_dict={
                    **self.globals,
                    **func._locals,
                },
            )

            func_cache = {}
            funcdef = funcdef.xreplace(
                {
                    f: func_cache.setdefault(f.name.lower(), sp.Function(f.name.lower()))(
                        *f.args
                    )
                    for f in funcdef.atoms(AppliedUndef)
                }
            )
        except (SyntaxError, TokenError, TypeError, sp.SympifyError) as e:
            raise ParserError(
                f"Invalid expression encountered: {e}",
                self.og_line,
                self.nline,
                self.file,
            )

        funcdef = float_piecewise_branches(funcdef)

        func.expr = funcdef
        # Locals are already substituted into the body; drop them.
        func._locals = {}
        self.current_func = ""
        self.__set_scope("global")

    def __resolve_func_deps(self):
        """Inline all inter-function calls using a DFS over the call graph.

        Each function whose body calls another ``.jfunc``-defined function
        has that call replaced by the callee's body with arguments
        substituted.  This produces fully expanded, self-contained SymPy
        expressions for each function.

        Raises
        ------
        ParserError
            If a circular dependency is detected between functions.
        """
        resolved_defs = {}
        visiting = set()

        def dfs_resolve_func(name: str):
            if name in resolved_defs:
                return resolved_defs[name]

            if name in visiting:
                raise ParserError(
                    f"Circular function dependency: {name}", fname=self.file
                )

            visiting.add(name)

            expr = self.func_dict[name].expr

            while True:
                calls = sorted(
                    (
                        f
                        for f in expr.atoms(AppliedUndef)
                        if f.func.__name__.lower() in self.func_dict
                    ),
                    key=sp.srepr,
                )
                if not calls:
                    break

                repl = {}
                for func in calls:
                    func_name = func.func.__name__.lower()
                    nested_f_def = dfs_resolve_func(func_name)
                    nested_f_args = self.func_dict[func_name].args

                    if len(func.args) != len(nested_f_args):
                        raise ParserError(
                            f"Function '{func_name}' called with "
                            f"{len(func.args)} argument(s) but expects "
                            f"{len(nested_f_args)}",
                            fname=self.file,
                        )

                    arg_map = dict(zip(nested_f_args, func.args))
                    repl[func] = nested_f_def.subs(arg_map, simultaneous=True)

                expr = expr.xreplace(repl)

            visiting.remove(name)
            resolved_defs[name] = expr
            return expr

        for func_name in self.func_dict:
            self.func_dict[func_name].expr = dfs_resolve_func(func_name)

    @staticmethod
    def __strip_trailing_comment(line: str) -> tuple[str, str]:
        """Split *line* at the first ``#`` and return ``(code, comment)``.

        Parameters
        ----------
        line : str
            Source line, possibly containing an inline comment.

        Returns
        -------
        tuple[str, str]
            ``(code_part, comment_part)``.  If no ``#`` is present,
            *comment_part* is an empty string.
        """
        parts: list[str] = line.split("#", maxsplit=1)
        if len(parts) == 1:
            parts.append("")

        return parts[0], parts[1]

    @cached_property
    def tokens(self) -> dict[str, Callable[..., Any]]:
        """Mapping from ``@token`` names to their handler methods.

        Returns
        -------
        dict[str, Callable]
            Keys are the supported token names (``"var"``, ``"function"``).
        """
        names: dict[str, Callable[..., Any]] = {
            "var": self.__handle_var,
            "function": self.__handle_function_decleration,
        }

        return names


if __name__ == "__main__":
    import time

    start_time = time.perf_counter()
    afp = AuxiliaryFunctionParser(Path("networks/GOW/GOW.jfunc"))

    end_time = time.perf_counter()
    print(f"Execution time: {end_time - start_time:.6f} seconds")
