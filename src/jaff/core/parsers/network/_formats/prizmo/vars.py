# ABOUTME: PRIZMO VARIABLES { } block handler: toggles block state and stores variables
# ABOUTME: Variable assignments become symbolic globals

"""PRIZMO ``VARIABLES { }`` block lines and variable assignments (directive)."""

import re

from sympy import parse_expr

from ......common import f90_convert
from ......errors import ParserError
from ..._typing import prizmoFormatProps


class PrizmoVars:
    """PRIZMO ``VARIABLES { }`` block / variable-assignment directive handler."""

    name = "prizmo_vars"
    priority = 30
    is_reaction = False

    global_re = re.compile(
        r"^\s*(?:"
        r"(?:(?i:variables)\s*\{|\})(?P<segment>.*?)"
        r"|"
        r"(?P<assignment>\w+\s*=.*?)"
        r")\s*$"
    )

    local_re = re.compile(
        r"^\s*(?P<begin>(?i:variables)\s*\{)\s*$"
        r"|"
        r"^\s*(?P<end>\}\s*)$"
        r"|"
        r"^\s*(?P<var>\w+)\s*=\s*\s*(?P<expr>.*?)\s*$"
    )

    def default_state(self) -> prizmoFormatProps:
        return {"parse_vars": False}

    def apply(self, line, nline, state, globals, file, logger) -> None:
        """Handle a PRIZMO ``VARIABLES { }`` block line or variable assignment.

        Toggles ``parse_vars`` on ``VARIABLES {`` and ``}`` tokens, and stores a
        parsed SymPy expression for any ``var = expr`` line encountered while
        inside the block.

        Raises
        ------
        ParserError
            If the line is malformed.
        """
        local = self.local_re.match(line)
        if not local:
            self._handle_errors(line, nline, state, file)

        assert local is not None

        if local.group("begin"):
            state["parse_vars"] = True
            return

        if local.group("end"):
            state["parse_vars"] = False
            return

        if local.group("var") and local.group("expr") and state["parse_vars"]:
            try:
                globals[local.group("var").lower()] = parse_expr(
                    f90_convert(local.group("expr").lower())
                )

            except (SyntaxError, NameError, TypeError):
                logger.warning(
                    f"Skipping variable: {local.group('var')}\n"
                    f"at line: {nline} since the expression is invalid sympy syntax"
                )

    def _handle_errors(self, line, nline, state, file) -> None:
        """Raise a descriptive error for a malformed PRIZMO variables section line."""
        match = self.global_re.match(line)
        assert match is not None

        segment = match.group("segment")
        assignment = match.group("assignment")

        if segment is None and assignment is None:
            raise ParserError("Invalid PRIZMO variable section", line, nline, file)

        if assignment is not None:
            if not state["parse_vars"]:
                raise ParserError(
                    "PRIZMO variable assignment found outside VARIABLES block",
                    line,
                    nline,
                    file,
                )

            var_name, expr = assignment.split("=", 1)
            var_name = var_name.strip()
            expr = expr.strip()

            if not var_name.isidentifier():
                raise ParserError(
                    f"Invalid variable name '{var_name}'", line, nline, file
                )

            if not expr:
                raise ParserError("Expression cannot be empty", line, nline, file)

        segment = segment.strip()
        if segment:
            raise ParserError(
                "Extra characters found after PRIZMO block declarative",
                line,
                nline,
                file,
            )
