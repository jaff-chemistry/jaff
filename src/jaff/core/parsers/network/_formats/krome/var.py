# ABOUTME: KROME @var: directive handler storing a symbolic global expression
# ABOUTME: Invalid SymPy expressions are skipped with a warning

"""KROME ``@var:`` directive — stores a symbolic global expression."""

import re

from sympy import parse_expr

from ......common import f90_convert
from ......errors import ParserError


class KromeVar:
    """KROME ``@var:`` directive handler — stores a symbolic global expression."""

    name = "krome_var"
    priority = 20
    is_reaction = False

    global_re = re.compile(r"^\s*@var\s*:(?P<segment>.*?)$")

    local_re = re.compile(r"^\s*@var\s*:\s*(?P<var>\w+)\s*=\s*\s*(?P<expr>.*?)\s*$")

    def apply(self, line, nline, state, globals, file, logger) -> None:
        """Parse a KROME ``@var:`` directive and store the symbolic expression.

        Logs a warning and skips the variable when the expression is not valid
        SymPy syntax.

        Raises
        ------
        ParserError
            If the line is not of the form ``@var: <name> = <expr>``.
        """
        local = self.local_re.match(line)
        if not local:
            raise ParserError(
                "Invalid KROME variable assignment detected", line, nline, file
            )

        try:
            globals[local.group("var").lower()] = parse_expr(
                f90_convert(local.group("expr").lower())
            )
        except (SyntaxError, NameError, TypeError):
            logger.warning(
                f"Skipping variable: {local.group('var')}\n"
                f"at line: {nline} since the expression is invalid sympy syntax"
            )
