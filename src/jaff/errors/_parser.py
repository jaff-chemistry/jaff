# ABOUTME: Exceptions for network-file parsing and SymPy JSON
# ABOUTME: ParserError carries file/line/function context in its message

"""
Custom exception classes for the JAFF parser and serialization subsystems.

Classes
-------
:class:`ParserError`
    Rich contextual error raised when the JAFF network file parser encounters
    a malformed input (bad syntax, cyclic dependencies, unknown symbols, …).
:class:`SympyJsonError`
    Raised by the :mod:`~jaff.common._sympy_json` module when a SymPy
    expression cannot be serialized or deserialized.
"""

from pathlib import Path


class ParserError(Exception):
    """
    Contextual parser error with optional file, line, and function information.

    When raised, the exception message is assembled from all available context
    fields so that the traceback is human-readable without needing to inspect
    the exception attributes separately.

    Parameters
    ----------
    message : str
        The primary error description.
    line : str or None, optional
        The raw source line that caused the error.
    nline : int or None, optional
        1-based line number within *fname*.  Requires *fname* to be set.
    fname : str or pathlib.Path or None, optional
        Path to the source file being parsed.
    funcname : str or None, optional
        Name of the function or parser routine that raised the error.

    Raises
    ------
    AttributeError
        If *nline* is provided but *fname* is ``None`` (a line number without
        a file name is meaningless in a formatted error message).
    """

    def __init__(
        self,
        message: str,
        line: str | None = None,
        nline: int | None = None,
        fname: str | Path | None = None,
        funcname: str | None = None,
    ):
        """Construct the error, assembling a formatted message from all context.

        Parameters
        ----------
        message : str
            Primary error description.
        line : str or None, optional
            Raw source line that caused the error.
        nline : int or None, optional
            1-based line number within *fname*.
        fname : str or Path or None, optional
            Path to the source file being parsed.
        funcname : str or None, optional
            Name of the parser routine that raised the error.

        Raises
        ------
        AttributeError
            If *nline* is provided without *fname*.
        """
        self.message = message
        self.line = line
        self.nline = nline
        self.fname = fname
        self.funcname = funcname

        if fname is None and isinstance(nline, int):
            raise AttributeError(
                "File name must be specified if line number is specified"
            )

        super().__init__(self.__format_message())

    def __format_message(self) -> str:
        """
        Build the full error message string from all available context fields.

        Returns
        -------
        str
            Multi-line error string with optional function name, file path,
            line number, raw line, and the primary message.
        """
        error_str: str = ""

        if self.funcname is not None:
            error_str += f"Function: {self.funcname}\n"

        if self.fname is not None:
            error_str += f"File: {self.fname}\n"

        if self.nline is not None:
            error_str += f"Line number: {self.nline}\n"

        if self.line is not None:
            error_str += f"Line: {self.line}\n"

        error_str += self.message

        return error_str


class SympyJsonError(ValueError):
    """
    Raised by :mod:`~jaff.common._sympy_json` for serialization / deserialization errors.

    Inherits from :exc:`ValueError` so callers can catch it alongside other
    value-related errors when processing untrusted JSON input.
    """

    pass
