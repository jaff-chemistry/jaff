# ABOUTME: KROME @format: header handler declaring the reaction-line column layout
# ABOUTME: Updates the parser state used by later reaction lines

"""KROME ``@format:`` header — declares the column layout for reaction lines."""

import re

from ......errors import ParserError
from ..._typing import kromeFormatProps


class KromeFormatHeader:
    """KROME ``@format:`` header directive — declares column layout."""

    name = "krome_format"
    priority = 10
    is_reaction = False

    global_re = re.compile(r"^\s*@format\s*:(?P<format>.*?)$")

    local_re = re.compile(
        r"^\s*@format\s*:\s*"
        r"(?P<idx>(?i:idx)\s*,\s*)?"
        r"(?P<reactants>(?:(?i:R)\s*,\s*)+)"
        r"(?P<products>(?:(?i:P)\s*,\s*)+)"
        r"(?P<tmin>(?i:tmin)\s*,?\s*)?"
        r"(?P<tmax>(?i:tmax)\s*,?\s*)?"
        r"(?P<rate>(?i:rate)\s*)?\s*$"
    )

    def default_state(self) -> kromeFormatProps:
        return {
            "format_nline": 0,  # line where @format was declared (0 = not yet seen)
            "idx": True,
            "nreact": 3,
            "nprod": 4,
            "tmin": True,
            "tmax": True,
            "rate": True,
        }

    def apply(self, line, nline, state, globals, file, logger) -> None:
        """Parse a KROME ``@format:`` header and update the shared column state.

        Updates *state* with the field counts and flags detected in the format
        declaration so subsequent reaction lines are matched with the correct
        column counts.

        Raises
        ------
        ParserError
            If the format line is malformed.
        """
        local = self.local_re.match(line)
        if local is None:
            self._handle_errors(line, nline, file)

        assert local is not None
        state.update(
            {
                "format_nline": nline,
                "idx": bool(local.group("idx")),
                "nreact": local.group("reactants").lower().count("r"),
                "nprod": local.group("products").lower().count("p"),
                "tmin": bool(local.group("tmin")),
                "tmax": bool(local.group("tmax")),
                "rate": bool(local.group("rate")),
            }
        )

    def _handle_errors(self, line, nline, file) -> None:
        """Raise a descriptive error for a malformed KROME ``@format:`` line."""
        match = self.global_re.match(line)
        assert match is not None

        format = match.group("format")
        if format is None:
            raise ParserError("Empty @format KROME declerative", line, nline, file)

        format = format.strip()
        if not format:
            raise ParserError("Empty @format KROME declerative", line, nline, file)

        if "," not in format:
            raise ParserError(
                "Invalid @format KROME declerative\n"
                "@format decelerative must be separated by ','",
                line,
                nline,
                file,
            )

        expected_tokens = {"idx", "R", "P", "tmin", "tmax", "rate"}
        tokens = [token.strip() for token in format.split(",")]
        for token in tokens:
            if token not in expected_tokens:
                raise ParserError(
                    f"Invalid token in krome format: {token}\n"
                    f"Supported tokens are {','.join(expected_tokens)}",
                    line,
                    nline,
                    file,
                )

        raise ParserError("Invalid @format KROME declerative", line, nline, file)
