# ABOUTME: ReactionSource: the network-file line a reaction or rate segment came from
# ABOUTME: Frozen dataclass of the raw line, its line number, source index and format

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ReactionSource:
    """The network-file line a reaction, or one of its rate segments, came from.

    File-side facts only, for error messages and provenance; nothing here
    changes how a reaction behaves.

    Attributes
    ----------
    line : str
        The raw network-file line.
    nline : int
        1-based line number of ``line`` in the network file.
    index : int
        File-side source number (``source_index``): the file-order row number
        of the reaction (0-based).
    format : str or None
        Name of the network-format parser that read the line (e.g.
        ``"krome"``); ``None`` when unknown.
    """

    line: str
    nline: int
    index: int
    format: str | None = None

    def __str__(self) -> str:
        """Return ``"line <nline>: <line>"``.

        Returns
        -------
        str
            The line number and the raw line, for error messages.
        """
        return f"line {self.nline}: {self.line.strip()}"
