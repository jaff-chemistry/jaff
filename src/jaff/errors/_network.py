# ABOUTME: NetworkCheckError: a strict network check (conservation, sinks/sources,
# ABOUTME: recombinations, isomers) failed; carries the check's report
"""Error raised when a strict network check finds a problem."""

from __future__ import annotations

from typing import Any

from ._parser import ParserError

__all__ = ["NetworkCheckError"]


class NetworkCheckError(ParserError):
    """A network check failed with ``strict`` enabled.

    Raised by :meth:`~jaff.core.Network.check_conservation`,
    :meth:`~jaff.core.Network.check_sink_sources`,
    :meth:`~jaff.core.Network.check_recombinations` and
    :meth:`~jaff.core.Network.check_isomers` (run on construction).  A
    :class:`ParserError`, so the CLI tools print it and exit with status 1.

    Attributes
    ----------
    report : ConservationReport or SinkSourceReport or RecombinationReport or IsomerReport
        The failing check's report, listing what it found.
    """

    def __init__(self, message: str, report: Any) -> None:
        """Construct the error.

        Parameters
        ----------
        message : str
            What the check found.
        report : Any
            The failing check's report.
        """
        super().__init__(message)
        self.report: Any = report
