# ABOUTME: TOML driver: Toml parses a TOML file once with tomllib and exposes
# ABOUTME: top-level key lookup; usable as a context manager.
"""
TOML configuration file reader.

This module provides :class:`Toml`, a lightweight wrapper around the standard
library :mod:`tomllib` module that parses a TOML file on construction and
exposes top-level key lookup.
"""

from pathlib import Path

import tomllib


class Toml:
    """
    Read and expose the contents of a TOML configuration file.

    Parses the file once at construction time and stores the result in
    :attr:`data`.  Supports use as a context manager (though no special
    cleanup beyond deletion of the :attr:`file` attribute is performed on
    exit).

    Parameters
    ----------
    file : str or Path
        Path to the TOML file to read.

    Attributes
    ----------
    file : Path
        Path to the TOML source file, as given (a ``str`` is converted to
        :class:`~pathlib.Path`; it is not resolved).
    data : dict
        Parsed contents of the TOML file as a nested dictionary.

    Examples
    --------
    Given ``config.toml``:

    ```toml
    temperature = 300.0
    species = ["H", "O", "H2O"]

    [network]
    host = "localhost"
    port = 8080
    ```

    >>> from jaff.drivers import Toml
    >>> t = Toml("config.toml")
    >>> t.get_key("temperature")
    300.0
    >>> t.get_key("network")
    {'host': 'localhost', 'port': 8080}
    >>> t.get_key("flux") is None
    True

    As a context manager:

    >>> with Toml("config.toml") as t:
    ...     t.get_key("species")
    ['H', 'O', 'H2O']
    """

    def __init__(self, file: str | Path):
        """Parse the TOML file and store its contents in :attr:`data`.

        Parameters
        ----------
        file : str or Path
            Path to the TOML configuration file to read.
        """
        if isinstance(file, str):
            file = Path(file)

        self.file = file
        self.data = self.__get_dict()

    def __get_dict(self) -> dict:
        """
        Parse the TOML file and return its contents as a dictionary.

        Returns
        -------
        dict
            Nested dictionary representation of the TOML file.
        """
        with open(self.file, "rb") as f:
            data = tomllib.load(f)

        return data

    def get_key(self, key: str) -> dict | None:
        """
        Retrieve the value for a top-level TOML key.

        A table (section) yields a ``dict``; a scalar or array yields the
        matching Python value (``str``, ``int``, ``float``, ``bool``,
        ``list``, ...).  A missing key returns ``None`` rather than raising.

        Parameters
        ----------
        key : str
            Top-level key to look up (e.g. ``"jaffgen"``, ``"network"``).

        Returns
        -------
        Any or None
            The value associated with *key*, or ``None`` if the key is
            absent from the parsed data.

        Examples
        --------
        With ``config.toml`` as in the :class:`Toml` examples:

        >>> Toml("config.toml").get_key("network")
        {'host': 'localhost', 'port': 8080}
        """
        return self.data.get(key, None)

    # ------------------------------------------------------------------
    # Context manager protocol
    # ------------------------------------------------------------------

    def __enter__(self) -> "Toml":
        """Return self when entering a ``with`` block."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Release the file reference when exiting a ``with`` block."""
        del self.file
