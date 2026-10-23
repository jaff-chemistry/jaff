# ABOUTME: NetworkSpec: a Network's resolved file path and its .jfunc file path
# ABOUTME: Resolves predefined network names and auxiliary .jfunc files

"""The :class:`NetworkSpec` source-file facts of a :class:`~jaff.Network`."""

from __future__ import annotations

from pathlib import Path
from ...config import NETWORKS_DIR, predefined_networks
from ...errors import ParserError


class NetworkSpec:
    """Facts about a :class:`~jaff.Network`'s source files.

    ``Network`` builds one from its file name and its ``funcfile`` setting
    (exposed as :attr:`Network.spec`).  It holds what was found on disk, not
    settings: those are :attr:`Network.config`.

    Attributes
    ----------
    fname : Path
        Resolved, absolute path to the network file (a filesystem path or a
        predefined network name; see :meth:`_resolve_network_path`).
    jfunc_path : Path or None
        The auxiliary ``.jfunc`` file actually used: the ``funcfile`` setting
        when it is a path, else the file found next to the network file
        (``<file name>.jfunc``, then ``<stem>.jfunc``); ``None`` when there is
        none or ``funcfile`` is ``False``.  Its functions are
        :attr:`Network.functions`.
    """

    def __init__(self, fname: str | Path, funcfile: bool | Path):
        """Resolve the network file and its ``.jfunc`` file.

        Parameters
        ----------
        fname : str | Path
            Network file path or predefined network name.
        funcfile : bool | Path
            The ``network.funcfile`` setting: a ``.jfunc`` path, ``True`` to
            look next to the network file, or ``False`` for none.

        Raises
        ------
        FileNotFoundError
            If the network file, or an explicit ``.jfunc`` path, does not exist.
        ParserError
            If a predefined network directory has more than one ``.jet`` file.
        """
        self.fname: Path = self._resolve_network_path(fname)
        self.jfunc_path: Path | None = self._find_jfunc(funcfile)

    def _find_jfunc(self, funcfile: bool | Path) -> Path | None:
        """Return the ``.jfunc`` file *funcfile* selects, or ``None``.

        Raises
        ------
        FileNotFoundError
            If *funcfile* is an explicit path that does not exist.
        """
        if funcfile is False:
            return None

        if funcfile is True:
            candidates = [Path(f"{self.fname}.jfunc"), self.fname.with_suffix(".jfunc")]

            return next((c for c in candidates if c.exists()), None)

        if not funcfile.exists():
            raise FileNotFoundError(funcfile)

        return funcfile

    @staticmethod
    def _resolve_network_path(fname: str | Path) -> Path:
        """Resolve *fname* to a network file path or a predefined network name.

        A predefined name is a sub-directory of :data:`~jaff.config.NETWORKS_DIR`
        that contains exactly one ``.jet`` file.  Only the final component of
        *fname* is compared with the predefined names, so any *fname* ending in
        one (``COthin``, but also ``/some/dir/COthin``) resolves to that
        predefined network, even when a file exists at that path.  Other names
        are treated as filesystem paths, with relative paths resolved against
        the current working directory.

        Parameters
        ----------
        fname : str | Path
            A filesystem path or a predefined network name.

        Returns
        -------
        Path
            The resolved, absolute network file path.

        Raises
        ------
        FileNotFoundError
            If *fname* is neither an existing file nor a predefined name, or a
            predefined network directory has no ``.jet`` file.
        ParserError
            If a predefined network directory has more than one ``.jet`` file.
        """
        p = Path(fname)
        abspath = p.resolve()

        names = predefined_networks()
        if p.name not in names:
            if not abspath.exists():
                raise FileNotFoundError(f"Network file '{fname}' not found")

            return abspath

        ndir = NETWORKS_DIR / p.name
        jets = sorted(f for f in ndir.iterdir() if f.suffix.lower() == ".jet")
        if not jets:
            raise FileNotFoundError(f"No .jet file in predefined network '{p.name}'")

        if len(jets) > 1:
            raise ParserError(
                f"Predefined network '{p.name}' has multiple .jet files: "
                f"{[j.name for j in jets]}"
            )

        return jets[0].resolve()
