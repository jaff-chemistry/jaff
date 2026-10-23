# ABOUTME: NetworkSpec: a Network's resolved file path, settings and .jfunc functions
# ABOUTME: Resolves predefined network names and auxiliary .jfunc files

"""The :class:`NetworkSpec` normalized parameters of a :class:`~jaff.Network`."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from ...config import NETWORKS_DIR, predefined_networks
from ...errors import ParserError
from ..parsers import AuxiliaryFunctionParser
from ..parsers.auxiliary_func import AuxiliaryFunctionProps

if TYPE_CHECKING:
    from .config import NetworkConfig


class NetworkSpec:
    """Normalized construction parameters of a :class:`~jaff.Network`.

    ``Network`` bundles its validated :class:`~jaff.core.network.config.NetworkConfig`
    and the network file into one of these (exposed as :attr:`Network.spec`):

    * ``fname`` is resolved to an absolute path (a filesystem path or a
      predefined network name; see :meth:`_resolve_network_path`);
    * ``funcfile`` is resolved to the actual ``.jfunc`` file (scanning next to
      the network file when ``True``) and its contents are parsed into
      :attr:`aux_funcs`.

    Attributes
    ----------
    fname : Path
        Resolved, absolute path to the network file.
    config : NetworkConfig
        The network's validated config (also :attr:`Network.config`).
    strict : bool
        If ``True``, conservation violations and validation problems are fatal.
    label : str or None
        Optional human-readable label for the network.
    funcfile : bool | Path
        ``True`` to scan the network directory, ``False`` to skip, or the
        resolved path of the ``.jfunc`` file actually used.
    aux_funcs : dict[str, AuxiliaryFunctionProps]
        Parsed auxiliary functions from the ``.jfunc`` file, keyed by lower-cased
        name (``{}`` when none).
        Exposed for inspection/debugging.
    duplicate_policy : str
        Global policy for two rate coefficients sharing a reaction, mechanism,
        and temperature range (``config.reactions.duplicate_policy``), kept for
        reference only: the loader reads each reaction's resolved policy from
        ``config.reactions`` instead, where a reaction's own policy may
        override the global one.
    expand_nuclei : bool
        Whether ``n_<element>_nuc`` symbols expand into element-nucleus sums.
    _from_cli : bool
        Internal flag marking construction from the command-line interface.
    """

    def __init__(
        self,
        fname: str | Path,
        config: NetworkConfig,
        _from_cli: bool,
    ):
        self.fname: Path = self._resolve_network_path(fname)
        self.config: NetworkConfig = config
        self.strict: bool = config.strict
        self.label: str | None = config.label
        self.funcfile: bool | Path = config.funcfile
        self.duplicate_policy: str = config.reactions.duplicate_policy
        # Resolves funcfile to the actual .jfunc path (when True) and parses it.
        self.aux_funcs: dict[str, AuxiliaryFunctionProps] = self._load_aux_funcs()
        self.expand_nuclei: bool = config.expand_nuclei
        self._from_cli: bool = _from_cli

    def _load_aux_funcs(self) -> dict[str, AuxiliaryFunctionProps]:
        """Detect and parse the auxiliary ``.jfunc`` file into a dict.

        When :attr:`funcfile` is ``True``, looks next to the network file for
        ``<file name>.jfunc`` (e.g. ``net.dat.jfunc``) and then for the network
        file with its suffix replaced by ``.jfunc`` (e.g. ``net.jfunc``), and
        rewrites :attr:`funcfile` to the first one found.  ``False`` skips
        loading.

        Returns
        -------
        dict[str, AuxiliaryFunctionProps]
            Parsed auxiliary functions, or ``{}`` when there is no ``.jfunc``.

        Raises
        ------
        FileNotFoundError
            If :attr:`funcfile` is an explicit path that does not exist.
        """
        if self.funcfile is False:
            return {}

        if self.funcfile is True:
            candidates = [
                Path(f"{self.fname}.jfunc"),
                self.fname.with_suffix(".jfunc"),
            ]
            for candidate in candidates:
                if candidate.exists():
                    self.funcfile = candidate
                    break
            else:
                return {}

        assert isinstance(self.funcfile, Path)
        if not self.funcfile.exists():
            raise FileNotFoundError(self.funcfile)

        with AuxiliaryFunctionParser(self.funcfile) as afp:
            return afp.get_dict()

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
