# ABOUTME: Package paths (data, db, templates, networks, plugins, shielding), listings of
# ABOUTME: the bundled networks and templates, and network-name resolution
from __future__ import annotations

from pathlib import Path

from ..errors import ParserError

CONFIG_DIR = Path(__file__).resolve().parent
JAFF_DIR = CONFIG_DIR.parent
SRC_DIR = JAFF_DIR.parent
NETWORKS_DIR = SRC_DIR.parent / "networks"
TEMPLATES_DIR = SRC_DIR.parent / "templates"
PLUGINS_DIR = SRC_DIR.parent / "plugins"
DATA_DIR = JAFF_DIR / "data"
XSECS_DATA_DIR = DATA_DIR / "xsecs"
SHIELDING_DATA_DIR = DATA_DIR / "shielding"
DB_DIR = JAFF_DIR / "db"


def list_subdirs(directory: Path) -> set[str]:
    """Return the names of the immediate sub-directories of *directory*."""
    return {f.name for f in directory.iterdir() if f.is_dir()}


def predefined_networks() -> set[str]:
    """Return the names of the built-in networks (sub-dirs of ``NETWORKS_DIR``)."""
    return list_subdirs(NETWORKS_DIR)


def predefined_templates() -> set[str]:
    """Return the names of the built-in templates (sub-dirs of ``TEMPLATES_DIR``)."""
    return list_subdirs(TEMPLATES_DIR)


def resolve_network_path(fname: str | Path) -> Path:
    """Resolve *fname* to a network file path or a predefined network name.

    A predefined name is a sub-directory of :data:`NETWORKS_DIR` that contains
    exactly one ``.jet`` file.  Only the final component of *fname* is compared
    with the predefined names, so any *fname* ending in one (``COthin``, but also
    ``/some/dir/COthin``) resolves to that predefined network, even when a file
    exists at that path.  Other names are treated as filesystem paths, with
    relative paths resolved against the current working directory.

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

    Examples
    --------
    >>> from jaff.config import resolve_network_path
    >>> resolve_network_path("GOW").name
    'GOW.jet'
    """
    p = Path(fname)
    if p.name not in predefined_networks():
        abspath = p.resolve()
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
