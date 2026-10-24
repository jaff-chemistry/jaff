# ABOUTME: Package paths (data, db, templates, networks, shielding) and listings of
# ABOUTME: the bundled networks and templates
from pathlib import Path

CONFIG_DIR = Path(__file__).resolve().parent
JAFF_DIR = CONFIG_DIR.parent
SRC_DIR = JAFF_DIR.parent
NETWORKS_DIR = SRC_DIR.parent / "networks"
DATA_DIR = JAFF_DIR / "data"
XSECS_DATA_DIR = DATA_DIR / "xsecs"
SHIELDING_DATA_DIR = DATA_DIR / "shielding"
TEMPLATES_DIR = JAFF_DIR / "templates"
DB_DIR = JAFF_DIR / "db"


def list_subdirs(directory: Path) -> set[str]:
    """Return the names of the immediate sub-directories of *directory*."""
    return {f.name for f in directory.iterdir() if f.is_dir()}


def predefined_networks() -> set[str]:
    """Return the names of the built-in networks (sub-dirs of ``NETWORKS_DIR``)."""
    return list_subdirs(NETWORKS_DIR)


def predefined_templates() -> set[str]:
    """Return the names of the built-in generator template collections."""
    return list_subdirs(TEMPLATES_DIR / "generator")
