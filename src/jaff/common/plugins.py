# ABOUTME: load_plugins(): imports the repo-root plugins folder as the package jaff.plugins
# ABOUTME: (every public subpackage, at any depth, once) so plugins register themselves
"""Loading of the repo-root ``plugins`` folder."""

from __future__ import annotations

import importlib.util
import sys
from importlib import import_module
from pathlib import Path

from ..config._config import PLUGINS_DIR

PACKAGE = "jaff.plugins"
"""The package name the ``plugins`` folder is imported as."""

_LOADED = False


def load_plugins() -> None:
    """Import the ``plugins`` folder (:data:`~jaff.config.PLUGINS_DIR`) once.

    The folder is imported as the package ``jaff.plugins``, then every subpackage
    whose name does not start with ``_``, at any depth, so the ``register``
    decorators of the plugins run (a package's ``__init__.py`` imports its modules).
    Private packages (``_utils``) are imported only by the plugins that use them.
    A missing folder loads nothing.  Called by every
    :class:`~jaff.common.registry.Registry` on its first lookup, so importing jaff
    loads no plugin.
    """
    global _LOADED
    if _LOADED:
        return

    _LOADED = True
    init = PLUGINS_DIR / "__init__.py"
    if not init.is_file():
        return

    spec = importlib.util.spec_from_file_location(
        PACKAGE, init, submodule_search_locations=[str(PLUGINS_DIR)]
    )
    assert spec is not None and spec.loader is not None
    package = importlib.util.module_from_spec(spec)
    sys.modules[PACKAGE] = package
    spec.loader.exec_module(package)
    _import_tree(PACKAGE, PLUGINS_DIR)


def _import_tree(name: str, path: Path) -> None:
    """Import every public subpackage of the package *name* at *path*, recursively."""
    for entry in sorted(path.iterdir()):
        public = not entry.name.startswith("_")
        if public and entry.is_dir() and (entry / "__init__.py").is_file():
            subpackage = f"{name}.{entry.name}"
            import_module(subpackage)
            _import_tree(subpackage, entry)
