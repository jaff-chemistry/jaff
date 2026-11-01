# ABOUTME: Config model of a code-generation run ([jaffgen]) and of a jaffgen.toml
# ABOUTME: file, assembled from the jaffgen.toml (or a template's) and the command line
"""Config of a code-generation run (``[jaffgen]``) and of a ``jaffgen.toml`` file.

:func:`load_jaffgen_config` assembles a run's config from the ``jaffgen.toml``
(or the one its template bundles) and the command line, lowest precedence first.
It validates ``[jaffgen]`` and gathers ``[network]`` for the network, which
merges it over its own ``jaff.toml`` and validates it.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Mapping, NamedTuple

from ...config.base import ConfigModel, PathLike
from ...config import TEMPLATES_DIR, predefined_templates
from ...config.loader import (
    atomic_tables,
    deep_merge,
    load_config,
    load_layer,
    resolve_paths,
)
from ...core.network.config import JaffTomlConfig, jaff_toml_for
from ...errors import ParserError

__all__ = [
    "JaffgenConfig",
    "JaffgenRunConfig",
    "JaffgenTomlConfig",
    "bundled_config",
    "load_jaffgen_config",
]


class JaffgenConfig(ConfigModel):
    """Code-generation run (``[jaffgen]``, ``jaffgen.toml`` only).

    Attributes
    ----------
    template : str or None
        Built-in template collection.
    input_dir : Path or None
        Directory of template files to process.
    input_files : list of Path
        Individual template files.
    output_dir : Path or None
        Output directory; without it, ``<repo_root>/generated``. Each run writes
        into a new folder in it, ``<network>_<template>_<YYYYmmdd-HHMMSS>``.
    network_file : Path or None
        Network file, or the name of a predefined network.
    network_config : Path or None
        Another ``jaff.toml``, merged over the one next to the network file.
    echo : bool
        Echo the template directives as comments in the output.
    table : list of dict
        ``[[jaffgen.table]]`` data-table conversions, passed on unchecked.
    """

    template: str | None = None
    input_dir: PathLike | None = None
    input_files: List[PathLike] = []
    output_dir: PathLike | None = None
    network_file: PathLike | None = None
    network_config: PathLike | None = None
    echo: bool = True
    table: List[Dict[str, Any]] = []


class JaffgenTomlConfig(JaffTomlConfig):
    """A ``jaffgen.toml`` file: a code-generation run and its network settings."""

    jaffgen: JaffgenConfig = JaffgenConfig()


class JaffgenRunConfig(NamedTuple):
    """A code-generation run's config, as :func:`load_jaffgen_config` assembles it.

    Attributes
    ----------
    jaffgen : JaffgenConfig
        The validated ``[jaffgen]`` section.
    network : dict
        The run's ``[network]`` settings, merged but not validated: a custom
        ``jaff.toml`` (``network_config``) < the ``jaffgen.toml`` < the command
        line.  :class:`~jaff.core.Network` merges it over the ``jaff.toml``
        next to the network and validates it (``Network(..., config=network)``).
    config_file : Path or None
        The ``jaffgen.toml`` used: the one given, else the one the template
        bundles, else ``None``.
    """

    jaffgen: JaffgenConfig
    network: Dict[str, Any]
    config_file: Path | None


def load_jaffgen_config(
    config_file: Path | None = None,
    cli_layer: Mapping[str, Any] | None = None,
    cwd: Path | None = None,
) -> JaffgenRunConfig:
    """Load a code-generation run's config: jaffgen.toml < command line.

    The ``jaffgen.toml`` is *config_file*, else (with a ``[jaffgen] template``)
    the one the template bundles.  ``[jaffgen]`` is validated here; ``[network]``
    is gathered for the network, which validates it over its own ``jaff.toml``.

    Parameters
    ----------
    config_file : Path or None, optional
        The ``jaffgen.toml``; its paths resolve against its directory.
    cli_layer : Mapping or None, optional
        Command-line settings in the shape of a ``jaffgen.toml``; its paths
        resolve against *cwd*.
    cwd : Path or None, optional
        Working directory; defaults to the current one.

    Returns
    -------
    JaffgenRunConfig
        The validated ``[jaffgen]`` section, the ``[network]`` settings and the
        ``jaffgen.toml`` used.

    Raises
    ------
    FileNotFoundError
        If a named config file does not exist.
    ParserError
        If a file is not valid TOML, a top-level table is unknown, or
        ``[jaffgen]`` is invalid.
    """
    atomic = atomic_tables(JaffgenTomlConfig)
    cwd = cwd or Path.cwd()
    cli = _resolve_table_sources(
        resolve_paths(cli_layer or {}, JaffgenTomlConfig, cwd), cwd
    )

    if config_file is None:
        run = cli.get("jaffgen", {})
        config_file = (
            bundled_config(run.get("template")) if isinstance(run, Mapping) else None
        )
    file_layer = (
        _resolve_table_sources(
            load_layer(config_file, JaffgenTomlConfig), config_file.parent
        )
        if config_file
        else {}
    )

    merged = deep_merge(file_layer, cli, atomic=atomic)
    network = merged.pop("network", {})
    jaffgen = load_config([merged], JaffgenTomlConfig).jaffgen

    custom = _custom_jaff_toml(jaffgen)
    if custom is not None:
        base = load_layer(custom, JaffTomlConfig)
        network = deep_merge(base, {"network": network}, atomic=atomic)["network"]

    return JaffgenRunConfig(jaffgen, network, config_file)


def _resolve_table_sources(layer: Dict[str, Any], base_dir: Path) -> Dict[str, Any]:
    """Return *layer* with each ``[[jaffgen.table]]`` source path made absolute.

    A relative ``source.path`` resolves against *base_dir* (the directory of the
    file that set it, or the working directory for the command line);
    ``"default"`` (the network's data file) is kept.
    """
    run = layer.get("jaffgen")
    if not isinstance(run, Mapping) or not isinstance(run.get("table"), list):
        return layer

    resolved = [_resolve_table_source(table, base_dir) for table in run["table"]]

    return {**layer, "jaffgen": {**run, "table": resolved}}


def _resolve_table_source(table: Any, base_dir: Path) -> Any:
    """*table* with a relative, non-``"default"`` ``source.path`` made absolute."""
    if not isinstance(table, Mapping):
        return table

    source = table.get("source")
    if not isinstance(source, Mapping):
        return table

    path = source.get("path")
    if not isinstance(path, str) or path == "default":
        return table

    return {**table, "source": {**source, "path": str((base_dir / path).resolve())}}


def bundled_config(template: str | None) -> Path | None:
    """The ``jaffgen.toml`` the built-in *template* bundles, if any.

    Looked for anywhere in the template's folder,
    ``jaff/templates/generator/<template>/``.

    Raises
    ------
    ParserError
        If the template holds more than one ``jaffgen.toml``.
    """
    if template is None or template not in predefined_templates():
        return None

    found = [
        p.resolve()
        for p in (TEMPLATES_DIR / "generator" / template).rglob("jaffgen.toml")
    ]

    if len(found) > 1:
        raise ParserError("More than one jaffgen.toml file found in template directory")

    return next(iter(found), None)


def _custom_jaff_toml(run: JaffgenConfig) -> Path | None:
    """The ``network_config`` to merge, unless it is the ``jaff.toml`` beside the network.

    That one is read by the network itself.
    """
    custom = run.network_config
    if custom is None:
        return None

    if not custom.exists():
        raise FileNotFoundError(custom)

    beside = jaff_toml_for(run.network_file) if run.network_file is not None else None
    if beside is not None and beside.resolve() == custom.resolve():
        return None

    return custom
