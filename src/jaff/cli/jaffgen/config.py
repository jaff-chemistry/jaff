"""Config of a code-generation run (``[jaffgen]``) and of a ``jaffgen.toml`` file.

:func:`load_jaffgen_config` assembles a run's config from the network's
``jaff.toml``, the ``jaffgen.toml`` and the command line, in that order of
precedence (lowest first).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Mapping

from ...config.base import ConfigModel, PathLike
from ...config.loader import deep_merge, load_config, load_layer, resolve_paths
from ...core.network.config import JaffTomlConfig, jaff_toml_for


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
        Output directory; ``None`` uses the default output directory.
    network_file : Path or None
        Network file, or the name of a predefined network.
    network_config : Path or None
        ``jaff.toml`` to load instead of the one next to the network file.
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


def load_jaffgen_config(
    config_file: Path | None = None,
    cli_layer: Mapping[str, Any] | None = None,
    cwd: Path | None = None,
) -> JaffgenTomlConfig:
    """Load a code-generation run's config: jaff.toml < jaffgen.toml < command line.

    Loads in two phases: the ``[jaffgen]`` sections of the ``jaffgen.toml`` and
    of the command line name the network (``network_file``) and its config
    (``network_config``, else the ``jaff.toml`` next to the network, if any);
    that ``jaff.toml`` is then loaded and all three are merged and validated.

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
    JaffgenTomlConfig
        The validated config.

    Raises
    ------
    FileNotFoundError
        If a named config file or the network does not exist.
    ParserError
        If a file is not valid TOML or the config is invalid.
    """
    file_layer = load_layer(config_file, JaffgenTomlConfig) if config_file else {}
    cli = resolve_paths(cli_layer or {}, JaffgenTomlConfig, cwd or Path.cwd())

    run = deep_merge(file_layer, cli).get("jaffgen")
    jaff_toml = _jaff_toml(run if isinstance(run, Mapping) else {})
    jaff_layer = load_layer(jaff_toml, JaffTomlConfig) if jaff_toml else {}

    return load_config([jaff_layer, file_layer, cli], JaffgenTomlConfig)


def _jaff_toml(run: Mapping[str, Any]) -> Path | None:
    """Return the ``jaff.toml`` a ``[jaffgen]`` section selects, or ``None``.

    An explicit ``network_config`` is returned even if missing, so loading it
    reports the missing file; the one next to the network only if it exists.
    """
    explicit = run.get("network_config")
    if isinstance(explicit, Path):
        return explicit

    network_file = run.get("network_file")
    if not isinstance(network_file, Path):
        return None

    candidate = jaff_toml_for(network_file)

    return candidate if candidate.exists() else None
