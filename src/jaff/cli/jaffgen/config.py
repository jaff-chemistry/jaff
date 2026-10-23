"""Config of a code-generation run (``[jaffgen]``) and of a ``jaffgen.toml`` file."""

from __future__ import annotations

from typing import Any, Dict, List


from ...config.base import ConfigModel, PathLike
from ...core.network.config import JaffTomlConfig


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
