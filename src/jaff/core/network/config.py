# ABOUTME: Config models of the network ([network]) and of a jaff.toml file;
# ABOUTME: NetworkConfig lists the per-module config models as its fields
"""Config of the network (``[network]``) and of a ``jaff.toml`` file.

:class:`NetworkConfig` lists the models of the modules it configures as its
fields.  A network's ``jaff.toml`` lives next to its network file
(:func:`jaff_toml_for`) and holds its defaults; :func:`load_network_config`
merges the caller's config and overrides over it.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Mapping, Union, cast

from ...config import resolve_network_path
from ...config.loader import load_config, load_layer, resolve_paths
from ...config.base import ConfigModel, PathLike, as_layer
from ...physics.dust.config import DustConfig
from ...physics.photo_reactions.config import RadiationConfig
from ...physics.thermodynamics.config import ThermodynamicsConfig
from ..reaction.config import ReactionsConfig

__all__ = [
    "JaffTomlConfig",
    "NetworkConfig",
    "NetworkConfigSource",
    "jaff_toml_for",
    "load_network_config",
]


class NetworkConfig(ConfigModel):
    """Network (``[network]``).

    Attributes
    ----------
    label : str or None
        Human-readable network name; without it, the network file stem.
    strict : bool
        Treat validation problems (conservation, sinks/sources, ...) as fatal.
    funcfile : bool or Path
        ``.jfunc`` auxiliary file: ``true`` looks next to the network file,
        ``false`` loads none, a path names it.
    expand_nuclei : bool
        Expand ``n_<element>_nuc`` into the element-nucleus sum.
    reactions : ReactionsConfig
        Settings of the reactions.
    radiation : RadiationConfig or None
        Radiation field; the table enables radiation.
    dust : DustConfig or None
        Dust; the table enables the dust module.
    thermodynamics : ThermodynamicsConfig
        Thermodynamics.
    """

    label: str | None = None
    strict: bool = False
    funcfile: bool | PathLike = True
    expand_nuclei: bool = True
    reactions: ReactionsConfig = ReactionsConfig()
    radiation: RadiationConfig | None = None
    dust: DustConfig | None = None
    thermodynamics: ThermodynamicsConfig = ThermodynamicsConfig()


class JaffTomlConfig(ConfigModel):
    """A ``jaff.toml`` file: the settings of the network next to it."""

    network: NetworkConfig = NetworkConfig()


def jaff_toml_for(network_file: Path) -> Path:
    """Return the ``jaff.toml`` path next to *network_file* (it may not exist).

    Parameters
    ----------
    network_file : Path
        Network file, or the name of a predefined network.

    Returns
    -------
    Path
        ``jaff.toml`` in the network file's directory.

    Raises
    ------
    FileNotFoundError
        If *network_file* is neither a file nor a predefined network.
    """
    return resolve_network_path(network_file).parent / "jaff.toml"


#: The ``config`` of a network: a ``jaff.toml``-style file (its ``[network]``
#: table), a dict of the ``[network]`` table's keys, or a ``NetworkConfig``.
NetworkConfigSource = Union[str, Path, Mapping[str, Any], NetworkConfig]


def load_network_config(
    network_file: Path,
    config: NetworkConfigSource | None = None,
    overrides: Mapping[str, Any] | None = None,
    cwd: Path | None = None,
) -> NetworkConfig:
    """Load a network's config: defaults < its ``jaff.toml`` < *config* < *overrides*.

    The layers merge key by key (a table given in a higher layer updates the
    lower one; an ``eos`` or a shielding table replaces it whole) and the
    result is validated once.

    Parameters
    ----------
    network_file : Path
        Network file, or the name of a predefined network; the ``jaff.toml``
        next to it, if any, holds the network's own defaults.
    config : NetworkConfigSource or None, optional
        One more source merged over that ``jaff.toml``: a TOML file (its
        ``[network]`` table), a dict of the ``[network]`` table's keys
        (``{"strict": True, "reactions": {...}}``), or a
        :class:`NetworkConfig` (only the keys explicitly set on it count).
        Relative paths resolve against *cwd*.
    overrides : Mapping or None, optional
        Top layer, a dict of the ``[network]`` table's keys (the keyword
        arguments of :class:`~jaff.core.Network`).
    cwd : Path or None, optional
        Working directory; defaults to the current one.

    Returns
    -------
    NetworkConfig
        The validated ``[network]`` config.

    Raises
    ------
    FileNotFoundError
        If a config file or the network does not exist.
    ParserError
        If a file is not valid TOML or the config is invalid.
    TypeError
        If *config* is not one of the accepted kinds.
    """
    cwd = cwd or Path.cwd()
    layers: List[Dict[str, Any]] = []

    jaff_toml = jaff_toml_for(network_file)
    if jaff_toml.exists():
        layers.append(load_layer(jaff_toml, JaffTomlConfig))

    if config is not None:
        layers.append(_source_layer(config, cwd))

    layers.append(_network_layer(overrides or {}, cwd))

    return load_config(layers, JaffTomlConfig).network


def _source_layer(config: NetworkConfigSource, cwd: Path) -> Dict[str, Any]:
    """Turn the *config* of :func:`load_network_config` into a ``jaff.toml`` layer."""
    if isinstance(config, (str, Path)):
        path = Path(config)

        return load_layer(path if path.is_absolute() else cwd / path, JaffTomlConfig)

    if isinstance(config, (NetworkConfig, Mapping)):
        return _network_layer(config, cwd)

    raise TypeError(
        f"Invalid network config {config!r}: expected a TOML file path, a dict of "
        "the [network] keys or a NetworkConfig"
    )


def _network_layer(table: NetworkConfig | Mapping[str, Any], cwd: Path) -> Dict[str, Any]:
    """Wrap a ``[network]`` table as a ``jaff.toml`` layer, paths resolved against *cwd*."""
    layer = {"network": as_layer(table)}

    return resolve_paths(cast(Mapping[str, Any], layer), JaffTomlConfig, cwd)
