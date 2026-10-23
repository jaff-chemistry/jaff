"""Config of the network (``[network]``) and of a ``jaff.toml`` file.

:class:`NetworkConfig` lists the models of the modules it configures as its
fields.  A network's ``jaff.toml`` lives next to its network file
(:func:`jaff_toml_for`); :func:`load_network_config` stacks it under the
caller's config layers.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence, Union, cast

from ...config.loader import load_config, load_layer, resolve_paths
from ...config.base import ConfigModel, PathLike
from ...physics.dust.config import DustConfig
from ...physics.photo_reactions.config import RadiationConfig
from ...physics.thermodynamics.config import ThermodynamicsConfig
from ..reaction.config import ReactionsConfig
from ._spec import NetworkSpec


class NetworkConfig(ConfigModel):
    """Network (``[network]``).

    Attributes
    ----------
    label : str or None
        Human-readable network name; ``None`` uses the network file stem.
    strict : bool
        Treat validation problems (conservation, sinks/sources, ...) as fatal.
    funcfile : bool or Path
        ``.jfunc`` auxiliary file: ``True`` looks next to the network file,
        ``False`` loads none, a path names it.
    expand_nuclei : bool
        Expand ``n_<element>_nuc`` into the element-nucleus sum.
    reactions : ReactionsConfig
        Settings of the reactions.
    radiation : RadiationConfig or None
        Radiation field; ``None`` disables radiation.
    dust : DustConfig or None
        Dust; ``None`` disables the dust module.
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
    return NetworkSpec._resolve_network_path(network_file).parent / "jaff.toml"


#: One config layer: a TOML file, a dict shaped like one, or a parsed file.
ConfigLayer = Union[str, Path, Mapping[str, Any], JaffTomlConfig]


def load_network_config(
    network_file: Path,
    config: ConfigLayer | Sequence[ConfigLayer] | None = None,
    use_jaff_toml: bool = True,
    overrides: Mapping[str, Any] | None = None,
    cwd: Path | None = None,
    defaults: Mapping[str, Any] | None = None,
) -> NetworkConfig:
    """Load a network's config: *defaults* < jaff.toml < *config* layers < *overrides*.

    Parameters
    ----------
    network_file : Path
        Network file, or the name of a predefined network.
    config : ConfigLayer or sequence of them, optional
        Layers stacked over the network's ``jaff.toml``, lowest first: a TOML
        file path, a ``dict`` shaped like a ``jaff.toml``
        (``{"network": {...}}``), or a :class:`JaffTomlConfig`, of which only
        the keys explicitly set count.  Relative paths resolve against *cwd*.
    use_jaff_toml : bool, optional
        Stack the ``jaff.toml`` next to the network file (if any) at the
        bottom (default ``True``).
    overrides : Mapping or None, optional
        Top layer, shaped like a ``jaff.toml``.
    cwd : Path or None, optional
        Working directory; defaults to the current one.
    defaults : Mapping or None, optional
        Bottom layer, below the ``jaff.toml``, shaped like a ``jaff.toml``.

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
    """
    cwd = cwd or Path.cwd()
    layers: List[Dict[str, Any]] = [resolve_paths(defaults or {}, JaffTomlConfig, cwd)]

    jaff_toml = jaff_toml_for(network_file)
    if use_jaff_toml and jaff_toml.exists():
        layers.append(load_layer(jaff_toml, JaffTomlConfig))

    layers += [_as_layer(item, cwd) for item in _config_items(config)]
    layers.append(resolve_paths(overrides or {}, JaffTomlConfig, cwd))

    return load_config(layers, JaffTomlConfig).network


def _config_items(
    config: ConfigLayer | Sequence[ConfigLayer] | None,
) -> List[ConfigLayer]:
    """Return the *config* argument of :func:`load_network_config` as a list."""
    if config is None:
        return []

    if isinstance(config, (str, Path, JaffTomlConfig, Mapping)):
        return [cast(ConfigLayer, config)]

    return list(config)


def _as_layer(item: ConfigLayer, cwd: Path) -> Dict[str, Any]:
    """Turn one *config* item of :func:`load_network_config` into a layer."""
    if isinstance(item, (str, Path)):
        path = Path(item)

        return load_layer(path if path.is_absolute() else cwd / path, JaffTomlConfig)

    if isinstance(item, JaffTomlConfig):
        return item.model_dump(exclude_unset=True)

    if isinstance(item, Mapping):
        return resolve_paths(cast(Mapping[str, Any], item), JaffTomlConfig, cwd)

    raise TypeError(
        f"Invalid config layer {item!r}: expected a path, dict or JaffTomlConfig"
    )
