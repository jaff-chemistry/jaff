"""Config of the network (``[network]``) and of a ``jaff.toml`` file.

:class:`NetworkConfig` lists the models of the modules it configures as its
fields.  A network's ``jaff.toml`` lives next to its network file
(:func:`jaff_toml_for`).
"""

from __future__ import annotations

from pathlib import Path

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
