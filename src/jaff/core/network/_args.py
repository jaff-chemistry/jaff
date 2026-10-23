"""The :class:`NetworkArgs` raw accumulator of :class:`~jaff.Network` arguments."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ...physics import constants


@dataclass
class NetworkArgs:
    """Raw, mutable accumulator of :class:`~jaff.Network` constructor arguments.

    The CLI tools fill this incrementally as they resolve command-line flags
    and ``jaffgen.toml`` values, then forward the fields to
    :class:`~jaff.Network` (which parses them into its own
    :class:`~jaff.NetworkSpec`).  Unlike ``NetworkSpec`` this performs no
    validation or normalization — it is pure CLI state.  Network settings
    default to ``None`` (not given), so they never override a ``jaff.toml``;
    ``_from_cli`` defaults to ``True`` (the CLI prints its own MOTD banner).
    """

    fname: str | Path = Path()
    config: str | Path | None = None
    strict: bool | None = None
    label: str | None = None
    funcfile: bool | str | Path | None = None
    duplicate_policy: str | None = None
    expand_nuclei: bool | None = None
    rad_bands: list = field(default_factory=list)
    rad_profile_index: int | float | list[int | float] = 0
    rad_mode: str = "nph"
    use_proxy_photoreaction: bool | None = None
    dust: bool = False
    dust_rv: float = 3.1
    dust_u_reduction: str | None = "absorption"
    dust_f_reduction: str | None = "transport"
    dust_pe_threshold_low: float = 6.0
    dust_pe_threshold_high: float = 13.6
    eos: dict[str, Any] | None = None
    background_field: str = "draine"
    c: float = constants.c.cgs.value
    pi_database: str = "norad"
    override_custom_rate: bool = False
    override_custom_cross_section: bool = False
    absorption_attenuates_flux: bool = True
    _from_cli: bool = True
    _metadata: dict[str, Any] = field(default_factory=dict)
