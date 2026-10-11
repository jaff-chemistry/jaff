# ABOUTME: EosProps: validated EOS configuration (type plus gamma parameters)
# ABOUTME: consumed by EosFactory

from numbers import Real
from typing import Any, Dict, Tuple


class EosProps:
    """EOS configuration passed to :class:`~jaff.core.network.Network`.

    ``type`` selects the EOS; the keyword arguments are the parameters that
    type requires (see :attr:`_REQUIRED`), with any omitted parameter taken
    from :attr:`_DEFAULTS`.  Every argument is validated on
    construction, so a bad configuration fails here rather than in codegen.

    Attributes
    ----------
    type : str
        EOS type, one of the keys of :attr:`_REQUIRED`.
    gamma : float
        Adiabatic index (``ideal``); must be > 1.  Default ``1.6666666666667``
        (monoatomic ideal gas).
    default_gamma : float
        Adiabatic index for species absent from ``gamma_map``
        (``multi_gamma``); must be > 1.
    gamma_map : dict[str, float]
        Per-species adiabatic index keyed by species name (``multi_gamma``);
        every value must be > 1.
    """

    _REQUIRED: Dict[str, Tuple[str, ...]] = {
        "ideal": ("gamma",),
        "multi_gamma": ("default_gamma", "gamma_map"),
        "fermi_degenerate": (),
        "relativistic_fermi_degenerate": (),
    }

    _DEFAULTS: Dict[str, Dict[str, Any]] = {
        "ideal": {"gamma": 1.6666666666667},
    }

    def __init__(self, type: str, **kwargs: Any) -> None:
        """Validate and store the EOS configuration.

        Parameters
        ----------
        type : str
            EOS type, one of the keys of :attr:`_REQUIRED`.
        **kwargs : Any
            Parameters listed for *type* in :attr:`_REQUIRED`; any omitted
            one is taken from :attr:`_DEFAULTS`.

        Raises
        ------
        ValueError
            If ``type`` is unknown, a required key is missing, an unexpected
            key is given, or an adiabatic index is not > 1.
        TypeError
            If an adiabatic index is not a real number or ``gamma_map`` is not
            a ``dict``.
        """
        if type not in self._REQUIRED:
            raise ValueError(
                f"Invalid eos: '{type}'. Valid eos types are: {', '.join(self._REQUIRED)}"
            )

        kwargs = {**self._DEFAULTS.get(type, {}), **kwargs}
        required = self._REQUIRED[type]
        missing = [k for k in required if k not in kwargs]
        if missing:
            raise ValueError(f"eos '{type}' requires: {', '.join(missing)}")

        unknown = [k for k in kwargs if k not in required]
        if unknown:
            raise ValueError(
                f"Unexpected parameters for eos '{type}': {', '.join(unknown)}"
            )

        self.type: str = type
        for key in ("gamma", "default_gamma"):
            if key in kwargs:
                setattr(self, key, self._validate_gamma(key, kwargs[key]))

        if "gamma_map" in kwargs:
            self.gamma_map: Dict[str, float] = self._validate_gamma_map(
                kwargs["gamma_map"]
            )

    @staticmethod
    def _validate_gamma(name: str, value: Any) -> float:
        """Return *value* as a float, requiring a real number > 1."""
        if isinstance(value, bool) or not isinstance(value, Real):
            raise TypeError(f"'{name}' must be a real number, got {value!r}")

        if value <= 1.0:
            raise ValueError(f"'{name}' must be > 1, got {value}")

        return float(value)

    @classmethod
    def _validate_gamma_map(cls, value: Any) -> Dict[str, float]:
        """Return a copy of *value* with every adiabatic index validated."""
        if not isinstance(value, dict):
            raise TypeError(f"'gamma_map' must be a dict, got {type(value).__name__}")

        return {
            name: cls._validate_gamma(f"gamma_map['{name}']", gamma)
            for name, gamma in value.items()
        }
