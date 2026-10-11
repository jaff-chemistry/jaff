import logging
from typing import cast

from sympy import Basic, oo

from ...errors import ParserError
from ...io import JaffLogger
from .. import constants

#: Photoionization cross-section databases, in fallback order after the user's choice.
PI_DATABASES: tuple[str, ...] = ("norad", "verner", "leiden")


def validate_pi_database(value: object) -> str:
    """Validate a photoionization database name and return it lower-cased.

    Raises
    ------
    ParserError
        If *value* is not a string naming one of :data:`PI_DATABASES`.
    """
    if not isinstance(value, str) or value.lower() not in PI_DATABASES:
        raise ParserError(
            f"Invalid pi_database {value!r}. Valid values are: {', '.join(PI_DATABASES)}"
        )

    return value.lower()


class RadiationProps:
    """Radiation-field configuration passed to :class:`Network` / :class:`Radiation`.

    Holds the settings that describe the discretised radiation field and its
    assumed spectrum, validating each on construction.

    The ``mode`` selects how the field is tracked: ``"nph"`` for photon number
    density (cm⁻³) or ``"u"`` for energy density (erg cm⁻³).  ``bands`` is the
    ordered list of photon-energy band edges (eV) defining the frequency bands.
    ``profile_index`` is the spectral index *α* of the assumed photon-number
    spectrum ``n(E) ∝ E^(α-2)``: a scalar applied to every band, or a list with
    one *α_i* per band (a piecewise power law).  ``c`` is the speed of light
    (a float in cm/s, or a string such as ``"c_hat"`` for a reduced speed of
    light that becomes a symbol downstream).  ``background_field`` names the background radiation
    field.

    Attributes
    ----------
    profile_index : float or list of float
        Validated spectral index *α*, or one index per band.
    mode : str
        Validated mode, ``"nph"`` or ``"u"`` (lower-cased).
    bands : list of (float or sympy.Basic)
        Validated band-edge list, with ``"inf"`` replaced by ``sympy.oo``.
    c : float or str
        Validated speed of light in cm/s, or a string to become a symbol.
    background_field : str
        Validated background-field name (lower-cased).
    pi_database : str
        Validated photoionization database name (lower-cased).
    """

    # Vaiid modes
    # nph: Photon number density
    # u: Enregy desnity
    _valid_modes: tuple[str, ...] = ("nph", "u")

    _valid_fields: tuple[str, ...] = (
        "bb_4000",
        "bb_10000",
        "bb_20000",
        "draine",
        "habing",
        "mathis",
        "solar",
        "tw_hydra",
    )

    def __init__(
        self,
        bands: list[str | float | Basic] = [],
        profile_index: float | list[float] = 0.0,
        mode: str = "nph",  # nph or u,
        c: float | str = constants.c.cgs.value,
        background_field: str = "draine",
        pi_database: str = "norad",
    ):
        """Validate and store the radiation-field configuration.

        Parameters
        ----------
        bands : list of (str, float, or sympy.Basic), optional
            Ordered photon-energy band edges in eV (default ``[]``).  The
            string ``"inf"`` is accepted in the last slot and replaced with
            ``sympy.oo``.
        profile_index : float or list of float, optional
            Spectral index *α* of the photon-number spectrum
            ``n(E) ∝ E^(α-2)`` (default ``0.0``).  Must be an int or float,
            or a list of ``len(bands) - 1`` ints/floats (one per band).
        mode : str, optional
            Radiation-tracking mode (default ``"nph"``); one of
            ``("nph", "u")`` -- photon number density or energy density.
        c : float or str, optional
            Speed of light in cm/s (default the CGS value); may also be a
            string such as ``"c_hat"`` for a reduced speed of light.
        background_field : str, optional
            Background radiation field (default ``"draine"``); one of
            ``("bb_4000", "bb_10000", "bb_20000", "draine", "habing",
            "mathis", "solar", "tw_hydra")`` (case-insensitive).
        pi_database : str, optional
            Photoionization cross-section database (default ``"norad"``); one of
            ``("norad", "verner", "leiden")`` (case-insensitive).  Photo-
            dissociation always uses Leiden.  See
            :meth:`~jaff.physics.photo_reactions._photochemistry.Photochemistry.get_xsec`
            for the fallback order.

        Raises
        ------
        ParserError
            If any argument fails validation (invalid type, mode, or field).
        """
        self.logger: logging.Logger = JaffLogger().get_logger()
        self.profile_index: float | list[float] = self._validate_profile_index(
            profile_index
        )
        self.mode: str = self._validate_mode(mode)
        self.bands: list[float | Basic] = self._validate_bands(bands)
        self.c: float | str = self._validate_c(c)
        self.background_field: str = self._validate_field(background_field)
        self.pi_database: str = validate_pi_database(pi_database)

    def _validate_c(self, c) -> float | str:
        if isinstance(c, (float, int, str)):
            return c

        raise ParserError(
            f"Speed of light must be of <float> or <string>. Found {type(c)}"
        )

    def _validate_field(self, field: str) -> str:
        if not isinstance(field, str):
            raise ParserError(
                f"Background radiation field must be a string: Supported fields are: {', '.join(self._valid_fields)}"
            )

        if field.lower() not in self._valid_fields:
            raise ParserError(
                f"Invalid background field {field}. Supported fields are {', '.join(self._valid_fields)}"
            )

        return field.lower()

    def _validate_mode(self, mode: str):
        if not isinstance(mode, str):
            raise ParserError(
                f"Radiation mode must be of type <str>. Valid modes are: {', '.join(self._valid_modes)}"
            )

        if mode.lower() not in self._valid_modes:
            raise ParserError(
                f"Invalid radiation mode {mode}. Valid modes are: {', '.join(self._valid_modes)}"
            )

        return mode.lower()

    @staticmethod
    def _is_number(value: object) -> bool:
        # bool is a subclass of int but is never a meaningful spectral index.
        return isinstance(value, (float, int)) and not isinstance(value, bool)

    def _validate_profile_index(self, index: float | list[float]) -> float | list[float]:
        if self._is_number(index):
            return index

        if isinstance(index, list) and index and all(self._is_number(i) for i in index):
            return list(index)

        raise ParserError(
            f"Invalid radiation profile index: {index!r} ({type(index)})\n"
            "Radiation profile index must be an integer, float or a non-empty list "
            "of int/float (one per band)"
        )

    def _band_profile_indices(self, nbands: int) -> list[float]:
        """Return one spectral index per band, validating a list's length.

        Parameters
        ----------
        nbands : int
            Number of bands (``len(bands) - 1``).

        Returns
        -------
        list of float
            The per-band spectral indices.

        Raises
        ------
        ParserError
            If ``profile_index`` is a list whose length is not ``nbands``.
        """
        if not isinstance(self.profile_index, list):
            return [float(self.profile_index)] * nbands

        if len(self.profile_index) != nbands:
            raise ParserError(
                f"profile_index has {len(self.profile_index)} entries but there are "
                f"{nbands} radiation bands; supply one index per band or a scalar"
            )

        return [float(i) for i in self.profile_index]

    def _validate_bands(self, bands: list[float | str | Basic]) -> list[float | Basic]:
        """
        Validate and store the band-edge list, replacing ``"inf"`` with ``sympy.oo``.

        Also checks that the power-law photon-number spectrum is integrable
        over the supplied band range when energy-density mode is active.
        The average energy ``<E>_i = ∫ E·n(E) dE / ∫ n(E) dE`` must
        converge; this requires:

        - The lower edge to be non-zero when the first band's spectral index
          is steep enough to cause a divergence at ``E → 0``.
        - The upper edge to be finite when the last band's spectral index is
          shallow enough to cause a divergence at ``E → ∞``.

        Interior bands have finite, non-zero edges and always converge.

        Parameters
        ----------
        bands : list of (float, int, str, or sympy.Basic)
            Mutable band-edge list; modified in-place to replace any
            ``"inf"`` string with ``sympy.oo``.

        Raises
        ------
        ParserError
            If ``bands`` contains a string entry other than ``"inf"`` in the
            last slot, has fewer than two edges, or a list ``profile_index``
            does not have one entry per band.
        RuntimeError
            If the average-energy integral would diverge given the supplied
            band edges and power-law index.
        """
        # Replace the sentinel string "inf" with SymPy's infinity symbol.
        if "inf" in bands:
            inf_index = bands.index("inf")
            bands[inf_index] = oo

        if any(isinstance(v, str) for v in bands):
            raise ParserError(
                f"Only 'inf' is supported as a string entry for radiation bands in the last slot. Found: {bands}"
            )

        # self.bands = cast(list[int | float | Basic], bands)

        if len(bands) < 2:
            raise ParserError(
                f"Radiation bands need at least two edges (one band). Found: {bands}"
            )

        alphas = self._band_profile_indices(len(bands) - 1)
        # Only the first band reaches E -> 0 and only the last reaches E -> inf,
        # so each edge is checked against its own band's index.
        alpha_lo, alpha_hi = alphas[0], alphas[-1]
        starts_at_zero = isinstance(bands[0], (float, int)) and float(bands[0]) == 0.0
        ends_at_inf = bands[-1] == oo

        if self.mode == "u":
            # The average-energy integral uses the *energy-density* spectrum
            # u(E) ∝ E^(α-1), so the integral ∫ E · u(E) dE ∝ ∫ E^α dE.
            # The effective power-law index for the ∫ E·n(E) dE integral is
            # pl_index = (α-2) + 1 = α - 1.
            # pl_index + 1 = α: α < 0 diverges at E → 0, α > 0 diverges at
            # E → ∞, and α == 0 (pl_index == -1) log-diverges at both ends.
            if starts_at_zero and alpha_lo <= 0.0:
                raise RuntimeError(
                    f"The integral for average energy will diverge since the radiation band starts from bands[0]: {bands[0]}\n"
                    "Please try a non-zero value"
                )
            if ends_at_inf and alpha_hi >= 0.0:
                raise RuntimeError(
                    f'The integral for average energy will diverge since the radiation band ends at bands[{len(bands) - 1}]: "inf"\n'
                    "Please try a non-infinite value or change the profile_index"
                )

        if (
            alpha_lo <= 1.0
            and isinstance(bands[0], (float, int))
            and float(bands[0]) < 1.0
        ):
            self.logger.warning(
                f"Radiation band starts at bands[0]={bands[0]} eV with "
                f"profile_index={alpha_lo} (first band): the photon-number "
                "normalisation integral ∫E^(α-2)dE is lower-edge divergent "
                "(exponent ≤ -1) and near E→0 becomes ill-conditioned, which "
                "can yield a negative/garbage photon density. Use a non-zero "
                "bands[0] ≳ 1 eV or a larger profile_index."
            )

        return cast(list[float | Basic], bands)
