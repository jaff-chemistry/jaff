"""UCLCHEM format: comma-delimited reactions with a ``NAN`` sentinel column."""

import math
import re

from ......errors import ParserError


class UclchemReaction:
    """UCLCHEM comma-delimited reaction line handler (``NAN``-sentinel format)."""

    name = "uclchem"
    priority = 70
    is_reaction = True

    global_re = re.compile(r"^(?!\s*[!]|(?:\s*#\s)).*,\s*(?i:NAN)\s*(?:,|$)")

    local_re = re.compile(
        r"^\s*"
        r"(?=.*,\s*(?i:NAN)\s*(?:,|$))"
        r"(?P<reactants>(?:[#@\w\d\+-]*\s*,\s*){3})"
        r"(?P<products>(?:[#@\w\d\+-]*\s*,\s*){4})"
        r"(?P<ka>[^,]*)\s*,\s*"
        r"(?P<kb>[^,]*)\s*,\s*"
        r"(?P<kc>[^,]*)\s*,\s*"
        r"(?P<tmin>[^,]*)\s*,\s*"
        r"(?P<tmax>[^,]*)\s*,\s*"
        r"(?P<extrapolate>.*?)"
        r"\s*$"
    )

    SPECIAL_MAP = {
        "CR": "_CR",
        "CRP": "_CRP",
        "CRPHOT": "_CRPHOT",
        "PHOTON": "_PHOTON",
    }

    ELEMENT_REPS = {"HE": "He", "SI": "Si", "CL": "Cl", "MG": "Mg"}

    IGNORE_SPECIES = {
        "NAN",
        "",
        "ER",
        "ERDES",
        "FREEZE",
        "H2FORM",
        "BULKSWAP",
        "DESCR",
        "DESOH2",
        "DEUVCR",
        "LH",
        "LHDES",
        "SURFSWAP",
        "THERM",
    }

    MECHANISM_TYPE = {
        "FREEZE": "freeze",
        "THERM": "desorption_thermal",
        "DESCR": "desorption_cr",
        "DEUVCR": "desorption_uvcr",
        "DESOH2": "desorption_h2",
        "ER": "eley_rideal",
        "ERDES": "eley_rideal_desorption",
        "LH": "langmuir_hinshelwood",
        "LHDES": "langmuir_hinshelwood_desorption",
        "H2FORM": "h2_formation",
        "BULKSWAP": "bulk_swap",
        "SURFSWAP": "surface_swap",
    }

    # Mechanisms whose rate coefficient JAFF cannot yet emit as a closed-form
    # expression (need dust temperature, a surface-diffusion competition
    # formula, or coupling to another reaction's rate).  Forced to ``0.0``.
    # See README.md for the per-mechanism rationale and source references.
    UNSUPPORTED_RATE = {
        "THERM",
        "DESOH2",
        "LH",
        "LHDES",
        "ER",
        "ERDES",
        "H2FORM",
        "BULKSWAP",
        "SURFSWAP",
    }

    def parse(self, line: str, nline: int, state: dict, file) -> dict:
        """Parse a UCLCHEM-format reaction line into its reaction fields.

        Extracts reactants, products, rate parameters, temperature bounds, and
        an extrapolation flag from the comma-delimited UCLCHEM format (identified
        by the ``NAN`` sentinel column).  Species names are normalised via
        :meth:`_normalize_species`.

        Raises
        ------
        ParserError
            If the line does not match the expected UCLCHEM format.
        """
        local = self.local_re.match(line)
        if not local:
            raise ParserError("Invalid UCLCHEM reaction detected", line, nline, file)

        reactants: str = local.group("reactants")
        products: str = local.group("products")
        ka: float = float(local.group("ka"))
        kb: float = float(local.group("kb"))
        kc: float = float(local.group("kc"))
        tmin: float = float(local.group("tmin"))
        tmax: float = float(local.group("tmax"))
        extrapolate: bool = local.group("extrapolate").strip().lower() == "true"

        # UCLCHEM only switches a rate off below T_min, and gas temperature is
        # never negative, so T_min <= 0 (makerates' placeholder 0.0) is unbounded.
        t_min: float | None = 3.0 if extrapolate else (tmin if tmin > 0 else None)
        t_max: float = 1e6 if extrapolate else tmax

        rr: list[str] = [self._normalize_species(r) for r in reactants.split(",")]
        pp: list[str] = [
            self._normalize_species(p)
            for p in products.split(",")
            if p.strip().upper() not in self.IGNORE_SPECIES
        ]

        mechanism_type: str | None = next(
            (
                self.MECHANISM_TYPE[r.strip().upper()]
                for r in rr
                if r.strip().upper() in self.MECHANISM_TYPE
            ),
            None,
        )

        omega = 0.5
        rate_dict = {
            # Cosmic-ray ionisation:  k = alpha * zeta
            "CRP": f"{ka:.2e} * crate",
            # CR-induced photoreaction (UCLCHEM):
            #   k = alpha * gama/(1-omega) * (T/300)**beta * zeta
            "CRPHOT": (f"{ka * kc / (1.0 - omega):.2e} * (tgas/3e2)**({kb:.2f}) * crate"),
            # UV photoreaction:  k = alpha * radfield * exp(-gama*Av) / 1.7
            # (the 1.7 converts the Draine field to Habing units).
            "PHOTON": f"{ka / 1.7:.2e} * chi * exp(-{kc:.2f} * av)",
            # Freeze-out:  k = freezeFactor * alpha * v_th * sqrt(T/m) * sigma_grain
            "FREEZE": f"(1e0 + {kb:.2e} * 1.671e-3/tgas/asize)*nuth*sigmah*sqrt(tgas/m)",
            # CR thermal desorption (Hasegawa & Herbst 1993):
            #   k = 4*pi*zeta*1.64e-4 * surfaceArea * phi * alpha,
            #   surfaceArea per H = 4 * sigma_grain.
            "DESCR": f"{ka * 4.0 * math.pi * 1.64e-4 * 4.0:.2e} * sigmah * phi * crate",
            # CR-induced UV photodesorption:
            #   k = sigma_grain * uv_yield * 4.875e3 * zeta
            #        * (1 + (radfield/uvcreff)/zeta * exp(-1.8*Av)) * alpha
            "DEUVCR": (
                f"{ka * 4.875e3:.2e} * sigmah * uv_yield * crate"
                f" * (1e0 + (chi/uvcreff)/crate * exp(-1.8 * av))"
            ),
        }

        # Two-body Kooij/Arrhenius default (zero exponents dropped).
        rate = f"{ka:.2e}"
        if kb != 0.0:
            rate += f" * (tgas/3e2)**({kb:.2f})"
        if kc != 0.0:
            rate += f" * exp(-{kc:.2f}/tgas)"

        for r in rr:
            token = r.strip().upper()
            if token in rate_dict:
                rate = rate_dict[token]
                break
            if token in self.UNSUPPORTED_RATE:
                rate = "0.0"
                break

        rr = [r for r in rr if r.strip().upper() not in self.IGNORE_SPECIES]

        # Normalise exotics after rate selection (which keys off raw tokens).
        rr = [self.SPECIAL_MAP.get(r, r) for r in rr]
        pp = [self.SPECIAL_MAP.get(p, p) for p in pp]

        return {
            "r": rr,
            "p": pp,
            "tmin": t_min,
            "tmax": t_max,
            "rate": rate,
            "type": mechanism_type
            if mechanism_type is not None
            else self._reaction_type(rate, rr),
            "string": line.strip(),
        }

    @staticmethod
    def _reaction_type(rate: str, rr: list[str]) -> str:
        """Conclude the reaction type from the reactants, falling back to rate.

        Structural signals are checked first so the result survives custom
        auxiliary-function rates: a ``_PHOTON`` reactant -> ``"photo"``, a
        cosmic-ray pseudo-species (``_CR``/``_CRP``/``_CRPHOT``) ->
        ``"cosmic_ray"``, three or more real reactants -> ``"3_body"``. Only
        then is the rate inspected (``photo``/``av`` -> ``"photo"``, ``crate``
        -> ``"cosmic_ray"``, ``ntot`` -> ``"3_body"``); otherwise ``"unknown"``.
        This fallback only runs for reactions without a mechanism keyword (the
        keyworded ones are classified via :data:`MECHANISM_TYPE`).
        """
        if "_PHOTON" in rr:
            return "photo"
        if any(c in rr for c in ("_CR", "_CRP", "_CRPHOT")):
            return "cosmic_ray"
        if sum(1 for r in rr if not r.startswith("_")) >= 3:
            return "3_body"

        r = rate.lower()
        if "photo" in r:
            return "photo"
        if "crate" in r:
            return "cosmic_ray"
        if "av" in r:
            return "photo"
        if "ntot" in r:
            return "3_body"

        return "unknown"

    @staticmethod
    def _normalize_species(s: str) -> str:
        """Normalise a UCLCHEM species token to the JAFF canonical form.

        Transformations applied:
        - ``#X`` → ``X_DUST`` (grain-surface species prefix)
        - ``@X`` → ``X_BULK`` (bulk ice species prefix)
        - ``E-`` → ``e-`` (electron lower-case)
        - ``HE`` → ``He``, ``SI`` → ``Si``, ``CL`` → ``Cl``, ``MG`` → ``Mg``

        Parameters
        ----------
        s : str
            Raw species token from the UCLCHEM file.

        Returns
        -------
        str
            Normalised species name.
        """
        s = s.strip()
        if s.startswith("#"):
            s = s[1:] + "_DUST"
        if s.startswith("@"):
            s = s[1:] + "_BULK"
        if s == "E-":
            s = "e-"

        for k, v in UclchemReaction.ELEMENT_REPS.items():
            s = s.replace(k, v)

        return s
