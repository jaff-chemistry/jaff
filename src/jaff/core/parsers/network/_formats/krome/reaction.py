# ABOUTME: KROME comma-delimited reaction handler
# ABOUTME: Parses a reaction line against the active @format column layout

"""KROME comma-delimited reaction handler."""

import re

from ......common import f90_convert
from ......errors import ParserError


class KromeReaction:
    """KROME comma-delimited reaction line handler."""

    name = "krome"
    priority = 60
    is_reaction = True

    global_re = re.compile(
        r"^(?!\s*[!#@])"
        r"(?!.*,\s*(?i:NAN)\s*(?:,|$))"
        r"(?=.*,)"
        r"(?P<segment>.*)$"
    )

    SPECIAL_MAP = {
        "CR": "_CR",
        "CRP": "_CRP",
        "CRPHOT": "_CRPHOT",
        "PHOTON": "_PHOTON",
    }

    WHOLE_TOKEN_MAP = {"E": "e-", "e": "e-", "g": "_GRAIN", **SPECIAL_MAP}

    SUBSTRING_MAP = {"HE": "He"}

    TMINMAX_REPS = {
        "d": "e",
        ".le.": "",
        ".ge.": "",
        ".lt.": "",
        ".gt.": "",
        ">": "",
        "<": "",
    }

    RATE_REPS = {
        "user_crflux": "crate",
        "user_crate": "crate",
        "user_av": "av",
    }

    def local_re(self, state: dict) -> re.Pattern:
        """Build the field-extraction regex for the active KROME column layout."""
        return re.compile(
            r"^\s*"
            r"(?!.*,\s*(?i:NAN)\s*(?:,|$))"
            + (r"(?P<idx>[^,]*)\s*,\s*" if state["idx"] else "")
            + rf"(?P<reactants>(?:[^,]*\s*,\s*){{{state['nreact']}}})"
            + rf"(?P<products>(?:[^,]*\s*,\s*){{{state['nprod']}}})"
            + (r"(?P<tmin>[^,]*)\s*,\s*" if state["tmin"] else "")
            + (r"(?P<tmax>[^,]*)\s*,\s*" if state["tmax"] else "")
            + (r"(?P<rate>.*)" if state["rate"] else "")
            + r"\s*$"
        )

    def parse(self, line: str, nline: int, state: dict, file) -> dict:
        """Parse a KROME-format reaction line into its reaction fields.

        Extracts the index, reactants, products, temperature bounds, and rate
        expression from the comma-delimited KROME format.  Applies species
        normalisation (``E``/``e`` → ``e-``, ``g`` → ``_GRAIN``,
        ``HE`` → ``He``), normalises exotic pseudo-species
        (``CR``/``CRP``/``CRPHOT``/``PHOTON`` → underscore form), and converts
        ``user_crflux``/``user_av`` aliases.  Fortran exponent
        notation is converted to Python notation via :func:`~jaff.common.f90_convert`.

        Raises
        ------
        ParserError
            If the line structure is inconsistent with the declared KROME format.
        """
        local = self.local_re(state).match(line)
        if not local:
            self._handle_errors(line, nline, state, file)

        reactants: str = local.group("reactants")
        products: str = local.group("products")
        tmin: str = local.groupdict().get("tmin", "").strip().lower()
        tmax: str = local.groupdict().get("tmax", "").strip().lower()
        rate: str = local.groupdict().get("rate", "").strip()

        rr: list[str] = [r.strip() for r in reactants.split(",")[:-1]]
        pp: list[str] = [p.strip() for p in products.split(",")[:-1]]

        if len(rr) != state["nreact"]:
            raise ParserError(
                "Invalid KROME line detected\n"
                f"Expected {state['nreact']} reactants\n"
                f"from line {state['format_nline']}.\n"
                f"Instead got {len(rr)} reactants",
                line,
                nline,
                file,
            )

        if len(pp) != state["nprod"]:
            raise ParserError(
                "Invalid KROME line detected\n"
                f"Expected {state['nprod']} products \n"
                f"from line {state['format_nline']}.\n"
                f"Instead got {len(pp)} products",
                line,
                nline,
                file,
            )

        t_min: None | float = None
        t_max: None | float = None

        rr = [self.WHOLE_TOKEN_MAP.get(r, r) for r in rr]
        pp = [self.WHOLE_TOKEN_MAP.get(p, p) for p in pp]

        for k, v in self.SUBSTRING_MAP.items():
            rr = [x.replace(k, v) for x in rr]
            pp = [x.replace(k, v) for x in pp]

        rr = [r for r in rr if r != ""]
        pp = [p for p in pp if p != ""]

        if tmin != "none" and tmin != "":
            for k, v in self.TMINMAX_REPS.items():
                tmin = tmin.replace(k, v)
            t_min = float(tmin)

        if tmax != "none" and tmax != "":
            for k, v in self.TMINMAX_REPS.items():
                tmax = tmax.replace(k, v)
            t_max = float(tmax)

        for k, v in self.RATE_REPS.items():
            rate = rate.replace(k, v)

        rate = f90_convert(rate)
        # ``KromeParser.post_process`` numbers the placeholder into photorates().
        if "auto" in rate:
            rate = rate.replace("auto", "krome_photo(1e99, 1e99)")

        if "photo" in rate.lower() and "_PHOTON" not in rr:
            rr.append("_PHOTON")

        return {
            "r": rr,
            "p": pp,
            "tmin": t_min,
            "tmax": t_max,
            "rate": rate,
            "type": self._reaction_type(rate, rr),
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

    def _handle_errors(self, line, nline, state, file) -> None:
        """Raise a descriptive error for a malformed KROME reaction line.

        Diagnoses the most likely cause (wrong field count, wrong reactant or
        product count) before falling back to a generic error message.
        """
        match = self.global_re.match(line)
        assert match is not None

        segment = match.group("segment").lower()
        props = state
        num_fields = (
            int(props["idx"])
            + props["nreact"]
            + props["nreact"]
            + int(props["tmin"])
            + int(props["tmax"])
            + int(props["rate"])
        )
        num_fields_detected: int = segment.count(",") + 1

        if num_fields != num_fields_detected:
            raise ParserError(
                "Number of fields in KROME reaction doesn't match\n"
                f"Number of fields detected: {num_fields_detected}\n"
                f"Number of fields expected: {num_fields}\n"
                + (
                    f"KROME format defined on line: {props['format_nline']}"
                    if props["format_nline"]
                    else ""
                ),
                line,
                nline,
                file,
            )

        if segment.count("r") != props["nreact"]:
            raise ParserError(
                "Expected number of reactants did not match krome format\n"
                f"Number of reactants expected: {props['nreact']}\n"
                f"Number of reactants detected: {segment.count('r')}\n"
                + (
                    f"KROME format defined on line: {props['format_nline']}"
                    if props["format_nline"]
                    else ""
                ),
                line,
                nline,
                file,
            )

        if segment.count("p") != props["nprod"]:
            raise ParserError(
                "Expected number of products did not match krome format\n"
                f"Number of products expected: {props['nprod']}\n"
                f"Number of products detected: {props['nprod']}\n"
                + (
                    f"KROME format defined on line: {props['format_nline']}"
                    if props["format_nline"]
                    else ""
                ),
                line,
                nline,
                file,
            )

        raise ParserError("Invalid KROME reaction detected", line, nline, file)
