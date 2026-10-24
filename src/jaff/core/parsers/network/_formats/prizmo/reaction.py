# ABOUTME: PRIZMO arrow-notation reaction handler
# ABOUTME: Parses reactants -> products [tmin, tmax] rate lines

"""PRIZMO arrow-notation reaction handler."""

import re

from ......errors import ParserError


class PrizmoReaction:
    """PRIZMO arrow-notation reaction line handler."""

    name = "prizmo"
    priority = 40
    is_reaction = True

    global_re = re.compile(r"^(?!\s*[!#]).*->.*$")

    #: A photo rate: ``photo(lo, hi)``, ``PHOTO, lo[, hi]`` or a bare ``PHOTO``.
    photo_re = re.compile(r"^photo\s*(\(|,|$)", re.IGNORECASE)

    local_re = re.compile(
        r"^\s*"
        r"(?P<reactants>[\w\+\-\s]+)"
        r"\s*->\s*"
        r"(?P<products>[\w\+\-\s]+)"
        r"\s*\[\s*"
        r"(?P<tmin>[^,\]]*)?"
        r"\s*,?\s*"
        r"(?P<tmax>[^,\]]*)?"
        r"\s*\]\s*"
        r"(?P<rate>.*)"
        r"\s*$"
    )

    SPECIAL_MAP = {
        "GRAIN0": "_GRAIN",
        "GRAIN": "_GRAIN",
        "CR": "_CR",
        "CRP": "_CRP",
        "CRPHOT": "_CRPHOT",
        "PHOTON": "_PHOTON",
        "dummy": "_DUMMY",
    }

    def parse(self, line: str, nline: int, state: dict, file) -> dict:
        """Parse a PRIZMO-format reaction line into its reaction fields.

        Extracts reactants, products, optional temperature bounds, and rate
        expression from the ``R1 + R2 -> P1 + P2 [tmin, tmax] rate`` pattern.
        Applies species-name normalisation (``HE`` → ``He``, ``E`` → ``e-``)
        and exotic pseudo-species normalisation (``GRAIN0``/``GRAIN`` →
        ``_GRAIN``, ``CR`` → ``_CR``, ``PHOTON`` → ``_PHOTON``, ``dummy`` →
        ``_DUMMY``, ...), and converts ``user_crflux``/``user_av`` aliases to
        canonical JAFF symbols.

        Raises
        ------
        ParserError
            If the line does not match the expected PRIZMO format.
        """
        local = self.local_re.match(line)
        if not local:
            raise ParserError("Invalid PRIZMO reaction detected", line, nline, file)

        reactants: str = local.group("reactants")
        products: str = local.group("products")
        tmin: str | None = local.group("tmin")
        tmax: str | None = local.group("tmax")
        rate: str = local.group("rate").strip()

        reactants = (
            reactants.replace("HE", "He").replace(" E", " e-").replace("E ", "e- ")
        )
        products = products.replace("HE", "He").replace(" E", " e-").replace("E ", "e- ")

        rr: list[str] = [
            self.SPECIAL_MAP.get(r.strip(), r.strip()) for r in reactants.split(" + ")
        ]
        pp: list[str] = [
            self.SPECIAL_MAP.get(p.strip(), p.strip()) for p in products.split(" + ")
        ]

        t_min: float | None = (
            float(tmin.strip().replace("d", "e")) if tmin and tmin.strip() else None
        )
        t_min = t_min if (t_min is not None and t_min > 0) else None

        t_max: float | None = (
            float(tmax.strip().replace("d", "e")) if tmax and tmax.strip() else None
        )
        t_max = t_max if (t_max is not None and t_max < 1e8) else None

        rate = rate.replace("user_crflux", "crate").replace("user_av", "av")

        if self.photo_re.match(rate):
            rate = self._photo_placeholder(rate, line, nline, file)
            if "_PHOTON" not in rr:
                rr.append("_PHOTON")

        return {
            "reactants": rr,
            "products": pp,
            "tmin": t_min,
            "tmax": t_max,
            "rate": rate,
            "type": self._reaction_type(rate, rr),
            "line": line.strip(),
        }

    @staticmethod
    def _photo_placeholder(rate: str, line: str, nline: int, file) -> str:
        """Turn a photo rate into the placeholder ``prizmo_photo(lo, hi)``.

        Accepts ``photo(lo, hi)``, ``PHOTO, lo[, hi]`` and a bare ``PHOTO``; a
        missing band edge defaults to ``1e99``.  ``PrizmoParser.post_process``
        later replaces the placeholder with ``photorates(n, lo, hi)``.

        Raises
        ------
        ParserError
            If a band edge is not a number.
        """
        match = re.match(r"(?i:photo)\s*\((.*?)\)", rate)
        if match:
            args = [a.strip() for a in match.group(1).split(",") if a.strip()]
        else:
            args = [a.strip() for a in rate.split(",")[1:] if a.strip()]

        band: list[float] = []
        for arg in (args + ["1e99", "1e99"])[:2]:
            try:
                band.append(float(arg.replace("d", "e")))
            except ValueError:
                raise ParserError(
                    f"Invalid photo band edge '{arg}'", line, nline, file
                ) from None

        return f"prizmo_photo({band[0]!r}, {band[1]!r})"

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
