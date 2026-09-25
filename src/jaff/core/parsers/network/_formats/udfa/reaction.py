"""UDFA (UMIST) format: colon-delimited fixed-column reaction handler."""

import re

from ......errors import ParserError


class UdfaReaction:
    """UDFA colon-delimited reaction line handler."""

    name = "udfa"
    priority = 50
    is_reaction = True

    global_re = re.compile(r"^(?!\s*[!#@]).*:.*$")

    local_re = re.compile(
        r"^\s*\d+\s*:"
        r"\s*(?P<rtype>[^:]*?)\s*:"
        r"\s*(?P<reactants>(?:[^:]*:){2})"
        r"\s*(?P<products>(?:[^:]*:){4})"
        r"\s*(?P<flag>[^:]*)\s*:"
        r"\s*(?P<blocks>.*)$"
    )

    _MIN_FIELDS_PER_RANGE = 5

    _ZETA0 = 1.36e-17
    _ALBEDO = 0.5

    SPECIAL_MAP = {
        "CR": "_CR",
        "CRP": "_CRP",
        "CRPHOT": "_CRPHOT",
        "PHOTON": "_PHOTON",
    }

    def parse(self, line: str, nline: int, state: dict, file) -> list[dict]:
        """Parse a UDFA (UMIST/RATE22)-format reaction line into segment records.

        Extracts the reaction type, reactants, products, and the declared
        number of temperature ranges.  One parameter block (rate coefficients
        ``ka``, ``kb``, ``kc`` and temperature bounds) follows per range; each
        becomes its own segment record with shared participants and provenance.
        A rate expression is built per block based on the reaction type:
        cosmic-ray (``"CR"``), photo-desorption (``"PH"``), or standard
        Arrhenius.

        Returns
        -------
        list[dict]
            One reaction dict per declared temperature range, in file order.

        Raises
        ------
        ParserError
            If the line does not match the expected UDFA format, if the
            range-count flag is not a positive integer, or if the number of
            parameter blocks does not match the declared range count.
        """
        local = self.local_re.match(line)
        if not local:
            raise ParserError("Invalid UDFA reaction detected", line, nline, file)

        rtype: str = local.group("rtype")
        reactants: str = local.group("reactants")
        products: str = local.group("products")

        try:
            nranges = int(local.group("flag"))
        except ValueError:
            raise ParserError(
                f"Invalid temperature-range count: {local.group('flag')!r}",
                line,
                nline,
                file,
            )

        if nranges < 1:
            raise ParserError(
                f"Temperature-range count must be positive, got {nranges}",
                line,
                nline,
                file,
            )

        tokens = self._tokenize(local.group("blocks"))

        if len(tokens) < nranges * self._MIN_FIELDS_PER_RANGE or len(tokens) % nranges:
            raise ParserError(
                f"Declared {nranges} temperature range(s) but parameter "
                f"blocks do not align ({len(tokens)} fields present)",
                line,
                nline,
                file,
            )
        fields_per_range = len(tokens) // nranges

        rr = [
            self.SPECIAL_MAP.get(r.strip(), r.strip())
            for r in reactants.split(":")[:-1]
            if r.strip() != ""
        ]
        pp = [
            self.SPECIAL_MAP.get(p.strip(), p.strip())
            for p in products.split(":")[:-1]
            if p.strip() != ""
        ]

        if rtype == "PH" and "_PHOTON" not in rr:
            rr.append("_PHOTON")
        elif rtype in ("CR", "CP") and not any(
            cr in rr for cr in ("_CR", "_CRP", "_CRPHOT")
        ):
            rr.append("_CR")

        rtype_str = self._reaction_type(rtype, rr)
        string = line.strip()

        segments: list[dict] = []
        for i in range(nranges):
            base = i * fields_per_range
            try:
                ka = float(tokens[base])
                kb = float(tokens[base + 1])
                kc = float(tokens[base + 2])
                tmin = float(tokens[base + 3])
                tmax = float(tokens[base + 4])
            except ValueError as e:
                raise ParserError(
                    f"Invalid rate parameter in temperature range {i + 1}: {e}",
                    line,
                    nline,
                    file,
                )

            segments.append(
                {
                    "r": rr,
                    "p": pp,
                    "tmin": tmin if tmin > 0 else None,
                    "tmax": tmax if tmax < 41000.0 else None,
                    "rate": self._build_rate(rtype, ka, kb, kc),
                    "type": rtype_str,
                    "string": string,
                }
            )

        return segments

    @classmethod
    def _build_rate(cls, rtype: str, ka: float, kb: float, kc: float) -> str:
        """Build the rate expression string for one parameter block.

        RATE22 cosmic-ray reactions (see the format specification, eqs. 2 and 4)
        are scaled by ``crate / zeta0`` so that they vanish when the ionisation
        rate ``crate`` is zero:

        - ``CP`` (direct cosmic-ray ionisation): ``k = alpha * zeta/zeta0``.
        - ``CR`` (cosmic-ray-induced photoreaction):
          ``k = alpha * (T/300)**beta * gamma/(1-omega) * zeta/zeta0``.

        ``PH`` is a UV photoreaction; anything else is Kooij/Arrhenius.
        """
        rate_dict = {
            "CP": f"{ka / cls._ZETA0:.2e} * crate",
            "CR": (
                f"{ka * kc / (1.0 - cls._ALBEDO) / cls._ZETA0:.2e}"
                f" * (tgas / 3e2)**({kb:.2f}) * crate"
            ),
            "PH": f"{ka:.2e} * exp(-{kc:.2f} * av)",
        }
        if rtype in rate_dict:
            return rate_dict[rtype]

        rate = f"{ka:.2e}"
        if kb:
            rate = f"{rate} * (tgas / 3e2)**({kb:.2f})"
        if kc:
            rate = f"{rate} * exp(-{kc:.2f} / tgas)"

        return rate

    @staticmethod
    def _tokenize(blocks: str) -> list[str]:
        """Split the parameter-block string on ``:`` outside quoted metadata.

        Reference and note fields are double-quoted and may themselves contain
        commas or colons, so a naive ``split(":")`` would corrupt the field
        alignment.  A trailing empty token from the line's final ``:`` is
        dropped.
        """
        tokens: list[str] = []
        current: list[str] = []
        in_quote = False
        for ch in blocks:
            if ch == '"':
                in_quote = not in_quote
                current.append(ch)
            elif ch == ":" and not in_quote:
                tokens.append("".join(current).strip())
                current = []
            else:
                current.append(ch)

        tail = "".join(current).strip()
        if tail:
            tokens.append(tail)

        return tokens

    @staticmethod
    def _reaction_type(rtype: str, rr: list[str]) -> str:
        """Conclude the reaction type from the UDFA code and reactants.

        ``"CR"``/``"CP"`` = cosmic-ray, ``"PH"`` = photoprocess. Otherwise a
        reaction with three or more real (non-pseudo) reactants is three-body;
        else ``"unknown"``. Reactant-count classification is rate-independent,
        so it survives custom auxiliary-function rates.
        """
        agent = {"CR": "cosmic_ray", "CP": "cosmic_ray", "PH": "photo"}.get(rtype)
        if agent:
            return agent
        if sum(1 for r in rr if not r.startswith("_")) >= 3:
            return "3_body"

        return "unknown"
