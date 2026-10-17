"""KIDA format: fixed-width column reaction handler."""

import re
from pathlib import Path

from ......errors import ParserError

KIDA_OPEN_TMAX = 9999.0


class KidaReaction:
    """KIDA fixed-width reaction line handler."""

    name = "kida"
    priority = 80
    is_reaction = True

    global_re = re.compile(r"^(?!\s*[!#@]).{34}.{57}")

    local_re = re.compile(
        r"^(?P<reactants>.{34})"
        r"(?P<products>.{57})"
        r"\s*(?P<ka>[^\s]+)"
        r"\s*(?P<kb>[^\s]+)"
        r"\s*(?P<kc>[^\s]+)"
        r"\s*[^\s]+\s*[^\s]+\s*[^\s]+"
        r"\s*(?P<itype>[^\s]+)"
        r"\s*(?P<tmin>[^\s]+)"
        r"\s*(?P<tmax>[^\s]+)"
        r"\s*(?P<formula>[^\s]+)"
        r".*$"
    )

    SPECIAL_MAP = {
        "CR": "_CR",
        "CRP": "_CRP",
        "CRPHOT": "_CRPHOT",
        "Photon": "_PHOTON",
        "PHOTON": "_PHOTON",
    }

    #: KIDA ``itype`` -> JAFF reaction type (see the network-file header).
    ITYPE_TYPE = {
        1: "cosmic_ray",
        2: "cosmic_ray",
        3: "photo",
        4: "two_body",
        5: "two_body",
        6: "radiative_association",
        7: "associative_detachment",
        8: "recombination",
        9: "grain",
        10: "special",
    }

    #: KIDA ``itype`` -> agent/catalyst pseudo-species injected as a reactant.
    ITYPE_AGENT = {1: "_CR", 2: "_CRP", 3: "_PHOTON", 9: "_GRAIN"}

    def parse(self, line: str, nline: int, state: dict, file: Path) -> dict:
        """Parse a KIDA-format reaction line into its reaction fields.

        Extracts reactants, products, rate parameters (``ka``, ``kb``, ``kc``),
        temperature bounds, and formula index from the fixed-width KIDA column
        format.  Rate expressions are selected from a formula dictionary keyed
        by the integer formula index (1-5).  A ``Tmax`` equal to
        :data:`KIDA_OPEN_TMAX` is returned unchanged;
        :meth:`KidaParser.process` decides whether it is open.

        Raises
        ------
        ParserError
            If the line does not match the expected KIDA format.
        """
        local = self.local_re.match(line)
        if not local:
            raise ParserError("Invalid KIDA reaction detected", line, nline, file)

        reactants: str = local.group("reactants")
        products: str = local.group("products")
        ka: float = float(local.group("ka"))
        kb: float = float(local.group("kb"))
        kc: float = float(local.group("kc"))
        tmin: float = float(local.group("tmin"))
        tmax: float = float(local.group("tmax"))
        formula: int = int(local.group("formula"))
        itype: int = int(local.group("itype"))

        t_min = tmin if tmin > 0 else None

        rr = [r.strip() for r in reactants.split() if r != "+"]
        pp = [p.strip() for p in products.split() if p != "+"]
        rates_dict = {
            1: (
                f"{ka:.2e} * crate"
                if "CRP" not in rr
                else f"{ka:.2e} * crate * 2.0 * n_H2 / n_H_nuc"
            ),
            2: f"{ka:.2e} * chi * exp(-{kc:.2e} * av)",
            3: f"{ka:.2e}"
            + (f" * (tgas / 300) ** ({kb:.2e})" if kb != 0.0 else "")
            + (f" * exp(-{kc:.2f} / tgas)" if kc != 0.0 else ""),
            4: f"{ka * kb:.2e} * (0.62 + 0.4767 * {kc:2e} * sqrt(300 / tgas))",
            5: f"{ka * kb:.2e} * (1 + 0.0967 * {kc:.2e} * sqrt(300 / tgas) + {kc**2:.2e} * 3e2 / 10.526 / tgas)",
        }
        rate = rates_dict.get(formula, "0.0")

        rr = [self.SPECIAL_MAP.get(r, r) for r in rr]
        pp = [self.SPECIAL_MAP.get(p, p) for p in pp]

        rtype = self._reaction_type(itype, rr)

        agent = self.ITYPE_AGENT.get(itype)
        if agent is not None and not (set(rr) & {agent, "_CR", "_CRP", "_CRPHOT"}):
            rr.append(agent)

        return {
            "r": rr,
            "p": pp,
            "tmin": t_min,
            "tmax": tmax,
            "rate": rate,
            "type": rtype,
            "string": line.strip(),
        }

    @staticmethod
    def _reaction_type(itype: int, rr: list[str]) -> str:
        """Conclude the reaction type from the KIDA ``itype`` column.

        The KIDA ``itype`` integer is the authoritative physical classification
        (defined in the network-file header):

        ===== =======================================================
        itype meaning
        ===== =======================================================
        1     direct cosmic-ray ionization/dissociation
        2     cosmic-ray-induced (secondary UV) ionization/dissociation
        3     FUV photo-dissociation/ionization
        4     neutral-neutral / ion-neutral
        5     exchange
        6     radiative association
        7     associative detachment
        8     recombination
        9     grain-assisted (H2 formation, recombination)
        10    other special
        ===== =======================================================

        A reaction with three or more real (non-pseudo) reactants is labelled
        ``"3_body"`` regardless of ``itype``, reflecting the flux arity.

        Parameters
        ----------
        itype : int
            KIDA reaction-type index (1-10).
        rr : list[str]
            Reactant names, after ``SPECIAL_MAP`` normalization but before agent
            pseudo-species injection.
        """
        if sum(1 for r in rr if not r.startswith("_")) >= 3:
            return "3_body"

        return KidaReaction.ITYPE_TYPE.get(itype, "unknown")
