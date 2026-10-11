"""Download NORAD (Nahar/OSU) ground-state photoionisation cross sections.

For every ion of elements Z = 1..26 this fetches the *ground-state*
photoionisation cross section from the NORAD-Atomic-Data archive at Ohio State
(S. N. Nahar), reformats it, and writes one text file per ion into
``data/xsecs/op/`` using the serialized reaction naming convention
``<ion>__<ion+>.e-.dat`` (e.g. ``H__H+.e-.dat``, ``He+__He++.e-.dat``).

Source URL scheme
-----------------
``<base>/<sym><stage>/<sym><stage>.px.gd[.fs|.ls].txt`` where ``sym`` is the
lower-case element symbol and ``stage`` is the spectroscopic ionisation stage
(1 = neutral).  Three accuracy variants may exist per ion; the most accurate
available is used, in priority order:

1. ``.px.gd.fs.txt`` -- fine-structure resolved (relativistic, radiation
   damped).  Most accurate.
2. ``.px.gd.ls.txt`` -- LS-coupling.
3. ``.px.gd.txt``    -- plain ground state (may itself be hydrogenic 3-column
   or LS-format 2-column).

File formats
------------
The three on-disk formats differ in both header and column layout; the format
is detected from the first record after the long dashed separator:

- **hydrogenic** (header ``Z nc isc ilc ipc``, 5 ints): per ``(2S+1) L pi n``
  block, ``BE N`` line, then ``N`` rows of ``Ephoton(Ry) Epe(Ry) PX(Mb)``.
  Cross section is column 3.
- **LS** (header ``zz nn P``): per ``2S+1 L ip ns`` block, ``nr ntot`` line,
  ``BE ac`` line, then rows of ``E(Ry) sig(Mb)``.  Cross section is column 2.
- **FS** (header ``zz nn ntg``, 3 ints): ``ntg`` core energies, then per
  ``is 2j ip ns`` block, ``BE ntot`` line, ``ac`` line, then ``ntot`` rows of
  ``E(Ry) Epe(Ry) signd(Mb) sig(Mb)``.  Cross section is column 4 (``sig``,
  radiation-damped final; column 3 ``signd`` is ignored per the data docs).

Ground-state combination
-------------------------
A ``px.gd`` file holds only ground-state blocks, but the ground term may be
split into several fine-structure (or LS) component blocks, each carrying its
own resonance structure.  To retain *every* resonance peak in a single curve,
all component blocks are combined by a statistical-weight (``g``) average onto
the sorted union of their photon-energy grids::

    sigma(E) = sum_i g_i * sigma_i(E) / sum_i g_i

with each component linearly interpolated onto the common grid and taken as
zero below its own threshold.

Output
------
Two columns: photon energy in eV (``E_Ry * 13.605693``) and cross section in
cm^2 (``sigma_Mb * 1e-18``), preceded by a ``#`` comment header recording
provenance (variant, Z, NE, charge, component terms, threshold).
"""

from __future__ import annotations

import urllib.error
import urllib.request
from pathlib import Path

import numpy as np

from jaff.config import XSECS_DATA_DIR
from jaff.io import JaffLogger

#: Rydberg energy in eV (CODATA), for Ry -> eV photon-energy conversion.
RY_TO_EV: float = 13.605693
#: Megabarn in cm^2, for Mb -> cm^2 cross-section conversion.
MB_TO_CM2: float = 1e-18

#: NORAD-Atomic-Data archive root (current OSU NORAD site).
BASE_URL: str = "https://norad.astronomy.osu.edu"

#: Atomic numbers to fetch: everything through Argon (1..18), plus Ca and Fe.
Z_SET: list[int] = list(range(1, 19)) + [20, 26]

#: Element symbols indexed by atomic number (index 0 unused).
ELEMENTS: list[str] = [
    "",
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
]

#: Accuracy variants in descending priority (URL infix, label).
VARIANTS: list[tuple[str, str]] = [
    ("px.gd.fs", "fs"),
    ("px.gd.ls", "ls"),
    ("px.gd", "plain"),
]

#: Per-ion non-standard filename infixes (tried before VARIANTS).  H I on the
#: current site has no ``px.gd`` file; its data lives in ``h1.px.1-10.txt``
#: (all nl levels n=1..10 -- the ground 1s block is selected at parse time).
SPECIAL_INFIXES: dict[str, list[str]] = {"h1": ["px.1-10"]}


def _is_float_token(tok: str) -> bool:
    """Return ``True`` if *tok* is a floating-point literal (has ``.``/``e``).

    NORAD records use plain integers for headers/counts and decimal or
    E-notation for energies and cross sections, so this cleanly separates the
    two token classes.
    """
    return ("." in tok) or ("e" in tok.lower())


def serialized_name(symbol: str, charge: int) -> str:
    """Build the serialized reaction key for photoionisation of an ion.

    Parameters
    ----------
    symbol : str
        Element symbol (e.g. ``"He"``).
    charge : int
        Charge of the ion being ionised (0 for neutral).

    Returns
    -------
    str
        e.g. ``charge=0`` -> ``"He__He+.e-"``, ``charge=1`` ->
        ``"He+__He++.e-"``.
    """
    ion = symbol + "+" * charge
    product = symbol + "+" * (charge + 1)
    return (
        f"{'.'.join(sorted([ion, '_PHOTON']))}__{'.'.join(sorted([product, 'e-']))}"
    )


def fetch(url: str) -> str | None:
    """Download *url*; return its text, or ``None`` on 404/other HTTP error."""
    req = urllib.request.Request(url, headers={"User-Agent": "jaff-norad-fetch"})
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            return resp.read().decode("utf-8", errors="replace")
    except urllib.error.HTTPError:
        return None
    except urllib.error.URLError:
        return None


def _data_rows(text: str) -> list[list[str]]:
    """Return tokenised non-empty lines of the data section.

    The NORAD files carry a free-text explanation terminated by a long dashed
    separator; everything after the final such line is numeric data.  Inline
    ``:label`` annotations (used by the hydrogenic files to tag header rows)
    are stripped.
    """
    lines = text.splitlines()
    sep_idx = max(
        (i for i, ln in enumerate(lines) if ln.strip().startswith("-" * 20)),
        default=-1,
    )
    rows: list[list[str]] = []
    for ln in lines[sep_idx + 1 :]:
        toks = ln.split(":", 1)[0].split()
        if toks:
            rows.append(toks)
    return rows


def _all_int(row: list[str]) -> bool:
    """Return ``True`` if every token in *row* is a plain integer literal."""
    return all(not _is_float_token(t) for t in row)


def _is_block_header(row: list[str]) -> bool:
    """Return ``True`` if *row* is a 4-integer level/term header line.

    Both LS (``2S+1 L ip ns``) and FS (``is 2j ip ns``) blocks begin with a
    line of exactly four integers; data rows always contain decimals/E-notation
    so they never match.
    """
    return len(row) == 4 and _all_int(row)


def _is_numeric(row: list[str]) -> bool:
    """Return ``True`` if every token in *row* parses as a float.

    Some files append a duplicate table whose rows carry spectroscopic labels
    (e.g. ``0.000000 1 F``); such a non-numeric row terminates the data block.
    """
    try:
        for t in row:
            float(t)
        return True
    except ValueError:
        return False


def parse(text: str) -> tuple[str, list[tuple[float, np.ndarray, np.ndarray]]]:
    """Detect the format and parse all ground-state component blocks.

    Line-oriented: one data row per line (photon energy is the first column,
    cross section the last in every NORAD variant); headers, counts and
    binding-energy lines occupy their own lines.  A ``0 0 0 0`` header line
    terminates the data.

    Format detection from the first record after the separator:

    - 5 integers (``Z nc isc ilc ipc``)        -> hydrogenic.
    - 3 tokens, last alphabetic (``zz nn P``)   -> LS.
    - 3 integers (``zz nn ntg``) + core energies-> FS.
    - 4 integers (block header, no global line) -> LS (e.g. Ar XIII).

    Returns
    -------
    fmt : str
        ``"hydrogenic"``, ``"ls"``, or ``"fs"``.
    blocks : list of (g, energy_Ry, xsec_Mb)
        One entry per ground-state component (fine-structure / term block).
    """
    rows = _data_rows(text)
    r0 = rows[0]

    if len(r0) == 5 and _all_int(r0):
        fmt, start = "hydrogenic", 1
    elif len(r0) == 3 and r0[2].isalpha():
        fmt, start = "ls", 1
    elif len(r0) == 3 and _all_int(r0):
        # FS: skip the ``zz nn ntg`` line and the following ntg core energies.
        fmt = "fs"
        ntg = int(r0[2])
        start, consumed = 1, 0
        while consumed < ntg:
            consumed += len(rows[start])
            start += 1
    elif _is_block_header(r0):
        fmt, start = "ls", 0  # LS data with no global header line
    else:
        raise ValueError(f"Unrecognised NORAD header: {r0}")

    blocks: list[tuple[float, np.ndarray, np.ndarray]] = []
    i = start
    n = len(rows)
    while i < n:
        if not _is_block_header(rows[i]):
            i += 1
            continue
        hdr = [int(x) for x in rows[i]]
        if hdr == [0, 0, 0, 0]:
            break  # end-of-data terminator
        # Statistical weight: FS header is (is, 2j, ip, ns) -> g = 2j + 1;
        # LS/hydrogenic header is (2S+1, L, ip, ns) -> g = (2S+1)(2L+1).
        g = (hdr[1] + 1) if fmt == "fs" else hdr[0] * (2 * hdr[1] + 1)
        i += 1

        # Skip per-block metadata lines preceding the data rows.
        if fmt == "ls":
            while i < n and _all_int(rows[i]):  # nr ntot [extra] count line(s)
                i += 1
            i += 1  # BE ac line
        elif fmt == "fs":
            i += 2  # BE ntot line, then ac line
        else:  # hydrogenic
            i += 1  # BE N line

        energies: list[float] = []
        xsecs: list[float] = []
        while i < n and not _is_block_header(rows[i]) and _is_numeric(rows[i]):
            row = rows[i]
            energies.append(float(row[0]))
            xsecs.append(float(row[-1]))  # cross section is always the last column
            i += 1
        # Keep only genuine photon-energy blocks (first column = photon energy,
        # starting at the threshold > 0).  Duplicate photoelectron-energy tables
        # start at E = 0 and are discarded.
        if energies and min(energies) > 1e-6:
            blocks.append((float(g), np.array(energies), np.array(xsecs)))

    return fmt, blocks


def keep_ground(
    blocks: list[tuple[float, np.ndarray, np.ndarray]],
) -> list[tuple[float, np.ndarray, np.ndarray]]:
    """Retain only the most-bound (ground-state) component blocks.

    Some files list excited levels too (e.g. ``h1.px.1-10.txt`` holds n=1..10).
    The ground state has the highest ionisation threshold; fine-structure
    components of a single ground term share that threshold to within a tiny
    splitting.  This keeps every block whose threshold is within 5% of the
    maximum and discards the (much lower-threshold) excited levels.
    """
    if len(blocks) <= 1:
        return blocks
    thresholds = [float(e.min()) for _, e, _ in blocks]
    emax = max(thresholds)
    return [b for b, t in zip(blocks, thresholds) if t >= 0.95 * emax]


def combine(
    blocks: list[tuple[float, np.ndarray, np.ndarray]],
) -> tuple[np.ndarray, np.ndarray]:
    """Statistical-weight average component blocks onto a unified grid.

    Parameters
    ----------
    blocks : list of (g, energy_Ry, xsec_Mb)

    Returns
    -------
    energy_Ry : numpy.ndarray
        Sorted union of all component photon-energy grids.
    xsec_Mb : numpy.ndarray
        ``sum_i g_i sigma_i(E) / sum_i g_i`` on that grid (zero below each
        component's threshold).

    Some raw files repeat a photon energy with a different cross section
    (e.g. ``fe21``); only the first occurrence in file order is kept.
    """
    if not blocks:
        return np.array([]), np.array([])

    blocks = [(g, *_dedupe_first(e, x)) for g, e, x in blocks]

    if len(blocks) == 1:
        _, e, x = blocks[0]
        return e, x

    grid = np.unique(np.concatenate([e for _, e, _ in blocks]))
    num = np.zeros_like(grid)
    gtot = 0.0
    for g, e, x in blocks:
        # left/right=0 so a component contributes nothing outside its range.
        num += g * np.interp(grid, e, x, left=0.0, right=0.0)
        gtot += g
    return grid, num / gtot


def _dedupe_first(e: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Sort by energy, keeping the first-listed point of each repeated energy."""
    e_unique, first = np.unique(e, return_index=True)
    return e_unique, x[first]


#: Element symbol -> atomic number, derived from ELEMENTS.
SYMBOL_TO_Z: dict[str, int] = {sym: z for z, sym in enumerate(ELEMENTS) if sym}


def download_raw(raw_dir: Path, logger) -> int:
    """Download the most-accurate raw NORAD file per ion into *raw_dir*.

    For each ion of the elements in :data:`Z_SET`, the variants are tried in
    descending accuracy (:data:`VARIANTS`); the first that exists is saved
    verbatim under its original name ``<stem>.<infix>.txt`` (e.g.
    ``o1.px.gd.txt``, ``c3.px.gd.fs.txt``).  Returns the number of files saved.
    """
    raw_dir.mkdir(parents=True, exist_ok=True)
    saved = 0
    for z in Z_SET:
        symbol = ELEMENTS[z]
        for stage in range(1, z + 1):
            stem = f"{symbol.lower()}{stage}"
            infixes = SPECIAL_INFIXES.get(stem, []) + [inf for inf, _ in VARIANTS]
            for infix in infixes:
                text = fetch(f"{BASE_URL}/{stem}/{stem}.{infix}.txt")
                if text is None:
                    continue
                (raw_dir / f"{stem}.{infix}.txt").write_text(text)
                saved += 1
                logger.info(f"downloaded {stem}.{infix}.txt")
                break  # most-accurate variant found; stop
    logger.info(f"Downloaded {saved} raw NORAD files to {raw_dir}")
    return saved


def _stem_variant(filename: str) -> tuple[str, str, str]:
    """Split a raw filename into ``(stem, infix, label)``.

    ``"c3.px.gd.fs.txt"`` -> ``("c3", "px.gd.fs", "fs")``.
    """
    base = filename[:-4] if filename.endswith(".txt") else filename
    stem, dot, infix = base.partition(".")
    label = next((lbl for inf, lbl in VARIANTS if inf == infix), infix)
    return stem, infix, label


def parse_local(raw_dir: Path, outdir: Path, logger) -> int:
    """Parse every raw file in *raw_dir* into a serialized ``.dat`` in *outdir*.

    Reading from the local raw store (rather than the network) makes parsing
    deterministic and reproducible.  Returns the number of files written.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    written = 0
    for raw_path in sorted(raw_dir.glob("*.txt")):
        stem, infix, label = _stem_variant(raw_path.name)
        symbol = "".join(c for c in stem if c.isalpha()).capitalize()
        stage = int("".join(c for c in stem if c.isdigit()))
        z = SYMBOL_TO_Z[symbol]
        charge = stage - 1
        ne = z - charge

        text = raw_path.read_text()
        try:
            fmt, blocks = parse(text)
            blocks = keep_ground(blocks)
            energy_ry, xsec_mb = combine(blocks)
        except Exception as exc:  # noqa: BLE001 - skip & report malformed files
            logger.warning(f"Failed to parse {raw_path.name}: {exc}")
            continue

        if energy_ry.size == 0:
            logger.warning(f"No data points for {raw_path.name}")
            continue

        energy_ev = energy_ry * RY_TO_EV
        xsec_cm2 = xsec_mb * MB_TO_CM2

        ser = serialized_name(symbol, charge)
        ion = symbol + "+" * charge
        with open(outdir / f"{ser}.dat", "w") as f:
            f.write("# NORAD (Nahar/OSU) ground-state photoionisation cross section\n")
            f.write(f"# reaction: {ion} -> {ion}+ + e-   serialized: {ser}\n")
            f.write(
                f"# source: {BASE_URL}/{stem}/{raw_path.name}\n"
                f"# variant={label} format={fmt} Z={z} NE={ne} charge={charge} "
                f"n_components={len(blocks)} n_points={energy_ry.size}\n"
            )
            f.write("# E(eV)         xsec(cm2)\n")
            for ev, cm2 in zip(energy_ev, xsec_cm2):
                # 10 significant figures: resonance points sit ~1e-6 Ry apart, so
                # 7 figures (.6E) merged neighbouring energies after the eV conversion.
                f.write(f"{ev:.9E}  {cm2:.6E}\n")
        written += 1
        logger.info(f"{ser}.dat  <- {raw_path.name} ({fmt}, {energy_ry.size} pts)")

    logger.info(f"Wrote {written} NORAD ground-state cross-section files to {outdir}")
    return written


def main(download: bool = True, do_parse: bool = True) -> None:
    """Download raw NORAD files and/or parse them into serialized ``.dat`` files.

    Parameters
    ----------
    download : bool
        Fetch raw files from the network into ``op/raw/``.
    do_parse : bool
        Parse the local raw store into ``op/<serialized>.dat`` files.
    """
    logger = JaffLogger().get_logger()
    op_dir = XSECS_DATA_DIR / "op"
    raw_dir = op_dir / "raw"

    if download:
        download_raw(raw_dir, logger)
    if do_parse:
        parse_local(raw_dir, op_dir, logger)


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--no-download",
        action="store_true",
        help="skip fetching; parse the existing op/raw/ store",
    )
    ap.add_argument(
        "--no-parse", action="store_true", help="only download raw files, do not parse"
    )
    args = ap.parse_args()
    main(download=not args.no_download, do_parse=not args.no_parse)
