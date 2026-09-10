"""Collapse the background-radiation ``.dat`` files into one combined HDF5 file.

``data/background_radiation/*.dat`` holds one tabulated radiation field per file
(blackbodies, ISRF variants, solar, TW-Hya, ...).  This utility merges them into
a single ``data/background_radiation.hdf5`` with one group per file (group name =
the file stem, e.g. ``"bb_10000"``, ``"solar"``, ``"tw_hydra"``).

Source ``.dat`` format
----------------------
A ``#`` header followed by two whitespace-separated columns::

    # Numerical radiation field:  blackbody_10000K
    # wavelength: nm
    # intensity: photons.s-1.cm-2.nm-1
    # wavelength intensity
    92.00000 2.88992e+04
    ...

Schema
------
Root attrs: ``description``, ``created``.

Each group has two datasets, co-sorted by **ascending wavelength**:

- ``wavelength`` -- carries a ``unit`` attr from the ``# wavelength:`` header.
- ``intensity``  -- carries a ``unit`` attr from the ``# intensity:`` header.

Group attrs:

- ``radiation_field`` -- the field name from the ``# Numerical radiation field:``
  header line.
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

import h5py
import numpy as np

from jaff.config import DATA_DIR
from jaff.io import JaffLogger

#: Lossless dataset compression (gzip level 4 + chunking).
COMPRESSION_KW: dict = {"compression": "gzip", "compression_opts": 4, "chunks": True}


def _parse_header(lines: list[str]) -> dict[str, str]:
    """Extract metadata from the ``#`` header lines of a radiation ``.dat`` file.

    Returns a dict with ``radiation_field``, ``wavelength_unit`` and
    ``intensity_unit`` keys (missing entries default to an empty string).
    """
    attrs = {"radiation_field": "", "wavelength_unit": "", "intensity_unit": ""}
    for ln in lines:
        body = ln.lstrip("#").strip()
        if body.lower().startswith("numerical radiation field:"):
            attrs["radiation_field"] = body.split(":", 1)[1].strip()
        elif body.lower().startswith("wavelength:"):
            attrs["wavelength_unit"] = body.split(":", 1)[1].strip()
        elif body.lower().startswith("intensity:"):
            attrs["intensity_unit"] = body.split(":", 1)[1].strip()
    return attrs


def collapse_background_radiation(rad_dir: Path, out_path: Path, logger) -> int:
    """Merge all background-radiation ``.dat`` files into one HDF5 file.

    One group per ``.dat`` file, keyed by the file stem; returns the number of
    groups written.
    """
    files = sorted(rad_dir.glob("*.dat"))
    with h5py.File(out_path, "w") as h5:
        h5.attrs["description"] = (
            "Background radiation fields, one group per field (group name = file "
            "stem); wavelength and intensity datasets carry their unit attrs."
        )
        h5.attrs["created"] = date.today().isoformat()

        for f in files:
            lines = f.read_text().splitlines()
            hdr = _parse_header([ln for ln in lines if ln.startswith("#")])
            # genfromtxt + invalid_raise=False drops malformed rows (e.g. a
            # truncated final line) instead of aborting the whole file.
            data = np.genfromtxt(f, comments="#", usecols=(0, 1), invalid_raise=False)
            data = data[np.all(np.isfinite(data), axis=1)]
            wavelength, intensity = data[:, 0], data[:, 1]
            order = np.argsort(wavelength)  # ascending wavelength

            grp = h5.create_group(f.stem)
            grp.attrs["radiation_field"] = hdr["radiation_field"]

            w_ds = grp.create_dataset(
                "wavelength", data=wavelength[order], **COMPRESSION_KW
            )
            w_ds.attrs["unit"] = hdr["wavelength_unit"]
            i_ds = grp.create_dataset(
                "intensity", data=intensity[order], **COMPRESSION_KW
            )
            i_ds.attrs["unit"] = hdr["intensity_unit"]

    logger.info(f"Wrote {len(files)} radiation-field groups to {out_path}")
    return len(files)


def main() -> None:
    """Build ``background_radiation.hdf5`` in ``data/``."""
    logger = JaffLogger().get_logger()
    collapse_background_radiation(
        DATA_DIR / "background_radiation",
        DATA_DIR / "background_radiation.hdf5",
        logger,
    )


if __name__ == "__main__":
    main()
