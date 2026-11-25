# ABOUTME: Collapses per-reaction Leiden .h5 and NORAD .dat cross-section files into
# ABOUTME: combined HDF5 files (one group per reaction channel, eV + cm^2).
"""Collapse the per-reaction cross-section files into two combined HDF5 files.

The Leiden (``data/xsecs/leiden/*.h5``) and NORAD/OP (``data/xsecs/op/*.dat``)
folders each hold one file per reaction.  This utility merges each folder into a
single HDF5 file -- ``data/xsecs/leiden.h5`` and ``data/xsecs/op.h5`` -- with one
group per reaction (group name = the source file's serialized stem, e.g.
``"CH__C.H"``).

This script no longer reproduces the shipped data as-is: JAFF reads
``data/xsecs/leiden.hdf5`` and ``data/xsecs/norad.hdf5`` (not ``leiden.h5`` /
``op.h5``), whose group keys and ``reactants`` attrs all carry ``_PHOTON`` on
the reactant side (e.g. ``"CH._PHOTON__C.H"``) and whose root attrs include
``jaff_source_sha256``.

Schema
------
Root attrs: ``database`` ("leiden"/"norad"), ``description``, ``created``.

Each group has datasets, all co-sorted by **ascending photon energy**:

- ``photon_energy``   -- eV   (Leiden wavelengths in nm are converted via
  ``E = h c / lambda``).
- ``photodecay``      -- cm^2, the reaction's single decay-channel cross section
  (dissociation xsec for a dissociation reaction, ionisation xsec for an
  ionisation reaction).
- ``photoabsorption`` -- cm^2, optional, shared between a species' channels.

Each source Leiden file bundles both decay channels, so it is split into one
group per channel: the dissociation reaction keeps the serialized stem, while
the ionisation reaction is keyed ``<R>._PHOTON__<R+>.e-`` (the ``_PHOTON``
agent is injected on the reactant side).  NORAD files are photoionisation only.

Every dataset has a ``unit`` attr ("eV" or "cm2").

Group attrs:

- ``reactants`` / ``products`` -- string lists for the channel reaction.
- ``decay_type`` -- ``"dissociation"`` or ``"ionization"``.
- NORAD groups additionally: ``Z``, ``n_electrons``, ``charge``,
  ``accuracy_variant``, ``data_format``, ``source_url`` (from the ``.dat``
  ``#`` header).

Excluded: ``cross_section_properties.csv``, the per-file ``README``, the stray
``HI.dat``/``HeI.dat`` spectroscopic-name files, and ``op/raw/``.
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

import h5py
import numpy as np

from jaff.config import XSECS_DATA_DIR
from jaff.io import JaffLogger
from jaff.physics import constants

#: eV per erg (1 eV = 1.602176634e-12 erg).
EV_PER_ERG: float = 1.0 / 1.602176634e-12
#: cm per nm.
CM_PER_NM: float = 1e-7

_STR_DT = h5py.string_dtype(encoding="utf-8")

#: Lossless dataset compression (gzip level 4 + chunking). ~45% smaller files;
#: decompression is cheap per small reaction, costlier only for the few very
#: large datasets (read once and cache).
COMPRESSION_KW: dict = {"compression": "gzip", "compression_opts": 4, "chunks": True}


def split_reaction(stem: str) -> tuple[list[str], list[str]]:
    """Split a serialized stem into ``(reactants, products)`` string lists.

    ``"CH__C.H"`` -> ``(["CH"], ["C", "H"])``.
    """
    react, _, prod = stem.partition("__")
    return react.split("."), prod.split(".")


def wavelength_nm_to_eV(wavelength_nm: np.ndarray) -> np.ndarray:
    """Convert wavelength (nm) to photon energy (eV) via ``E = h c / lambda``."""
    lam_cm = wavelength_nm * CM_PER_NM
    energy_erg = constants.h.cgs.value * constants.c.cgs.value / lam_cm
    return energy_erg * EV_PER_ERG


def _ionize(species: str) -> str:
    """Singly-ionise a species name: ``X`` -> ``X+``, anion ``X-`` -> ``X``."""
    return species[:-1] if species.endswith("-") else species + "+"


def _has_signal(arr: np.ndarray) -> bool:
    """True if the array has at least one non-zero, all-finite value."""
    return bool(np.any(arr != 0)) and bool(np.all(np.isfinite(arr)))


def _write_channel(
    h5: h5py.File,
    key: str,
    reactants: list[str],
    products: list[str],
    decay_type: str,
    energy: np.ndarray,
    photodecay: np.ndarray,
    photoabsorption: np.ndarray | None,
) -> None:
    """Create one per-channel reaction group with a ``photodecay`` dataset."""
    grp = h5.create_group(key)
    grp.attrs["reactants"] = np.array(reactants, dtype=_STR_DT)
    grp.attrs["products"] = np.array(products, dtype=_STR_DT)
    grp.attrs["decay_type"] = decay_type
    e = grp.create_dataset("photon_energy", data=energy, **COMPRESSION_KW)
    e.attrs["unit"] = "eV"
    d = grp.create_dataset("photodecay", data=photodecay, **COMPRESSION_KW)
    d.attrs["unit"] = "cm2"
    if photoabsorption is not None:
        a = grp.create_dataset("photoabsorption", data=photoabsorption, **COMPRESSION_KW)
        a.attrs["unit"] = "cm2"


def collapse_leiden(leiden_dir: Path, out_path: Path, logger) -> int:
    """Merge Leiden per-reaction ``.h5`` files, split by decay channel."""
    files = sorted(leiden_dir.glob("*.h5"))
    emitted = 0
    ionis_seen: set[str] = set()
    with h5py.File(out_path, "w") as h5:
        h5.attrs["database"] = "leiden"
        h5.attrs["description"] = (
            "Leiden photo cross sections, one group per reaction (split by decay "
            "channel); photon_energy in eV, photodecay/photoabsorption in cm^2."
        )
        h5.attrs["created"] = date.today().isoformat()

        for f in files:
            stem = f.stem
            reactants, products = split_reaction(stem)
            reactant = reactants[0]
            with h5py.File(f, "r") as src:
                energy_ev = wavelength_nm_to_eV(src["wavelength"][:].astype(float))
                order = np.argsort(energy_ev)  # ascending energy
                energy = energy_ev[order]

                photoabs = None
                if "photoabsorption" in src:
                    pa = src["photoabsorption"][:].astype(float)[order]
                    photoabs = pa if _has_signal(pa) else None

                # Dissociation channel -> serialized stem.
                if "photodissociation" in src:
                    pd_xs = src["photodissociation"][:].astype(float)[order]
                    if _has_signal(pd_xs):
                        _write_channel(
                            h5,
                            stem,
                            reactants,
                            products,
                            "dissociation",
                            energy,
                            pd_xs,
                            photoabs,
                        )
                        emitted += 1

                # Ionisation channel -> constructed key (deduped per reactant).
                if "photoionisation" in src:
                    pi_xs = src["photoionisation"][:].astype(float)[order]
                    ion_products = sorted([_ionize(reactant), "e-"])
                    ion_key = (
                        f"{'.'.join(sorted([reactant, '_PHOTON']))}"
                        f"__{'.'.join(ion_products)}"
                    )
                    if _has_signal(pi_xs) and ion_key not in ionis_seen:
                        ionis_seen.add(ion_key)
                        _write_channel(
                            h5,
                            ion_key,
                            reactants,
                            ion_products,
                            "ionization",
                            energy,
                            pi_xs,
                            photoabs,
                        )
                        emitted += 1
    logger.info(f"Wrote {emitted} Leiden channel groups to {out_path}")
    return emitted


def _parse_op_header(lines: list[str]) -> dict:
    """Extract NORAD provenance attrs from the ``#`` header of a ``.dat`` file."""
    attrs: dict = {}
    for ln in lines:
        if ln.startswith("# source:"):
            attrs["source_url"] = ln.split(":", 1)[1].strip()
        elif ln.startswith("# variant="):
            for tok in ln[1:].split():
                key, _, val = tok.partition("=")
                if key in ("variant", "format", "Z", "NE", "charge"):
                    attrs[key] = val
    return attrs


def collapse_op(op_dir: Path, out_path: Path, logger) -> int:
    """Merge all NORAD/OP per-reaction serialized ``.dat`` files into one HDF5."""
    files = sorted(p for p in op_dir.glob("*.dat") if "__" in p.stem)
    with h5py.File(out_path, "w") as h5:
        h5.attrs["database"] = "norad"
        h5.attrs["description"] = (
            "NORAD (Nahar/OSU) ground-state photoionisation cross sections, one "
            "group per ion; photon_energy in eV, cross section in cm^2."
        )
        h5.attrs["created"] = date.today().isoformat()

        for f in files:
            stem = f.stem
            lines = f.read_text(encoding="utf-8").splitlines()
            hdr = _parse_op_header([ln for ln in lines if ln.startswith("#")])
            data = np.loadtxt(f, comments="#")
            energy_ev, xsec_cm2 = data[:, 0], data[:, 1]
            order = np.argsort(energy_ev)  # ascending energy

            grp = h5.create_group(stem)
            reactants, products = split_reaction(stem)
            grp.attrs["reactants"] = np.array(reactants, dtype=_STR_DT)
            grp.attrs["products"] = np.array(products, dtype=_STR_DT)
            grp.attrs["decay_type"] = "ionization"  # NORAD is photoionisation only
            if "Z" in hdr:
                grp.attrs["Z"] = int(hdr["Z"])
                grp.attrs["n_electrons"] = int(hdr["NE"])
                grp.attrs["charge"] = int(hdr["charge"])
                grp.attrs["accuracy_variant"] = hdr["variant"]
                grp.attrs["data_format"] = hdr["format"]
                grp.attrs["source_url"] = hdr.get("source_url", "")

            e_ds = grp.create_dataset(
                "photon_energy", data=energy_ev[order], **COMPRESSION_KW
            )
            e_ds.attrs["unit"] = "eV"
            x_ds = grp.create_dataset(
                "photodecay", data=xsec_cm2[order], **COMPRESSION_KW
            )
            x_ds.attrs["unit"] = "cm2"
    logger.info(f"Wrote {len(files)} NORAD groups to {out_path}")
    return len(files)


def main() -> None:
    """Build ``leiden.h5`` and ``op.h5`` in ``data/xsecs/`` (see the module notes)."""
    logger = JaffLogger().get_logger()
    xsecs = XSECS_DATA_DIR
    collapse_leiden(xsecs / "leiden", xsecs / "leiden.h5", logger)
    collapse_op(xsecs / "op", xsecs / "op.h5", logger)


if __name__ == "__main__":
    main()
