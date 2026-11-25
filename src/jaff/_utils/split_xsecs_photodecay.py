# ABOUTME: One-shot migration splitting bundled Leiden/NORAD cross-section HDF5
# ABOUTME: groups into per-decay-channel groups with a single photodecay dataset.
"""Split the bundled Leiden/NORAD cross-section HDF5 files into per-channel
groups carrying a single ``photodecay`` dataset.

Background
----------
The collapsed ``data/xsecs/leiden.hdf5`` bundled **both** decay channels
(``photodissociation`` + ``photoionization``) under one dissociation-keyed
group (e.g. ``CH__C.H``).  Dissociation and ionisation are physically distinct
reactions with different products and rates, so a single bundled row made the
molecular-ionisation reaction invisible to the network.  This migration splits
each bundled group into one reaction per decay channel:

- ``<R>__<dissoc products>``  -- ``photodecay`` = photodissociation xsec
- ``<R>._PHOTON__<R+>.e-``    -- ``photodecay`` = photoionisation xsec

emitted only for channels that carry signal.  This was a one-shot migration of
the older bundled layout; the shipped files are already split, and their
dissociation keys have since gained ``_PHOTON`` too
(``<R>._PHOTON__<dissoc products>``).  ``photoabsorption`` (shared
between a species' channels) and ``photon_energy`` are copied into each group.

``norad.hdf5`` is single-channel already (photoionisation); its
``photoionization`` dataset is renamed to ``photodecay``.

New schema
----------
::

    <reaction_key>/
        attrs: reactants, products, decay_type ("dissociation" | "ionization")
        photon_energy    (eV)
        photoabsorption  (cm^2, optional)
        photodecay       (cm^2)

The operation rewrites the files in place and is guarded: if a file already
carries ``photodecay`` datasets it is left untouched.
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

import h5py
import numpy as np

from jaff.config import XSECS_DATA_DIR
from jaff.io import JaffLogger

_STR_DT = h5py.string_dtype(encoding="utf-8")
COMPRESSION_KW: dict = {"compression": "gzip", "compression_opts": 4, "chunks": True}

#: Package data directory holding the collapsed xsec files.
XSECS_DIR: Path = XSECS_DATA_DIR


def _ionize(species: str) -> str:
    """Singly-ionise a species name: ``X`` -> ``X+``, ``X+`` -> ``X++``,
    anion ``X-`` -> neutral ``X`` (photodetachment)."""
    if species.endswith("-"):
        return species[:-1]
    return species + "+"


def _ionis_key(reactant: str) -> tuple[str, list[str]]:
    """Return ``(serialized_key, products)`` for the photoionisation channel."""
    products = sorted([_ionize(reactant), "e-"])
    return f"{'.'.join(sorted([reactant, '_PHOTON']))}__{'.'.join(products)}", products


def _has_signal(dataset: h5py.Dataset) -> bool:
    """True if the dataset exists with at least one non-zero, all-finite value."""
    arr = dataset[:]
    return bool(np.any(arr != 0)) and bool(np.all(np.isfinite(arr)))


def _already_split(h5_path: Path) -> bool:
    """True if every group already carries a ``photodecay`` dataset."""
    with h5py.File(h5_path, "r") as f:
        return all("photodecay" in g for g in f.values())


def _write_group(
    parent: h5py.File,
    key: str,
    reactants: list[str],
    products: list[str],
    decay_type: str,
    energy: np.ndarray,
    photodecay: np.ndarray,
    photoabsorption: np.ndarray | None,
) -> None:
    """Create one per-channel reaction group with a ``photodecay`` dataset."""
    if key in parent:
        raise ValueError(f"duplicate reaction key on split: {key!r}")
    grp = parent.create_group(key)
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


def split_leiden(h5_path: Path, logger) -> int:
    """Split the bundled Leiden file in place into per-channel reactions."""
    if _already_split(h5_path):
        logger.info(f"{h5_path} already split; skipping")
        return 0

    # Read every bundled group into memory before rewriting the file.
    emitted: list[dict] = []
    ionis_seen: set[str] = set()
    with h5py.File(h5_path, "r") as f:
        root_attrs = dict(f.attrs)
        for name, g in f.items():
            reactant = g.attrs["reactants"][0]
            reactant = reactant.decode() if isinstance(reactant, bytes) else str(reactant)
            stem_products = [
                p.decode() if isinstance(p, bytes) else str(p)
                for p in g.attrs["products"]
            ]
            energy = g["photon_energy"][:]
            photoabs = (
                g["photoabsorption"][:]
                if "photoabsorption" in g and _has_signal(g["photoabsorption"])
                else None
            )

            # Dissociation channel -> stem key (products as stored).
            if "photodissociation" in g and _has_signal(g["photodissociation"]):
                emitted.append(
                    {
                        "key": name,
                        "reactants": [reactant],
                        "products": stem_products,
                        "decay_type": "dissociation",
                        "energy": energy,
                        "photodecay": g["photodissociation"][:],
                        "photoabsorption": photoabs,
                    }
                )

            # Ionisation channel -> constructed key.  A reactant with several
            # dissociation channels (e.g. C3H3) appears in several bundled
            # groups but has a single ionisation reaction, so dedupe by key.
            if "photoionization" in g and _has_signal(g["photoionization"]):
                ion_key, ion_products = _ionis_key(reactant)
                if ion_key not in ionis_seen:
                    ionis_seen.add(ion_key)
                    emitted.append(
                        {
                            "key": ion_key,
                            "reactants": [reactant],
                            "products": ion_products,
                            "decay_type": "ionization",
                            "energy": energy,
                            "photodecay": g["photoionization"][:],
                            "photoabsorption": photoabs,
                        }
                    )

    with h5py.File(h5_path, "w") as f:
        for k, v in root_attrs.items():
            f.attrs[k] = v
        f.attrs["created"] = date.today().isoformat()
        f.attrs["description"] = (
            "Leiden photo cross sections, one group per reaction (split by decay "
            "channel); photon_energy in eV, photodecay/photoabsorption in cm^2."
        )
        for e in emitted:
            _write_group(
                f,
                e["key"],
                e["reactants"],
                e["products"],
                e["decay_type"],
                e["energy"],
                e["photodecay"],
                e["photoabsorption"],
            )
    logger.info(f"Split {h5_path} into {len(emitted)} per-channel reactions")
    return len(emitted)


def split_norad(h5_path: Path, logger) -> int:
    """Rename the NORAD ``photoionization`` dataset to ``photodecay`` in place."""
    if _already_split(h5_path):
        logger.info(f"{h5_path} already split; skipping")
        return 0

    groups: list[dict] = []
    with h5py.File(h5_path, "r") as f:
        root_attrs = dict(f.attrs)
        for name, g in f.items():
            reactants = [
                r.decode() if isinstance(r, bytes) else str(r)
                for r in g.attrs["reactants"]
            ]
            products = [
                p.decode() if isinstance(p, bytes) else str(p)
                for p in g.attrs["products"]
            ]
            groups.append(
                {
                    "key": name,
                    "reactants": reactants,
                    "products": products,
                    "energy": g["photon_energy"][:],
                    "photodecay": g["photoionization"][:],
                    "attrs": {
                        k: v
                        for k, v in g.attrs.items()
                        if k not in ("reactants", "products")
                    },
                }
            )

    with h5py.File(h5_path, "w") as f:
        for k, v in root_attrs.items():
            f.attrs[k] = v
        for grp in groups:
            _write_group(
                f,
                grp["key"],
                grp["reactants"],
                grp["products"],
                "ionization",
                grp["energy"],
                grp["photodecay"],
                None,
            )
            for k, v in grp["attrs"].items():
                f[grp["key"]].attrs[k] = v
    logger.info(f"Renamed photoionization->photodecay in {len(groups)} NORAD groups")
    return len(groups)


def main() -> None:
    """Split the bundled Leiden and NORAD xsec files in place."""
    logger = JaffLogger().get_logger()
    split_leiden(XSECS_DIR / "leiden.hdf5", logger)
    split_norad(XSECS_DIR / "norad.hdf5", logger)


if __name__ == "__main__":
    main()
