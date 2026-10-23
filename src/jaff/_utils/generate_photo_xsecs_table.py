# ABOUTME: Builds the photo_reaction_cross_sections table in jaff.db from the groups
# ABOUTME: of the shipped leiden.hdf5 / norad.hdf5 cross-section files.
"""Build the ``photo_reaction_cross_sections`` table in ``jaff.db``.

One row per serialized reaction found in the collapsed cross-section HDF5
files (``data/xsecs/leiden.hdf5`` and ``data/xsecs/norad.hdf5``), keyed by the
HDF5 group name (the serialized reaction stem, e.g. ``"CO._PHOTON__C.O"``).

Each HDF5 group is a single decay channel carrying one ``photodecay`` dataset
(plus an optional ``photoabsorption``) and a ``decay_type`` attribute, so a
reaction is unambiguously either a dissociation or an ionisation -- never both.

Columns
-------
- ``reaction``          -- PK, serialized reaction / HDF5 group name.
- ``photo_absorption``  -- 1 only for the H2 dissociation ``H2._PHOTON__H.H``, else 0.
- ``decay_type``        -- ``"dissociation"`` or ``"ionization"`` (from the
  group's ``decay_type`` attribute).
- ``leiden``            -- ``data/xsecs/leiden.hdf5::<group>`` if present, else NULL.
- ``norad``             -- ``data/xsecs/norad.hdf5::<group>`` if present, else NULL.
"""

from pathlib import Path

import h5py
import pandas as pd

from jaff.config import JAFF_DIR
from jaff.drivers.sqlite import JaffDb
from jaff.io import JaffLogger

#: Reaction stem flagged as photo-absorption (H2 dissociation).
H2_DISSOCIATION: str = "H2._PHOTON__H.H"

#: Package root (``src/jaff``) -- HDF5 paths are stored relative to this.
PKG_ROOT: Path = JAFF_DIR

#: source name -> collapsed HDF5 file (relative to ``PKG_ROOT``).
XSEC_FILES: dict[str, Path] = {
    "leiden": Path("data") / "xsecs" / "leiden.hdf5",
    "norad": Path("data") / "xsecs" / "norad.hdf5",
}


def read_decay_types(h5_path: Path) -> dict[str, str]:
    """Map each group name to its ``decay_type`` attribute."""
    with h5py.File(h5_path, "r") as h5:
        return {name: grp.attrs["decay_type"] for name, grp in h5.items()}


def main() -> None:
    sources = {src: read_decay_types(PKG_ROOT / rel) for src, rel in XSEC_FILES.items()}
    reactions = sorted(set().union(*(g.keys() for g in sources.values())))

    rows = []
    for reaction in reactions:
        present = {src for src, groups in sources.items() if reaction in groups}
        # decay_type is identical across sources for a shared reaction.
        decay_type = next(sources[src][reaction] for src in present)
        rows.append(
            {
                "reaction": reaction,
                "photo_absorption": int(reaction == H2_DISSOCIATION),
                "decay_type": decay_type,
                **{
                    src: (
                        f"{XSEC_FILES[src].as_posix()}::{reaction}"
                        if src in present
                        else None
                    )
                    for src in XSEC_FILES
                },
            }
        )

    df = pd.DataFrame(rows).set_index("reaction")
    table_name = "photo_reaction_cross_sections"

    with JaffDb() as jdb:
        table = jdb.table_from_dataframe(table_name, df)
        logger = JaffLogger().get_logger()
        logger.info(f"'{table_name}' table created in {jdb.db_path}\n")
        print(pd.DataFrame(table.all_rows()))


if __name__ == "__main__":
    main()
