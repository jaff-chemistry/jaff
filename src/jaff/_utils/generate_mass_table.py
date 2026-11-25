# ABOUTME: Builds the atomic-mass table in jaff.db from data/atom_mass.csv.
# ABOUTME: Run as a one-off maintenance script.
import pandas as pd

from jaff.config import DATA_DIR
from jaff.drivers.sqlite import JaffDb
from jaff.io import JaffLogger


def main():
    masses = DATA_DIR / "atom_mass.csv"
    df = pd.read_csv(
        masses,
        sep=r"\s*,\s*",
        index_col=0,
        dtype={"Protons": "Int64", "Neutrons": "Int64", "Electrons": "Int64"},
    )
    masses_df = df.rename(
        columns={
            "Name": "name",
            "Mass": "mass",
            "AtomicMass": "atomic_mass",
            "Protons": "protons",
            "Neutrons": "neutrons",
            "Electrons": "electrons",
        }
    )
    masses_df.index.name = "element"
    table_name = "atomic_masses"

    with JaffDb() as jdb:
        table = jdb.table_from_dataframe(table_name, masses_df)
        logger = JaffLogger().get_logger()
        logger.info(f"'{table_name}' table created in {jdb.db_path}\n")

        print(pd.DataFrame(table.all_rows()))


if __name__ == "__main__":
    main()
