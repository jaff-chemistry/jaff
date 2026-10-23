# ABOUTME: Builds the Verner (1996) photoionisation cross-section table in jaff.db from
# ABOUTME: verner_1996.csv as piecewise SymPy expressions.
import pandas as pd
from sympy import Expr, Piecewise, srepr

from jaff.config import XSECS_DATA_DIR
from jaff.core.network import NetworkSymbols
from jaff.drivers.sqlite import JaffDb
from jaff.io import JaffLogger


def verner_xsecs(
    sigma0: float,
    E0: float,
    y0: float,
    y1: float,
    yw: float,
    ya: float,
    P: float,
    Emin: float,
    Emax: float,
) -> Expr:
    def F(x, y):
        p1 = (x - 1) * (x - 1) + yw * yw
        p2 = y ** (0.5 * P - 5.5)
        p3 = (1 + (y / ya) ** 0.5) ** (-P)

        return p1 * p2 * p3

    E = NetworkSymbols.ph_E
    x = E / E0 - y0
    y = (x * x + y1 * y1) ** 0.5
    mb = 1e-18  # Units of cm^2

    return Piecewise((sigma0 * mb * F(x, y), (E >= Emin) & (E <= Emax)), (0, True))


def main():
    verner_data = XSECS_DATA_DIR / "verner_1996.csv"
    df = pd.read_csv(verner_data, sep=r"\s+", index_col=0)
    rows = [
        {
            "reaction": (
                f"{'.'.join(sorted([ion, '_PHOTON']))}"
                f"__{'.'.join(sorted([f'{ion}+', 'e-']))}"
            ),
            "Z": row["Z"],
            "N": row["N"],
            "xsecs": srepr(
                verner_xsecs(
                    sigma0=row["sigma_0(Mb)"],
                    E0=row["E_0(eV)"],
                    y0=row["y_0"],
                    y1=row["y_1"],
                    yw=row["y_w"],
                    ya=row["y_a"],
                    P=row["P"],
                    Emin=row["E_th(eV)"],
                    Emax=row["E_max(eV)"],
                )
            ),
        }
        for ion, row in df.iterrows()
    ]
    del df
    xsecs_df = pd.DataFrame(rows).set_index("reaction")
    table_name: str = "verner_cross_sections"

    with JaffDb() as jdb:
        table = jdb.table_from_dataframe(table_name, xsecs_df)
        logger = JaffLogger().get_logger()
        logger.info(f"'{table_name}' table created in {jdb.db_path}\n")

        print(pd.DataFrame(table.all_rows()))


if __name__ == "__main__":
    main()
