# ABOUTME: I/O drivers for the file formats JAFF reads and writes: HDF5, SQLite,
# ABOUTME: TOML and CSV; re-exports the driver classes and their modules.
"""I/O drivers for the file formats used by JAFF.

``jaff.drivers`` wraps the libraries JAFF uses to read and write data:

- :class:`HDF5` converts between HDF5 files and the nested
  :class:`~jaff.types.HDF5Dict` representation, with optional compression,
  and exports HDF5 data to CSV.
- :class:`Db` is a low-level SQLite connection and query wrapper;
  :class:`Table` gives row- and cell-level access to one of its tables.
- :class:`JaffDb` is a :class:`Db` pre-configured for the database bundled
  with JAFF.
- :class:`Toml` reads a TOML configuration file.
- The :mod:`jaff.drivers.csv` module loads CSV files into
  :class:`pandas.DataFrame` objects (``csv_to_df``) or converts them to HDF5
  (``csv_to_hdf5``).

Examples
--------
Write a dictionary to an HDF5 file and read it back:

>>> import numpy as np
>>> from jaff.drivers import HDF5
>>> h5 = HDF5(compression="gzip")
>>> h5.from_dict(
...     "rates.hdf5",
...     {"T": {"_kind": "linear", "_data": np.array([100.0, 200.0]), "_dtype": "f64"}},
... )
>>> h5.to_dict("rates.hdf5")["T"]["_data"]
array([100., 200.])

Read a TOML file:

>>> from jaff.drivers import Toml
>>> Toml("config.toml").get_key("network")
{'host': 'localhost', 'port': 8080}

Query the bundled SQLite database:

>>> from jaff.drivers import JaffDb
>>> with JaffDb() as db:
...     db.get_table_names()
['photo_reaction_cross_sections', 'atomic_masses', 'verner_cross_sections']

Load a CSV file:

>>> from pathlib import Path
>>> from jaff.drivers.csv import csv_to_df
>>> df = csv_to_df(Path("grid.csv"))
"""

from . import csv, hdf5, sqlite, toml
from .hdf5 import HDF5
from .sqlite import Db, JaffDb, Table
from .toml import Toml

__all__ = ["sqlite", "toml", "csv", "Db", "JaffDb", "Table", "Toml", "hdf5", "HDF5"]
