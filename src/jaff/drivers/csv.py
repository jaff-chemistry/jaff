# ABOUTME: CSV driver: read CSV files into pandas DataFrames and convert them to
# ABOUTME: HDF5 (per-column datasets, flattened coordinate axes or PyTables tables).
"""
CSV file driver utilities.

This module provides thin convenience wrappers around :func:`pandas.read_csv`
and :mod:`h5py` for reading CSV files and converting their contents to HDF5
format:

- ``csv_to_df(file, *args, **kwargs)`` reads a CSV file and returns it as a
  :class:`pandas.DataFrame`; extra arguments go to :func:`pandas.read_csv`.
  Raises :class:`FileNotFoundError` if *file* (a :class:`~pathlib.Path`)
  does not exist.
- ``csv_to_hdf5(file, outfile, ...)`` reads a CSV file and writes it to an
  HDF5 file, adding to *outfile* if it exists (its suffix becomes ``.hdf5``
  unless it is ``.hdf5``, ``.hdf`` or ``.h5``); writing a column into a
  group that already holds it raises :class:`ValueError`.  By default each column is
  stored as a 1-D dataset in the group ``<group_path>/<hdf5_key>``
  (*hdf5_key* defaults to the CSV file stem).  Columns named in
  *flatten_columns* are stored as their sorted unique values instead, which
  suits coordinate axes that repeat across rows.  ``as_table=True`` instead
  writes a PyTables-compatible table with :meth:`pandas.DataFrame.to_hdf`
  (needs the optional ``tables`` package; incompatible with
  *flatten_columns*).

Examples
--------
Given ``grid.csv``:

```text
T,n,k
10,1,0.1
10,2,0.2
20,1,0.3
20,2,0.4
```

>>> from pathlib import Path
>>> from jaff.drivers.csv import csv_to_df, csv_to_hdf5
>>> csv_to_df(Path("grid.csv"))
    T  n    k
0  10  1  0.1
1  10  2  0.2
2  20  1  0.3
3  20  2  0.4

Convert it to HDF5, storing the ``T`` and ``n`` axes once each:

>>> import h5py
>>> csv_to_hdf5("grid.csv", "grid.h5", flatten_columns=["T", "n"])
>>> with h5py.File("grid.h5") as f:
...     {name: f["grid"][name][()].tolist() for name in f["grid"]}
{'T': [10, 20], 'k': [0.1, 0.2, 0.3, 0.4], 'n': [1, 2]}
"""

from __future__ import annotations

from os import PathLike
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from ..common import HDF_EXTENSIONS


def csv_to_df(file: str | PathLike[str], *args, **kwargs) -> pd.DataFrame:
    """
    Read a CSV file and return it as a :class:`pandas.DataFrame`.

    A thin wrapper around :func:`pandas.read_csv` that validates existence
    before delegating.

    Parameters
    ----------
    file : str or PathLike
        Path to the CSV file to read.
    *args
        Positional arguments forwarded to :func:`pandas.read_csv`.
    **kwargs
        Keyword arguments forwarded to :func:`pandas.read_csv`.

    Returns
    -------
    pandas.DataFrame
        Contents of the CSV file, with columns inferred from the file header
        (unless overridden through *kwargs*).

    Raises
    ------
    FileNotFoundError
        If *file* does not exist on the filesystem.

    Examples
    --------
    Given ``grid.csv``:

    ```text
    T,n,k
    10,1,0.1
    10,2,0.2
    20,1,0.3
    20,2,0.4
    ```

    >>> from pathlib import Path
    >>> from jaff.drivers.csv import csv_to_df
    >>> csv_to_df(Path("grid.csv"))
        T  n    k
    0  10  1  0.1
    1  10  2  0.2
    2  20  1  0.3
    3  20  2  0.4

    Extra arguments go to :func:`pandas.read_csv`:

    >>> csv_to_df(Path("grid.csv"), index_col=0).index.tolist()
    [10, 10, 20, 20]
    """
    file = Path(file)
    if not file.exists():
        raise FileNotFoundError(file)

    return pd.read_csv(file, *args, **kwargs)


def csv_to_hdf5(
    file: str | Path,
    outfile: str | Path,
    group_path: str = "",
    hdf5_key: str = "",
    flatten_columns: list[str] = [],
    as_table: bool = False,
    *args,
    **kwargs,
) -> None:
    """
    Convert a CSV file to HDF5 format.

    Reads *file* with :func:`csv_to_df` and writes the result to *outfile*
    (adding to the HDF5 file if it already exists; other groups are kept, but
    an existing dataset of the same name is not overwritten, see Raises).
    Supports three storage
    modes, chosen by the arguments passed:

    * **Default** – Each DataFrame column is stored as a 1-D dataset inside
      an HDF5 group at ``<group_path>/<hdf5_key>``.
    * **Table** (``as_table=True``) – The DataFrame is stored as a
      PyTables-compatible HDF5 table via :meth:`pandas.DataFrame.to_hdf`,
      keeping the column metadata so the data can be read back with
      :func:`pandas.read_hdf`.  Incompatible with *flatten_columns*.
    * **Flatten** (``flatten_columns=[...]``) – Like the default mode, but the
      listed columns are deduplicated with :func:`numpy.unique` (sorted unique
      values) before storage; useful for coordinate axes that repeat across
      rows, such as a temperature column of a gridded table.

    Parameters
    ----------
    file : str or Path
        Path to the source CSV file.
    outfile : str or Path
        Path to the output HDF5 file.  If the suffix is not a recognised HDF5
        extension (``.hdf5``, ``.hdf``, ``.h5``) it is replaced with
        ``.hdf5``.
    group_path : str, optional
        Parent HDF5 group path under which *hdf5_key* is created.
        Defaults to ``""`` (the root group).
    hdf5_key : str, optional
        Name of the HDF5 group/dataset key.  Defaults to the stem of *file*.
    flatten_columns : list of str, optional
        Column names whose values should be stored as the unique-sorted set of
        values (via :func:`numpy.unique`) rather than the raw array.  The
        remaining columns are stored in full, after sorting the rows by these
        columns.
    as_table : bool, optional
        If ``True``, write the DataFrame as a PyTables-compatible HDF5 table.
        Requires the optional ``tables`` (PyTables) package, which is not a
        JAFF dependency.  Cannot be combined with *flatten_columns*.
    *args
        Additional positional arguments forwarded to :func:`csv_to_df` (and
        on to :func:`pandas.read_csv`).
    **kwargs
        Additional keyword arguments forwarded to :func:`csv_to_df` (and on
        to :func:`pandas.read_csv`).

    Returns
    -------
    None

    Raises
    ------
    FileNotFoundError
        If *file* does not exist.
    ValueError
        If both *as_table* and *flatten_columns* are specified, or (default
        and flatten modes) if the target group already holds a dataset named
        after one of the columns, e.g. when converting the same CSV into the
        same file twice.
    RuntimeError
        If any name in *flatten_columns* is not a column in the DataFrame.
        The message lists every name in *flatten_columns*, not only the
        missing ones, and the target group has already been created in
        *outfile* (left empty if it was new) when this is raised.

    Examples
    --------
    With ``grid.csv`` holding columns ``T``, ``n`` and ``k`` on a 2x2 grid
    (see :func:`csv_to_df`), flatten the two coordinate columns:

    >>> import h5py
    >>> from jaff.drivers.csv import csv_to_hdf5
    >>> csv_to_hdf5("grid.csv", "grid.h5", flatten_columns=["T", "n"])
    >>> with h5py.File("grid.h5") as f:
    ...     {name: f["grid"][name][()].tolist() for name in f["grid"]}
    {'T': [10, 20], 'k': [0.1, 0.2, 0.3, 0.4], 'n': [1, 2]}

    Store every column in full under ``/tables/raw``; the missing extension
    becomes ``.hdf5``:

    >>> csv_to_hdf5("grid.csv", "grid", group_path="tables", hdf5_key="raw")
    >>> with h5py.File("grid.hdf5") as f:
    ...     list(f["tables/raw"])
    ['T', 'k', 'n']
    """
    if isinstance(file, str):
        file = Path(file)

    if not group_path:
        group_path = ""

    if as_table and flatten_columns:
        raise ValueError("Cannot flatten arrays if data is to be saved as table")

    # Default the HDF5 key to the CSV file stem.
    if not hdf5_key:
        hdf5_key = file.stem

    if not file.exists():
        raise FileNotFoundError(file)

    outfile = Path(outfile) if isinstance(outfile, str) else outfile
    # Ensure the output has a valid HDF5 file extension.
    if outfile.suffix.lower() not in HDF_EXTENSIONS:
        outfile = outfile.with_suffix(".hdf5")

    df = csv_to_df(file, *args, **kwargs)

    if as_table:
        # PyTables table format — lets pandas store with full column metadata.
        df.to_hdf(outfile, key=f"{group_path}/{hdf5_key}", format="table")
        return

    # Open in append mode so existing groups/datasets are preserved.
    with h5py.File(outfile, "a") as f:
        grp = f.require_group(f"{group_path}/{hdf5_key}")
        if flatten_columns:
            # Validate that every requested flatten column actually exists.
            if any(col not in df.columns for col in flatten_columns):
                raise RuntimeError(
                    f"Invalid columns found in flatten_columns argument: {', '.join(flatten_columns)}"
                )
            # Sort the DataFrame by the flattened coordinate columns so that
            # numpy.unique returns a contiguous, ordered coordinate axis.
            df = df.sort_values(by=flatten_columns)

        for col in df.columns:
            grp.create_dataset(
                col,
                # For flatten columns store only unique values; all others get
                # the full array.
                data=np.unique(df[col].to_numpy())
                if col in flatten_columns
                else df[col].to_numpy(),
            )
