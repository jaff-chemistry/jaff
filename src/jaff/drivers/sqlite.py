# ABOUTME: SQLite driver: Db (connection wrapper), Table (single-table queries) and
# ABOUTME: JaffDb (handle on the bundled jaff.db photo/Verner/atomic-mass database).
"""
SQLite database driver.

This module provides two general-purpose SQLite wrappers (:class:`Db` and
:class:`Table`) and one JAFF-specific convenience subclass (:class:`JaffDb`)
that points at the database bundled with JAFF (photo-reaction cross-section
sources, Verner photoionization fits and atomic masses).

Classes
-------
Db
    Context-manager wrapper around a :class:`sqlite3.Connection`.
Table
    Thin query interface for a single table inside a :class:`Db`.
JaffDb
    Subclass of :class:`Db` pre-configured for the built-in ``jaff.db``
    database.

Examples
--------
>>> from jaff.drivers import JaffDb
>>> with JaffDb() as db:
...     table = db.table("verner_cross_sections")
...     rows = table.rows(cols=["reaction", "xsecs"], conditions="Z = 1")
"""

from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Any

import pandas as pd

from ..config import DB_DIR
from .csv import csv_to_df


class Db:
    """
    Context-manager wrapper around a SQLite database connection.

    Manages connection lifecycle (open/close), provides table enumeration,
    and exposes helpers for constructing :class:`Table` objects and running
    raw SQL queries.  The constructor only validates and stores the path;
    the connection is opened by :meth:`connect`, or on entering a ``with``
    block (which also calls :meth:`close` on exit).

    Parameters
    ----------
    db_path : Path or str
        Path to the ``.db`` SQLite file.  The file must already exist; the
        driver does not create new databases.

    Attributes
    ----------
    connection : sqlite3.Connection or None
        Active database connection; ``None`` before :meth:`connect` is called.
    cursor : sqlite3.Cursor or None
        Active cursor; ``None`` before :meth:`connect` is called.
    db_path : Path
        Path to the database file, as given (a ``str`` is converted to
        :class:`~pathlib.Path`; it is not resolved).

    Examples
    --------
    Build a table from a DataFrame in an existing (here empty) database file:

    >>> import sqlite3
    >>> import pandas as pd
    >>> from jaff.drivers import Db
    >>> sqlite3.connect("demo.db").close()
    >>> df = pd.DataFrame(
    ...     {"mass": [1.008, 4.003], "charge": [0, 0]},
    ...     index=pd.Index(["H", "He"], name="species"),
    ... )
    >>> with Db("demo.db") as db:
    ...     table = db.table_from_dataframe("species", df)
    ...     db.get_table_names()
    ['species']

    Without a ``with`` block, open and close the connection explicitly:

    >>> db = Db("demo.db")
    >>> db.connect()
    >>> [tuple(row) for row in db.query("SELECT species, mass FROM species")]
    [('H', 1.008), ('He', 4.003)]
    >>> db.close()
    """

    def __init__(self, db_path: Path | str):
        """Validate and store the database path without opening a connection.

        Parameters
        ----------
        db_path : Path or str
            Path to the ``.db`` SQLite file.  The file must already exist; the
            driver does not create new databases.

        Raises
        ------
        ValueError
            If *db_path* is not a :class:`str` or :class:`~pathlib.Path`.
        FileNotFoundError
            If *db_path* does not exist on the filesystem.
        """
        if not isinstance(db_path, (str, Path)):
            raise ValueError("Database path must be a string or pathlib.Path object")

        self.connection: sqlite3.Connection | None = None
        self.cursor: sqlite3.Cursor | None = None
        self.table_obj: Table

        if isinstance(db_path, str):
            db_path = Path(db_path)

        if not db_path.exists():
            raise FileNotFoundError(db_path)

        self.db_path = db_path

    # ------------------------------------------------------------------
    # Connection lifecycle
    # ------------------------------------------------------------------

    def connect(self) -> None:
        """
        Open the database connection and create a cursor.

        Uses :class:`sqlite3.Row` as the row factory so that column values
        can be accessed by name (``row["column_name"]``) as well as by index.
        Called automatically on entering a ``with`` block.

        Returns
        -------
        None
        """
        self.connection = sqlite3.connect(self.db_path)
        # Row factory enables dict-like access: row["column_name"]
        self.connection.row_factory = sqlite3.Row
        self.cursor = self.connection.cursor()

    def close(self) -> None:
        """
        Close the cursor and the underlying database connection.

        Releases the database file handle.  Called automatically on exiting a
        ``with`` block.

        Returns
        -------
        None

        Raises
        ------
        RuntimeError
            If called before :meth:`connect`.
        """
        self.__verify_cursor_instance()
        assert self.cursor is not None
        assert self.connection is not None

        self.cursor.close()
        self.connection.close()

    def delete(self) -> None:
        """
        Delete the database file from disk.

        Does not close an open connection first; call :meth:`close` before
        deleting if one is active.  This cannot be undone.

        Returns
        -------
        None
        """
        self.db_path.unlink()

    # ------------------------------------------------------------------
    # Table enumeration
    # ------------------------------------------------------------------

    def get_tables(self) -> list[Any]:
        """
        Return all table records from the SQLite master table.

        Prefer :meth:`get_table_names` for plain strings.

        Returns
        -------
        list[sqlite3.Row]
            Each element is a row from ``sqlite_master`` with ``type='table'``,
            holding only its ``name`` column.

        Raises
        ------
        RuntimeError
            If called before :meth:`connect`.
        """
        self.__verify_cursor_instance()
        assert self.cursor is not None

        self.cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")

        return self.cursor.fetchall()

    def get_table_names(self) -> list[str]:
        """
        Return the names of all tables in the database.

        Returns
        -------
        list[str]
            Table names, in the (unsorted) order of the SQLite schema.

        Raises
        ------
        RuntimeError
            If called before :meth:`connect`.
        """
        self.__verify_cursor_instance()
        assert self.cursor is not None

        self.cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")

        return [row[0] for row in self.cursor.fetchall()]

    # ------------------------------------------------------------------
    # Table access
    # ------------------------------------------------------------------

    def table(self, name: str) -> "Table":
        """
        Return a :class:`Table` object for the named table.

        The wrapper is bound to the current connection and cursor.  The table
        must already exist; use :meth:`table_from_dataframe` or
        :meth:`table_from_csv` to create one.

        Parameters
        ----------
        name : str
            Name of the table to wrap.

        Returns
        -------
        Table
            Query interface for the named table.

        Raises
        ------
        RuntimeError
            If called before :meth:`connect`.
        ValueError
            If *name* does not correspond to an existing table (delegated
            to the :class:`Table` constructor).
        """
        self.__verify_cursor_instance()
        assert self.cursor is not None
        assert self.connection is not None

        self.table_obj = Table(name, self.connection, self.cursor)

        return self.table_obj

    def table_from_dataframe(self, name: str, df: pd.DataFrame) -> "Table":
        """
        Create a new table from a :class:`pandas.DataFrame` and return it.

        Drops any existing table with the same name, infers SQLite column
        types from the DataFrame dtypes (INTEGER / REAL / TEXT), creates the
        table with a primary-key index column, and populates it via
        :meth:`pandas.DataFrame.to_sql`.

        Parameters
        ----------
        name : str
            Name of the table to create.  Any existing table with this name is
            dropped first.
        df : pandas.DataFrame
            Source data.  The DataFrame's index becomes a ``TEXT`` primary key
            (named after ``df.index.name`` or ``"index"`` if unnamed).  Integer
            columns become ``INTEGER``, float columns ``REAL`` and all others
            ``TEXT``.

        Returns
        -------
        Table
            Query interface for the newly created table.

        Raises
        ------
        RuntimeError
            If called before :meth:`connect`.

        Examples
        --------
        With *db* a connected :class:`Db`:

        >>> import pandas as pd
        >>> df = pd.DataFrame(
        ...     {"mass": [1.008, 4.003], "charge": [0, 0]},
        ...     index=pd.Index(["H", "He"], name="species"),
        ... )
        >>> table = db.table_from_dataframe("species", df)
        >>> [tuple(row) for row in table.all_rows()]
        [('H', 1.008, 0), ('He', 4.003, 0)]
        """
        self.__verify_cursor_instance()
        assert self.cursor is not None
        assert self.connection is not None

        # Drop any pre-existing table of the same name.
        self.cursor.execute(f"DROP TABLE IF EXISTS {name}")

        index_col = df.index.name or "index"

        # Map pandas dtypes to the closest SQLite affinity type.
        col_defs = []
        for col_name, dtype in df.dtypes.items():
            if pd.api.types.is_integer_dtype(dtype):
                sql_type = "INTEGER"
            elif pd.api.types.is_float_dtype(dtype):
                sql_type = "REAL"
            else:
                sql_type = "TEXT"
            col_defs.append(f"{col_name} {sql_type}")

        columns_str = ", ".join(col_defs)

        self.connection.execute(f"""
            CREATE TABLE {name} (
                {index_col} TEXT PRIMARY KEY,
                {columns_str}
            )
        """)

        df.to_sql(name=name, con=self.connection, if_exists="append", index=True)

        return Table(name, self.connection, self.cursor)

    def table_from_csv(self, name: str, file: Path, delimiter: str) -> "Table":
        """
        Create a new table from a CSV file and return it.

        Delegates to :meth:`table_from_dataframe` after reading the CSV with
        :func:`~jaff.drivers.csv.csv_to_df`.  The first column of the CSV is
        used as the DataFrame index (and therefore the primary key).  Column
        types are inferred as in :meth:`table_from_dataframe`, and any
        existing table with the same name is replaced.

        Parameters
        ----------
        name : str
            Name of the table to create.
        file : Path
            Path to the source CSV file.
        delimiter : str
            Field delimiter character for the CSV file, e.g. ``","``.

        Returns
        -------
        Table
            Query interface for the newly created table.

        Raises
        ------
        RuntimeError
            If called before :meth:`connect`.
        FileNotFoundError
            If *file* does not exist.
        """
        return self.table_from_dataframe(
            name, csv_to_df(file, sep=delimiter, index_col=0)
        )

    # ------------------------------------------------------------------
    # Raw query
    # ------------------------------------------------------------------

    def query(self, query: str) -> list[Any]:
        """
        Execute a raw SQL query and return all result rows.

        Parameters
        ----------
        query : str
            Valid SQL query string, e.g. ``"SELECT * FROM rates WHERE T > 100"``.

        Returns
        -------
        list[sqlite3.Row]
            All rows returned by the query (empty if there are none).  Rows
            support both index and column-name access (``row["mass"]``).

        Raises
        ------
        RuntimeError
            If called before :meth:`connect`, or if the query is invalid.

        Examples
        --------
        >>> from jaff.drivers import JaffDb
        >>> with JaffDb() as db:
        ...     [
        ...         tuple(row)
        ...         for row in db.query(
        ...             "SELECT reaction, Z, N FROM verner_cross_sections WHERE Z = 1"
        ...         )
        ...     ]
        [('H._PHOTON__H+.e-', 1, 1)]
        """
        try:
            self.__verify_cursor_instance()
            assert self.cursor is not None

            self.cursor.execute(query)
        except sqlite3.Error:
            raise RuntimeError(f"Invalid query: {query}")

        return self.cursor.fetchall()

    # ------------------------------------------------------------------
    # Context manager protocol
    # ------------------------------------------------------------------

    def __enter__(self) -> "Db":
        """Open the connection and return self for use in a ``with`` block."""
        self.connect()

        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Close the connection when exiting a ``with`` block."""
        self.close()

    # ------------------------------------------------------------------
    # Dunder helpers
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        """Return detailed string representation of this database wrapper.

        Returns
        -------
        str
            Multi-line string showing path, connection, and cursor state.
        """
        info = (
            "JAFF SQLite database object\n"
            f"Path: {self.db_path}\n"
            f"Connection: {self.connection}\n"
            f"Cursor: {self.cursor}"
        )

        return info

    def __str__(self) -> str:
        """Return a short human-readable description of this database wrapper.

        Returns
        -------
        str
        """
        return "JAFF database object"

    def __verify_cursor_instance(self) -> None:
        """
        Raise if the cursor has not been initialised yet.

        Raises
        ------
        RuntimeError
            If ``self.cursor`` is not a :class:`sqlite3.Cursor` instance.
        """
        if not isinstance(self.cursor, sqlite3.Cursor):
            raise RuntimeError(
                "No cursor object found.\nPlease instantiate the database first."
            )


class Table:
    """
    Query interface for a single table inside a :class:`Db` connection.

    Provides methods to select and insert rows, add columns, read and update
    single cells, and drop the whole table (:meth:`delete`) without writing
    raw SQL for each operation.  There is no method to delete rows.
    Instances are usually obtained from :meth:`Db.table`,
    :meth:`Db.table_from_dataframe` or :meth:`Db.table_from_csv`.

    Parameters
    ----------
    name : str
        Name of the table.  Must already exist in the database.
    conn : sqlite3.Connection
        Active database connection (typically obtained from :class:`Db`).
    cur : sqlite3.Cursor
        Active cursor for the connection.

    Attributes
    ----------
    name : str
        Name of the table in the database.
    conn : sqlite3.Connection
        Database connection, shared with the parent :class:`Db`.
    cur : sqlite3.Cursor
        Cursor used to execute queries on the table.

    Examples
    --------
    With *db* a connected :class:`Db` holding the ``species`` table built in
    the :class:`Db` examples:

    >>> table = db.table("species")
    >>> table.get_cell("species", "He", "mass")
    4.003
    >>> table.set_cell("species", "He", "charge", 1)
    >>> table.add_column("note", default="none")
    >>> table.insert_row(["C", 12.011, 0, "atom"])
    >>> [
    ...     tuple(row)
    ...     for row in table.rows(cols=["species", "note"], conditions="mass > 2")
    ... ]
    [('He', 'none'), ('C', 'atom')]
    """

    def __init__(self, name: str, conn: sqlite3.Connection, cur: sqlite3.Cursor):
        """Validate the connection and cursor, then bind to *name* in the database.

        Parameters
        ----------
        name : str
            Name of the table; must already exist in the database.
        conn : sqlite3.Connection
            Active database connection.
        cur : sqlite3.Cursor
            Active cursor for the connection.

        Raises
        ------
        ValueError
            If *conn* is not a :class:`sqlite3.Connection`, *cur* is not a
            :class:`sqlite3.Cursor`, or *name* does not match any table in the
            database.
        """
        if not isinstance(conn, sqlite3.Connection):
            raise ValueError(
                f"Invalid connection passed: {conn}\n"
                "Connection must be an object of sqlite3.Connection"
            )
        if not isinstance(cur, sqlite3.Cursor):
            raise ValueError(
                f"Invalid cursor passed: {cur}\n"
                "Cursor must be an object of sqlite3.Cursor"
            )
        # Verify the table actually exists in the database.
        comm = f"SELECT * FROM sqlite_master WHERE type='table' and name='{name}'"
        if not cur.execute(comm).fetchall():
            raise ValueError(f"Invalid table name: {name}")

        self.name = name
        self.cur = cur
        self.conn = conn

    # ------------------------------------------------------------------
    # Read operations
    # ------------------------------------------------------------------

    def all_rows(self, cols: list[str] = ["*"]) -> list[Any]:
        """
        Fetch every row in the table.

        Parameters
        ----------
        cols : list of str, optional
            Columns to retrieve.  Defaults to all columns (``["*"]``).

        Returns
        -------
        list[sqlite3.Row]
            All rows in the table for the selected columns.  Rows support both
            index and column-name access.
        """
        comm = f"SELECT {','.join(cols) if len(cols) > 1 else cols[0]} FROM {self.name}"
        self.cur.execute(comm)

        return self.cur.fetchall()

    def rows(
        self,
        cols: list[str] = ["*"],
        conditions: str = "",
    ) -> list[Any]:
        """
        Fetch rows that match an optional SQL ``WHERE`` clause.

        Parameters
        ----------
        cols : list of str, optional
            Columns to retrieve.  Defaults to all columns (``["*"]``).
        conditions : str, optional
            SQL ``WHERE`` clause body (without the ``WHERE`` keyword), e.g.
            ``"nion = 1 AND ne < 5"``.  If empty, all rows are returned.

        Returns
        -------
        list[sqlite3.Row]
            Matching rows for the selected columns.
        """
        columns = ",".join(cols) if len(cols) > 1 else cols[0]

        comm = (
            f"SELECT {columns} FROM {self.name} WHERE {conditions}"
            if conditions
            else f"SELECT {columns} FROM {self.name}"
        )
        self.cur.execute(comm)

        return self.cur.fetchall()

    # ------------------------------------------------------------------
    # Write operations
    # ------------------------------------------------------------------

    def insert_row(self, values: list[str | float | int]) -> None:
        """
        Insert a single row into the table.

        Uses parameterised placeholders to guard against SQL injection, and
        commits immediately.

        Parameters
        ----------
        values : list of str, float, or int
            One value per column (including the primary key), in schema
            order.

        Returns
        -------
        None
        """
        placeholders = ",".join(["?"] * len(values))
        comm = f"INSERT INTO {self.name} VALUES ({placeholders})"
        self.cur.execute(comm, values)
        self.conn.commit()

    def insert_rows(self, rows: list[list[str | float | int]]) -> None:
        """
        Insert multiple rows into the table.

        Calls :meth:`insert_row` for each row, committing after every
        individual insertion.  For bulk imports prefer
        :meth:`Db.table_from_dataframe`, which uses pandas' faster path.

        Parameters
        ----------
        rows : list of list
            Each inner list must be a valid argument to :meth:`insert_row`.

        Returns
        -------
        None
        """
        for row in rows:
            self.insert_row(row)

    def add_column(
        self,
        col_name: str,
        col_type: str = "TEXT",
        default: str | float | int | None = None,
    ) -> None:
        """
        Add a new column to the table.

        Issues an ``ALTER TABLE ... ADD COLUMN`` statement and commits.
        Existing rows receive *default* (or ``NULL`` if *default* is not
        given) for the new column.

        Parameters
        ----------
        col_name : str
            Name of the column to add.
        col_type : str, optional
            SQLite column type/affinity (e.g. ``"TEXT"``, ``"INTEGER"``,
            ``"REAL"``, ``"BLOB"``).  Defaults to ``"TEXT"``.
        default : str, float, int, or None, optional
            Default value applied to existing rows and used when no value is
            supplied on insert.  If ``None`` (the default), no ``DEFAULT``
            clause is emitted and existing rows receive ``NULL``.

        Returns
        -------
        None
        """
        comm = f"ALTER TABLE {self.name} ADD COLUMN {col_name} {col_type}"
        if default is not None:
            literal = f"'{default}'" if isinstance(default, str) else default
            comm += f" DEFAULT {literal}"
        self.cur.execute(comm)
        self.conn.commit()

    def get_cell(
        self,
        index_col: str,
        index_value: str | float | int,
        col_name: str,
    ) -> Any:
        """
        Fetch the value of a single cell, located by row and column.

        The target row is identified by the value of its index column
        (typically the primary key), and the target column by name.

        Parameters
        ----------
        index_col : str
            Name of the index/primary-key column used to locate the row.
        index_value : str, float, or int
            Value of *index_col* identifying the target row.
        col_name : str
            Name of the column whose cell is fetched.

        Returns
        -------
        Any
            The cell value if exactly one row matches; ``None`` if no row
            matches; a list of values if multiple rows match *index_value*.

        Raises
        ------
        ValueError
            If *col_name* does not exist in the table.
        """
        # Validate the target column exists (consistent ValueError on typo).
        self.cur.execute(
            f"SELECT type FROM pragma_table_info('{self.name}') WHERE name = ?",
            [col_name],
        )
        if self.cur.fetchone() is None:
            raise ValueError(f"Invalid column name: {col_name}")

        comm = f"SELECT {col_name} FROM {self.name} WHERE {index_col} = ?"
        self.cur.execute(comm, [index_value])
        rows = self.cur.fetchall()

        if not rows:
            return None
        if len(rows) == 1:
            return rows[0][col_name]

        return [row[col_name] for row in rows]

    def set_cell(
        self,
        index_col: str,
        index_value: str | float | int,
        col_name: str,
        value: str | float | int | bytes,
    ) -> None:
        """
        Set the value of a single cell, located by row and column.

        The target row is identified by the value of its index column
        (typically the primary key), and the target column by name.  Before
        writing, the Python type of *value* is checked against the declared
        SQLite type of *col_name*; the change is then committed.

        Parameters
        ----------
        index_col : str
            Name of the index/primary-key column used to locate the row.
        index_value : str, float, or int
            Value of *index_col* identifying the target row.
        col_name : str
            Name of the column whose cell is being set.
        value : str, float, int, or bytes
            New value for the cell.  Its type must match the declared type of
            *col_name* (TEXT->str, INTEGER->int, REAL->int/float, BLOB->bytes);
            ``bool`` is rejected for INTEGER and REAL columns.  Columns of any
            other declared type accept any value.

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If *col_name* does not exist in the table.
        TypeError
            If the type of *value* does not match the column's declared type.
        """
        # Map declared SQLite type -> acceptable Python type(s).
        type_map: dict[str, tuple[type, ...]] = {
            "TEXT": (str,),
            "INTEGER": (int,),
            "REAL": (int, float),
            "BLOB": (bytes,),
        }

        # Look up the declared type of the target column.
        self.cur.execute(
            f"SELECT type FROM pragma_table_info('{self.name}') WHERE name = ?",
            [col_name],
        )
        col = self.cur.fetchone()
        if col is None:
            raise ValueError(f"Invalid column name: {col_name}")
        col_type = col["type"].upper()

        # Validate the value's Python type against the column affinity.
        # bool is a subclass of int but is rejected for INTEGER/REAL columns.
        allowed = type_map.get(col_type)
        if allowed is not None and (
            isinstance(value, bool) or not isinstance(value, allowed)
        ):
            raise TypeError(
                f"Value {value!r} ({type(value).__name__}) does not match "
                f"column '{col_name}' type {col_type}"
            )

        comm = f"UPDATE {self.name} SET {col_name} = ? WHERE {index_col} = ?"
        self.cur.execute(comm, [value, index_value])
        self.conn.commit()

    def delete(self) -> None:
        """
        Drop this table from the database.

        Issues ``DROP TABLE`` and commits; this cannot be undone, and the
        :class:`Table` object is no longer usable afterwards.

        Returns
        -------
        None
        """
        comm = f"DROP TABLE {self.name}"
        self.cur.execute(comm)
        self.conn.commit()


class JaffDb(Db):
    """
    Pre-configured database handle for the bundled JAFF database.

    Points at ``jaff/db/jaff.db``, which holds the data used by JAFF's
    photo-reaction processing and species handling, one table each:

    - ``photo_reaction_cross_sections``: per photo-reaction, its decay type,
      photo-absorption flag and Leiden/NORAD cross-section sources;
    - ``verner_cross_sections``: the Verner (1996) photoionization fits;
    - ``atomic_masses``: per element, its name, ``mass``, ``atomic_mass`` and
      proton, neutron and electron counts.

    The database path is resolved automatically from the package directory,
    so the constructor takes no arguments.  Access the data via
    :meth:`~Db.table`.

    Examples
    --------
    >>> from jaff.drivers import JaffDb
    >>> with JaffDb() as db:
    ...     table = db.table("verner_cross_sections")
    ...     [
    ...         tuple(row)
    ...         for row in table.rows(cols=["reaction", "Z", "N"], conditions="Z = 1")
    ...     ]
    [('H._PHOTON__H+.e-', 1, 1)]
    """

    def __init__(self):
        """Initialise the database wrapper pointing at the bundled ``jaff.db`` file.

        The path is resolved automatically relative to the JAFF package
        directory; no arguments are required.
        """
        jaff_db_path = DB_DIR / "jaff.db"
        super().__init__(jaff_db_path)
