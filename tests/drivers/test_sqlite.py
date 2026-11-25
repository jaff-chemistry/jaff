# ABOUTME: Tests for drivers.sqlite: table and column names are quoted (spaces, brackets,
# ABOUTME: reserved words work) and values are bound, never pasted into the SQL

import inspect
import sqlite3

import pandas as pd
import pytest

from jaff.drivers import Db
from jaff.drivers.sqlite import Table

_NAME = "odd table"


@pytest.fixture
def db(tmp_path):
    path = tmp_path / "t.db"
    sqlite3.connect(path).close()
    df = pd.DataFrame(
        {"E (eV)": [13.6, 24.6], "order": [1, 2]},
        index=pd.Index(["H", "it's"], name="species"),
    )
    with Db(path) as handle:
        handle.table_from_dataframe(_NAME, df)
        yield handle


def _rows(rows):
    return [tuple(row) for row in rows]


class TestNames:
    def test_all_rows(self, db):
        assert _rows(db.table(_NAME).all_rows()) == [("H", 13.6, 1), ("it's", 24.6, 2)]

    def test_selected_columns(self, db):
        table = db.table(_NAME)

        assert _rows(table.all_rows(cols=["E (eV)", "order"])) == [(13.6, 1), (24.6, 2)]
        assert _rows(table.rows(cols=["order"])) == [(1,), (2,)]

    def test_get_and_set_cell(self, db):
        table = db.table(_NAME)
        table.set_cell("species", "H", "E (eV)", 13.598)

        assert table.get_cell("species", "H", "E (eV)") == pytest.approx(13.598)
        assert table.get_cell("species", "it's", "order") == 2

    def test_add_column_with_a_quoted_text_default(self, db):
        table = db.table(_NAME)
        table.add_column("note (text)", default="it's")

        assert _rows(table.all_rows(cols=["note (text)"])) == [("it's",), ("it's",)]

    def test_insert_and_delete(self, db):
        table = db.table(_NAME)
        table.insert_row(["C", 11.3, 3])

        assert len(table.all_rows()) == 3

        table.delete()

        assert db.get_table_names() == []

    def test_a_missing_table_is_a_value_error(self, db):
        with pytest.raises(ValueError, match="Invalid table name"):
            Table("x' OR '1'='1", db.connection, db.cursor)


class TestValues:
    def test_conditions_take_bound_params(self, db):
        rows = db.table(_NAME).rows(conditions="species = ?", params=["it's"])

        assert _rows(rows) == [("it's", 24.6, 2)]

    @pytest.mark.parametrize("method", ["all_rows", "rows"])
    def test_cols_has_no_mutable_default(self, method):
        default = inspect.signature(getattr(Table, method)).parameters["cols"].default

        assert default is None
