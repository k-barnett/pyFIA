"""Tests for :func:`fia.data_io.read_fia` with SQLite FIADB files."""

from __future__ import annotations

import sqlite3
import tempfile
from pathlib import Path

import pandas as pd
import pytest

from fia.data_io import FiaDatabase, read_fia


def test_read_fia_sqlite_loads_table_with_uppercase_key() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        dbpath = Path(tmp) / "fiadb.sqlite"
        conn = sqlite3.connect(dbpath)
        conn.execute("CREATE TABLE PLOT (PLT_CN INTEGER, INVYR INTEGER)")
        conn.execute("INSERT INTO PLOT VALUES (1, 2020)")
        conn.commit()
        conn.close()

        db = read_fia(str(dbpath), tables=["PLOT"], common=False)
        assert isinstance(db, FiaDatabase)
        assert list(db.keys()) == ["PLOT"]
        assert len(db["PLOT"]) == 1
        assert int(db["PLOT"]["PLT_CN"].iloc[0]) == 1


def test_read_fia_sqlite_case_insensitive_table_name() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        dbpath = Path(tmp) / "test.db"
        conn = sqlite3.connect(dbpath)
        conn.execute('CREATE TABLE "plot" (x REAL)')
        conn.execute("INSERT INTO plot VALUES (3.14)")
        conn.commit()
        conn.close()

        db = read_fia(str(dbpath), tables=["PLOT"], common=False)
        assert "PLOT" in db.tables
        assert len(db["PLOT"]) == 1


def test_read_fia_rejects_non_sqlite_file_with_clear_error(tmp_path: Path) -> None:
    bad = tmp_path / "not_sqlite.txt"
    bad.write_text("hello", encoding="utf-8")
    with pytest.raises(ValueError, match="SQLite FIADB"):
        read_fia(str(bad), tables=["PLOT"], common=False)
