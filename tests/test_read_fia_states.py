"""Tests for :func:`fia.data_io.read_fia` optional ``states`` filter."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import pandas as pd
import pytest

from fia.data_io import read_fia


def test_read_fia_csv_states_filters_prefix(tmp_path: Path) -> None:
    (tmp_path / "MT_PLOT.csv").write_text("CN,STATECD\n1,30\n", encoding="utf-8")
    (tmp_path / "AZ_PLOT.csv").write_text("CN,STATECD\n2,4\n", encoding="utf-8")

    db = read_fia(str(tmp_path), tables=["PLOT"], common=False, states=["MT"])
    assert len(db["PLOT"]) == 1
    assert int(db["PLOT"]["CN"].iloc[0]) == 1


def test_read_fia_csv_unprefixed_requires_entire(tmp_path: Path) -> None:
    (tmp_path / "PLOT.csv").write_text("CN,STATECD\n1,30\n", encoding="utf-8")
    (tmp_path / "MT_TREE.csv").write_text("PLT_CN,SUBP,TREE\n1,1,1\n", encoding="utf-8")

    db_no = read_fia(str(tmp_path), tables=["PLOT"], common=False, states=["MT"])
    assert "PLOT" not in db_no.tables or len(db_no["PLOT"]) == 0

    db_yes = read_fia(str(tmp_path), tables=["PLOT"], common=False, states=["MT", "ENTIRE"])
    assert len(db_yes["PLOT"]) == 1


def test_read_fia_sqlite_states_filters_by_statecd_and_pltcn(tmp_path: Path) -> None:
    dbpath = tmp_path / "fiadb.sqlite"
    conn = sqlite3.connect(dbpath)
    conn.execute("CREATE TABLE PLOT (CN INTEGER, STATECD INTEGER, INVYR INTEGER)")
    conn.execute("INSERT INTO PLOT VALUES (100, 30, 2020)")
    conn.execute("INSERT INTO PLOT VALUES (200, 4, 2020)")
    conn.execute("CREATE TABLE TREE (PLT_CN INTEGER, SUBP INTEGER, TREE INTEGER)")
    conn.execute("INSERT INTO TREE VALUES (100, 1, 1)")
    conn.execute("INSERT INTO TREE VALUES (200, 1, 1)")
    conn.commit()
    conn.close()

    db = read_fia(
        str(dbpath),
        tables=["PLOT", "TREE"],
        common=False,
        states="MT",
    )
    assert len(db["PLOT"]) == 1
    assert int(db["PLOT"]["CN"].iloc[0]) == 100
    assert len(db["TREE"]) == 1
    assert int(db["TREE"]["PLT_CN"].iloc[0]) == 100


def test_read_fia_sqlite_states_requires_plot_in_tables(tmp_path: Path) -> None:
    dbpath = tmp_path / "fiadb.sqlite"
    conn = sqlite3.connect(dbpath)
    conn.execute("CREATE TABLE TREE (PLT_CN INTEGER)")
    conn.execute("INSERT INTO TREE VALUES (1)")
    conn.commit()
    conn.close()

    with pytest.raises(ValueError, match="requires the PLOT table"):
        read_fia(str(dbpath), tables=["TREE"], common=False, states=["MT"])


def test_read_fia_states_invalid_code(tmp_path: Path) -> None:
    (tmp_path / "MT_PLOT.csv").write_text("CN\n1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="unrecognized"):
        read_fia(str(tmp_path), states=["XX"], tables=["PLOT"], common=False)
