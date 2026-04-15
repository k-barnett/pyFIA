"""Tests for FIADB column normalization (SQLite / CN vs PLT_CN) used by estimators."""

from __future__ import annotations

import pandas as pd

from fia.data_io import FiaDatabase
from fia.estimators import tpa


def test_tpa_without_clip_with_cn_only_plot_and_lowercase_columns() -> None:
    """Unclipped FIADB-style frames: PLOT may use CN only; SQLite may lowercase names."""
    plot = pd.DataFrame(
        {
            "cn": [42],
            "invyr": [2020],
            "plot_status_cd": [1],
            "measyear": [2020],
        }
    )
    cond = pd.DataFrame(
        {
            "plt_cn": [42],
            "condid": [1],
            "cond_status_cd": [1],
            "condprop_unadj": [1.0],
            "siteclcd": [1],
            "reservcd": [0],
        }
    )
    tree = pd.DataFrame(
        {
            "plt_cn": [42],
            "condid": [1],
            "dia": [10.0],
            "tpa_unadj": [1.0],
            "statuscd": [1],
            "treeclcd": [1],
        }
    )
    db = FiaDatabase({"PLOT": plot, "COND": cond, "TREE": tree})
    out = tpa(db, land_type="forest", tree_type="live", tree_list=True)
    assert not out.empty
    assert "PLT_CN" in out.columns
    assert int(out["PLT_CN"].iloc[0]) == 42
