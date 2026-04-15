"""Tests for :func:`fia.estimators.tpa` with ``by_cond=True``."""

from __future__ import annotations

import pandas as pd
import pytest

from fia.data_io import FiaDatabase
from fia.estimators import tpa


def _two_cond_db() -> FiaDatabase:
    plot = pd.DataFrame(
        {
            "PLT_CN": [1],
            "INVYR": [2020],
            "PLOT_STATUS_CD": [1],
            "MEASYEAR": [2020],
        }
    )
    cond = pd.DataFrame(
        {
            "PLT_CN": [1, 1],
            "CONDID": [1, 2],
            "COND_STATUS_CD": [1, 1],
            "CONDPROP_UNADJ": [0.5, 0.5],
            "SITECLCD": [1, 1],
            "RESERVCD": [0, 0],
        }
    )
    tree = pd.DataFrame(
        {
            "PLT_CN": [1, 1],
            "CONDID": [1, 2],
            "SUBP": [1, 1],
            "TREE": [1, 1],
            "DIA": [10.0, 20.0],
            "TPA_UNADJ": [1.0, 2.0],
            "STATUSCD": [1, 1],
            "TREECLCD": [1, 1],
        }
    )
    return FiaDatabase({"PLOT": plot, "COND": cond, "TREE": tree})


def test_tpa_by_cond_two_rows_distinct_condid() -> None:
    db = _two_cond_db()
    out = tpa(db, land_type="forest", tree_type="live", by_cond=True)
    assert len(out) == 2
    assert set(out["CONDID"]) == {1, 2}
    assert "PROP_FOREST" in out.columns
    c1 = out.loc[out["CONDID"].eq(1)].iloc[0]
    c2 = out.loc[out["CONDID"].eq(2)].iloc[0]
    assert float(c1["TPA"]) == pytest.approx(1.0)
    assert float(c2["TPA"]) == pytest.approx(2.0)
    assert float(c1["PROP_FOREST"]) == pytest.approx(0.5)
    assert float(c2["PROP_FOREST"]) == pytest.approx(0.5)


def test_tpa_by_plot_and_by_cond_rejected() -> None:
    db = _two_cond_db()
    with pytest.raises(ValueError, match="by_plot"):
        tpa(db, by_plot=True, by_cond=True)
