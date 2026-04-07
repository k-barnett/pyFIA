"""Unit tests for :func:`fia.estimators.cond_mean_crown_ratio`."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from fia.data_io import FiaDatabase
from fia.estimators import cond_mean_crown_ratio


def _minimal_db(tree: pd.DataFrame) -> FiaDatabase:
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
            "PLT_CN": [1],
            "CONDID": [1],
            "COND_STATUS_CD": [1],
            "CONDPROP_UNADJ": [1.0],
            "SITECLCD": [1],
            "RESERVCD": [0],
        }
    )
    return FiaDatabase({"PLOT": plot, "COND": cond, "TREE": tree})


def test_weighted_mean_matches_numpy_average_over_100() -> None:
    tree = pd.DataFrame(
        {
            "PLT_CN": [1, 1],
            "CONDID": [1, 1],
            "DIA": [10.0, 20.0],
            "TPA_UNADJ": [1.0, 3.0],
            "STATUSCD": [1, 1],
            "TREECLCD": [1, 1],
            "CR": [20.0, 40.0],
        }
    )
    db = _minimal_db(tree)
    out = cond_mean_crown_ratio(
        db,
        tree_type="live",
        weight_by="tpa",
        weighted=True,
    )
    expected_pct = float(np.average([20.0, 40.0], weights=[1.0, 3.0]))
    assert out["mean_crown_ratio"].iloc[0] == pytest.approx(expected_pct / 100.0)
    assert out["n_trees"].iloc[0] == 2


def test_unweighted_mean_over_100() -> None:
    tree = pd.DataFrame(
        {
            "PLT_CN": [1, 1],
            "CONDID": [1, 1],
            "DIA": [10.0, 30.0],
            "TPA_UNADJ": [1.0, 1.0],
            "STATUSCD": [1, 1],
            "TREECLCD": [1, 1],
            "CR": [10.0, 50.0],
        }
    )
    db = _minimal_db(tree)
    out = cond_mean_crown_ratio(
        db, tree_type="live", weighted=False,
    )
    assert out["mean_crown_ratio"].iloc[0] == pytest.approx(0.3)


def test_live_excludes_dead_from_mean() -> None:
    tree = pd.DataFrame(
        {
            "PLT_CN": [1, 1],
            "CONDID": [1, 1],
            "DIA": [10.0, 10.0],
            "TPA_UNADJ": [1.0, 1.0],
            "STATUSCD": [1, 2],
            "TREECLCD": [1, 1],
            "CR": [50.0, 0.0],
        }
    )
    db = _minimal_db(tree)
    out = cond_mean_crown_ratio(db, tree_type="live", weighted=False)
    assert out["n_trees"].iloc[0] == 1
    assert out["mean_crown_ratio"].iloc[0] == pytest.approx(0.5)


def test_baa_weights_differ_from_tpa_when_diameters_differ() -> None:
    tree = pd.DataFrame(
        {
            "PLT_CN": [1, 1],
            "CONDID": [1, 1],
            "DIA": [5.0, 15.0],
            "TPA_UNADJ": [1.0, 1.0],
            "STATUSCD": [1, 1],
            "TREECLCD": [1, 1],
            "CR": [0.0, 100.0],
        }
    )
    db = _minimal_db(tree)
    tpa_m = cond_mean_crown_ratio(db, weight_by="tpa", weighted=True)["mean_crown_ratio"].iloc[0]
    baa_m = cond_mean_crown_ratio(db, weight_by="baa", weighted=True)["mean_crown_ratio"].iloc[0]
    assert tpa_m == pytest.approx(0.5)
    assert baa_m != pytest.approx(tpa_m, abs=1e-6)
    assert baa_m > tpa_m
