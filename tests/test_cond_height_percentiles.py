"""Unit tests for :func:`fia.estimators.cond_height_percentiles`."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from fia.data_io import FiaDatabase
from fia.estimators import cond_height_percentiles, _weighted_percentile


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


def test_weighted_percentile_matches_reference_formula() -> None:
    values = np.array([1.0, 2.0, 3.0, 4.0])
    weights = np.array([1.0, 1.0, 1.0, 1.0])
    out = _weighted_percentile(values, weights, [50.0])
    assert out.shape == (1,)
    assert 2.0 < out[0] < 4.0


def test_cond_height_percentiles_raw_median_three_trees() -> None:
    tree = pd.DataFrame(
        {
            "PLT_CN": [1, 1, 1],
            "CONDID": [1, 1, 1],
            "DIA": [10.0, 20.0, 20.0],
            "TPA_UNADJ": [1.0, 1.0, 1.0],
            "STATUSCD": [1, 1, 1],
            "TREECLCD": [1, 1, 1],
            "HT": [10.0, 20.0, 30.0],
        }
    )
    db = _minimal_db(tree)
    out = cond_height_percentiles(
        db,
        tree_type="live",
        percentile_method="raw",
        percentiles=[50.0],
    )
    assert len(out) == 1
    assert out["n_trees"].iloc[0] == 3
    assert out["HT_P50"].iloc[0] == pytest.approx(20.0)


def test_cond_height_percentiles_live_excludes_dead() -> None:
    tree = pd.DataFrame(
        {
            "PLT_CN": [1, 1, 1],
            "CONDID": [1, 1, 1],
            "DIA": [10.0, 20.0, 20.0],
            "TPA_UNADJ": [1.0, 1.0, 1.0],
            "STATUSCD": [1, 1, 2],
            "TREECLCD": [1, 1, 1],
            "HT": [10.0, 20.0, 99.0],
        }
    )
    db = _minimal_db(tree)
    out = cond_height_percentiles(
        db,
        tree_type="live",
        percentile_method="raw",
        percentiles=[50.0],
    )
    assert out["n_trees"].iloc[0] == 2
    assert out["HT_P50"].iloc[0] == pytest.approx(15.0)


def test_cond_height_percentiles_weighted_equal_weights_matches_raw() -> None:
    tree = pd.DataFrame(
        {
            "PLT_CN": [1, 1, 1],
            "CONDID": [1, 1, 1],
            "DIA": [10.0, 20.0, 20.0],
            "TPA_UNADJ": [1.0, 1.0, 1.0],
            "STATUSCD": [1, 1, 1],
            "TREECLCD": [1, 1, 1],
            "HT": [10.0, 20.0, 30.0],
        }
    )
    db = _minimal_db(tree)
    raw = cond_height_percentiles(
        db, tree_type="live", percentile_method="raw", percentiles=[50.0]
    )
    wtp = cond_height_percentiles(
        db,
        tree_type="live",
        percentile_method="weighted",
        weight_by="tpa",
        percentiles=[50.0],
    )
    assert raw["HT_P50"].iloc[0] == pytest.approx(wtp["HT_P50"].iloc[0], rel=0, abs=1e-6)


def test_tree_type_gs_rejected() -> None:
    tree = pd.DataFrame(
        {
            "PLT_CN": [1],
            "CONDID": [1],
            "DIA": [10.0],
            "TPA_UNADJ": [1.0],
            "STATUSCD": [1],
            "TREECLCD": [2],
            "HT": [40.0],
        }
    )
    db = _minimal_db(tree)
    with pytest.raises(ValueError, match="tree_type"):
        cond_height_percentiles(db, tree_type="gs", percentiles=[50.0])
