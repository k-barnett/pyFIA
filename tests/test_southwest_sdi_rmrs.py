"""Southwest Table 9 SDI% = 100 * SDI_OFE / SDI_RMRS (rFIA-ofe-region3.R style)."""

from __future__ import annotations

import math

import pandas as pd
import pytest

from fia_mog.engine import ConditionContext
from fia_mog.southwest.evaluate import _southwest_relative_sdi_and_qmd, _southwest_sdi_ofe_live_18


def _ctx(trees: pd.DataFrame, rmrs: float | None) -> ConditionContext:
    return ConditionContext(
        region="southwest",
        forest_type=201,
        condition_area_acres=1.0,
        stand_age=80.0,
        trees=trees,
        condition_sdi_rmrs=rmrs,
    )


def test_sdi_ofe_single_live_large_tree():
    trees = pd.DataFrame(
        {
            "DIA": [20.0],
            "STATUSCD": [1],
            "TPA_UNADJ": [10.0],
        }
    )
    ofe = _southwest_sdi_ofe_live_18(trees)
    expected = 10.0 * (2.0**1.6)
    assert abs(ofe - expected) < 1e-6


def test_sdi_ofe_excludes_small_and_dead():
    trees = pd.DataFrame(
        {
            "DIA": [10.0, 20.0, 22.0],
            "STATUSCD": [1, 1, 2],
            "TPA_UNADJ": [100.0, 6.0, 6.0],
        }
    )
    ofe = _southwest_sdi_ofe_live_18(trees)
    assert abs(ofe - 6.0 * (2.0**1.6)) < 1e-6


def test_relative_sdi_percent_matches_r_formula():
    trees = pd.DataFrame(
        {
            "DIA": [20.0],
            "STATUSCD": [1],
            "TPA_UNADJ": [10.0],
        }
    )
    rmrs = 250.0
    ctx = _ctx(trees, rmrs)
    pct, _qmd, ofe = _southwest_relative_sdi_and_qmd(ctx)
    assert abs(ofe - 10.0 * (2.0**1.6)) < 1e-6
    assert abs(pct - 100.0 * ofe / rmrs) < 1e-6


@pytest.mark.parametrize("bad_rmrs", [None, 0.0, -1.0])
def test_relative_sdi_nan_without_valid_rmrs(bad_rmrs: float | None):
    trees = pd.DataFrame(
        {
            "DIA": [20.0],
            "STATUSCD": [1],
            "TPA_UNADJ": [10.0],
        }
    )
    ctx = _ctx(trees, bad_rmrs)
    pct, _qmd, ofe = _southwest_relative_sdi_and_qmd(ctx)
    assert ofe > 0
    assert isinstance(pct, float)
    assert math.isnan(pct)
