"""OG_FLAG must follow regional OG rules only, never Table 19 maturity alone."""

import numpy as np
import pandas as pd

from fia_mog.engine import ConditionContext, MOGEngine, compute_tree_metrics
from fia_mog.estimators import _region_is_southwest


def _minimal_live_tree_row() -> dict[str, object]:
    return {
        "PLT_CN": 1,
        "CONDID": 1,
        "SUBP": 1,
        "TREE": 1,
        "DIA": 22.0,
        "STATUSCD": 1,
        "CCLCD": 1,
        "TPA_UNADJ": 50.0,
        "HT": 72.0,
    }


def test_southwest_og_flag_matches_first_mog_component():
    trees = pd.DataFrame([_minimal_live_tree_row()])
    ctx = ConditionContext(
        region="southwest",
        forest_type=201,
        condition_area_acres=1.0,
        stand_age=120.0,
        trees=trees,
        plot_statecd=4,
        condition_habtypcd1=11030.0,
        condition_sdi_rmrs=800.0,
    )
    eng = MOGEngine()
    vec = eng.mog_vector(ctx)
    assert len(vec) >= 1
    assert eng.old_growth_flag(ctx) == (1.0 if float(vec[0]) >= 1.0 else 0.0)


def test_eastern_og_flag_zero_when_og_rules_fail_even_if_mog_score_high():
    trees = pd.DataFrame([_minimal_live_tree_row()])
    ctx = ConditionContext(
        region="eastern",
        forest_type=805,
        condition_area_acres=1.0,
        stand_age=50.0,
        trees=trees,
    )
    eng = MOGEngine()
    metrics = compute_tree_metrics(ctx)
    assert metrics is not None
    vec = eng.mog_vector(ctx)
    ogf = eng.old_growth_flag(ctx)
    assert ogf == 0.0
    assert vec[0] == 0.0
    if max(vec, default=0.0) >= 1.0:
        assert ogf == 0.0


def test_region_is_southwest_na_safe():
    assert not _region_is_southwest(None)
    assert not _region_is_southwest(np.nan)
    assert not _region_is_southwest(pd.NA)
