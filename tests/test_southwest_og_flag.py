"""Southwest OG_FLAG must follow Table 9 (vector[0]), not max(MOG) (maturity can be 1)."""

import numpy as np
import pandas as pd

from fia_mog.estimators import _og_flag_for_condition, _region_is_southwest


def test_southwest_og_flag_uses_first_element_only():
    assert _og_flag_for_condition("southwest", [0.0, 1.0]) == 0.0
    assert _og_flag_for_condition("Southwest ", [0.0, 0.45, 1.0]) == 0.0
    assert _og_flag_for_condition("southwest", [1.0]) == 1.0
    assert _og_flag_for_condition("southwest", [1.0, 0.2]) == 1.0


def test_non_southwest_legacy_max_equals_one():
    assert _og_flag_for_condition("eastern", [0.0, 1.0]) == 1.0
    assert _og_flag_for_condition("", []) == 0.0


def test_region_is_southwest_na_safe():
    assert not _region_is_southwest(None)
    assert not _region_is_southwest(np.nan)
    assert not _region_is_southwest(pd.NA)
