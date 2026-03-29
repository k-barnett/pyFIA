"""Regression: custom_pse (TI, AREA_BASIS) should match area() totals."""

from __future__ import annotations

from pathlib import Path

import pytest

from fia import area, clip_fia, custom_pse, read_fia

_TESTS = Path(__file__).resolve().parent


@pytest.mark.parametrize(
    "fiadb_subdir",
    ["MT"],
)
def test_custom_pse_prop_forest_matches_area_total_by_adforcd(fiadb_subdir: str) -> None:
    db = clip_fia(read_fia(str(_TESTS / "fiadb" / fiadb_subdir)))

    res_area = area(db, land_type="forest", grp_by=["ADFORCD"])
    cond_list = area(
        db,
        land_type="forest",
        cond_list=True,
        grp_by=["ADFORCD"],
    )
    cond_meets = cond_list.assign(EVAL_TYP="CURR")

    cp = custom_pse(
        db=db,
        x=cond_meets,
        x_vars=["PROP_FOREST"],
        x_grp_by=["ADFORCD"],
        method="TI",
        totals=True,
        variance=True,
    )

    merged = res_area.merge(cp, on=["YEAR", "ADFORCD"], how="inner")
    assert not merged.empty
    diff = (merged["AREA_TOTAL"] - merged["PROP_FOREST_TOTAL"]).abs()
    assert (diff < 1e-3).all(), diff.max()
