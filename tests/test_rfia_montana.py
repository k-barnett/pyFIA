import os
from pathlib import Path

import pandas as pd
import pytest
import requests

from fia import FiaDatabase, area, clip_fia, custom_pse, read_fia, tpa

_TESTS = Path(__file__).resolve().parent

# Same logical tables as the previous Datamart subset; CSVs live under
# tests/fiadb/mt or tests/fiadb/MT (case varies by checkout / OS).
_MT_TABLES = [
    "COND",
    "COND_DWM_CALC",
    "INVASIVE_SUBPLOT_SPP",
    "PLOT",
    "POP_ESTN_UNIT",
    "POP_EVAL",
    "POP_EVAL_GRP",
    "POP_EVAL_TYP",
    "POP_PLOT_STRATUM_ASSGN",
    "POP_STRATUM",
    "SUBPLOT",
    "TREE",
    "TREE_GRM_COMPONENT",
    "TREE_GRM_MIDPT",
    "TREE_GRM_BEGIN",
    "SUBP_COND_CHNG_MTRX",
    "SEEDLING",
    "SURVEY",
    "SUBP_COND",
    "P2VEG_SUBP_STRUCTURE",
    "PLOTGEOM",
]


def _montana_fiadb_dir() -> Path:
    base = _TESTS / "fiadb"
    for name in ("mt", "MT"):
        p = base / name
        if p.is_dir():
            return p
    raise FileNotFoundError(
        "Montana FIADB not found. Add state CSVs under "
        f"{base / 'mt'} or {base / 'MT'} (see repo layout next to fiadb/AZ)."
    )


FIADB_FULLREPORT = "https://apps.fs.usda.gov/fiadb-api/fullreport"
FIADB_WC_PARAMS = "https://apps.fs.usda.gov/fiadb-api/fullreport/parameters/wc"

# FIADB-API /fullreport/parameters/snum (EXPVOL, forest land, tree basal area)
FIADB_SNUM_LIVE_BASAL_AREA_SF_FOREST_LAND = 1004
# "Basal area of live trees (at least 1 inch d.b.h./d.r.c.), in square feet, on forest land"


def fiadb_api_get(params: dict) -> dict:
    """Call FIADB-API /fullreport and return parsed pieces."""
    resp = requests.get(FIADB_FULLREPORT, params=params, timeout=120)
    resp.raise_for_status()
    data = resp.json()

    out = {"raw": data}
    out["estimates"] = pd.DataFrame(data.get("estimates", []))
    out["subtotals"] = {
        k: pd.DataFrame(v) for k, v in data.get("subtotals", {}).items()
    }
    out["totals"] = pd.DataFrame(data.get("totals", []))
    out["metadata"] = data.get("metadata", {})
    return out


def fiadb_grand_total_estimate(api_result: dict, value_col: str = "ESTIMATE") -> float:
    """
    Prefer the ``totals`` table from fullreport (NJSON) for the domain estimate;
    ``estimates.iloc[0]`` is often a subgroup row, not the state/grand total.
    """
    totals = api_result.get("totals")
    if isinstance(totals, pd.DataFrame) and not totals.empty and value_col in totals.columns:
        return float(pd.to_numeric(totals[value_col], errors="coerce").iloc[0])
    est = api_result["estimates"]
    assert not est.empty and value_col in est.columns
    return float(est.iloc[0][value_col])


def latest_montana_wc() -> int:
    """
    Query FIADB-API wc parameter dictionary and return the most recent
    Montana (FIPS 30) wc code.
    """
    # Allow manual override so tests can always be run, even if the
    # parameter endpoint is misbehaving or its format changes.
    env_wc = os.getenv("FIA_MT_WC")
    if env_wc is not None:
        return int(env_wc)

    resp = requests.get(FIADB_WC_PARAMS, timeout=60)
    resp.raise_for_status()
    try:
        entries = resp.json()
    except requests.exceptions.JSONDecodeError:
        # Fallback: use a known-good Montana wc code. This should be kept
        # up to date if FIADB publishes new inventories, or overridden via
        # the FIA_MT_WC environment variable above.
        return 302020

    mt_entries = []
    for e in entries:
        label = str(e.get("label", "")).lower()
        value = int(e.get("value"))
        if "montana" in label or str(value).startswith("30"):
            mt_entries.append(value)

    if not mt_entries:
        raise RuntimeError("Could not locate Montana wc codes from FIADB-API.")

    return max(mt_entries)


@pytest.fixture(scope="session")
def mt_db() -> FiaDatabase:
    """
    Load bundled Montana FIADB CSVs from ``tests/fiadb/mt`` or ``tests/fiadb/MT``,
    then clip to the most recent inventories (``clip_fia(..., most_recent=True)``).
    """
    db = read_fia(str(_montana_fiadb_dir()), tables=_MT_TABLES)
    return clip_fia(db, most_recent=True)


@pytest.mark.integration
def test_tpa_tree_list_integration(mt_db: FiaDatabase) -> None:
    """
    - ``tree_list=True``: required columns on per-tree output (live trees, forest land).
    - Design-based ``tpa`` (``tree_list=False``): ``BA_TOTAL`` (ft²) vs FIADB-API
      ``snum=1004`` (same Montana ``wc``). Agreement is within a **relative tolerance**
      (not exact integer match): the Python TI path uses a simplified EU / variance stack
      compared to the full FIADB estimator, so a small residual bias is expected.
    """
    res_list = tpa(
        mt_db,
        land_type="forest",
        tree_type="live",
        tree_list=True,
    )
    assert not res_list.empty
    for col in ["PLT_CN", "CONDID", "SUBP", "TREE", "YEAR", "TPA", "BAA"]:
        assert col in res_list.columns

    wc_code = latest_montana_wc()
    base_params = {
        "wc": wc_code,
        "rselected": "None",
        "cselected": "None",
        "outputFormat": "NJSON",
    }

    api_ba = fiadb_api_get(
        {**base_params, "snum": FIADB_SNUM_LIVE_BASAL_AREA_SF_FOREST_LAND}
    )
    api_ba_sf = fiadb_grand_total_estimate(api_ba)

    res_agg = tpa(
        mt_db,
        land_type="forest",
        tree_type="live",
        tree_list=False,
    )
    res_latest = res_agg.sort_values("YEAR").iloc[-1]
    our_ba_total = float(res_latest["BA_TOTAL"])

    denom = max(abs(api_ba_sf), 1.0)
    rel_err = abs(our_ba_total - api_ba_sf) / denom
    assert rel_err < 0.01, (
        f"BA_TOTAL within 1% of FIADB-API: rel_err={rel_err:.5f}, "
        f"ours={our_ba_total:.0f}, api={api_ba_sf:.0f}"
    )


@pytest.mark.integration
def test_area_cond_list_integration(mt_db: FiaDatabase) -> None:
    """Integration test: area with cond_list=True returns expected structure."""
    res = area(
        mt_db,
        land_type="forest",
        cond_list=True,
    )
    assert not res.empty
    for col in ["PLT_CN", "CONDID", "YEAR", "AREA_BASIS", "PROP_FOREST"]:
        assert col in res.columns
    assert (res["PROP_FOREST"] >= 0).all()


@pytest.mark.integration
def test_custom_pse_integration(mt_db: FiaDatabase) -> None:
    """Integration test: custom_pse on a simple forest/non-forest indicator."""
    cond_list = area(
        mt_db,
        land_type="forest",
        cond_list=True,
    )

    cond_meets = cond_list.assign(
        EVAL_TYP="CURR",
        forest_flag=lambda d: (d["PROP_FOREST"] > 0).astype(float),
    )

    res = custom_pse(
        db=mt_db,
        x=cond_meets,
        x_vars=["forest_flag"],
        x_grp_by=["YEAR"],
        method="TI",
        totals=True,
        variance=True,
    )
    assert not res.empty
    assert "YEAR" in res.columns
    assert any(c.startswith("forest_flag") for c in res.columns)


@pytest.mark.integration
def test_area_mt_matches_fiadb_api_forest_total(mt_db: FiaDatabase) -> None:
    """
    Compare our TI design-based forest area for Montana with FIADB-API
    for the most recent inventory (wc determined dynamically).
    """
    wc_code = latest_montana_wc()

    # Forest land area in acres (see snum table: ATTRIBUTE_NBR = 2)
    snum_forest_area = 2

    api = fiadb_api_get(
        {
            "wc": wc_code,
            "snum": snum_forest_area,
            "rselected": "None",
            "cselected": "None",
            "outputFormat": "NJSON",
        }
    )
    api_total = fiadb_grand_total_estimate(api)

    # Our estimate: forest land area for MT, no grouping; take latest YEAR
    res = area(
        mt_db,
        land_type="forest",
        grp_by=None,
    )
    res_latest = res.sort_values("YEAR").iloc[-1]
    our_total = float(res_latest["AREA_TOTAL"])

    # Compare with integer precision: rounded values must match exactly
    assert int(round(our_total)) == int(round(api_total))

