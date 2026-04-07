"""
Run from the ``fia_py`` repo root, e.g. ``python scripts/test_mog.py``,
or install the project in editable mode (``pip install -e .``).
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

# Allow ``import fia`` / ``import fia_mog`` when the script is executed directly.
_repo_root = Path(__file__).resolve().parent.parent
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

import pandas as pd

from fia import FiaDatabase, clip_fia, read_fia
from fia_mog import ConditionContext, MOGEngine, classify_region
from fia_mog.estimators import _mog_forest_type_from_cond_row


def compute_mog_for_arizona(db: FiaDatabase) -> pd.DataFrame:
    """
    One row per Arizona (PLT_CN, CONDID) with a MOG score:

        0   = not mature / OG
        ~0.5= mature (index)
        1   = old-growth
    """

    conds = db["COND"].copy()
    trees = db["TREE"].copy()

    fld_az = pd.to_numeric(conds["FLDTYPCD"], errors="coerce")
    fort_az = (
        pd.to_numeric(conds["FORTYPCD"], errors="coerce")
        if "FORTYPCD" in conds.columns
        else pd.Series(pd.NA, index=conds.index)
    )
    conds_az = conds[
        (conds["STATECD"] == 4)
        & (conds["COND_STATUS_CD"] == 1)
        & (fld_az.notna() | fort_az.notna())
    ].copy()

    # Join condition-level fields we need onto the tree table:
    # - HABTYPCD1 (for Southwest ERU)
    # - PROP_BASIS, CONDPROP_UNADJ, FLDTYPCD, STDAGE, FLDAGE, SITECLCD, ADFORCD, STATECD
    cond_cols_for_trees = [
        "PLT_CN",
        "CONDID",
        "HABTYPCD1",
        "PROP_BASIS",
        "CONDPROP_UNADJ",
        "FLDTYPCD",
        "FORTYPCD",
        "STDAGE",
        "FLDAGE",
        "SITECLCD",
        "ADFORCD",
        "STATECD",
        "SDI_RMRS",
    ]
    cond_cols_for_trees = [c for c in cond_cols_for_trees if c in conds_az.columns]
    trees_az = trees.merge(
        conds_az[cond_cols_for_trees],
        on=["PLT_CN", "CONDID"],
        how="inner",
        suffixes=("", "_COND"),
    )

    engine = MOGEngine()
    out_rows = []

    for (plt_cn, condid), cond_group in trees_az.groupby(["PLT_CN", "CONDID"]):
        c = cond_group.iloc[0]

        # Region (for AZ this should come out as "southwest")
        region = classify_region(
            adforcd=c["ADFORCD"],
            state_abbrev="AZ",
        )
        if region != "southwest":
            continue  # skip if somehow not classified as SW

        # Condition area (acres), per the R code
        basis = c["PROP_BASIS"]
        sub_plot_area = 0.2460 if (not pd.isna(basis) and basis == "MACR") else 0.0417
        condition_area = sub_plot_area * (4 * float(c["CONDPROP_UNADJ"]))

        # Stand age (max of STDAGE, FLDAGE)
        ages = [c["STDAGE"], c["FLDAGE"]]
        ages = [float(a) for a in ages if pd.notna(a)]
        stand_age = max(ages) if ages else 0.0

        forest_type = _mog_forest_type_from_cond_row(c, region=region)

        # Safely coerce ADFORCD and SITECLCD
        adforcd_val = c["ADFORCD"]
        adforcd_int = None if pd.isna(adforcd_val) else int(adforcd_val)
        
        sitecl_val = c["SITECLCD"]
        sitecl_float = None if pd.isna(sitecl_val) else float(sitecl_val)

        sdi_rmrs_ctx = None
        if "SDI_RMRS" in c.index:
            rmrs = c.get("SDI_RMRS")
            if rmrs is not None and not pd.isna(rmrs):
                try:
                    rmrs_f = float(rmrs)
                    if math.isfinite(rmrs_f) and rmrs_f > 0:
                        sdi_rmrs_ctx = rmrs_f
                except (TypeError, ValueError):
                    pass

        ctx = ConditionContext(
            region=region,
            forest_type=forest_type,
            condition_area_acres=condition_area,
            stand_age=stand_age,
            trees=cond_group,
            plot_statecd=int(c["STATECD"]),
            condition_fortypcd=None,
            condition_physclcd=None,
            condition_siteclcd=sitecl_float,
            condition_adforcd=adforcd_int,
            condition_habtypcd1=c["HABTYPCD1"],  # <-- single value from COND
            condition_sdi_rmrs=sdi_rmrs_ctx,
            ecosubcd=None,
        )

        mog_vec = engine.mog_vector(ctx)
        mog_score = max(mog_vec, default=0.0)

        out_rows.append(
            {
                "PLT_CN": plt_cn,
                "CONDID": condid,
                "MOG_SCORE": mog_score,
            }
        )

    return pd.DataFrame(out_rows)


def compute_mog_for_montana(db: FiaDatabase) -> pd.DataFrame:
    """
    One row per Montana (PLT_CN, CONDID) with a MOG score.

    Same score semantics as :func:`compute_mog_for_arizona` (0 / ~0.5 / 1), using
    ``classify_region(..., state_abbrev="MT")`` (typically ``"northern"``) and
    the shared :class:`~fia_mog.engine.MOGEngine`.

    Expects ``db`` to contain at least ``COND`` and ``TREE`` with Montana
    ``STATECD`` (30).
    """
    conds = db["COND"].copy()
    trees = db["TREE"].copy()

    fld_mt = pd.to_numeric(conds["FLDTYPCD"], errors="coerce")
    fort_mt = (
        pd.to_numeric(conds["FORTYPCD"], errors="coerce")
        if "FORTYPCD" in conds.columns
        else pd.Series(pd.NA, index=conds.index)
    )
    conds_mt = conds[
        (conds["STATECD"] == 30)
        & (conds["COND_STATUS_CD"] == 1)
        & (fld_mt.notna() | fort_mt.notna())
    ].copy()

    cond_cols_for_trees = [
        "PLT_CN",
        "CONDID",
        "HABTYPCD1",
        "PROP_BASIS",
        "CONDPROP_UNADJ",
        "FLDTYPCD",
        "FORTYPCD",
        "STDAGE",
        "FLDAGE",
        "SITECLCD",
        "ADFORCD",
        "STATECD",
    ]
    cond_cols_for_trees = [c for c in cond_cols_for_trees if c in conds_mt.columns]
    trees_mt = trees.merge(
        conds_mt[cond_cols_for_trees],
        on=["PLT_CN", "CONDID"],
        how="inner",
        suffixes=("", "_COND"),
    )

    engine = MOGEngine()
    out_rows = []

    for (plt_cn, condid), cond_group in trees_mt.groupby(["PLT_CN", "CONDID"]):
        c = cond_group.iloc[0]

        region = classify_region(
            c["ADFORCD"],
            state_abbrev="MT",
        )
        if region is None:
            continue

        basis = c["PROP_BASIS"]
        sub_plot_area = 0.2460 if (not pd.isna(basis) and str(basis).upper() == "MACR") else 0.0417
        condition_area = sub_plot_area * (4 * float(c["CONDPROP_UNADJ"]))

        ages = [c["STDAGE"], c["FLDAGE"]]
        ages = [float(a) for a in ages if pd.notna(a)]
        stand_age = max(ages) if ages else 0.0

        forest_type = _mog_forest_type_from_cond_row(c, region=region)

        adforcd_val = c["ADFORCD"]
        adforcd_int = None if pd.isna(adforcd_val) else int(adforcd_val)

        sitecl_val = c["SITECLCD"]
        sitecl_float = None if pd.isna(sitecl_val) else float(sitecl_val)

        ctx = ConditionContext(
            region=region,
            forest_type=forest_type,
            condition_area_acres=condition_area,
            stand_age=stand_age,
            trees=cond_group,
            plot_statecd=int(c["STATECD"]),
            condition_fortypcd=None,
            condition_physclcd=None,
            condition_siteclcd=sitecl_float,
            condition_adforcd=adforcd_int,
            condition_habtypcd1=c["HABTYPCD1"],
            ecosubcd=None,
        )

        mog_vec = engine.mog_vector(ctx)
        mog_score = max(mog_vec, default=0.0)

        out_rows.append(
            {
                "PLT_CN": plt_cn,
                "CONDID": condid,
                "MOG_SCORE": mog_score,
            }
        )

    return pd.DataFrame(out_rows)

if __name__ == "__main__":

    # db = read_fia("./tests/fiadb/AZ")

    # # Suppose `db` is your FiaDatabase containing at least PLOT, COND, TREE for AZ
    # mog_az = compute_mog_for_arizona(db)

    # # Merge back onto the condition table
    # mog_az.to_csv("./tests/az_mog.csv", index=False)

    # Montana (bundled FIADB: tests/fiadb/MT or tests/fiadb/mt), optional:
    db_mt = clip_fia(read_fia("./tests/fiadb/MT"))
    mog_mt = compute_mog_for_montana(db_mt)
    mog_mt.to_csv("./tests/mt_mog_v2.csv", index=False)
