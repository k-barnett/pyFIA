from __future__ import annotations

from typing import Iterable, List, Sequence

import pandas as pd

from .data_io import FiaDatabase, normalize_fiadb_dataframe_columns


def combine_mr(df: pd.DataFrame, year_col: str = "YEAR") -> pd.DataFrame:
    """
    Combine most‑recent population estimates across states with potentially
    different reporting schedules.

    Direct analogue of `combineMR()` in rFIA::util.R: it replaces the YEAR
    column with the maximum YEAR found in the data frame.
    """
    if year_col not in df.columns:
        raise ValueError(f"'{year_col}' column not found in DataFrame.")
    df = df.copy()
    df[year_col] = df[year_col].max()
    return df


def ratio_var(
    x: pd.Series,
    y: pd.Series,
    x_var: pd.Series,
    y_var: pd.Series,
    cv: pd.Series,
) -> pd.Series:
    """
    Variance of a ratio estimator, ported from `ratioVar()` in rFIA.

    In R:

        r.var <- (1 / (y^2)) * (x.var + ((x/y)^2 * y.var) - (2 * (x/y) * cv))
        r.var[r.var < 0] <- 0
    """
    r_var = (1.0 / (y**2)) * (x_var + ((x / y) ** 2) * y_var - 2.0 * (x / y) * cv)
    # Guard against small negative values from numerical rounding
    r_var = r_var.where(r_var >= 0, 0.0)
    return r_var


def get_design_info(
    db: FiaDatabase,
    eval_types: Sequence[str] = ("VOL",),
    most_recent: bool = True,
) -> pd.DataFrame:
    """
    Python port of `getDesignInfo()` for in‑memory FIA databases.

    This focuses on the common use case needed by estimators like `tpa`
    under the temporally‑indifferent (TI) method. It expects that `db`
    contains the POP_* tables used by rFIA:

        POP_EVAL, POP_EVAL_TYP, POP_ESTN_UNIT, POP_STRATUM, POP_PLOT_STRATUM_ASSGN
    """

    required = {
        "PLOT",
        "POP_EVAL",
        "POP_EVAL_TYP",
        "POP_ESTN_UNIT",
        "POP_STRATUM",
        "POP_PLOT_STRATUM_ASSGN",
    }
    missing = required.difference(set(db.keys()))
    if missing:
        raise ValueError(
            "Missing required tables in `db` for design info: "
            + ", ".join(sorted(missing))
        )

    pop_eval = normalize_fiadb_dataframe_columns(db["POP_EVAL"].copy())
    pop_eval_typ = normalize_fiadb_dataframe_columns(db["POP_EVAL_TYP"].copy())
    pop_estn_unit = normalize_fiadb_dataframe_columns(db["POP_ESTN_UNIT"].copy())
    pop_stratum = normalize_fiadb_dataframe_columns(db["POP_STRATUM"].copy())
    pop_plot_stratum = normalize_fiadb_dataframe_columns(
        db["POP_PLOT_STRATUM_ASSGN"].copy()
    )

    # Basic column checks (abbreviated)
    for col in ["CN", "STATECD", "END_INVYR", "EVALID", "ESTN_METHOD"]:
        if col not in pop_eval.columns:
            raise ValueError(f"`POP_EVAL` must contain column '{col}'.")
    for col in ["EVAL_CN", "EVAL_TYP"]:
        if col not in pop_eval_typ.columns:
            raise ValueError(f"`POP_EVAL_TYP` must contain column '{col}'.")

    # EVAL table with type attached
    evals = (
        pop_eval.loc[:, ["CN", "STATECD", "END_INVYR", "EVALID", "ESTN_METHOD"]]
        .rename(columns={"CN": "EVAL_CN", "END_INVYR": "YEAR"})
        .merge(pop_eval_typ.loc[:, ["EVAL_CN", "EVAL_TYP"]], on="EVAL_CN", how="left")
    )
    evals = evals[evals["EVAL_TYP"].notna()].copy()

    # Filter by evaluation type(s). Accept both EXPCURR and CURR (etc.) so we
    # match evals regardless of DB naming (rFIA uses type = 'CURR' → EXPCURR).
    eval_types = list({t.upper() for t in eval_types})
    wanted = set()
    for t in eval_types:
        wanted.add(f"EXP{t}")
        wanted.add(t)
    evals["_EVAL_TYP_upper"] = evals["EVAL_TYP"].astype(str).str.upper()
    evals = evals[evals["_EVAL_TYP_upper"].isin(wanted)].copy()
    evals = evals.drop(columns=["_EVAL_TYP_upper"], errors="ignore")

    # Drop early periodic inventories (pre‑2003)
    evals = evals[evals["YEAR"] >= 2003].copy()

    # Optionally keep only most‑recent inventory per state
    if most_recent:
        evals = (
            evals.sort_values("YEAR")
            .groupby(["STATECD", "EVAL_TYP"], as_index=False)
            .tail(1)
        )

    # Join estimation units
    for col in ["CN", "P1PNTCNT_EU", "AREA_USED", "EVAL_CN"]:
        if col not in pop_estn_unit.columns:
            raise ValueError(f"`POP_ESTN_UNIT` must contain column '{col}'.")
    strata = evals.merge(
        pop_estn_unit.loc[:, ["CN", "P1PNTCNT_EU", "AREA_USED", "EVAL_CN"]],
        on="EVAL_CN",
        how="left",
        suffixes=("", "_EU"),
    ).rename(columns={"CN": "ESTN_UNIT_CN"})

    # Join stratum info
    needed_stratum_cols = [
        "CN",
        "ESTN_UNIT_CN",
        "P1POINTCNT",
        "P2POINTCNT",
        "ADJ_FACTOR_MICR",
        "ADJ_FACTOR_SUBP",
        "ADJ_FACTOR_MACR",
    ]
    for col in needed_stratum_cols:
        if col not in pop_stratum.columns:
            raise ValueError(f"`POP_STRATUM` must contain column '{col}'.")
    strata = strata.merge(
        pop_stratum.loc[:, needed_stratum_cols],
        on="ESTN_UNIT_CN",
        how="left",
        suffixes=("", "_STR"),
    ).rename(columns={"CN": "STRATUM_CN"})

    # Stratum weights within estimation units
    strata["STRATUM_WGT"] = strata["P1POINTCNT"] / strata["P1PNTCNT_EU"]

    # Join plots to stratum
    needed_ppsa_cols = ["STRATUM_CN", "PLT_CN", "UNITCD", "COUNTYCD", "PLOT"]
    for col in needed_ppsa_cols:
        if col not in pop_plot_stratum.columns:
            raise ValueError(
                f"`POP_PLOT_STRATUM_ASSGN` must contain column '{col}'."
            )

    strata = strata.merge(
        pop_plot_stratum.loc[:, needed_ppsa_cols],
        on="STRATUM_CN",
        how="left",
    )

    # Add plot identifier (STATECD is in evals/strata already)
    strata["pltID"] = (
        strata["UNITCD"].astype("Int64").astype("string")
        + "_"
        + strata["STATECD"].astype("Int64").astype("string")
        + "_"
        + strata["COUNTYCD"].astype("Int64").astype("string")
        + "_"
        + strata["PLOT"].astype("Int64").astype("string")
    )

    # Minimal column subset approximating rFIA::getDesignInfo output
    # P2POINTCNT: second-phase plot count in stratum (stratum mean uses sum(y)/P2, not mean(y))
    keep_cols = [
        "STATECD",
        "YEAR",
        "EVAL_TYP",
        "EVALID",
        "ESTN_METHOD",
        "ESTN_UNIT_CN",
        "P1PNTCNT_EU",
        "AREA_USED",
        "STRATUM_CN",
        "P2POINTCNT",
        "P1POINTCNT",
        "ADJ_FACTOR_MICR",
        "ADJ_FACTOR_SUBP",
        "ADJ_FACTOR_MACR",
        "STRATUM_WGT",
        "pltID",
        "PLT_CN",
    ]
    strata = strata[keep_cols].drop_duplicates().reset_index(drop=True)

    # Sum number of plots per estimation unit
    p2eu = (
        strata[["ESTN_UNIT_CN", "STRATUM_CN", "P2POINTCNT"]]
        .drop_duplicates()
        .groupby("ESTN_UNIT_CN", as_index=False)["P2POINTCNT"]
        .sum()
        .rename(columns={"P2POINTCNT": "P2PNTCNT_EU"})
    )
    strata = strata.merge(p2eu, on="ESTN_UNIT_CN", how="left")

    return strata


def handle_pops(
    db: FiaDatabase,
    eval_types: Sequence[str],
    method: str = "TI",
    most_recent: bool = True,
) -> pd.DataFrame:
    """
    Simplified Python port of `handlePops()` for the temporally‑indifferent (TI)
    case.

    For non‑TI methods (SMA, EMA, etc.), the R code adds additional INVYR‑level
    columns that are not yet implemented here.
    """
    # For the TI case we can delegate "most‑recent" handling to get_design_info,
    # which already mirrors rFIA::getDesignInfo for EXP* eval types.
    pops = get_design_info(db, eval_types=eval_types, most_recent=most_recent)

    # Recode ESTN_METHOD to match rFIA's simplified labels
    method_map = {
        "Post-Stratification": "strat",
        "Stratified random sampling": "strat",
        "Simple random sampling": "simple",
        "Subsampling units of unequal size": "simple",
        "Double sampling for stratification": "double",
    }
    pops = pops.copy()
    pops["ESTN_METHOD"] = pops["ESTN_METHOD"].map(
        lambda m: method_map.get(m, "double")
    )

    # For TI, we do not modify P2POINTCNT / P2PNTCNT_EU further here.
    return pops


__all__ = ["combine_mr", "ratio_var", "get_design_info", "handle_pops"]


