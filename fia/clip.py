from __future__ import annotations

from typing import Optional

import pandas as pd

from .data_io import FiaDatabase


def clip_fia(db: FiaDatabase, most_recent: bool = True) -> FiaDatabase:
    """
    Lightweight Python analogue of rFIA::clipFIA for in-memory databases.

    Current behavior (TI-focused):
    - If `most_recent=True` and POP design tables are present, restricts the
      database to plots belonging to the most recent inventory year in each
      state (based on POP_EVAL.END_INVYR).
    - Filters:
        - POP_EVAL, POP_EVAL_TYP, POP_ESTN_UNIT, POP_STRATUM,
          POP_PLOT_STRATUM_ASSGN
        - PLOT, TREE, COND, PLOTGEOM
      to those EVALIDs / plots.
    - If required tables are missing, returns `db` unchanged.

    This is sufficient for workflows like the `fia_fse.R` script, where
    `clipFIA` is used with default arguments (`mostRecent = TRUE`) to obtain
    a “current” subset (`az_current`) before computing plot-level summaries
    and calling estimators.
    """

    # If no clipping requested, or db lacks design tables, return unchanged
    if not most_recent:
        return db

    table_names = set(db.keys())
    required_design = {
        "POP_EVAL",
        "POP_EVAL_TYP",
        "POP_ESTN_UNIT",
        "POP_STRATUM",
        "POP_PLOT_STRATUM_ASSGN",
    }
    if not required_design.issubset(table_names):
        return db

    # Work on copies of all tables
    tables = {name: df.copy() for name, df in db.items()}

    pop_eval = tables["POP_EVAL"]
    ppsa = tables["POP_PLOT_STRATUM_ASSGN"]

    # Basic column checks
    for col in ["EVALID", "STATECD", "END_INVYR"]:
        if col not in pop_eval.columns:
            return db
    if "PLT_CN" not in ppsa.columns:
        return db

    # Find most recent EVALID per state (END_INVYR max)
    latest_years = (
        pop_eval.groupby("STATECD", as_index=False)["END_INVYR"].max()
    )
    latest = pop_eval.merge(
        latest_years, on=["STATECD", "END_INVYR"], how="inner"
    )
    latest_evalids = latest["EVALID"].unique()

    # Filter POP tables to these EVALIDs
    tables["POP_EVAL"] = pop_eval[pop_eval["EVALID"].isin(latest_evalids)]

    pop_estn_unit = tables["POP_ESTN_UNIT"]
    if "EVALID" in pop_estn_unit.columns:
        tables["POP_ESTN_UNIT"] = pop_estn_unit[
            pop_estn_unit["EVALID"].isin(latest_evalids)
        ]

    pop_stratum = tables["POP_STRATUM"]
    if "EVALID" in pop_stratum.columns:
        tables["POP_STRATUM"] = pop_stratum[
            pop_stratum["EVALID"].isin(latest_evalids)
        ]

    pop_eval_typ = tables["POP_EVAL_TYP"]
    if "EVAL_CN" in pop_eval_typ.columns and "CN" in pop_eval.columns:
        valid_eval_cns = tables["POP_EVAL"]["CN"].unique()
        tables["POP_EVAL_TYP"] = pop_eval_typ[
            pop_eval_typ["EVAL_CN"].isin(valid_eval_cns)
        ]

    # Filter POP_PLOT_STRATUM_ASSGN to these EVALIDs and derive plot list
    ppsa_filtered = ppsa[ppsa["EVALID"].isin(latest_evalids)].copy()
    tables["POP_PLOT_STRATUM_ASSGN"] = ppsa_filtered
    valid_plots = ppsa_filtered["PLT_CN"].unique()

    # Helper to ensure PLT_CN exists in PLOT
    def _ensure_plt_cn(plot_df: pd.DataFrame) -> pd.DataFrame:
        if "PLT_CN" not in plot_df.columns and "CN" in plot_df.columns:
            plot_df = plot_df.copy()
            plot_df["PLT_CN"] = plot_df["CN"]
        return plot_df

    # Filter PLOT, TREE, COND, PLOTGEOM
    if "PLOT" in tables:
        plot_df = _ensure_plt_cn(tables["PLOT"])
        tables["PLOT"] = plot_df[plot_df["PLT_CN"].isin(valid_plots)].copy()

    if "TREE" in tables and "PLT_CN" in tables["TREE"].columns:
        tables["TREE"] = tables["TREE"][
            tables["TREE"]["PLT_CN"].isin(valid_plots)
        ].copy()

    if "COND" in tables and "PLT_CN" in tables["COND"].columns:
        tables["COND"] = tables["COND"][
            tables["COND"]["PLT_CN"].isin(valid_plots)
        ].copy()

    if "PLOTGEOM" in tables and "CN" in tables["PLOTGEOM"].columns:
        tables["PLOTGEOM"] = tables["PLOTGEOM"][
            tables["PLOTGEOM"]["CN"].isin(valid_plots)
        ].copy()

    return FiaDatabase(tables)


__all__ = ["clip_fia"]

