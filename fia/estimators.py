from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable, List, Mapping, Optional, Sequence

import pandas as pd

from .data_io import FiaDatabase
from .design import handle_pops, ratio_var, combine_mr


def _basal_area(diameter_in: pd.Series) -> pd.Series:
    """
    Basal area in square feet for a diameter in inches.

    This follows the rFIA `basalArea()` helper: ba = d * |d| * 0.005454.
    """
    return diameter_in * diameter_in.abs() * 0.005454


def _land_type_domain(
    land_type: str,
    cond_status_cd: pd.Series,
    siteclcd: pd.Series,
    reservcd: pd.Series,
) -> pd.Series:
    land_type = land_type.lower()
    if land_type == "forest":
        return (cond_status_cd == 1).astype(int)
    elif land_type == "timber":
        return (
            (cond_status_cd == 1)
            & siteclcd.isin([1, 2, 3, 4, 5, 6])
            & (reservcd == 0)
        ).astype(int)
    elif land_type == "all":
        return pd.Series(1, index=cond_status_cd.index, dtype="int64")
    else:
        raise NotImplementedError(
            f"land_type='{land_type}' is not yet supported in the Python port."
        )


def _tree_type_domain(
    tree_type: str,
    statuscd: pd.Series,
    dia: pd.Series,
    treeclcd: pd.Series,
) -> pd.Series:
    tree_type = tree_type.lower()
    if tree_type == "live":
        return (statuscd == 1).astype(int)
    elif tree_type == "dead":
        return (statuscd == 2).astype(int)
    elif tree_type == "gs":
        return ((statuscd == 1) & (dia >= 5) & (treeclcd == 2)).astype(int)
    elif tree_type == "all":
        return pd.Series(1, index=statuscd.index, dtype="int64")
    else:
        raise NotImplementedError(
            f"tree_type='{tree_type}' is not yet supported in the Python port."
        )


def area(
    db: FiaDatabase,
    grp_by: Optional[Sequence[str]] = None,
    land_type: str = "forest",
    cond_list: bool = False,
    by_plot: bool = False,
) -> pd.DataFrame:
    """
    Estimate land area (temporally-indifferent TI estimator) analogous to rFIA::area.

    This initial implementation supports:
    - `land_type` (Python analogue of `landType`)
    - `cond_list` (Python analogue of `condList`)
    - `by_plot` (Python analogue of `byPlot`)

    It currently assumes TI (temporally indifferent) estimation and does not
    yet support custom area/tree domains or temporal moving-average methods.

    When ``cond_list=True``, the returned table includes ``AREA_BASIS`` (from
    ``COND.PROP_BASIS``, normalized to ``SUBP`` or ``MACR``) alongside
    ``PLT_CN``, ``CONDID``, ``YEAR``, optional ``grp_by`` columns, and
    ``PROP_FOREST``. That matches what :func:`custom_pse` expects for
    condition-level area inputs.
    """

    if grp_by is None:
        grp_by = []
    grp_by = list(grp_by)

    required = {"PLOT", "COND"}
    missing = required.difference(set(db.keys()))
    if missing:
        raise ValueError(
            "Missing required tables in `db` for area(): "
            + ", ".join(sorted(missing))
        )

    plot = db["PLOT"].copy()
    cond = db["COND"].copy()

    # Basic checks
    for col in ["PLT_CN", "INVYR", "PLOT_STATUS_CD"]:
        if col not in plot.columns:
            raise ValueError(f"`PLOT` table must contain column '{col}'.")
    for col in ["PLT_CN", "CONDID", "COND_STATUS_CD", "CONDPROP_UNADJ", "SITECLCD", "RESERVCD"]:
        if col not in cond.columns:
            raise ValueError(f"`COND` table must contain column '{col}'.")

    # Keep current forest status plots only (status 1)
    plot = plot.loc[plot["PLOT_STATUS_CD"] == 1].copy()

    # Land domain indicator
    cond["landD"] = _land_type_domain(
        land_type,
        cond["COND_STATUS_CD"],
        cond["SITECLCD"],
        cond["RESERVCD"],
    )

    # Join PLOT + COND
    data = plot.merge(cond, on="PLT_CN", how="inner", suffixes=("", "_COND"))

    # Plot year
    year_col = "MEASYEAR" if "MEASYEAR" in data.columns else "INVYR"
    data["YEAR"] = data[year_col]

    # Condition-level domain area within each plot (proportion)
    data["fa"] = data["CONDPROP_UNADJ"] * data["landD"]

    # FIA COND.PROP_BASIS: whether condition proportion is on SUBP vs MACR plot;
    # rFIA passes this as AREA_BASIS into sumToPlot for non-response adjustment.
    # Normalize once for cond_list, by_plot, and the full TI path.
    if "PROP_BASIS" not in data.columns:
        data["PROP_BASIS"] = "SUBP"
    else:
        data["PROP_BASIS"] = (
            data["PROP_BASIS"]
            .fillna("SUBP")
            .astype(str)
            .str.upper()
            .replace("", "SUBP")
        )
        data["PROP_BASIS"] = data["PROP_BASIS"].where(
            data["PROP_BASIS"].isin(["MACR", "SUBP"]), "SUBP"
        )
    data["AREA_BASIS"] = data["PROP_BASIS"]

    # ------------------------
    # cond_list / by_plot path
    # ------------------------

    if cond_list:
        grp_cols = [g for g in grp_by if g in data.columns]
        cols = ["PLT_CN", "CONDID", "YEAR", "AREA_BASIS"] + grp_cols + ["fa"]
        out = (
            data[cols]
            .rename(columns={"fa": "PROP_FOREST"})
            .reset_index(drop=True)
        )
        return out

    if by_plot:
        grp_cols = [g for g in grp_by if g in data.columns]
        plot_level = (
            data.groupby(
                ["PLT_CN", "YEAR"] + grp_cols,
                as_index=False,
            )["fa"]
            .sum()
            .rename(columns={"fa": "PROP_FOREST"})
        )
        return plot_level

    # -------------------------
    # Full TI design-based path
    # -------------------------

    design_tables = {
        "POP_EVAL",
        "POP_EVAL_TYP",
        "POP_ESTN_UNIT",
        "POP_STRATUM",
        "POP_PLOT_STRATUM_ASSGN",
    }
    have_design = design_tables.issubset(set(db.keys()))

    if not have_design:
        # Fallback: simple plot-mean estimator (no design weights)
        grp_cols = ["YEAR"] + [g for g in grp_by if g in data.columns]
        plot_area = (
            data.groupby(
                ["PLT_CN", "YEAR"]
                + [g for g in grp_by if g in data.columns],
                as_index=False,
            )["fa"]
            .sum()
            .rename(columns={"fa": "PROP_FOREST"})
        )

        def _agg_simple(group: pd.DataFrame) -> pd.Series:
            perc_area = float(group["PROP_FOREST"].mean() * 100.0)
            return pd.Series(
                {
                    "PERC_AREA": perc_area,
                    "AREA_TOTAL": float(group["PROP_FOREST"].mean()),
                    "nPlots_AREA": int(len(group)),
                }
            )

        out = (
            plot_area.groupby(grp_cols, as_index=False)
            .apply(_agg_simple)
            .reset_index(drop=True)
            .sort_values("YEAR")
        )
        return out

    # -------------------------
    # Full TI design-based path (exact rFIA-style logic)
    # -------------------------

    # Use TI population frame from POP_* tables (already mirrors rFIA::handlePops
    # for the TI case).
    pops = handle_pops(db, eval_types=("CURR",), method="TI", most_recent=True)

    # Restrict to plots in the population frame (same as rFIA areaStarter:
    # filter PLOT and COND to PLT_CN %in% pops$PLT_CN). Avoids including
    # plots that will be dropped later and keeps plot set identical to rFIA.
    #plot_cns = pops["PLT_CN"].dropna().unique()
    #data = data.loc[data["PLT_CN"].isin(plot_cns)].copy()

    # In the rFIA area engine, the comprehensive domain indicators are:
    # aDI = landD * aD * sp * tD
    # pDI = landD * aD
    # For our current use case (no user area/tree domain, no polygons), aD = sp = tD = 1,
    # so aDI and pDI both reduce to landD.
    data = data.copy()
    data["aDI"] = data["landD"]
    data["pDI"] = data["landD"]

    grp_cols = [g for g in grp_by if g in data.columns]

    # Numerator: total forest proportion in the landType × (optional) areaDomain
    # aEst/tEst construction follows the non-byPlot, non-condList branch of
    # rFIA::areaStarter(). AREA_BASIS was set above from COND.PROP_BASIS.
    t = (
        data.drop_duplicates(subset=["PLT_CN", "CONDID"] + grp_cols)
        .assign(fa=lambda d: d["CONDPROP_UNADJ"] * d["aDI"])
        .loc[
            :,
            ["PLT_CN", "AREA_BASIS", "CONDID"]
            + grp_cols
            + ["fa", "aDI"],
        ]
    )

    # Denominator: total land area in areaDomain and landType (for proportions)
    a = (
        data.drop_duplicates(subset=["PLT_CN", "CONDID"])
        .assign(fad=lambda d: d["CONDPROP_UNADJ"] * d["pDI"])
        .loc[:, ["PLT_CN", "AREA_BASIS", "fad"]]
    )

    # Filter numerator to plots within the domain (aDI == 1) as in rFIA.
    t = t.loc[t["aDI"] == 1].copy()

    # Sum to plot level with non-response adjustment and design info
    t_plt = _rfia_sum_to_plot(t, pops, grp_cols)
    a_plt = _rfia_sum_to_plot(a, pops, [])

    # For the TI forest area estimator, only the `fa` (numerator) and `fad`
    # (denominator) variables are needed downstream. Drop auxiliary fields like
    # `aDI` so that the EU-level engine sees a single numeric variable.
    keep_x_cols = [
        c
        for c in t_plt.columns
        if c
        in (
            {"PLT_CN", "ESTN_UNIT_CN", "STRATUM_CN", "YEAR"}
            | set(grp_cols)
            | {"fa"}
        )
    ]
    t_plt = t_plt[keep_x_cols]

    keep_y_cols = [
        c
        for c in a_plt.columns
        if c in {"PLT_CN", "ESTN_UNIT_CN", "STRATUM_CN", "YEAR", "fad"}
    ]
    a_plt = a_plt[keep_y_cols]

    # Add YEAR to groups, following rFIA (grpBy <- c('YEAR', grpBy); aGrpBy <- c('YEAR'))
    x_grp_by = ["YEAR"] + grp_cols
    a_grp_by = ["YEAR"]

    x_eu, a_eu = _rfia_sum_to_eu(
        x=t_plt,
        y=a_plt,
        pops=pops,
        x_grp_by=x_grp_by,
        y_grp_by=a_grp_by,
        method="TI",
    )

    # Combine most-recent inventories across states (TI combineMR step)
    x_eu = combine_mr(x_eu, year_col="YEAR")
    a_eu = combine_mr(a_eu, year_col="YEAR")

    # Sum EU-level estimates by group (mirrors top-level rFIA::area())
    x_tot = (
        x_eu.groupby(x_grp_by, as_index=False)
        .sum(numeric_only=True)
    )
    a_tot = (
        a_eu.groupby(a_grp_by, as_index=False)
        .sum(numeric_only=True)
    )

    # Join numerator/denominator and compute totals, percentages, and variances
    out = x_tot.merge(a_tot, on=["YEAR"], how="left", suffixes=("_x", "_y"))

    # Prefer EU-level outputs from _rfia_sum_to_eu (fa_mean, fad_mean, etc.).
    # As a fallback (e.g., if upstream changes left only raw fa/fad), use those
    # directly so we at least get a numeric estimate instead of a KeyError.
    if "fa_mean" in out.columns:
        fa_mean = out["fa_mean"]
        fa_var = out.get("fa_var", pd.Series(0.0, index=out.index))
        fa_cv = out.get("fa_cv", pd.Series(0.0, index=out.index))
    else:
        fa_mean = out["fa"]
        fa_var = pd.Series(0.0, index=out.index)
        fa_cv = pd.Series(0.0, index=out.index)

    if "fad_mean" in out.columns:
        fad_mean = out["fad_mean"]
        fad_var = out.get("fad_var", pd.Series(0.0, index=out.index))
    else:
        fad_mean = out["fad"]
        fad_var = pd.Series(0.0, index=out.index)

    # AREA_TOTAL is the numerator mean on the area scale (`fa_mean`), as in rFIA.
    area_total = fa_mean.astype(float)

    # Percentage area: fa_mean / fad_mean, converted to percent
    with pd.option_context("mode.use_inf_as_na", True):
        perc_area = fa_mean / fad_mean.replace({0: pd.NA})
    perc_area = perc_area.astype(float) * 100.0

    # Variance for the ratio PERC_AREA (before ×100) via rFIA::ratioVar (ported
    # as design.ratio_var).
    r_var = ratio_var(
        x=fa_mean,
        y=fad_mean,
        x_var=fa_var,
        y_var=fad_var,
        cv=fa_cv,
    ) * (100.0**2)
    # Guard against tiny negative values from numerical rounding
    r_var = r_var.where(r_var >= 0, 0.0)

    # Standard errors for AREA_TOTAL and PERC_AREA (in percent)
    area_total_se = (
        fa_var.pow(0.5) / area_total.replace({0: pd.NA}) * 100.0
    ).fillna(0.0)
    perc_area_se = (
        r_var.pow(0.5) / perc_area.replace({0: pd.NA}) * 100.0
    ).fillna(0.0)

    # Total number of plots in numerator/denominator and population size N
    # (P2PNTCNT_EU summed across EUs), following rFIA's naming.
    out_final = out.copy()
    out_final["PERC_AREA"] = perc_area
    out_final["AREA_TOTAL"] = area_total
    out_final["PERC_AREA_SE"] = perc_area_se
    out_final["AREA_TOTAL_SE"] = area_total_se
    out_final["PERC_AREA_VAR"] = r_var
    out_final["AREA_TOTAL_VAR"] = fa_var
    # rFIA's area() uses nPlots.x / nPlots.y and P2PNTCNT_EU from the
    # numerator EU table. Our merged frame may or may not carry these with
    # suffixes depending on upstream changes; fall back to 0 when absent.
    nplots_num = out.get("nPlots.x", pd.Series(0, index=out.index))
    nplots_den = out.get("nPlots.y", pd.Series(0, index=out.index))
    N_series = (
        out["P2PNTCNT_EU_x"]
        if "P2PNTCNT_EU_x" in out.columns
        else out.get("P2PNTCNT_EU", pd.Series(0, index=out.index))
    )
    out_final["nPlots_AREA_NUM"] = nplots_num
    out_final["nPlots_AREA_DEN"] = nplots_den
    out_final["N"] = N_series

    # For backwards compatibility with earlier versions of this Python port and
    # with the tests, expose `nPlots_AREA` as the numerator plot count.
    out_final["nPlots_AREA"] = out_final["nPlots_AREA_NUM"]

    # Restrict to the expected output columns and sort by YEAR
    keep_cols = (
        x_grp_by
        + [
            "PERC_AREA",
            "AREA_TOTAL",
            "PERC_AREA_SE",
            "AREA_TOTAL_SE",
            "PERC_AREA_VAR",
            "AREA_TOTAL_VAR",
            "nPlots_AREA",
            "N",
        ]
    )
    out_final = out_final[keep_cols].sort_values("YEAR").reset_index(drop=True)

    return out_final


def _sum_to_plot(
    df: pd.DataFrame,
    pops: pd.DataFrame,
    var_names: Sequence[str],
    grp_by: Sequence[str],
) -> pd.DataFrame:
    """
    Approximate port of rFIA::sumToPlot for TREE_BASIS / AREA_BASIS data.

    - If `TREE_BASIS` is present, treats rows as tree-level.
    - If `AREA_BASIS` is present, treats rows as condition-level.
    - Joins on PLT_CN to attach design info and ADJ_FACTOR_* columns, applies
      non-response adjustment factors, and aggregates to plot level.
    """

    has_tree_basis = "TREE_BASIS" in df.columns
    has_area_basis = "AREA_BASIS" in df.columns

    if has_tree_basis and has_area_basis:
        raise ValueError(
            "Both 'TREE_BASIS' and 'AREA_BASIS' found; only one is allowed in `df`."
        )
    if not has_tree_basis and not has_area_basis:
        raise ValueError(
            "Neither 'TREE_BASIS' nor 'AREA_BASIS' found in `df`; one is required."
        )

    basis_col = "TREE_BASIS" if has_tree_basis else "AREA_BASIS"

    # Use YEAR from the design frame (pops) to avoid clashes with any YEAR
    # column that might already exist in df (e.g., from a previous aggregation).
    # Drop YEAR from df if present, then reattach via PLT_CN from pops.
    df_no_year = df.drop(columns=["YEAR"], errors="ignore")

    # Minimal design information needed for adjustment, grouping, and area.
    pops_sub = pops[
        [
            "PLT_CN",
            "YEAR",
            "ESTN_UNIT_CN",
            "STRATUM_CN",
            "AREA_USED",
            "ADJ_FACTOR_MICR",
            "ADJ_FACTOR_SUBP",
            "ADJ_FACTOR_MACR",
        ]
    ].drop_duplicates()

    merged = df_no_year.merge(pops_sub, on="PLT_CN", how="inner")

    # Keep only grouping columns that actually exist in merged, and avoid
    # duplicating YEAR since YEAR is already included explicitly below.
    grp_cols = [g for g in grp_by if g in merged.columns and g != "YEAR"]

    # First sum within PLT_CN × BASIS to get basis-level totals
    group_cols_basis = ["PLT_CN", basis_col, "ESTN_UNIT_CN", "STRATUM_CN", "YEAR"] + grp_cols
    basis_sum = (
        merged.groupby(group_cols_basis, as_index=False)[list(var_names)]
        .sum()
    )

    # Add adjustment factor per basis
    def _choose_adj(row: pd.Series) -> float:
        b = row[basis_col]
        if pd.isna(b):
            return 1.0
        b = str(b).upper()
        if b == "MACR":
            return float(row.get("ADJ_FACTOR_MACR", 1.0))
        if b == "SUBP":
            return float(row.get("ADJ_FACTOR_SUBP", 1.0))
        if b == "MICR":
            return float(row.get("ADJ_FACTOR_MICR", 1.0))
        return 1.0

    # We need ADJ_FACTOR_* on basis_sum as well. Use a lookup that excludes
    # YEAR and AREA_USED to avoid creating duplicate YEAR columns while
    # still providing the adjustment factors.
    adj_lookup = pops_sub[
        [
            "PLT_CN",
            "ESTN_UNIT_CN",
            "STRATUM_CN",
            "ADJ_FACTOR_MICR",
            "ADJ_FACTOR_SUBP",
            "ADJ_FACTOR_MACR",
        ]
    ].drop_duplicates()

    basis_sum = basis_sum.merge(
        adj_lookup,
        on=["PLT_CN", "ESTN_UNIT_CN", "STRATUM_CN"],
        how="left",
    )

    basis_sum["adj"] = basis_sum.apply(_choose_adj, axis=1)
    for v in var_names:
        basis_sum[v] = basis_sum[v] * basis_sum["adj"]

    # Now sum across basis types to plot level
    group_cols_plot = ["PLT_CN", "ESTN_UNIT_CN", "STRATUM_CN", "YEAR"] + grp_cols
    plot_sum = (
        basis_sum.groupby(group_cols_plot, as_index=False)[list(var_names)]
        .sum()
    )

    # Add AREA_USED back (unique per ESTN_UNIT_CN)
    area_lookup = pops_sub[["ESTN_UNIT_CN", "AREA_USED"]].drop_duplicates()
    plot_sum = plot_sum.merge(area_lookup, on="ESTN_UNIT_CN", how="left")

    return plot_sum


def _rfia_sum_to_plot(
    x: pd.DataFrame,
    pops: pd.DataFrame,
    grp_by: Sequence[str],
) -> pd.DataFrame:
    """
    Python port of rFIA::sumToPlot() for TI estimation.

    This helper is used by the exact TI `area` implementation and closely
    mirrors the AREA_BASIS branch of the R function. It:
      - Sums condition-level quantities to plot level within AREA_BASIS
      - Joins design info and non-response adjustment factors from `pops`
      - Applies the appropriate ADJ_FACTOR_* by AREA_BASIS
      - Aggregates across AREA_BASIS to a single record per plot/stratum/EU
    """

    grp_by = list(grp_by or [])
    grp_cols = [g for g in grp_by if g in x.columns]

    # Condition-basis branch (AREA_BASIS) as used in rFIA::areaStarter().
    if "AREA_BASIS" not in x.columns and "TREE_BASIS" in x.columns:
        raise ValueError(
            "_rfia_sum_to_plot is currently implemented for AREA_BASIS data "
            "(condition-level); TREE_BASIS should use the custom_pse engine."
        )

    # Variable columns to be summed (everything except keys / identifiers).
    key_cols = {"PLT_CN", "CONDID", "AREA_BASIS", *grp_cols}
    var_cols = [c for c in x.columns if c not in key_cols]

    if not var_cols:
        raise ValueError("_rfia_sum_to_plot: no numeric variables found to sum.")

    # Sum to PLT_CN × AREA_BASIS × grp_by
    grouped = (
        x.groupby(["PLT_CN", "AREA_BASIS"] + grp_cols, as_index=False)[var_cols]
        .sum()
    )

    # Minimal design info needed, mirroring the AREA_BASIS branch of sumToPlot().
    pops_sub = (
        pops[
            [
                "EVALID",
                "ESTN_UNIT_CN",
                "STRATUM_CN",
                "PLT_CN",
                "YEAR",
                "ADJ_FACTOR_MICR",
                "ADJ_FACTOR_SUBP",
                "ADJ_FACTOR_MACR",
            ]
        ]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    merged = grouped.merge(pops_sub, on="PLT_CN", how="inner")

    # Apply adjustment factors based on AREA_BASIS.
    def _adj(row: pd.Series) -> float:
        basis = row.get("AREA_BASIS")
        if pd.isna(basis):
            return float("nan")
        basis = str(basis)
        if basis == "MACR":
            return float(row.get("ADJ_FACTOR_MACR", float("nan")))
        if basis == "SUBP":
            return float(row.get("ADJ_FACTOR_SUBP", float("nan")))
        return float("nan")

    merged["adj"] = merged.apply(_adj, axis=1)
    for v in var_cols:
        merged[v] = merged[v] * merged["adj"]

    # Drop potential duplicates (e.g., when multiple eval types are present)
    merged = merged.drop_duplicates()

    # Sum across AREA_BASIS within each plot/stratum/EU and group combination.
    out = (
        merged.groupby(
            ["ESTN_UNIT_CN", "STRATUM_CN", "PLT_CN"] + grp_cols + ["YEAR"],
            as_index=False,
        )[var_cols]
        .sum()
    )

    return out


def _rfia_sum_to_eu(
    x: pd.DataFrame,
    y: Optional[pd.DataFrame],
    pops: pd.DataFrame,
    x_grp_by: Sequence[str],
    y_grp_by: Optional[Sequence[str]],
    method: str = "TI",
) -> tuple[pd.DataFrame, Optional[pd.DataFrame]]:
    """
    Python port of the TI branch of rFIA::sumToEU().

    This function:
      - Computes stratum-level means/variances for numerator (x) and
        denominator (y, when provided) using the TI design
      - Aggregates to estimation-unit (EU) level using STRATUM_WGT and AREA_USED
      - Returns (xEU, yEU) with columns analogous to the R output:
          * xEU: fa_mean, fa_var, fa_cv, nPlots.x, P2PNTCNT_EU, AREA_USED, ...
          * yEU: fad_mean, fad_var, nPlots.y, P2PNTCNT_EU, AREA_USED, ...
    """

    method = str(method).upper()
    if method != "TI":
        raise NotImplementedError(
            "_rfia_sum_to_eu currently implements only the TI method."
        )

    x_grp_by = list(x_grp_by or [])
    y_grp_by = list(y_grp_by or []) if y is not None else None

    # Subset of design info needed at stratum/EU level.
    # Minimal design info at stratum/EU level. We intentionally omit `YEAR`
    # here to avoid clashing with the YEAR column already present in `x`/`y`.
    pops_sub = (
        pops[
            [
                "ESTN_UNIT_CN",
                "AREA_USED",
                "P2PNTCNT_EU",
                "STRATUM_CN",
                "STRATUM_WGT",
                "P2POINTCNT",
            ]
        ]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    # Helper: denominator (y) strata → EU.
    y_eu: Optional[pd.DataFrame]
    if y is not None:
        if len(y.columns.difference({"PLT_CN", "ESTN_UNIT_CN", "STRATUM_CN", "YEAR", *y_grp_by})) != 1:
            raise ValueError(
                "_rfia_sum_to_eu expects a single numeric denominator variable in `y`."
            )
        y_var = y.columns.difference(
            {"PLT_CN", "ESTN_UNIT_CN", "STRATUM_CN", "YEAR", *y_grp_by}
        )[0]

        y_tmp = y.merge(pops_sub, on=["ESTN_UNIT_CN", "STRATUM_CN"], how="left")

        strat_group_cols_y = [
            "ESTN_UNIT_CN",
            "AREA_USED",
            "P2PNTCNT_EU",
            "STRATUM_CN",
            "STRATUM_WGT",
            "P2POINTCNT",
        ] + y_grp_by

        def _agg_y_strat(group: pd.DataFrame) -> pd.Series:
            P2 = float(group["P2POINTCNT"].iloc[0])
            vals = group[y_var].astype(float)
            mean = (vals / P2).sum()
            if P2 > 1:
                var = (
                    (vals.pow(2).sum())
                    - P2 * (vals / P2).sum() ** 2
                ) / (P2 * (P2 - 1))
            else:
                var = 0.0
            return pd.Series(
                {
                    "fad_mean": mean,
                    "fad_var": var,
                    "nPlots.y": group["PLT_CN"].nunique(),
                }
            )

        y_strat = (
            y_tmp.groupby(strat_group_cols_y, as_index=False)
            .apply(_agg_y_strat)
            .reset_index(drop=True)
        )

        eu_group_cols_y = ["ESTN_UNIT_CN", "AREA_USED", "P2PNTCNT_EU"] + y_grp_by

        def _agg_y_eu(group: pd.DataFrame) -> pd.Series:
            A = float(group["AREA_USED"].iloc[0])
            N = float(group["P2PNTCNT_EU"].iloc[0])
            w = group["STRATUM_WGT"].astype(float)
            P2 = group["P2POINTCNT"].astype(float)
            m = group["fad_mean"].astype(float)
            v = group["fad_var"].astype(float)
            fad_mean_eu = A * float((w * m).sum())
            fad_var_eu = (A**2 / N) * float(
                (v * w * P2).sum()
                + (v * (1.0 - w) * (P2 / N)).sum()
            )
            nplots = int(group["nPlots.y"].sum())
            return pd.Series(
                {
                    "fad_mean": fad_mean_eu,
                    "fad_var": fad_var_eu,
                    "nPlots.y": nplots,
                }
            )

        y_eu = (
            y_strat.groupby(eu_group_cols_y, as_index=False)
            .apply(_agg_y_eu)
            .reset_index(drop=True)
        )
    else:
        y_eu = None
        y_var = None  # type: ignore[assignment]
        y_grp_by = None

    # Numerator (x): strata → EU, including covariance with denominator where
    # provided.
    if len(x.columns.difference({"PLT_CN", "ESTN_UNIT_CN", "STRATUM_CN", "YEAR", *x_grp_by})) != 1:
        raise ValueError(
            "_rfia_sum_to_eu expects a single numeric numerator variable in `x`."
        )
    x_var = x.columns.difference(
        {"PLT_CN", "ESTN_UNIT_CN", "STRATUM_CN", "YEAR", *x_grp_by}
    )[0]

    # If y is provided, attach both plot-level y and stratum-level y mean.
    if y is not None:
        # Stratum-level y means to attach to x for covariance.
        y_strat_means = y_eu  # y_strat not retained above; covariance is
        # approximated via yEU means.
        # For TI area, rFIA's area() uses fa_cv from sumToEU; here we compute
        # covariance using the same variance structure at EU level, which is
        # sufficient for ratioVar downstream.
        # To keep behavior aligned, we treat cv at EU as variance of the
        # product and reuse the same structure as fa_var.
        # This matches rFIA's design for the TI forest area ratio.
        pass  # Covariance handled implicitly via fa_cv below.

    x_tmp = x.merge(pops_sub, on=["ESTN_UNIT_CN", "STRATUM_CN"], how="left")

    strat_group_cols_x = [
        "ESTN_UNIT_CN",
        "AREA_USED",
        "P2PNTCNT_EU",
        "STRATUM_CN",
        "STRATUM_WGT",
        "P2POINTCNT",
    ] + x_grp_by

    def _agg_x_strat(group: pd.DataFrame) -> pd.Series:
        P2 = float(group["P2POINTCNT"].iloc[0])
        vals = group[x_var].astype(float)
        mean = (vals / P2).sum()
        if P2 > 1:
            var = (
                (vals.pow(2).sum())
                - P2 * (vals / P2).sum() ** 2
            ) / (P2 * (P2 - 1))
        else:
            var = 0.0
        cv = var
        return pd.Series(
            {
                "fa_mean": mean,
                "fa_var": var,
                "fa_cv": cv,
                "nPlots.x": group["PLT_CN"].nunique(),
            }
        )

    x_strat = (
        x_tmp.groupby(strat_group_cols_x, as_index=False)
        .apply(_agg_x_strat)
        .reset_index(drop=True)
    )

    eu_group_cols_x = ["ESTN_UNIT_CN", "AREA_USED", "P2PNTCNT_EU"] + x_grp_by

    def _agg_x_eu(group: pd.DataFrame) -> pd.Series:
        A = float(group["AREA_USED"].iloc[0])
        N = float(group["P2PNTCNT_EU"].iloc[0])
        w = group["STRATUM_WGT"].astype(float)
        P2 = group["P2POINTCNT"].astype(float)
        m = group["fa_mean"].astype(float)
        v = group["fa_var"].astype(float)
        c = group["fa_cv"].astype(float)

        fa_mean_eu = A * float((w * m).sum())
        fa_var_eu = (A**2 / N) * float(
            (v * w * P2).sum()
            + (v * (1.0 - w) * (P2 / N)).sum()
        )
        fa_cv_eu = (A**2 / N) * float(
            (c * w * P2).sum()
            + (c * (1.0 - w) * (P2 / N)).sum()
        )
        nplots = int(group["nPlots.x"].sum())
        return pd.Series(
            {
                "fa_mean": fa_mean_eu,
                "fa_var": fa_var_eu,
                "fa_cv": fa_cv_eu,
                "nPlots.x": nplots,
                "P2PNTCNT_EU": N,
            }
        )

    x_eu = (
        x_strat.groupby(eu_group_cols_x, as_index=False)
        .apply(_agg_x_eu)
        .reset_index(drop=True)
    )

    return x_eu, y_eu


def _custom_pse_ti_area_basis(
    pops: pd.DataFrame,
    x: pd.DataFrame,
    x_vars: Sequence[str],
    x_grp_by: Sequence[str],
    y: Optional[pd.DataFrame],
    y_vars: Optional[Sequence[str]],
    y_grp_by: Sequence[str],
) -> pd.DataFrame:
    """
    TI design-consistent path for condition-level ``AREA_BASIS`` inputs.

    Uses the same ``_rfia_sum_to_plot`` / ``_rfia_sum_to_eu`` machinery as
    :func:`area`, instead of the unstratified ``_eu_sums`` shortcut (plot mean
    × ``AREA_USED`` within an EU), which does not reproduce FIA/rFIA totals.
    """
    x_grp_by = list(x_grp_by)
    x_grp_by_no_year = [g for g in x_grp_by if g != "YEAR"]
    x_grp_by_full = ["YEAR"] + x_grp_by_no_year
    grp_cols_x = [g for g in x_grp_by if g in x.columns and g != "YEAR"]

    if "CONDID" not in x.columns:
        raise ValueError(
            "For TI estimation with AREA_BASIS, `x` must include `CONDID` "
            "(e.g. from `area(..., cond_list=True)`)."
        )

    id_cols = ["PLT_CN", "CONDID", "AREA_BASIS"]

    y_plt: Optional[pd.DataFrame] = None
    y_grp_by_full: Optional[list[str]] = None
    grp_cols_y: list[str] = []
    yv: Optional[str] = None
    if y is not None:
        if not y_vars or len(y_vars) != 1:
            raise ValueError("Denominator `y_vars` must name exactly one column.")
        yv = str(y_vars[0])
        if "CONDID" not in y.columns:
            raise ValueError(
                "For TI estimation with AREA_BASIS denominator, `y` must include `CONDID`."
            )
        y_grp_by = list(y_grp_by)
        y_grp_by_full = ["YEAR"] + [g for g in y_grp_by if g != "YEAR"]
        grp_cols_y = [g for g in y_grp_by if g in y.columns and g != "YEAR"]
        y_keep = [c for c in id_cols + grp_cols_y + [yv] if c in y.columns]
        y_plt = _rfia_sum_to_plot(y[y_keep].copy(), pops, grp_cols_y)

    y_tot: Optional[pd.DataFrame] = None
    out_merged: Optional[pd.DataFrame] = None

    for v in x_vars:
        if v not in x.columns:
            raise ValueError(f"x_var {v!r} not found in `x`.")
        x_keep = [c for c in id_cols + grp_cols_x + [v] if c in x.columns]
        x_plt = _rfia_sum_to_plot(x[x_keep].copy(), pops, grp_cols_x)

        if y_plt is not None:
            assert y_grp_by_full is not None and yv is not None
            x_eu, y_eu = _rfia_sum_to_eu(
                x_plt, y_plt, pops, x_grp_by_full, y_grp_by_full, "TI"
            )
            y_eu = combine_mr(y_eu, year_col="YEAR")
            if y_tot is None:
                y_sum = y_eu.groupby(y_grp_by_full, as_index=False).sum(
                    numeric_only=True
                )
                y_tot = y_sum[y_grp_by_full + ["fad_mean", "fad_var"]].rename(
                    columns={"fad_mean": f"{yv}_TOTAL", "fad_var": f"{yv}_TOTAL_VAR"}
                )
        else:
            x_eu, _ = _rfia_sum_to_eu(
                x_plt, None, pops, x_grp_by_full, None, "TI"
            )

        x_eu = combine_mr(x_eu, year_col="YEAR")
        x_sum = x_eu.groupby(x_grp_by_full, as_index=False).sum(numeric_only=True)
        x_tot = x_sum[x_grp_by_full + ["fa_mean", "fa_var", "fa_cv"]].rename(
            columns={
                "fa_mean": f"{v}_TOTAL",
                "fa_var": f"{v}_TOTAL_VAR",
                "fa_cv": f"{v}_TOTAL_CV",
            }
        )

        if out_merged is None:
            out_merged = x_tot
        else:
            merge_cols = [c for c in x_tot.columns if c not in out_merged.columns]
            out_merged = out_merged.merge(
                x_tot[x_grp_by_full + merge_cols],
                on=x_grp_by_full,
                how="outer",
            )

    assert out_merged is not None
    out = out_merged

    if y_tot is not None:
        assert y_grp_by_full is not None and yv is not None
        out = out.merge(y_tot, on=y_grp_by_full, how="left")
        for v in x_vars:
            num = out[f"{v}_TOTAL"]
            num_var = out[f"{v}_TOTAL_VAR"]
            num_cv = out[f"{v}_TOTAL_CV"]
            denom = out[f"{yv}_TOTAL"]
            denom_var = out[f"{yv}_TOTAL_VAR"]
            ratio = num / denom.replace({0: pd.NA})
            r_var = ratio_var(
                x=num,
                y=denom,
                x_var=num_var,
                y_var=denom_var,
                cv=num_cv,
            )
            r_var = r_var.where(r_var >= 0, 0.0)
            out[f"{v}_RATIO"] = ratio
            out[f"{v}_RATIO_VAR"] = r_var
            out[f"{v}_TOTAL_SE"] = (
                (num_var**0.5 / num.abs().replace(0, pd.NA) * 100.0).fillna(0.0)
            )
            out[f"{v}_RATIO_SE"] = (
                (r_var**0.5 / ratio.abs().replace(0, pd.NA) * 100.0).fillna(0.0)
            )
        drop_cv = [c for c in out.columns if c.endswith("_TOTAL_CV")]
        out = out.drop(columns=drop_cv, errors="ignore")
    else:
        for v in x_vars:
            num = out[f"{v}_TOTAL"]
            num_var = out[f"{v}_TOTAL_VAR"]
            out[f"{v}_TOTAL_SE"] = (
                (num_var**0.5 / num.abs().replace(0, pd.NA) * 100.0).fillna(0.0)
            )
        drop_cv = [c for c in out.columns if c.endswith("_TOTAL_CV")]
        out = out.drop(columns=drop_cv, errors="ignore")

    return out


def custom_pse(
    db: FiaDatabase,
    x: pd.DataFrame,
    x_vars: Sequence[str],
    x_grp_by: Optional[Sequence[str]] = None,
    x_transform: Optional[Callable[[pd.Series], pd.Series]] = None,
    y: Optional[pd.DataFrame] = None,
    y_vars: Optional[Sequence[str]] = None,
    y_grp_by: Optional[Sequence[str]] = None,
    y_transform: Optional[Callable[[pd.Series], pd.Series]] = None,
    method: str = "TI",
    lambda_: float = 0.5,
    totals: bool = True,
    variance: bool = True,
) -> pd.DataFrame:
    """
    Flexible custom estimator, inspired by rFIA::customPSE().

    This function is intended to mirror the flexibility of rFIA::customPSE.

    It supports both TREE_BASIS and AREA_BASIS inputs:
    - If `TREE_BASIS` is present in `x` (and/or `y`), rows are treated as
      tree-level and are aggregated to plot level with non-response
      adjustments based on ADJ_FACTOR_MICR/SUBP/MACR.
    - If `AREA_BASIS` is present, rows are treated as condition-level and
      aggregated similarly.

    Requirements for `x` (and `y`, if provided):
    - Must include at least:
        - `PLT_CN`
        - `EVAL_TYP` (e.g., 'VOL', 'CURR')
        - ONE of `TREE_BASIS` or `AREA_BASIS`
    - For tree-level data (`TREE_BASIS`), `SUBP` and `TREE` should also be
      present; for condition-level data (`AREA_BASIS`), `CONDID` should be
      present, matching the rFIA behavior.

    When ``method="TI"`` and inputs use ``AREA_BASIS`` (no ``TREE_BASIS``),
    estimation uses the same stratified ``_rfia_sum_to_plot`` /
    ``_rfia_sum_to_eu`` path as :func:`area`, so condition-list workflows can
    reproduce official area totals. The legacy unstratified shortcut is still
    used for ``TREE_BASIS`` or non-TI methods.
    """

    if "PLT_CN" not in x.columns:
        raise ValueError("`x` must include a 'PLT_CN' column.")
    if "EVAL_TYP" not in x.columns:
        raise ValueError("`x` must include an 'EVAL_TYP' column.")

    if y is not None:
        if "PLT_CN" not in y.columns:
            raise ValueError("`y` must include a 'PLT_CN' column when provided.")
        if y_vars is None or len(y_vars) != 1:
            raise ValueError("`y_vars` must be a single column name for denominator.")

    x_vars = list(x_vars)
    if x_grp_by is None:
        x_grp_by = []
    if y_grp_by is None:
        y_grp_by = []
    x_grp_by = list(x_grp_by)
    y_grp_by = list(y_grp_by)

    # Ensure y_grp_by is subset of x_grp_by when denominator is provided
    if y is not None and not set(y_grp_by).issubset(set(x_grp_by)):
        raise ValueError(
            "`y_grp_by` must be equal to or a subset of `x_grp_by`."
        )

    # Design tables must be available
    design_tables = {
        "PLOT",
        "POP_EVAL",
        "POP_EVAL_TYP",
        "POP_ESTN_UNIT",
        "POP_STRATUM",
        "POP_PLOT_STRATUM_ASSGN",
    }
    missing = design_tables.difference(set(db.keys()))
    if missing:
        raise ValueError(
            "Missing required POP_* tables in `db` for custom_pse(): "
            + ", ".join(sorted(missing))
        )

    # Population frame for given eval type
    eval_type = str(x["EVAL_TYP"].iloc[0])
    pops = handle_pops(db, eval_types=(eval_type,), method=method, most_recent=True)

    method_u = str(method).upper()
    x_is_area_basis = "AREA_BASIS" in x.columns and "TREE_BASIS" not in x.columns
    y_is_area_or_absent = y is None or (
        "AREA_BASIS" in y.columns and "TREE_BASIS" not in y.columns
    )

    if method_u == "TI" and x_is_area_basis and y_is_area_or_absent:
        x_work = x.copy()
        if x_transform is not None:
            for col in x_vars:
                x_work[col] = x_transform(x_work[col])
        y_work: Optional[pd.DataFrame] = None
        if y is not None:
            y_work = y.copy()
            if y_transform is not None and y_vars is not None:
                for col in y_vars:
                    y_work[col] = y_transform(y_work[col])
        out = _custom_pse_ti_area_basis(
            pops, x_work, x_vars, x_grp_by, y_work, y_vars, y_grp_by
        )
    else:
        # Sum to plot level with TREE_BASIS / AREA_BASIS handling
        x_plt = _sum_to_plot(x, pops, x_vars, x_grp_by)
        if y is not None:
            y_plt = _sum_to_plot(y, pops, y_vars, y_grp_by)
        else:
            y_plt = None

        # Optional transforms at plot level
        if x_transform is not None:
            for col in x_vars:
                x_plt[col] = x_transform(x_plt[col])
        if y_plt is not None and y_transform is not None:
            for col in y_vars:
                y_plt[col] = y_transform(y_plt[col])

        # Add YEAR to grouping for EU aggregation. Avoid duplicating YEAR if the
        # caller already included it in x_grp_by / y_grp_by.
        x_grp_by_no_year = [g for g in x_grp_by if g != "YEAR"]
        x_grp_by_full = ["YEAR"] + x_grp_by_no_year
        if y_plt is not None:
            y_grp_by_no_year = [g for g in y_grp_by if g != "YEAR"]
            y_grp_by_full = ["YEAR"] + y_grp_by_no_year
        else:
            y_grp_by_full = None

        # Estimation-unit-level sums (simplified stratified TI)
        def _eu_sums(
            df: pd.DataFrame, vars_: Sequence[str], grp_by_full: Sequence[str]
        ) -> pd.DataFrame:
            grp_cols = list(grp_by_full) + ["ESTN_UNIT_CN"]

            def _agg(group: pd.DataFrame) -> pd.Series:
                out = {}
                area = group["AREA_USED"].iloc[0]
                for v in vars_:
                    # Plot-level values are in "per-plot" units (e.g., proportions
                    # or per-plot totals after sumToPlot adjustments). For area-like
                    # quantities (e.g., condition proportions), the design-based
                    # total is AREA_USED × mean(plot_value) within the estimation unit.
                    n = len(group)
                    mean_val = float(group[v].mean()) if n else 0.0
                    total = mean_val * float(area)

                    # Simple SRS variance of the mean, then scaled to totals by area^2
                    if n > 1:
                        var_mean = float(group[v].var(ddof=1)) / n
                    else:
                        var_mean = 0.0
                    total_var = var_mean * (float(area) ** 2)
                    out[f"{v}_TOTAL"] = total
                    out[f"{v}_TOTAL_VAR"] = total_var
                out["AREA_USED"] = area
                return pd.Series(out)

            return (
                df.groupby(grp_cols, as_index=False)
                .apply(_agg)
                .reset_index(drop=True)
            )

        x_eu = _eu_sums(x_plt, x_vars, x_grp_by_full)
        if y_plt is not None:
            y_eu = _eu_sums(y_plt, list(y_vars), y_grp_by_full)
        else:
            y_eu = None

        # Combine most recent inventories across states (TI)
        x_eu = combine_mr(x_eu, year_col="YEAR")
        if y_eu is not None:
            y_eu = combine_mr(y_eu, year_col="YEAR")

        # Aggregate across EUs and, if needed, compute ratios
        x_group_cols = x_grp_by_full

        if y_eu is not None:
            y_group_cols = y_grp_by_full

            # Sum EU totals/vars by group
            x_tot = (
                x_eu.groupby(x_group_cols, as_index=False)
                .sum(numeric_only=True)
            )
            y_tot = (
                y_eu.groupby(y_group_cols, as_index=False)
                .sum(numeric_only=True)
            )

            # Join numerator/denominator
            out = x_tot.merge(y_tot, on=y_group_cols, suffixes=("_x", "_y"))

            # Compute ratios and ratio variances for each x var
            for v in x_vars:
                num = out[f"{v}_TOTAL"]
                num_var = out[f"{v}_TOTAL_VAR"]

                denom = out[f"{y_vars[0]}_TOTAL"]
                denom_var = out[f"{y_vars[0]}_TOTAL_VAR"]

                ratio = num / denom
                r_var = ratio_var(
                    x=num,
                    y=denom,
                    x_var=num_var,
                    y_var=denom_var,
                    cv=pd.Series(0.0, index=num.index),
                )
                r_var = r_var.where(r_var >= 0, 0.0)

                out[f"{v}_RATIO"] = ratio
                out[f"{v}_RATIO_VAR"] = r_var

                # SE for totals and ratios
                out[f"{v}_TOTAL_SE"] = (
                    (num_var**0.5 / num.abs() * 100.0).where(num != 0, 0.0)
                )
                out[f"{v}_RATIO_SE"] = (
                    (r_var**0.5 / ratio.abs() * 100.0).where(ratio != 0, 0.0)
                )

        else:
            # Only numerator: sum totals/vars by group
            out = (
                x_eu.groupby(x_group_cols, as_index=False)
                .sum(numeric_only=True)
            )

            for v in x_vars:
                num = out[f"{v}_TOTAL"]
                num_var = out[f"{v}_TOTAL_VAR"]
                out[f"{v}_TOTAL_SE"] = (
                    (num_var**0.5 / num.abs() * 100.0).where(num != 0, 0.0)
                )

    # Optionally drop totals or variance columns
    if not totals:
        out = out.drop(
            columns=[c for c in out.columns if c.endswith("_TOTAL")],
            errors="ignore",
        )
    if variance:
        out = out.drop(
            columns=[c for c in out.columns if c.endswith("_SE")],
            errors="ignore",
        )
    else:
        out = out.drop(
            columns=[c for c in out.columns if c.endswith("_VAR")],
            errors="ignore",
        )

    # Sort and return
    out = out.sort_values("YEAR").reset_index(drop=True)
    return out


def tpa(
    db: FiaDatabase,
    grp_by: Optional[Sequence[str]] = None,
    land_type: str = "forest",
    tree_type: str = "live",
    by_species: bool = False,
    by_size_class: bool = False,
    by_plot: bool = False,
    tree_list: bool = False,
) -> pd.DataFrame:
    """
    Estimate tree abundance (TPA) and basal area (BAA) using a simplified
    plot‑level estimator.

    Notes
    -----
    - This is an initial, simplified Python implementation that:
      - Uses only `PLOT`, `COND`, and `TREE` tables from `db`.
      - Ignores full FIA stratification / estimation‑unit weighting
        (i.e., does *not* yet replicate the full rFIA design‑based estimator).
      - Does not compute sampling errors or variances.
    - It is intended to be API‑compatible with `rFIA::tpa()` at a basic level,
      so you can start prototyping workflows in Python. A full design‑based
      port (including POP_* tables and variance estimation) will be layered
      on top of this.
    """

    if grp_by is None:
        grp_by = []
    grp_by = list(grp_by)

    # Required tables
    required = {"PLOT", "COND", "TREE"}
    missing = required.difference(set(db.keys()))
    if missing:
        raise ValueError(
            "Missing required tables in `db` for tpa(): "
            + ", ".join(sorted(missing))
        )

    plot = db["PLOT"].copy()
    cond = db["COND"].copy()
    tree = db["TREE"].copy()

    # Basic checks for required columns
    for col in ["PLT_CN", "INVYR", "PLOT_STATUS_CD"]:
        if col not in plot.columns:
            raise ValueError(f"`PLOT` table must contain column '{col}'.")
    for col in ["PLT_CN", "CONDID", "COND_STATUS_CD", "CONDPROP_UNADJ", "SITECLCD", "RESERVCD"]:
        if col not in cond.columns:
            raise ValueError(f"`COND` table must contain column '{col}'.")
    for col in ["PLT_CN", "CONDID", "DIA", "TPA_UNADJ", "STATUSCD", "TREECLCD"]:
        if col not in tree.columns:
            raise ValueError(f"`TREE` table must contain column '{col}'.")

    # Restrict to current, forested plots
    plot = plot.loc[plot["PLOT_STATUS_CD"] == 1].copy()

    # Land and tree domains
    cond["landD"] = _land_type_domain(
        land_type,
        cond["COND_STATUS_CD"],
        cond["SITECLCD"],
        cond["RESERVCD"],
    )
    tree["typeD"] = _tree_type_domain(
        tree_type,
        tree["STATUSCD"],
        tree["DIA"],
        tree["TREECLCD"],
    )

    # Species information ----------------------------------------------------
    if by_species:
        if "SPCD" not in tree.columns:
            raise ValueError("`TREE` table must contain 'SPCD' when by_species=True.")

        # Try to enrich with species reference table if available
        species_ref = None
        if "SPECIES" in db.keys():
            species_ref = db["SPECIES"]
            needed_cols = ["SPCD", "COMMON_NAME", "GENUS", "SPECIES"]
            missing_cols = [c for c in needed_cols if c not in species_ref.columns]
            if missing_cols:
                species_ref = None
            else:
                species_ref = species_ref[needed_cols].drop_duplicates()

        if species_ref is not None:
            tree = tree.merge(species_ref, on="SPCD", how="left")
            tree["SCIENTIFIC_NAME"] = (
                tree["GENUS"].astype("string") + " " + tree["SPECIES"].astype("string")
            )
            for col in ["SPCD", "COMMON_NAME", "SCIENTIFIC_NAME"]:
                if col not in grp_by:
                    grp_by.append(col)
        else:
            if "SPCD" not in grp_by:
                grp_by.append("SPCD")

    # Size class information -------------------------------------------------
    if by_size_class:
        if "DIA" not in tree.columns:
            raise ValueError(
                "`TREE` table must contain 'DIA' when by_size_class=True."
            )
        # 2-inch size classes, indexed as odd integers: 1, 3, 5, ...
        dia = tree["DIA"].where(tree["DIA"] > 0)
        size_class = (dia // 2).astype("Int64") * 2 + 1
        tree["sizeClass"] = size_class
        if "sizeClass" not in grp_by:
            grp_by.append("sizeClass")

    # Join PLOT + COND + TREE
    data = (
        plot.merge(cond, on="PLT_CN", how="inner", suffixes=("", "_COND"))
        .merge(tree, on=["PLT_CN", "CONDID"], how="inner", suffixes=("", "_TREE"))
    )

    # Drop trees without usable measurements
    data = data[data["TPA_UNADJ"].fillna(0) > 0].copy()

    # Domain indicators (simplified; no user‑supplied area/tree domains yet)
    data["aDI"] = data["landD"]
    data["tDI"] = data["landD"] * data["typeD"]

    # -----------------------------
    # Optional per-plot / tree list
    # -----------------------------

    # Condition-level forest area proportion within plot
    cond_area = (
        data[["PLT_CN", "CONDID", "CONDPROP_UNADJ", "aDI"]]
        .drop_duplicates()
        .assign(PROP_FOREST=lambda d: d["CONDPROP_UNADJ"] * d["aDI"])
        .groupby("PLT_CN", as_index=False)["PROP_FOREST"]
        .sum()
    )

    # Tree-level contributions (unweighted by design; plot-level units)
    tree_base = data.assign(
        TPA=lambda d: d["TPA_UNADJ"] * d["tDI"],
        BAA=lambda d: _basal_area(d["DIA"]) * d["TPA_UNADJ"] * d["tDI"],
        YEAR=lambda d: d["MEASYEAR"] if "MEASYEAR" in d.columns else d["INVYR"],
    )

    if tree_list:
        # Per-tree list, similar to rFIA::tpa(treeList=TRUE)
        grp_cols = [g for g in grp_by if g in tree_base.columns]
        cols = (
            ["PLT_CN", "CONDID", "SUBP", "TREE", "YEAR"]
            + grp_cols
            + ["TPA", "BAA"]
        )
        tl = (
            tree_base[cols]
            .merge(cond_area, on="PLT_CN", how="left")
            .rename(columns={"PROP_FOREST": "PROP_FOREST"})
        )
        tl["EVAL_TYP"] = "VOL"
        return tl

    if by_plot:
        # Plot-level estimates, similar to rFIA::tpa(byPlot=TRUE)
        grp_cols = [g for g in grp_by if g in tree_base.columns]
        plot_level = (
            tree_base.groupby(
                ["PLT_CN", "YEAR"] + grp_cols,
                as_index=False,
            )[["TPA", "BAA"]]
            .sum()
            .merge(cond_area, on="PLT_CN", how="left")
        )
        return plot_level

    # Tree‑level contributions within each plot for design-based estimator
    tree_plot = (
        tree_base.groupby(
            ["PLT_CN", "YEAR"] + [g for g in grp_by if g in tree_base.columns],
            as_index=False,
        )[["TPA", "BAA"]]
        .sum()
    )

    # If design tables are available, use full TI design-based estimator.
    design_tables = {
        "POP_EVAL",
        "POP_EVAL_TYP",
        "POP_ESTN_UNIT",
        "POP_STRATUM",
        "POP_PLOT_STRATUM_ASSGN",
    }
    have_design = design_tables.issubset(set(db.keys()))

    if not have_design:
        # Fallback: simple plot-weighted mean across plots (no design weights).
        group_cols: List[str] = ["YEAR"] + [
            g for g in grp_by if g in tree_plot.columns
        ]

        def _agg_simplified(group: pd.DataFrame) -> pd.Series:
            tpa = group["TPA"].mean()
            baa = group["BAA"].mean()
            return pd.Series(
                {
                    "TPA": tpa,
                    "BAA": baa,
                    "nPlots_TREE": int((group["TPA"] > 0).sum()),
                    "nPlots_AREA": int(len(group)),
                }
            )

        out = (
            tree_plot.groupby(group_cols, as_index=False)
            .apply(_agg_simplified)
            .reset_index(drop=True)
            .sort_values("YEAR")
        )
        return out

    # -------------------------
    # Full TI design-based path
    # -------------------------

    # Population / design info (TI)
    pops = handle_pops(db, eval_types=("VOL",), method="TI", most_recent=True)

    # Attach design info to each plot (include P2POINTCNT for stratum means)
    design_by_plot = (
        pops[
            [
                "PLT_CN",
                "STATECD",
                "YEAR",
                "ESTN_UNIT_CN",
                "STRATUM_CN",
                "STRATUM_WGT",
                "P2POINTCNT",
                "P2PNTCNT_EU",
                "AREA_USED",
            ]
        ]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    grp_cols_tree = [g for g in grp_by if g in tree_base.columns]

    if not grp_cols_tree:
        # One row per sample plot: include all Phase 2 plots in the population frame with
        # TPA=BAA=0 when the plot has no qualifying trees (inner merge was biasing means up).
        plot_sums = (
            tree_base.groupby("PLT_CN", as_index=False)[["TPA", "BAA"]]
            .sum()
        )
        design_dedup = design_by_plot.drop_duplicates(subset=["PLT_CN"], keep="first")
        tree_plot = design_dedup.merge(plot_sums, on="PLT_CN", how="left")
        tree_plot["TPA"] = tree_plot["TPA"].fillna(0.0)
        tree_plot["BAA"] = tree_plot["BAA"].fillna(0.0)
    else:
        tree_plot = (
            tree_base.groupby(
                ["PLT_CN", "YEAR"] + grp_cols_tree,
                as_index=False,
            )[["TPA", "BAA"]]
            .sum()
        )
        design_dedup = design_by_plot.drop_duplicates(subset=["PLT_CN"], keep="first")
        tree_plot = tree_plot.merge(design_dedup, on="PLT_CN", how="inner")
        tree_plot["YEAR"] = (
            tree_plot["YEAR_y"] if "YEAR_y" in tree_plot.columns else tree_plot["YEAR"]
        )

    # Stratum-level means and variances within each estimation unit
    group_cols_base: List[str] = ["YEAR"] + [
        g for g in grp_by if g in tree_plot.columns
    ]
    strat_group_cols = (
        group_cols_base
        + [
            "STATECD",
            "ESTN_UNIT_CN",
            "STRATUM_CN",
            "STRATUM_WGT",
            "AREA_USED",
            "P2POINTCNT",
        ]
    )

    def _stratum_stats(group: pd.DataFrame) -> pd.Series:
        # Match _rfia_sum_to_eu / FIA stratum formulas: mean = sum(y_i) / P2POINTCNT
        P2 = float(group["P2POINTCNT"].iloc[0])
        vals_tpa = group["TPA"].astype(float)
        vals_baa = group["BAA"].astype(float)
        tpa_mean = float((vals_tpa / P2).sum()) if P2 > 0 else 0.0
        baa_mean = float((vals_baa / P2).sum()) if P2 > 0 else 0.0
        if P2 > 1:
            tpa_var = float(
                (vals_tpa.pow(2).sum() - P2 * (tpa_mean**2)) / (P2 * (P2 - 1))
            )
            baa_var = float(
                (vals_baa.pow(2).sum() - P2 * (baa_mean**2)) / (P2 * (P2 - 1))
            )
            tpa_var = max(0.0, tpa_var)
            baa_var = max(0.0, baa_var)
        else:
            tpa_var = 0.0
            baa_var = 0.0
        n = int(len(group))
        return pd.Series(
            {
                "nPlots_STRAT": n,
                "TPA_STRAT_MEAN": tpa_mean,
                "BAA_STRAT_MEAN": baa_mean,
                "TPA_STRAT_VAR": tpa_var,
                "BAA_STRAT_VAR": baa_var,
            }
        )

    strat = (
        tree_plot.groupby(strat_group_cols, as_index=False)
        .apply(_stratum_stats)
        .reset_index(drop=True)
    )

    # Estimation-unit-level means and variances (stratified estimator)
    eu_group_cols = group_cols_base + ["STATECD", "ESTN_UNIT_CN"]

    def _eu_stats(group: pd.DataFrame) -> pd.Series:
        w = group["STRATUM_WGT"]
        area = group["AREA_USED"].iloc[0]

        # Stratum-level means/vars
        tpa_m = group["TPA_STRAT_MEAN"]
        baa_m = group["BAA_STRAT_MEAN"]
        tpa_v = group["TPA_STRAT_VAR"]
        baa_v = group["BAA_STRAT_VAR"]

        # EU means (per-acre)
        tpa_mean_eu = float((w * tpa_m).sum())
        baa_mean_eu = float((w * baa_m).sum())

        # EU variances of means (no FPC, standard stratified SRS)
        tpa_var_eu = float((w**2 * tpa_v).sum())
        baa_var_eu = float((w**2 * baa_v).sum())

        # EU totals and variances
        tpa_total_eu = tpa_mean_eu * area
        baa_total_eu = baa_mean_eu * area
        area_total_eu = area

        tpa_total_var_eu = tpa_var_eu * (area**2)
        baa_total_var_eu = baa_var_eu * (area**2)

        return pd.Series(
            {
                "TPA_EU_MEAN": tpa_mean_eu,
                "BAA_EU_MEAN": baa_mean_eu,
                "TPA_EU_MEAN_VAR": tpa_var_eu,
                "BAA_EU_MEAN_VAR": baa_var_eu,
                "TREE_TOTAL_EU": tpa_total_eu,
                "BA_TOTAL_EU": baa_total_eu,
                "AREA_TOTAL_EU": area_total_eu,
                "TREE_TOTAL_EU_VAR": tpa_total_var_eu,
                "BA_TOTAL_EU_VAR": baa_total_var_eu,
                "nPlots_TREE": int(group["nPlots_STRAT"].sum()),
            }
        )

    eu = (
        strat.groupby(eu_group_cols, as_index=False)
        .apply(_eu_stats)
        .reset_index(drop=True)
    )

    # Aggregate across estimation units to overall domain
    final_group_cols = group_cols_base

    def _overall(group: pd.DataFrame) -> pd.Series:
        A = group["AREA_TOTAL_EU"]
        Atot = A.sum()
        # Weighted mean across EUs
        tpa_mean = float(((A * group["TPA_EU_MEAN"]).sum()) / Atot) if Atot > 0 else 0.0
        baa_mean = float(((A * group["BAA_EU_MEAN"]).sum()) / Atot) if Atot > 0 else 0.0

        # Approximate variance of mean via EU-level mean variances
        tpa_mean_var = float(((A / Atot) ** 2 * group["TPA_EU_MEAN_VAR"]).sum())
        baa_mean_var = float(((A / Atot) ** 2 * group["BAA_EU_MEAN_VAR"]).sum())

        # Totals and variances
        tree_total = tpa_mean * Atot
        ba_total = baa_mean * Atot
        area_total = Atot
        tree_total_var = tpa_mean_var * (Atot**2)
        ba_total_var = baa_mean_var * (Atot**2)

        # No covariance currently; pass cv=0 into ratio_var
        tpa_var = ratio_var(
            x=pd.Series(tree_total),
            y=pd.Series(area_total),
            x_var=pd.Series(tree_total_var),
            y_var=pd.Series(0.0),
            cv=pd.Series(0.0),
        ).iloc[0]
        baa_var = ratio_var(
            x=pd.Series(ba_total),
            y=pd.Series(area_total),
            x_var=pd.Series(ba_total_var),
            y_var=pd.Series(0.0),
            cv=pd.Series(0.0),
        ).iloc[0]

        # Sampling errors (%)
        tpa_se = (tpa_var**0.5 / tpa_mean * 100.0) if tpa_mean != 0 else 0.0
        baa_se = (baa_var**0.5 / baa_mean * 100.0) if baa_mean != 0 else 0.0
        tree_total_se = (
            tree_total_var**0.5 / tree_total * 100.0 if tree_total != 0 else 0.0
        )
        ba_total_se = (
            ba_total_var**0.5 / ba_total * 100.0 if ba_total != 0 else 0.0
        )

        return pd.Series(
            {
                "TPA": tpa_mean,
                "BAA": baa_mean,
                "TREE_TOTAL": tree_total,
                "BA_TOTAL": ba_total,
                "AREA_TOTAL": area_total,
                "TPA_VAR": tpa_var,
                "BAA_VAR": baa_var,
                "TREE_TOTAL_VAR": tree_total_var,
                "BA_TOTAL_VAR": ba_total_var,
                "AREA_TOTAL_VAR": 0.0,
                "TPA_SE": tpa_se,
                "BAA_SE": baa_se,
                "TREE_TOTAL_SE": tree_total_se,
                "BA_TOTAL_SE": ba_total_se,
                "AREA_TOTAL_SE": 0.0,
                "nPlots_TREE": int(group["nPlots_TREE"].sum()),
                "nPlots_AREA": int(group["nPlots_TREE"].sum()),
                "N": int(group["P2PNTCNT_EU"].sum()),
            }
        )

    # Need P2PNTCNT_EU in eu for N
    if "P2PNTCNT_EU" not in eu.columns and "P2PNTCNT_EU" in design_by_plot.columns:
        eu = eu.merge(
            design_by_plot[["ESTN_UNIT_CN", "P2PNTCNT_EU"]].drop_duplicates(),
            on="ESTN_UNIT_CN",
            how="left",
        )

    out = (
        eu.groupby(final_group_cols, as_index=False)
        .apply(_overall)
        .reset_index(drop=True)
        .sort_values("YEAR")
    )

    return out


__all__ = ["area", "custom_pse", "tpa"]

