from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional, Sequence

import pandas as pd

# Bump together with ``scripts/test_mog_area_mt.py`` ``MOG_MT_SCRIPT_BUILD_STAMP`` when
# changing condition filters or forest-type logic; ``--verify-mog-code`` asserts they match.
MOG_ESTIMATORS_REVISION = "radical-verify-2026-03-29-t5"

from fia.data_io import FiaDatabase
from fia.estimators import _land_type_domain, custom_pse

from .auxiliary import (
    load_usfs_master_tree_species_list,
    montana_plot_in_cont_divide_east,
    resolve_mog_auxiliary_dir,
    species_lookup_from_master_list,
)
from .crosswalk import classify_region, statecd_to_abbrev
from .engine import ConditionContext, MOGEngine
from .northern.diagnostics import MOG_CONDITION_NR_COLUMNS, northern_veg_diagnostic_row
from .southwest.diagnostics import MOG_CONDITION_SW_COLUMNS
from .southwest.evaluate import southwest_og_diagnostic_row
from .paths import (
    master_tree_species_csv,
    oregon_counties_shp,
    pacific_northwest_nwfp_boundary_shp,
    pacific_northwest_paz_raster,
)
from .pnw.spatial import (
    oregon_plot_in_counties_layer,
    plot_inside_nwfp,
    sample_r5_paz_value,
)


def _ref_species_lookup(ref: pd.DataFrame | None) -> pd.DataFrame | None:
    if ref is None or ref.empty or "SPCD" not in ref.columns:
        return None
    for col in ("PLANTS_SYMBOL", "PLANTS_SYMBOL_NRCS", "PLANTS_CODE"):
        if col in ref.columns:
            out = ref[["SPCD", col]].drop_duplicates(subset=["SPCD"])
            return out.rename(columns={col: "PLANTS_Code"})
    return None


def _to_int_or_na(val: object) -> object:
    """Coerce FIADB codes to int for CSV output; missing/invalid → pandas NA."""
    if val is None:
        return pd.NA
    num = pd.to_numeric(val, errors="coerce")
    if pd.isna(num):
        return pd.NA
    try:
        return int(num)
    except (TypeError, ValueError, OverflowError):
        return pd.NA


def _mog_forest_type_from_cond_row(c: pd.Series, *, region: str | None = None) -> int:
    """
    Forest-type integer for ``ConditionContext.forest_type`` (R uses ``FLDTYPCD``).

    **Northern Region:** If ``FLDTYPCD`` does not map to a Green et al. OG label but
    ``FORTYPCD`` does (e.g. coarse group vs 221 ponderosa), use the code that maps so
    habitat OG vectors are not empty while letters still crosswalk from composition.

    Elsewhere, keep mapMOG order: ``FLDTYPCD`` if present, else ``FORTYPCD``.
    """

    from .northern.core import northern_og_forest_type

    fld = pd.to_numeric(c.get("FLDTYPCD"), errors="coerce")
    fort = pd.to_numeric(c.get("FORTYPCD"), errors="coerce")

    if _region_is_northern(region):
        for cand in (fld, fort):
            if pd.notna(cand):
                v = int(cand)
                if northern_og_forest_type(v) is not None:
                    return v

    if pd.notna(fld):
        return int(fld)
    if pd.notna(fort):
        return int(fort)
    raise RuntimeError("MOG forest type: FLDTYPCD and FORTYPCD are both missing")


@dataclass(frozen=True)
class MOGAreaResult:
    """
    Convenience wrapper for old-growth area estimation outputs.

    - ``df``: old-growth acres from ``custom_pse`` (``OG_PROP`` totals).
    - ``forest_area_df``: total forested acres on the same domain, from a
      parallel ``custom_pse`` run with ``FOREST_PROP = PROP_FOREST`` (no OG filter).
      Use ``FOREST_PROP_TOTAL`` for acres by group (same structure as ``df``).
    - ``cond_mog``: condition-level MOG scores and ``PROP_FOREST``.
    """

    df: pd.DataFrame
    cond_mog: pd.DataFrame
    forest_area_df: pd.DataFrame = field(default_factory=pd.DataFrame)


def _region_is_southwest(region: object) -> bool:
    """True for Region 3 MOG (``classify_region`` → ``southwest``), robust to NA-like values."""

    if region is None:
        return False
    try:
        if pd.isna(region):
            return False
    except (TypeError, ValueError):
        pass
    return str(region).strip().lower() == "southwest"


def _region_is_northern(region: object) -> bool:
    """True for Region 1 MOG (``classify_region`` → ``northern``)."""

    if region is None:
        return False
    try:
        if pd.isna(region):
            return False
    except (TypeError, ValueError):
        pass
    return str(region).strip().lower() == "northern"


def _safe_max_age(stdage: object, fldage: object) -> float:
    ages = []
    for a in (stdage, fldage):
        if a is None or pd.isna(a):
            continue
        try:
            ages.append(float(a))
        except (TypeError, ValueError):
            continue
    return max(ages) if ages else 0.0


def _condition_area_acres(prop_basis: object, condprop_unadj: object) -> float:
    basis = None if prop_basis is None or pd.isna(prop_basis) else str(prop_basis).upper()
    sub_plot_area = 0.2460 if basis == "MACR" else 0.0417
    cpu = 0.0 if condprop_unadj is None or pd.isna(condprop_unadj) else float(condprop_unadj)
    return float(sub_plot_area * (4.0 * cpu))


def _pnw_site_class_max(siteclcd: object, siteclcdest: object) -> float | None:
    """R ``max(SITECLCD, SITECLCDEST, na.rm=TRUE)`` for PNW Table 14."""

    vals: list[float] = []
    for x in (siteclcd, siteclcdest):
        if x is None or pd.isna(x):
            continue
        try:
            vals.append(float(x))
        except (TypeError, ValueError):
            continue
    return max(vals) if vals else None


def mog_condition_scores(
    db: FiaDatabase,
    *,
    states: Optional[Sequence[str]] = None,
    eval_typ: str = "CURR",
    land_type: str = "forest",
    mog_auxiliary_dir: Path | str | None = None,
    use_montana_divide_shapefile: bool = True,
) -> pd.DataFrame:
    """
    Compute MOG scores per (PLT_CN, CONDID) and return a condition-level table.

    ``mog_auxiliary_dir`` defaults to ``fia_py/mog_auxillary`` (R ``source.path``).
    When the National Master Tree Species List CSV is present there, species
    lookups follow ``FUNCTION_mapMOG.R`` (``USFStrees`` / ``FIA.Code`` → ``PLANTS.Code``).
    For Montana, ``geopandas`` + ``utility_MTcontDivide`` assign the eastern vs western
    northern subregion when ``use_montana_divide_shapefile`` is True.
    For Pacific Northwest plots, ``rasterio`` and ``geopandas`` sample ``utility_R5_PAZ.tif``,
    test NWFP boundary overlap, and (for white fir) Oregon county overlay when those
    files exist under the auxiliary directory.
    For Pacific Southwest (California / Region 5), the same NWFP shapefile is used to
    test plot overlap for **white fir** Table 12 logic (R loads PAZ+NWFP for CA/NV; MOG uses
    NWFP for PSW white fir only, not the PAZ raster).

    For **southwest** Table 9 SDI rules, **COND** should include **``SDI_RMRS``** (FIADB field)
    so **SW_REL_SDI** can be computed as ``100 × SDI_OFE / SDI_RMRS``; it is merged from
    **COND** when that column exists.

    Output columns:
    - PLT_CN, CONDID
    - MOG_SCORE (max of mog vector; can be 1 from Table 19 maturity weights alone)
    - OG_FLAG (binary **old growth** for area expansion: regional OG tables only—
      Table 9 / 11 / 12 / 13–14 style rules—not Table 19 maturity, even if
      ``MOG_SCORE`` is 1 from maturity weights)
    - For **southwest** rows only, Table 9 / Table 19 diagnostics (``SW_*`` columns;
      see :func:`fia_mog.southwest.diagnostics.summarize_southwest_og_by_eru`); ``pd.NA``
      elsewhere
    - For **northern** rows only, vegetation / habitat crosswalk diagnostics (``NR_*``;
      see :func:`fia_mog.northern.diagnostics.summarize_mt_northern_veg_gaps` for Montana
      summaries); ``pd.NA`` elsewhere
    - **STATECD** on each row (for filtering to a single state in post-processing)
    - AREA_BASIS (from COND.PROP_BASIS)
    - CONDPROP_UNADJ
    - PROP_FOREST (CONDPROP_UNADJ × land domain; same definition as ``area`` / ``tpa``)
    - EVAL_TYP (constant, for custom_pse)
    """

    cond = db["COND"].copy()
    tree = db["TREE"].copy()

    if states is not None:
        # Use STATECD mapping if present; otherwise rely on COND.STATECD directly.
        # States are expected as abbreviations like ["AZ"].
        # Minimal mapping for CONUS use; extend if needed.
        st_to_fips = {
            "AZ": 4,
            "CA": 6,
            "MT": 30,
            "ID": 16,
            "WA": 53,
            "OR": 41,
            "WY": 56,
            "SD": 46,
            "ND": 38,
        }
        wanted = {st_to_fips[s.upper()] for s in states if s.upper() in st_to_fips}
        if "STATECD" in cond.columns and wanted:
            cond = cond[cond["STATECD"].isin(wanted)].copy()

    # Match mapMOG intent: accessible forest conditions with a type code (FLDTYPCD and/or FORTYPCD).
    fld = pd.to_numeric(cond["FLDTYPCD"], errors="coerce") if "FLDTYPCD" in cond.columns else None
    fort = (
        pd.to_numeric(cond["FORTYPCD"], errors="coerce") if "FORTYPCD" in cond.columns else None
    )
    if fld is None:
        has_ft = fort.notna() if fort is not None else pd.Series(False, index=cond.index)
    elif fort is None:
        has_ft = fld.notna()
    else:
        has_ft = fld.notna() | fort.notna()
    cond = cond[(cond["COND_STATUS_CD"] == 1) & has_ft].copy()

    # Forest / timber land domain (matches fia.estimators.area and tpa).
    for col in ("SITECLCD", "RESERVCD"):
        if col not in cond.columns:
            cond[col] = pd.NA
    cond["landD"] = _land_type_domain(
        land_type,
        cond["COND_STATUS_CD"],
        cond["SITECLCD"],
        cond["RESERVCD"],
    )
    cond["PROP_FOREST"] = pd.to_numeric(cond["CONDPROP_UNADJ"], errors="coerce").fillna(0.0) * cond[
        "landD"
    ]

    # Join condition attributes onto trees so we can group once.
    join_cols = [
        "PLT_CN",
        "CONDID",
        "COND_STATUS_CD",
        "STATECD",
        "ADFORCD",
        "FLDTYPCD",
        "FORTYPCD",
        "PHYSCLCD",
        "SITECLCD",
        "PROP_BASIS",
        "CONDPROP_UNADJ",
        "STDAGE",
        "FLDAGE",
        "HABTYPCD1",
        "SITECLCDEST",
        "SDI_RMRS",
        "landD",
        "PROP_FOREST",
    ]
    join_cols = [c for c in join_cols if c in cond.columns]

    data = tree.merge(cond[join_cols], on=["PLT_CN", "CONDID"], how="inner")

    engine = MOGEngine()
    out = []
    _mog_diag_na = {
        **{k: pd.NA for k in MOG_CONDITION_SW_COLUMNS},
        **{k: pd.NA for k in MOG_CONDITION_NR_COLUMNS},
    }
    aux_dir = resolve_mog_auxiliary_dir(mog_auxiliary_dir)
    master_csv = master_tree_species_csv(aux_dir)
    if master_csv.is_file():
        try:
            ref_lookup = species_lookup_from_master_list(load_usfs_master_tree_species_list(aux_dir))
        except (OSError, ValueError, KeyError):
            ref_lookup = _ref_species_lookup(db.tables.get("REF_SPECIES"))
    else:
        ref_lookup = _ref_species_lookup(db.tables.get("REF_SPECIES"))
    p2veg = db.tables.get("P2VEG_SUBPLOT_SPP")

    plot_ll_by_cn: pd.DataFrame | None = None
    plo = db.tables.get("PLOT")
    if (
        plo is not None
        and not plo.empty
        and {"CN", "LON", "LAT"}.issubset(plo.columns)
    ):
        plot_ll_by_cn = plo[["CN", "LON", "LAT"]].drop_duplicates(subset=["CN"]).set_index("CN")

    plot_elev_by_cn: dict[object, float] = {}
    if (
        plo is not None
        and not plo.empty
        and {"CN", "ELEV"}.issubset(plo.columns)
    ):
        pe = plo[["CN", "ELEV"]].drop_duplicates(subset=["CN"])
        for _, r in pe.iterrows():
            try:
                if pd.notna(r["ELEV"]):
                    ev = float(r["ELEV"])
                    if math.isfinite(ev):
                        plot_elev_by_cn[r["CN"]] = ev
            except (TypeError, ValueError):
                pass

    dwm_by_pc: dict[tuple[object, object], pd.DataFrame] = {}
    dwm_tbl = db.tables.get("DWM_COARSE_WOODY_DEBRIS")
    if dwm_tbl is not None and not dwm_tbl.empty and {"PLT_CN", "CONDID"}.issubset(dwm_tbl.columns):
        for key, sub in dwm_tbl.groupby(["PLT_CN", "CONDID"]):
            dwm_by_pc[key] = sub

    ecosub_by_plt: dict[object, object] = {}
    plotgeom = db.tables.get("PLOTGEOM")
    if (
        plotgeom is not None
        and not plotgeom.empty
        and "CN" in plotgeom.columns
        and "ECOSUBCD" in plotgeom.columns
    ):
        for _, row in plotgeom.drop_duplicates(subset=["CN"]).iterrows():
            ecosub_by_plt[row["CN"]] = row["ECOSUBCD"]

    paz_tif = pacific_northwest_paz_raster(aux_dir)
    nwfp_shp = pacific_northwest_nwfp_boundary_shp(aux_dir)
    or_shp = oregon_counties_shp(aux_dir)

    for (plt_cn, condid), g in data.groupby(["PLT_CN", "CONDID"]):
        c = g.iloc[0]

        st_abbrev = statecd_to_abbrev(c.get("STATECD"))
        region = classify_region(c.get("ADFORCD"), state_abbrev=st_abbrev)
        stand_age = _safe_max_age(c.get("STDAGE"), c.get("FLDAGE"))
        condition_area = _condition_area_acres(c.get("PROP_BASIS"), c.get("CONDPROP_UNADJ"))

        plot_lon_f: float | None = None
        plot_lat_f: float | None = None
        if plot_ll_by_cn is not None:
            try:
                _row_ll = plot_ll_by_cn.loc[plt_cn]
                plot_lon_f = float(_row_ll["LON"])
                plot_lat_f = float(_row_ll["LAT"])
            except (KeyError, TypeError, ValueError):
                pass

        mt_div: bool | None = None
        if st_abbrev == "MT" and use_montana_divide_shapefile and plot_lon_f is not None and plot_lat_f is not None:
            try:
                mt_div = montana_plot_in_cont_divide_east(plot_lon_f, plot_lat_f, aux_dir)
            except (KeyError, TypeError, ValueError):
                mt_div = None

        pnw_paz: float | int | None = None
        pnw_inside: bool | None = None
        pnw_or_counties: bool | None = None
        pnw_woody = dwm_by_pc.get((plt_cn, condid))
        raw_eco = ecosub_by_plt.get(plt_cn)
        eco_str = None if raw_eco is None or pd.isna(raw_eco) else str(raw_eco).strip()
        pnw_site_max = _pnw_site_class_max(c.get("SITECLCD"), c.get("SITECLCDEST"))

        if plot_lon_f is not None and plot_lat_f is not None:
            reg = region or ""
            if reg == "pacific northwest":
                pv = sample_r5_paz_value(plot_lon_f, plot_lat_f, paz_tif)
                if pv is not None:
                    pnw_paz = int(round(float(pv)))
                pnw_or_counties = oregon_plot_in_counties_layer(plot_lon_f, plot_lat_f, or_shp)
            if reg in ("pacific northwest", "pacific southwest"):
                pnw_inside = plot_inside_nwfp(plot_lon_f, plot_lat_f, nwfp_shp)

        sdi_rmrs: float | None = None
        if "SDI_RMRS" in g.columns:
            raw_rmrs = c.get("SDI_RMRS")
            if raw_rmrs is not None and not pd.isna(raw_rmrs):
                try:
                    sdi_rmrs_f = float(raw_rmrs)
                    sdi_rmrs = sdi_rmrs_f if math.isfinite(sdi_rmrs_f) and sdi_rmrs_f > 0 else None
                except (TypeError, ValueError):
                    sdi_rmrs = None

        plot_elev_ft = plot_elev_by_cn.get(plt_cn)

        ctx = ConditionContext(
            region=region or "",
            forest_type=_mog_forest_type_from_cond_row(c, region=region),
            condition_area_acres=condition_area,
            stand_age=stand_age,
            trees=g,
            plot_statecd=None if pd.isna(c.get("STATECD")) else int(c.get("STATECD")),
            condition_fortypcd=None if pd.isna(c.get("FORTYPCD")) else int(c.get("FORTYPCD")),
            condition_physclcd=None if pd.isna(c.get("PHYSCLCD")) else float(c.get("PHYSCLCD")),
            condition_siteclcd=None if pd.isna(c.get("SITECLCD")) else float(c.get("SITECLCD")),
            condition_siteclcdest=None if pd.isna(c.get("SITECLCDEST")) else float(c.get("SITECLCDEST")),
            condition_adforcd=None if pd.isna(c.get("ADFORCD")) else int(c.get("ADFORCD")),
            condition_habtypcd1=None if pd.isna(c.get("HABTYPCD1")) else float(c.get("HABTYPCD1")),
            condition_sdi_rmrs=sdi_rmrs,
            ecosubcd=eco_str,
            northern_species_lookup=ref_lookup,
            northern_veg_subplot=p2veg,
            northern_mt_east_of_divide=mt_div,
            pnw_paz_raster_value=pnw_paz,
            pnw_inside_nwfp=pnw_inside,
            pnw_woody_debris=pnw_woody,
            pnw_site_class_max=pnw_site_max,
            pnw_plot_in_or_counties_layer=pnw_or_counties,
            plot_elev_ft=plot_elev_ft,
            plot_lon=plot_lon_f,
        )

        mog_vec = engine.mog_vector(ctx)
        mog_score = max(mog_vec, default=0.0)
        og_flag = float(engine.old_growth_flag(ctx))

        if _region_is_southwest(region):
            diag = dict(_mog_diag_na)
            diag.update(southwest_og_diagnostic_row(ctx))
            if len(mog_vec) > 1:
                diag["SW_MATURITY_SCORE"] = float(max(mog_vec[1:]))
            else:
                diag["SW_MATURITY_SCORE"] = pd.NA
        elif _region_is_northern(region):
            diag = dict(_mog_diag_na)
            diag.update(northern_veg_diagnostic_row(ctx, mog_vec))
        else:
            diag = dict(_mog_diag_na)

        statecd_out = c.get("STATECD")
        if statecd_out is not None and not pd.isna(statecd_out):
            try:
                statecd_out = int(statecd_out)
            except (TypeError, ValueError):
                pass

        out.append(
            {
                "PLT_CN": plt_cn,
                "CONDID": condid,
                "STATECD": statecd_out,
                "COND_STATUS_CD": _to_int_or_na(c.get("COND_STATUS_CD")),
                "FORTYPCD": _to_int_or_na(c.get("FORTYPCD")),
                "MOG_SCORE": float(mog_score),
                "OG_FLAG": og_flag,
                "AREA_BASIS": c.get("PROP_BASIS"),
                "CONDPROP_UNADJ": float(c.get("CONDPROP_UNADJ", 0.0) or 0.0),
                "PROP_FOREST": float(c.get("PROP_FOREST", 0.0) or 0.0),
                "ADFORCD": c.get("ADFORCD"),
                "EVAL_TYP": eval_typ,
                **diag,
            }
        )

    return pd.DataFrame(out)


def old_growth_area(
    db: FiaDatabase,
    *,
    states: Optional[Sequence[str]] = None,
    grp_by: Optional[Sequence[str]] = None,
    eval_typ: str = "CURR",
    land_type: str = "forest",
    method: str = "TI",
    totals: bool = True,
    variance: bool = True,
    mog_auxiliary_dir: Path | str | None = None,
    use_montana_divide_shapefile: bool = True,
) -> MOGAreaResult:
    """
    Estimate old-growth area using `custom_pse`.

    We treat each FIA condition as a record with AREA_BASIS and a numeric
    variable:

        OG_PROP = PROP_FOREST * OG_FLAG

    where ``PROP_FOREST`` is ``CONDPROP_UNADJ × landD`` (same as ``area`` / ``tpa``),
    and ``OG_FLAG`` is **old-growth only** (regional OG criteria via
    :meth:`fia_mog.engine.MOGEngine.old_growth_flag`, not ``max(MOG.vector)`` from
    maturity alone).

    This sums to plot level (with non-response adjustment by AREA_BASIS),
    then is expanded to area totals by the design frame inside `custom_pse`.
    """

    cond_mog = mog_condition_scores(
        db,
        states=states,
        eval_typ=eval_typ,
        land_type=land_type,
        mog_auxiliary_dir=mog_auxiliary_dir,
        use_montana_divide_shapefile=use_montana_divide_shapefile,
    )
    if cond_mog.empty:
        return MOGAreaResult(
            df=pd.DataFrame(),
            cond_mog=cond_mog,
            forest_area_df=pd.DataFrame(),
        )

    grp = list(grp_by or [])

    x = cond_mog.assign(
        AREA_BASIS=lambda d: d["AREA_BASIS"],
        OG_PROP=lambda d: d["PROP_FOREST"].fillna(0.0) * d["OG_FLAG"].fillna(0.0),
    )[["PLT_CN", "CONDID", "EVAL_TYP", "AREA_BASIS", "ADFORCD", "OG_PROP"]]

    res = custom_pse(
        db=db,
        x=x,
        x_vars=["OG_PROP"],
        x_grp_by=grp,
        method=method,
        totals=totals,
        variance=variance,
    )

    # Same design path as OG, but numerator is full condition forest proportion
    # (PROP_FOREST), for diagnosing OG as a share of estimated forest acres.
    x_forest = cond_mog.assign(
        AREA_BASIS=lambda d: d["AREA_BASIS"],
        FOREST_PROP=lambda d: d["PROP_FOREST"].fillna(0.0),
    )[["PLT_CN", "CONDID", "EVAL_TYP", "AREA_BASIS", "ADFORCD", "FOREST_PROP"]]

    forest_res = custom_pse(
        db=db,
        x=x_forest,
        x_vars=["FOREST_PROP"],
        x_grp_by=grp,
        method=method,
        totals=totals,
        variance=variance,
    )

    return MOGAreaResult(df=res, cond_mog=cond_mog, forest_area_df=forest_res)


__all__ = ["MOGAreaResult", "mog_condition_scores", "old_growth_area"]

