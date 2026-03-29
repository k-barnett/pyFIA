"""
Northern Region (Green et al. 1992 / R ``region %in% "northern"``) crosswalk helpers.

VEG prefix rules and habitat letter tables live alongside this module under
``fia_mog.northern``.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from fia.estimators import _basal_area

# First four characters of P2VEG `VEG_SPCD` must be in this set (R `keep.veg`).
NORTHERN_KEEP_VEG_PREFIXES: frozenset[str] = frozenset(
    """
    AGSP FEID FESC PURT SYAL BERE PRIV SHCA PHMA SMST CARU VAGL ARUV XETE VACA CAGE SPBE
    JUCO ARCO SYOR CLUN ASCA MEFE ARNU SETR TABR LIBO OPHO ATFI ADPE COOC GYDR EQAR GATR
    STAM LICA CACA LEGL LUHI VASC CLPS PSME RIMO ABLA FIED PUTR PRVI ASCR SEST ALSI ALFI
    UBO THOC GAGE JUSC MARE SYOC PSSP PIPO PIAB PIBR PIEN PIGL PIMA PIPU PIRU PISI PIOM
    """.split()
)


def northern_og_forest_type(fldtypcd: int | float) -> str | None:
    """
    Map FLDTYPCD to the Northern old-growth forest type label used in the R script
    (`northern.OG.forest.type`).
    """

    try:
        ft = int(fldtypcd)
    except (TypeError, ValueError):
        return None
    og: str | None = None
    if 200 <= ft <= 203:
        og = "DF"
    if 220 <= ft <= 226:
        og = "PP"
    if ft in (320, 321):
        og = "L"
    if ft in (280, 281):
        og = "LP"
    if ft == 368:
        og = "Y"
    if ft == 267:
        og = "GF"
    if ft in (265, 266, 268):
        og = "SAF"
    if ft in (240, 241):
        og = "WP"
    if ft == 301:
        og = "WH"
    if ft == 270:
        og = "MAF"
    if ft == 367:
        og = "WBP"
    if ft == 304:
        og = "C"
    return og


def northern_subregion(
    state_abbrev: str | None,
    *,
    mt_east_of_continental_divide: bool | None,
) -> str | None:
    """
    Northern MOG sub-zone: Idaho / western MT / eastern MT (+ plains states east of divide).

    For Montana, R uses ``st_filter(local.plot, ContDivideEast)`` and assigns
    eastern vs western using ``nrow(...) > 1`` and ``< 1``. For a single-geometry
    plot, :func:`fia_mog.auxiliary.montana_plot_in_cont_divide_east` uses
    **intersects** and sets this flag to ``True`` when the plot lies in that
    polygon (eastern Montana zone). If ``None`` for MT, habitat-based OG is
    skipped; maturity indices still run.
    """

    if not state_abbrev:
        return None
    st = state_abbrev.strip().upper()
    if st in {"ID", "WA"}:
        return "northern Idaho zone"
    if st in {"WY", "SD", "ND"}:
        return "eastern Montana zone"
    if st == "MT":
        if mt_east_of_continental_divide is True:
            return "eastern Montana zone"
        if mt_east_of_continental_divide is False:
            return "western Montana zone"
        return None
    return None


def northern_basal_area_per_acre(trees: pd.DataFrame, condition_area_acres: float) -> float:
    """
    Living trees DIA >= 5 in: unweighted sum of per-record basal area / condition acres.

    Matches ``FUNCTION_mapMOG.R`` (no ``TPA_UNADJ`` on each stem). This is **not** the
    same as :func:`fia.estimators.tpa` plot totals, which use ``BAA = basal_area(DIA)
    * TPA_UNADJ * domain``. We intentionally reuse :func:`fia.estimators._basal_area`
    so the inch→ft² factor matches ``tpa`` for each tree.
    """

    if condition_area_acres <= 0 or trees is None or len(trees) == 0:
        return 0.0
    t = trees.copy()
    t["DIA"] = pd.to_numeric(t.get("DIA"), errors="coerce")
    t["STATUSCD"] = pd.to_numeric(t.get("STATUSCD"), errors="coerce")
    sub = t[(t["DIA"] >= 5) & (t["STATUSCD"] == 1)]
    if sub.empty:
        return 0.0
    ba = _basal_area(sub["DIA"])
    return float(ba.sum() / condition_area_acres)


def _normalize_species_lookup(species_lookup: pd.DataFrame) -> tuple[pd.Series, pd.Series] | None:
    """Return (FIA species codes numeric series, PLANTS codes) for merge like R USFStrees (FIA.Code → PLANTS.Code)."""

    if species_lookup is None or len(species_lookup) == 0:
        return None
    df = species_lookup.copy()
    col_fia = None
    for c in df.columns:
        cu = str(c).upper().replace(" ", "_").replace(".", "_")
        if cu in {"SPCD", "FIA_CODE", "FIA_CD"}:
            col_fia = c
            break
    if col_fia is None:
        return None
    col_plants = None
    for c in df.columns:
        cu = str(c).upper().replace(" ", "_").replace(".", "_")
        if cu in {"PLANTS_CODE", "PLANTS_SYMBOL", "PLANTS_SYMBOL_NRCS"}:
            col_plants = c
            break
        if "PLANTS" in cu and "CODE" in cu:
            col_plants = c
            break
    if col_plants is None:
        return None
    fia = pd.to_numeric(df[col_fia], errors="coerce")
    plants = df[col_plants].astype(str)
    return fia, plants


def northern_dominant_tree_plants_prefix(
    trees: pd.DataFrame,
    species_lookup: pd.DataFrame | None,
) -> str | None:
    """
    Mean CCLCD by SPCD → species with minimum mean → PLANTS symbol first 4 chars (R northern).
    """

    if trees is None or len(trees) == 0 or "SPCD" not in trees.columns or "CCLCD" not in trees.columns:
        return None
    norm = _normalize_species_lookup(species_lookup)
    if norm is None:
        return None
    fia_lu, plants_lu = norm
    t = trees.copy()
    t["SPCD"] = pd.to_numeric(t["SPCD"], errors="coerce")
    t["CCLCD"] = pd.to_numeric(t["CCLCD"], errors="coerce")
    grp = t.groupby("SPCD", dropna=True)["CCLCD"].mean()
    if grp.empty:
        return None
    dom_spcd = float(grp.idxmin())
    m = pd.to_numeric(fia_lu, errors="coerce") == dom_spcd
    if not m.any():
        return None
    raw = plants_lu[m].iloc[0]
    if raw is None or pd.isna(raw) or str(raw).lower() in {"nan", "none"}:
        return None
    s = str(raw).strip().upper()
    return s[:4] if len(s) >= 4 else s


def _veg_spcd_prefix_r(s: object) -> str:
    """First four characters of VEG_SPCD, matching R substr(..., 1, 4) on trimmed codes."""

    t = str(s).strip().upper()
    return t[:4] if t else ""


def northern_dominant_understory_plants_prefix(
    veg_df: pd.DataFrame | None,
    plt_cn: Any,
    condid: Any,
    *,
    keep_prefixes: frozenset[str] | None = None,
) -> str | None:
    """
    Dominant understory 4-letter code (FUNCTION_mapMOG.R northern block).

    Mirrors the R loop over ``unique(local.veg$VEG_SPCD)``, mean ``COVER_PCT`` per
    distinct code, then ``substr(...,0,4)`` into ``keep.veg``, then max dominance
    with first tie-breaker.
    """

    if veg_df is None or len(veg_df) == 0:
        return None
    k = keep_prefixes or NORTHERN_KEEP_VEG_PREFIXES
    v = veg_df.copy()
    if "PLT_CN" in v.columns:
        v = v[v["PLT_CN"] == plt_cn]
    if "CONDID" in v.columns:
        v = v[v["CONDID"] == condid]
    if v.empty or "VEG_SPCD" not in v.columns or "COVER_PCT" not in v.columns:
        return None
    v["COVER_PCT"] = pd.to_numeric(v["COVER_PCT"], errors="coerce")
    v = v[v["VEG_SPCD"].map(_veg_spcd_prefix_r).isin(k)]
    if v.empty:
        return None

    dom_rows: list[tuple[str, float]] = []
    for u in v["VEG_SPCD"].drop_duplicates():
        sp4 = _veg_spcd_prefix_r(u)
        if sp4 not in k:
            continue
        sub = v[v["VEG_SPCD"] == u]
        mean_cov = float(sub["COVER_PCT"].mean(skipna=True))
        if np.isnan(mean_cov):
            continue
        dom_rows.append((sp4, mean_cov))

    if not dom_rows:
        return None
    max_cov = max(m for _, m in dom_rows)
    for sp4, m in dom_rows:
        if m == max_cov:
            return sp4
    return None


def northern_veg_code(
    trees: pd.DataFrame,
    species_lookup: pd.DataFrame | None,
    veg_df: pd.DataFrame | None,
    plt_cn: Any,
    condid: Any,
) -> str | None:
    """Dominant overstory PLANTS prefix [ / dominant understory prefix]."""

    base = northern_dominant_tree_plants_prefix(trees, species_lookup)
    if base is None:
        return None
    und = northern_dominant_understory_plants_prefix(veg_df, plt_cn, condid)
    if und:
        return f"{base}/{und}"
    return base


def northern_habitat_letters(subregion: str | None, veg_code: str | None) -> list[str]:
    if not subregion or not veg_code:
        return []
    from .veg_east_mt import east_mt_habitat_letters
    from .veg_idaho import idaho_habitat_letters
    from .veg_west_mt import west_mt_habitat_letters

    if subregion == "northern Idaho zone":
        return idaho_habitat_letters(veg_code)
    if subregion == "western Montana zone":
        return west_mt_habitat_letters(veg_code)
    if subregion == "eastern Montana zone":
        return east_mt_habitat_letters(veg_code)
    return []
