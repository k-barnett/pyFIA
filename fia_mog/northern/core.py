"""
Northern Region (Green et al. 1992 / R ``region %in% "northern"``) crosswalk helpers.

VEG prefix rules and habitat letter tables live alongside this module under
``fia_mog.northern``.
"""

from __future__ import annotations

import math
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

# If GIS divide polygons are unavailable: plots with longitude >= this value (less negative
# than the Rockies front) are treated as the eastern Montana MOG zone; west → western zone.
# Approximate coarse split (~110.5°W); prefer :func:`montana_plot_in_cont_divide_east` when possible.
_MT_LON_HEURISTIC_EASTERN_IF_GTE: float = -110.5


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
    plot_lon: float | None = None,
) -> str | None:
    """
    Northern MOG sub-zone: Idaho / western MT / eastern MT (+ plains states east of divide).

    For Montana, R uses ``st_filter(local.plot, ContDivideEast)`` and assigns
    eastern vs western using ``nrow(...) > 1`` and ``< 1``. For a single-geometry
    plot, :func:`fia_mog.auxiliary.montana_plot_in_cont_divide_east` uses
    **intersects** and sets this flag to ``True`` when the plot lies in that
    polygon (eastern Montana zone).

    If ``mt_east_of_continental_divide`` is ``None`` (shapefile missing, geopandas
    unavailable, or divide disabled) but ``plot_lon`` is set, a **longitude
    heuristic** (:data:`_MT_LON_HEURISTIC_EASTERN_IF_GTE`) assigns east vs west so
    habitat OG rules can still run. Without either, Montana returns ``None`` and
    habitat OG vectors stay empty.
    """

    if not state_abbrev:
        return None
    st = state_abbrev.strip().upper()
    if st in {"ID", "WA"}:
        return "northern Idaho zone"
    if st in {"WY", "SD", "ND"}:
        return "eastern Montana zone"
    if st == "MT":
        mt_east = mt_east_of_continental_divide
        if mt_east is None and plot_lon is not None:
            try:
                lon = float(plot_lon)
                if math.isfinite(lon):
                    mt_east = lon >= _MT_LON_HEURISTIC_EASTERN_IF_GTE
            except (TypeError, ValueError):
                pass
        if mt_east is True:
            return "eastern Montana zone"
        if mt_east is False:
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


# PLANTS 4-letter prefix (dominant live-tree BA) → Green et al. northern OG forest-type label.
# Only species that plausibly map to a published habitat OG rule set; unlisted / aspen → no inference.
_NORTHERN_PLANTS_PREFIX_OG: dict[str, str] = {
    "PSME": "DF",
    "TSME": "DF",
    "TSHE": "DF",
    "THPL": "DF",
    "ABGR": "DF",
    "PIPO": "PP",
    "PIED": "PP",
    "PIPU": "PP",
    "PIEN": "LP",
    "PIAB": "LP",
    "PIBR": "LP",
    "PIGL": "LP",
    "PIMA": "LP",
    "PIRU": "LP",
    "PISI": "LP",
    "PIOM": "LP",
    "PIFL": "LP",
    "LALA": "LP",
    "LAOC": "LP",
    "LALY": "LP",
    "ABLA": "LP",
    "ABLO": "LP",
    "PICO": "SAF",
    "ABCO": "GF",
    "ABMA": "GF",
    "PIAL": "WBP",
    "PIAR": "WBP",
    "THOC": "C",
    "JUOC": "C",
    "JUSC": "C",
}


def _code_matches_tree_base(code: str, base: str) -> bool:
    """True if a zone veg.code string refers to the given 4-letter overstory prefix."""

    b = base.strip().upper()
    if not b:
        return False
    cu = str(code).strip().upper()
    if cu == b or cu.startswith(b + "/") or cu.startswith(b + "-"):
        return True
    if cu.endswith("-" + b):
        return True
    first = cu.split("-", 1)[0]
    return first == b


def northern_ba_dominant_plants_prefix(
    trees: pd.DataFrame,
    species_lookup: pd.DataFrame | None,
) -> str | None:
    """
    Dominant overstory PLANTS prefix by **sum of basal area** on live trees (DIA ≥ 1 in, STATUSCD 1).

    Used when ``FLDTYPCD`` does not map to a northern OG forest type but composition still
    matches a Green et al. habitat rule group.
    """

    if trees is None or len(trees) == 0:
        return None
    norm = _normalize_species_lookup(species_lookup)
    if norm is None:
        return None
    fia_lu, plants_lu = norm
    t = trees.copy()
    t["SPCD"] = pd.to_numeric(t.get("SPCD"), errors="coerce")
    t["STATUSCD"] = pd.to_numeric(t.get("STATUSCD"), errors="coerce")
    t["DIA"] = pd.to_numeric(t.get("DIA"), errors="coerce")
    sub = t[(t["STATUSCD"] == 1) & (t["DIA"] >= 1)]
    if sub.empty or "SPCD" not in sub.columns:
        return None
    ba = _basal_area(sub["DIA"])
    sub = sub.assign(_ba=ba)
    by_sp = sub.groupby("SPCD", dropna=True)["_ba"].sum()
    if by_sp.empty:
        return None
    dom_spcd = float(by_sp.idxmax())
    m = pd.to_numeric(fia_lu, errors="coerce") == dom_spcd
    if not m.any():
        return None
    raw = plants_lu[m].iloc[0]
    if raw is None or pd.isna(raw) or str(raw).lower() in {"nan", "none"}:
        return None
    s = str(raw).strip().upper()
    return s[:4] if len(s) >= 4 else s


def infer_northern_og_type_from_species(
    trees: pd.DataFrame,
    species_lookup: pd.DataFrame | None,
) -> str | None:
    """Infer ``northern_og_forest_type`` label from BA-dominant species when ``FLDTYPCD`` is unmapped."""

    pfx = northern_ba_dominant_plants_prefix(trees, species_lookup)
    if not pfx:
        return None
    key = pfx[:4] if len(pfx) >= 4 else pfx
    return _NORTHERN_PLANTS_PREFIX_OG.get(key)


def _northern_habitat_rules_for_subregion(
    subregion: str | None,
) -> list[tuple[str, frozenset[str]]] | None:
    if subregion == "northern Idaho zone":
        from .veg_idaho import IDAHO_RULES

        return IDAHO_RULES
    if subregion == "western Montana zone":
        from .veg_west_mt import WEST_MT_RULES

        return WEST_MT_RULES
    if subregion == "eastern Montana zone":
        from .veg_east_mt import EAST_MT_RULES

        return EAST_MT_RULES
    return None


def fallback_northern_habitat_letters(
    subregion: str | None,
    veg_code: str | None,
    trees: pd.DataFrame,
    species_lookup: pd.DataFrame | None,
) -> list[str]:
    """
    When the exact ``veg_code`` is absent from the zone table, match **overstory prefix**
    (and understory prefix when present) against rule codes so habitat OG blocks can still run.

    Tries understory-constrained matches first, then relaxes to overstory-only.
    """

    rules = _northern_habitat_rules_for_subregion(subregion)
    if not rules:
        return []

    base = ""
    und_u = ""
    if veg_code and str(veg_code).strip():
        parts = str(veg_code).strip().upper().split("/", 1)
        base = parts[0].strip()[:4]
        und_u = (parts[1].strip()[:4] if len(parts) > 1 else "") or ""
    if not base:
        dp = northern_dominant_tree_plants_prefix(trees, species_lookup)
        base = (dp or "").strip().upper()[:4]
    if not base:
        ba_p = northern_ba_dominant_plants_prefix(trees, species_lookup)
        base = (ba_p or "").strip().upper()[:4]
    if not base:
        return []

    def collect(require_understory_match: bool) -> set[str]:
        out: set[str] = set()
        for letter, codeset in rules:
            for code in codeset:
                cu = str(code).strip().upper()
                if require_understory_match and und_u:
                    if "/" not in cu:
                        continue
                    pref, _, suff = cu.partition("/")
                    pref = pref.strip()
                    suff = suff.strip()
                    if pref[:4] != base:
                        continue
                    if suff.startswith(und_u) or suff[:4] == und_u:
                        out.add(letter)
                        break
                elif not require_understory_match:
                    if _code_matches_tree_base(cu, base):
                        out.add(letter)
                        break
        return out

    if und_u:
        s1 = collect(require_understory_match=True)
        if s1:
            return sorted(s1)
    s2 = collect(require_understory_match=False)
    return sorted(s2)


# Western MT: approximate subalpine break (feet) to narrow broad prefix matches.
_NORTHERN_WEST_MT_UPPER_ELEV_FT = 6500.0


def refine_habitat_letters_environmental(
    subregion: str | None,
    letters: list[str],
    *,
    siteclcd: float | None,
    elev_ft: float | None,
) -> list[str]:
    """
    When prefix fallback returns many letters, narrow using coarse FIA site class and elevation.

    Site class follows FIA ``SITECLCD`` (1 = lowest productivity, 6 = highest). Elevation uses
    plot feet above sea level (``PLOT.ELEV``). Heuristics apply only to Montana zones.
    """

    if len(letters) <= 2:
        return letters
    uniq = sorted(set(letters))
    out = list(uniq)

    if subregion == "western Montana zone" and elev_ft is not None and not pd.isna(elev_ft):
        try:
            e = float(elev_ft)
        except (TypeError, ValueError):
            e = None
        if e is not None and e >= _NORTHERN_WEST_MT_UPPER_ELEV_FT:
            upper = {"E", "F", "G", "H", "I"}
            hit = [L for L in out if L in upper]
            if hit:
                out = hit

    if siteclcd is not None and not pd.isna(siteclcd):
        try:
            sc = int(float(siteclcd))
        except (TypeError, ValueError):
            return sorted(set(out))
        if subregion in ("western Montana zone", "eastern Montana zone"):
            if sc <= 2:
                dry = {"A", "B", "C"}
                hit = [L for L in out if L in dry]
                if hit:
                    out = hit
            elif sc >= 5:
                mesic = {"D", "E", "F", "G", "H", "I", "J"}
                hit = [L for L in out if L in mesic]
                if hit:
                    out = hit

    return sorted(set(out))
