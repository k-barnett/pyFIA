from __future__ import annotations

"""
Crosswalking helpers for the MOG engine.

This module centralizes the code/label mappings that in the original R
implementation live inside ``FUNCTION_mapMOG.R`` / ``mapMOG``. It is part of the
stand-alone ``fia_mog`` package (see ``fia_py/fia_mog``).
"""

from typing import Any, Iterable

import numpy as np
import pandas as pd

from .northern.core import (
    NORTHERN_KEEP_VEG_PREFIXES,
    northern_basal_area_per_acre,
    northern_dominant_tree_plants_prefix,
    northern_dominant_understory_plants_prefix,
    northern_habitat_letters,
    northern_og_forest_type,
    northern_subregion,
    northern_veg_code,
)
from .pnw.paz import paz_value_to_group
from .southwest.core import southwest_eru


def classify_region(adforcd: Any, state_abbrev: str | None = None) -> str | None:
    """
    Map FIA `ADFORCD` (admin forest code) / state to a broad MOG region.

    This mirrors the region assignment logic at lines 176–203 in the R
    script, but only uses the leading digit of ADFORCD when available.

    `adforcd` may be a float, int, string, or pandas.NA; this helper
    normalizes those cases so callers don't have to pre-clean.

    State-only fallback (no ``ADFORCD``) is ordered: first match wins. Washington is
    **not** grouped with the Northern MOG region here—Region 6 uses the Pacific
    Northwest rules. (``northern_subregion`` still maps WA with ID for *internal*
    Northern-habitat labels when a caller explicitly uses the northern path.)
    """

    # Use pandas-friendly NA check to avoid ambiguous boolean for pd.NA
    if adforcd is None or pd.isna(adforcd):
        if not state_abbrev:
            return None
        st = state_abbrev.upper()
        if st in {"MT", "ID", "WY", "SD", "ND"}:
            return "northern"
        if st in {"CO", "KA", "NE"}:
            return "rocky"
        if st in {"AZ", "NM"}:
            return "southwest"
        if st in {"NV", "UT"}:
            return "intermountain"
        if st in {"TX", "OK", "AR", "LA", "MS", "KY", "TN", "AL", "GA", "FL", "SC", "NC", "VA"}:
            return "southern"
        if st == "CA":
            return "pacific southwest"
        if st in {"WA", "OR"}:
            return "pacific northwest"
        if st in {
            "MN",
            "IA",
            "MO",
            "IL",
            "WI",
            "IN",
            "MI",
            "OH",
            "WV",
            "DC",
            "PA",
            "MD",
            "DE",
            "NJ",
            "CT",
            "RI",
            "MA",
            "ME",
            "NH",
            "VT",
            "NY",
        }:
            return "eastern"
        return None

    # Normalize adforcd to an int for the first-digit dispatch
    try:
        adfor_int = int(adforcd)
    except (TypeError, ValueError):
        # Fallback to state-only logic if adforcd is not interpretable
        return classify_region(None, state_abbrev=state_abbrev)

    first = str(adfor_int)[0]
    return {
        "1": "northern",
        "2": "rocky",
        "3": "southwest",
        "4": "intermountain",
        "5": "pacific southwest",
        "6": "pacific northwest",
        "8": "southern",
        "9": "eastern",
    }.get(first)


def intermountain_type(
    *,
    statecd: int,
    fortypcd: int,
    physclcd: float | int | None,
    siteclcd: float | int | None,
    adforcd: int | None,
    ecosubcd: str | None,
    tree_species_codes: Iterable[int],
) -> str | None:
    """
    Port of the `intermountain.type` assignment logic (lines 1155–1198).

    Returns the string label used in the R code (e.g.
    "engelmann spruce - subalpine fir - warm - UT").
    """

    t: str | None = None
    phys = None if physclcd is None or np.isnan(physclcd) else float(physclcd)
    site = None if siteclcd is None or np.isnan(siteclcd) else float(siteclcd)
    ecos = ecosubcd

    def any_sp(codes: Iterable[int]) -> bool:
        s = set(tree_species_codes)
        return any(c in s for c in codes)

    # Warm/cold Englemann spruce-subalpine fir, UT vs ID
    if statecd != 16 and fortypcd in {265, 261}:
        t = "engelmann spruce - subalpine fir - warm - UT"
    elif statecd != 16 and fortypcd == 266 and phys is not None and phys > 20:
        t = "engelmann spruce - subalpine fir - warm - UT"
    elif statecd != 16 and fortypcd == 266 and not any_sp({113, 101, 72}):
        t = "engelmann spruce - subalpine fir - warm - UT"
    elif statecd == 16 and fortypcd in {265, 261}:
        t = "engelmann spruce - subalpine fir - warm - ID"
    elif statecd == 16 and fortypcd == 266 and phys is not None and phys > 20:
        t = "engelmann spruce - subalpine fir - warm - ID"
    elif statecd == 16 and fortypcd == 266 and not any_sp({113, 101, 72}):
        t = "engelmann spruce - subalpine fir - warm - ID"
    elif fortypcd == 266 and phys is not None and phys < 20:
        t = "engelmann spruce - subalpine fir - cold"
    elif fortypcd == 266 and any_sp({113, 101, 72}):
        t = "engelmann spruce - subalpine fir - cold"
    elif fortypcd == 268 and site is not None and site < 7:
        t = "engelmann spruce - subalpine fir - cold"
    elif fortypcd == 268 and site == 7:
        t = "engelmann spruce - subalpine fir - alpine"
    elif fortypcd == 367:
        t = "whitebark pine"
    elif fortypcd == 365:
        t = "bristlecone pine"
    elif fortypcd == 201 and site is not None and site < 6:
        t = "douglas - fir - high"
    elif fortypcd == 201 and site is not None and site >= 6:
        t = "douglas - fir - low"
    elif fortypcd == 267:
        t = "grand fir"
    elif fortypcd == 269:
        t = "blue spruce"
    elif fortypcd in {371, 262} and phys is not None and phys < 20:
        t = "conifer mixed forest - low"
    elif fortypcd in {371, 262} and phys is not None and phys >= 20:
        t = "conifer mixed forest - productive"
    elif fortypcd == 901 and phys is not None and phys < 20:
        t = "aspen - dry"
    elif fortypcd == 901 and phys is not None and phys > 20:
        t = "aspen - mesic"
    elif fortypcd == 281:
        t = "lodgepole pine"
    elif fortypcd == 366 and site is not None and site > 6:
        t = "limber pine - lower"
    elif fortypcd == 366 and site is not None and site <= 6:
        t = "limber pine - montane"
    elif fortypcd in {220, 221, 222, 225} and adforcd in {402, 412, 413, 414} and site is not None and site > 5:
        t = "ponderosa pine - n - seral"
    elif fortypcd in {220, 221, 222, 225} and adforcd in {402, 412, 413, 414} and site is not None and site <= 5:
        t = "ponderosa pine - n - climax"
    elif fortypcd in {220, 221, 222, 225} and (adforcd not in {402, 412, 413, 414}) and site is not None and site > 5:
        t = "ponderosa pine - rm - seral"
    elif fortypcd in {220, 221, 222, 225} and (adforcd not in {402, 412, 413, 414}) and site is not None and site <= 5:
        t = "ponderosa pine - rm - climax"
    elif fortypcd in {182, 184, 185} and adforcd in {402, 403, 412, 413, 414, 415, 417, 420} and phys is not None and phys < 20:
        t = "pinyon - juniper - nw - low"
    elif fortypcd in {182, 184, 185} and adforcd in {402, 403, 412, 413, 414, 415, 417, 420} and phys is not None and phys > 20:
        t = "pinyon - juniper - nw - high"
    elif fortypcd in {182, 184, 185} and adforcd in {401, 407, 408, 410} and phys is not None and phys < 20:
        t = "pinyon - juniper - se - low"
    elif fortypcd in {182, 184, 185} and adforcd in {401, 407, 408, 410} and phys is not None and phys > 20:
        t = "pinyon - juniper - se - high"

    # If still unassigned, fall back to ECOSUBCD rules (lines 1191–1196)
    if t is None and ecos is not None:
        if fortypcd in {182, 184, 185} and adforcd in {418, 419} and ecos in {"M331Dn", "M331Do", "M331Dv", "M331Di"}:
            t = "pinyon - juniper - nw - low" if phys is not None and phys < 20 else "pinyon - juniper - nw - high"
        elif fortypcd in {182, 184, 185} and adforcd in {418, 419}:
            t = "pinyon - juniper - se - low" if phys is not None and phys < 20 else "pinyon - juniper - se - high"

    return t


def pacific_southwest_veg_type(forest_type: int | float | None) -> str | None:
    """
    Map FLDTYPCD to `veg.type` for the Pacific Southwest region
    (lines 2419–2433).
    """

    if forest_type is None or np.isnan(forest_type):
        return None
    ft = int(forest_type)
    if ft == 341:
        return "coast redwood"
    if ft in {371, 226, 361}:
        return "conifer mixed forests"
    if ft == 261:
        return "white fir"
    if ft in {201, 202}:
        return "pacific douglas-fir"
    if ft == 941:
        return "douglas-fir/tanoak/madrone"
    if ft in {241, 342, 365, 366, 367}:
        return "mixed subalpine (western white pine assc)"
    if ft == 270:
        return "mixed subalpine (mountain hemlock assc)"
    if ft == 369:
        return "mixed subalpine (western juniper assc)"
    if ft == 901:
        return "mixed subalpine (quaking aspen assc)"
    if ft == 262:
        return "red fir"
    if ft == 225:
        return "jeffrey pine"
    if ft == 281:
        return "lodgepole pine"
    if ft == 221:
        return "ponderosa pine"
    return None


def pacific_southwest_site_index_class(
    *,
    trees: pd.DataFrame,
    stdage: float | int | None,
    fldage: float | int | None,
    siteclcd: float | int | None,
    siteclcdest: float | int | None = None,
) -> str | None:
    """
    Compute Dunning-style site productivity class for Pacific Southwest.

    This ports the logic around lines 2435–2449:
    - Prefer per-tree BHAGE/TOTAGE if present to compute a mean site
      index from heights.
    - Fallback to ``max(SITECLCD, SITECLCDEST, na.rm=TRUE)`` (R lines 2444–2447).
    """

    # Determine raw stand age (max of STDAGE/FLDAGE) as a fallback.
    ages = [a for a in [stdage, fldage] if a is not None and not np.isnan(a)]
    raw_stand_age = max(ages) if ages else None

    # Compute per-tree age vector
    if "BHAGE" in trees.columns or "TOTAGE" in trees.columns:
        age_series = trees.get("BHAGE").copy()
        if age_series is None:
            age_series = pd.Series(np.nan, index=trees.index)
        if "TOTAGE" in trees.columns:
            age_series = age_series.fillna(trees["TOTAGE"])
        if raw_stand_age is not None:
            age_series = age_series.fillna(raw_stand_age)
    else:
        age_series = pd.Series(raw_stand_age, index=trees.index)

    if "HT" in trees.columns and not age_series.isna().all():
        ht = trees["HT"]
        valid = (~ht.isna()) & (~age_series.isna()) & (age_series > 0)
        if valid.any():
            site_index = (ht[valid] * (0.25489 + (29.377 / age_series[valid]))).mean()
        else:
            site_index = np.nan
    else:
        site_index = np.nan

    if not np.isnan(site_index):
        return "low" if site_index < 45 else "productive"

    # Fallback: max of measured / estimated site class (R ``site.clcd``)
    sc_vals: list[float] = []
    for x in (siteclcd, siteclcdest):
        if x is None or (isinstance(x, float) and np.isnan(x)):
            continue
        try:
            sc_vals.append(float(x))
        except (TypeError, ValueError):
            continue
    if not sc_vals:
        return None
    site_clcd_max = max(sc_vals)
    return "low" if site_clcd_max >= 5 else "productive"


def pacific_northwest_paz_group(paz_val: int | float | None) -> str | None:
    """
    Map numeric PAZ raster cell value to plant association zone group
    (``FUNCTION_mapMOG.R`` lines 2987–3003). Implementation: ``fia_mog.pnw.paz``.
    """

    return paz_value_to_group(paz_val)


# FIA STATECD → USPS (CONUS + common); extend if needed.
STATECD_TO_ST_ABBREV: dict[int, str] = {
    1: "AL",
    2: "AK",
    4: "AZ",
    5: "AR",
    6: "CA",
    8: "CO",
    9: "CT",
    10: "DE",
    11: "DC",
    12: "FL",
    13: "GA",
    16: "ID",
    17: "IL",
    18: "IN",
    19: "IA",
    20: "KS",
    21: "KY",
    22: "LA",
    23: "ME",
    24: "MD",
    25: "MA",
    26: "MI",
    27: "MN",
    28: "MS",
    29: "MO",
    30: "MT",
    31: "NE",
    32: "NV",
    33: "NH",
    34: "NJ",
    35: "NM",
    36: "NY",
    37: "NC",
    38: "ND",
    39: "OH",
    40: "OK",
    41: "OR",
    42: "PA",
    44: "RI",
    45: "SC",
    46: "SD",
    47: "TN",
    48: "TX",
    49: "UT",
    50: "VT",
    51: "VA",
    53: "WA",
    54: "WV",
    55: "WI",
    56: "WY",
}


def statecd_to_abbrev(statecd: int | float | str | None) -> str | None:
    """Map FIA ``STATECD`` (numeric or string from CSV) to a two-letter abbreviation."""

    if statecd is None:
        return None
    try:
        if pd.isna(statecd):
            return None
    except (TypeError, ValueError):
        pass
    if isinstance(statecd, float) and np.isnan(statecd):
        return None
    num = pd.to_numeric(statecd, errors="coerce")
    if pd.isna(num):
        return None
    try:
        k = int(num)
    except (TypeError, ValueError, OverflowError):
        return None
    return STATECD_TO_ST_ABBREV.get(k)
