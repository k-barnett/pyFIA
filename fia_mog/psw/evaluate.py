"""
Pacific Southwest MOG (R ``FUNCTION_mapMOG.R`` ~2417–2976, Tables 12 & 19).
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, List

import pandas as pd

from ..crosswalk import pacific_southwest_site_index_class, pacific_southwest_veg_type
from ..engine import _weighted_index
from ..pnw.helpers import max_tree_and_stand_age

if TYPE_CHECKING:
    from ..engine import ConditionContext, TreeMetrics


def pacific_southwest_ponderosa_ecosub_class(ecosubcd: str | None) -> str | None:
    """
    R ponderosa ``pond.class`` (lines 2810–2817, footnote A Table 12).

    ``ECOSUBCD`` exact match ``M261Di,M`` is ported literally from R.
    """

    if ecosubcd is None:
        return None
    s = str(ecosubcd).strip()
    if s == "M261Di,M":
        return "pacific"
    if s in {"M261Ea", "M261Eb", "M261Ec", "M261Ei", "M261Ej"}:
        return "interior"
    if len(s) >= 5 and s[:5] in {"M261G", "M261D"}:
        return "interior"
    if len(s) >= 4 and s[:4] == "342B":
        return "interior"
    return "pacific"


def _pair_og(
    trees: pd.DataFrame,
    acres: float,
    stand_age: float,
    *,
    dia_min: float,
    trees_per_acre_min: float,
    age_min: float,
) -> float:
    if acres <= 0 or trees is None or len(trees) == 0:
        return 0.0
    dia = pd.to_numeric(trees["DIA"], errors="coerce")
    n = int((dia >= dia_min).sum())
    dens_ok = 1.0 if (n / acres) >= trees_per_acre_min else 0.0
    mx = max_tree_and_stand_age(stand_age, trees)
    age_ok = 0.0 if math.isnan(mx) else (1.0 if mx >= age_min else 0.0)
    return 1.0 if (dens_ok + age_ok) == 2.0 else 0.0


def _pair_og_large_subset(
    large_trees: pd.DataFrame,
    all_trees: pd.DataFrame,
    acres: float,
    stand_age: float,
    *,
    trees_per_acre_min: float,
    age_min: float,
) -> float:
    """Stem density / acre from ``large_trees`` (already diameter-filtered); age from full stand."""

    if acres <= 0 or large_trees is None:
        return 0.0
    n = len(large_trees)
    dens_ok = 1.0 if (n / acres) >= trees_per_acre_min else 0.0
    mx = max_tree_and_stand_age(stand_age, all_trees)
    age_ok = 0.0 if math.isnan(mx) else (1.0 if mx >= age_min else 0.0)
    return 1.0 if (dens_ok + age_ok) == 2.0 else 0.0


def _quaking_aspen_large(trees: pd.DataFrame) -> pd.DataFrame:
    """R SPGRPCD 44 & DIA>=18 plus conifers SPGRPCD 1:24 & DIA>=30."""

    sp = pd.to_numeric(trees.get("SPGRPCD"), errors="coerce")
    dia = pd.to_numeric(trees["DIA"], errors="coerce")
    asp = trees.loc[(sp == 44) & (dia >= 18.0)]
    con = trees.loc[sp.isin(range(1, 25)) & (dia >= 30.0)]
    if len(asp) == 0:
        return con
    if len(con) == 0:
        return asp
    return pd.concat([asp, con]).drop_duplicates()


def psw_mog_vector(ctx: "ConditionContext", metrics: "TreeMetrics") -> List[float]:
    trees = ctx.trees
    acres = float(ctx.condition_area_acres)
    stand_age = float(ctx.stand_age)
    ft = int(ctx.forest_type)

    veg = pacific_southwest_veg_type(ft)
    site = (
        pacific_southwest_site_index_class(
            trees=trees,
            stdage=ctx.stand_age,
            fldage=None,
            siteclcd=ctx.condition_siteclcd,
            siteclcdest=ctx.condition_siteclcdest,
        )
        if veg
        else None
    )

    vec: List[float] = []
    pond: str | None = None
    if veg == "ponderosa pine":
        pond = pacific_southwest_ponderosa_ecosub_class(ctx.ecosubcd)

    nwfp_inside = ctx.pnw_inside_nwfp

    # --- Table 12 old growth (requires mapped ``veg.type``; R leaves veg NA otherwise) ---
    if not veg:
        return list(_psw_mature_vector(None, None, ft, metrics))

    if veg == "coast redwood":
        dia = pd.to_numeric(trees["DIA"], errors="coerce")
        n = int((dia >= 40.0).sum())
        dens = n / acres if acres > 0 else 0.0
        vec.append(1.0 if dens >= 15.0 else 0.0)

    if veg == "conifer mixed forests" and site == "productive":
        vec.append(_pair_og(trees, acres, stand_age, dia_min=39.0, trees_per_acre_min=6.0, age_min=188.0))
    if veg == "conifer mixed forests" and site == "low":
        vec.append(_pair_og(trees, acres, stand_age, dia_min=29.0, trees_per_acre_min=5.0, age_min=256.0))

    if veg == "white fir":
        if nwfp_inside is True:
            if site == "productive":
                vec.append(_pair_og(trees, acres, stand_age, dia_min=30.0, trees_per_acre_min=5.0, age_min=160.0))
            if site == "low":
                vec.append(_pair_og(trees, acres, stand_age, dia_min=25.0, trees_per_acre_min=23.0, age_min=303.0))
        elif nwfp_inside is False:
            if site == "productive":
                vec.append(_pair_og(trees, acres, stand_age, dia_min=39.0, trees_per_acre_min=6.0, age_min=143.0))
            if site == "low":
                vec.append(_pair_og(trees, acres, stand_age, dia_min=29.0, trees_per_acre_min=8.0, age_min=239.0))

    if veg == "pacific douglas-fir" and site == "productive":
        vec.append(_pair_og(trees, acres, stand_age, dia_min=40.0, trees_per_acre_min=12.0, age_min=180.0))
    if veg == "pacific douglas-fir" and site == "low":
        vec.append(_pair_og(trees, acres, stand_age, dia_min=30.0, trees_per_acre_min=18.0, age_min=260.0))

    if veg == "douglas-fir/tanoak/madrone" and site == "productive":
        vec.append(_pair_og(trees, acres, stand_age, dia_min=30.0, trees_per_acre_min=10.0, age_min=180.0))
    if veg == "douglas-fir/tanoak/madrone" and site == "low":
        vec.append(_pair_og(trees, acres, stand_age, dia_min=30.0, trees_per_acre_min=8.0, age_min=300.0))

    if veg == "mixed subalpine (western white pine assc)" and site == "productive":
        vec.append(_pair_og(trees, acres, stand_age, dia_min=30.0, trees_per_acre_min=9.0, age_min=150.0))
    if veg == "mixed subalpine (western white pine assc)" and site == "low":
        vec.append(_pair_og(trees, acres, stand_age, dia_min=30.0, trees_per_acre_min=10.0, age_min=200.0))

    if veg == "mixed subalpine (mountain hemlock assc)" and site == "productive":
        vec.append(_pair_og(trees, acres, stand_age, dia_min=30.0, trees_per_acre_min=12.0, age_min=150.0))
    if veg == "mixed subalpine (mountain hemlock assc)" and site == "low":
        vec.append(_pair_og(trees, acres, stand_age, dia_min=30.0, trees_per_acre_min=6.0, age_min=200.0))

    if veg == "mixed subalpine (western juniper assc)":
        vec.append(_pair_og(trees, acres, stand_age, dia_min=30.0, trees_per_acre_min=5.0, age_min=200.0))

    if veg == "mixed subalpine (quaking aspen assc)" and site == "productive":
        lt = _quaking_aspen_large(trees)
        vec.append(
            _pair_og_large_subset(
                lt, trees, acres, stand_age, trees_per_acre_min=5.0, age_min=80.0
            )
        )

    if veg == "mixed subalpine (quaking aspen assc)" and site == "low":
        lt = _quaking_aspen_large(trees)
        vec.append(
            _pair_og_large_subset(
                lt, trees, acres, stand_age, trees_per_acre_min=1.0, age_min=80.0
            )
        )

    if veg == "red fir" and site == "productive":
        vec.append(_pair_og(trees, acres, stand_age, dia_min=30.0, trees_per_acre_min=8.0, age_min=150.0))
    if veg == "red fir" and site == "low":
        vec.append(_pair_og(trees, acres, stand_age, dia_min=36.0, trees_per_acre_min=5.0, age_min=200.0))

    if veg == "jeffrey pine" and site == "productive":
        vec.append(_pair_og(trees, acres, stand_age, dia_min=30.0, trees_per_acre_min=3.0, age_min=150.0))
    if veg == "jeffrey pine" and site == "low":
        vec.append(_pair_og(trees, acres, stand_age, dia_min=30.0, trees_per_acre_min=1.0, age_min=200.0))

    if veg == "lodgepole pine" and site == "productive":
        vec.append(_pair_og(trees, acres, stand_age, dia_min=36.0, trees_per_acre_min=7.0, age_min=150.0))
    if veg == "lodgepole pine" and site == "low":
        vec.append(_pair_og(trees, acres, stand_age, dia_min=36.0, trees_per_acre_min=4.0, age_min=200.0))

    if veg == "ponderosa pine" and pond == "interior" and site == "productive":
        vec.append(_pair_og(trees, acres, stand_age, dia_min=21.0, trees_per_acre_min=19.0, age_min=150.0))
    if veg == "ponderosa pine" and pond == "interior" and site == "low":
        vec.append(_pair_og(trees, acres, stand_age, dia_min=21.0, trees_per_acre_min=16.0, age_min=200.0))
    if veg == "ponderosa pine" and pond == "pacific":
        vec.append(_pair_og(trees, acres, stand_age, dia_min=30.0, trees_per_acre_min=9.0, age_min=125.0))

    last_og = vec[-1] if vec else 0.0
    out = list(vec)

    if last_og == 0.0:
        out.extend(_psw_mature_vector(veg, pond, ft, metrics))

    return out


def _ft_in_ranges(ft: int, ranges: list[tuple[int, int]]) -> bool:
    for a, b in ranges:
        if a <= ft <= b:
            return True
    return False


def _psw_mature_vector(
    veg: str | None, pond: str | None, ft: int, metrics: "TreeMetrics"
) -> List[float]:
    """R lines 2859–2970 (``local.MOG.status == 0`` gate handled by caller)."""

    m: List[float] = []

    if veg == "douglas-fir/tanoak/madrone":
        m.append(
            _weighted_index(
                metrics,
                [
                    ("ddi_score", ">=", 53.3, 0.45),
                    ("qmd_dom", ">=", 14.8, 0.29),
                    ("tpadom", "<=", 76.6, 0.25),
                ],
            )
        )

    if veg == "jeffrey pine":
        m.append(
            _weighted_index(
                metrics,
                [
                    ("qmd_dom", ">=", 10.3, 0.52),
                    ("ddi_score", ">=", 30.8, 0.25),
                    ("ht_sd", ">=", 31.5, 0.23),
                ],
            )
        )

    mixed_conifer_veg = {
        "conifer mixed forests",
        "lodgepole pine",
        "mixed subalpine (western white pine assc)",
        "mixed subalpine (mountain hemlock assc)",
    }
    if veg in mixed_conifer_veg:
        m.append(
            _weighted_index(
                metrics,
                [
                    ("qmd_dom", ">=", 13.1, 0.6),
                    ("ddi_score", ">=", 42.1, 0.4),
                ],
            )
        )

    if veg == "ponderosa pine" and pond == "interior":
        m.append(
            _weighted_index(
                metrics,
                [
                    ("qmd_dom", ">=", 13.1, 0.6),
                    ("ddi_score", ">=", 42.1, 0.4),
                ],
            )
        )

    if veg in {"coast redwood", "pacific douglas-fir"}:
        m.append(
            _weighted_index(
                metrics,
                [
                    ("ddi_score", ">=", 52.6, 0.4),
                    ("qmd_dom", ">=", 25.3, 0.35),
                    ("snag_ba_tot", ">=", 2.7, 0.26),
                ],
            )
        )

    if veg == "ponderosa pine" and pond == "pacific":
        m.append(
            _weighted_index(
                metrics,
                [
                    ("ddi_score", ">=", 52.6, 0.4),
                    ("qmd_dom", ">=", 25.3, 0.35),
                    ("snag_ba_tot", ">=", 2.7, 0.26),
                ],
            )
        )

    if veg == "red fir":
        m.append(
            _weighted_index(
                metrics,
                [
                    ("ddi_score", ">=", 48.3, 0.32),
                    ("qmd_dom", ">=", 18.1, 0.28),
                    ("ht_quart", ">=", 66.2, 0.23),
                    ("ht_sd", ">=", 43.6, 0.17),
                ],
            )
        )

    if veg == "white fir":
        m.append(
            _weighted_index(
                metrics,
                [
                    ("ddi_score", ">=", 47.5, 0.31),
                    ("ht_quart", ">=", 68.5, 0.31),
                    ("badom", ">=", 150.0, 0.21),
                    ("snag_ba_tot", ">=", 24.9, 0.16),
                ],
            )
        )

    hardwood_extra = _ft_in_ranges(
        ft,
        [(910, 912), (940, 943), (700, 722), (920, 935), (960, 976), (900, 905)],
    )
    if veg == "mixed subalpine (quaking aspen assc)" or hardwood_extra:
        m.append(
            _weighted_index(
                metrics,
                [
                    ("ddi_score", ">=", 47.5, 0.31),
                    ("ht_quart", ">=", 68.5, 0.31),
                    ("badom", ">=", 150.0, 0.21),
                    ("snag_ba_tot", ">=", 24.9, 0.16),
                ],
            )
        )

    softwood_extra = _ft_in_ranges(ft, [(180, 185), (360, 369)])
    if veg == "mixed subalpine (western juniper assc)" or softwood_extra:
        m.append(
            _weighted_index(
                metrics,
                [
                    ("qmd_dom", ">=", 14.2, 0.54),
                    ("badom", ">=", 30.9, 0.46),
                ],
            )
        )

    return m
