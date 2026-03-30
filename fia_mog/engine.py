from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, List, Mapping, Sequence

import numpy as np
import pandas as pd

from .crosswalk import (
    intermountain_type,
    northern_basal_area_per_acre,
    northern_habitat_letters,
    northern_og_forest_type,
    northern_subregion,
    northern_veg_code,
    statecd_to_abbrev,
)
from .northern.og_dispatch import (
    northern_east_mt_og_vector,
    northern_idaho_og_vector,
    northern_west_mt_og_vector,
)


@dataclass(frozen=True)
class ConditionContext:
    """
    Inputs needed to evaluate MOG status for a single (plot, condition).

    This module is intentionally limited to the *MOG vector* logic and expects the
    caller to provide already-filtered FIA rows for one condition and its trees.
    """

    region: str  # "eastern", "southern", "rocky", "intermountain", ...
    forest_type: int  # FLDTYPCD in the R code
    condition_area_acres: float
    stand_age: float  # max(STDAGE, FLDAGE) (raw, not binarized)
    trees: pd.DataFrame
    # optional columns used by some regions
    plot_statecd: int | None = None
    condition_fortypcd: int | None = None
    condition_physclcd: float | None = None
    condition_siteclcd: float | None = None
    condition_siteclcdest: float | None = None
    condition_adforcd: int | None = None
    condition_habtypcd1: float | int | None = None
    # Region 3 Table 9 SDI% denominator: COND.SDI_RMRS (Reinecke-style max SDI).
    condition_sdi_rmrs: float | None = None
    ecosubcd: str | None = None
    # Northern Region: optional FIA REF_SPECIES-style lookup (SPCD → PLANTS symbol)
    northern_species_lookup: pd.DataFrame | None = None
    # Full or subset of P2VEG_SUBPLOT_SPP (filtered by PLT_CN / CONDID inside helpers)
    northern_veg_subplot: pd.DataFrame | None = None
    # Montana: True if plot is east of the Continental Divide (R `st_filter` vs ContDivideEast)
    northern_mt_east_of_divide: bool | None = None
    # Pacific Northwest: PAZ raster value, NWFP membership, DWM rows, site class (Table 14), OR counties overlay.
    # ``pnw_inside_nwfp`` is also used for Pacific Southwest **white fir** (same NWFP polygon as R).
    pnw_paz_raster_value: float | int | None = None
    pnw_inside_nwfp: bool | None = None
    pnw_woody_debris: pd.DataFrame | None = None
    pnw_site_class_max: float | None = None
    pnw_plot_in_or_counties_layer: bool | None = None


@dataclass(frozen=True)
class TreeMetrics:
    tpadom: float
    badom: float
    qmd_dom: float
    ddi_score: float
    ht_quart: float
    ht_sd: float
    snag_ba_tot: float
    # southern-only (but harmless elsewhere)
    stand_basal_area_5in: float


def _basal_area_sqft(diameter_in: pd.Series) -> pd.Series:
    """
    Basal area in square feet for a diameter in inches.

    Matches `fia.estimators._basal_area`: ba = d * |d| * 0.005454.
    """

    d = pd.to_numeric(diameter_in, errors="coerce")
    return d * d.abs() * 0.005454


def _basal_term_r(dia_in: pd.Series) -> pd.Series:
    """
    Basal area term exactly as used in the R script for `badom`/`snagbatot`.

    Note: this is *not* the standard forestry basal area formula; it is ported
    to match the R implementation.
    """

    return math.pi * (dia_in / 24.0) * 2.0


def compute_tree_metrics(ctx: ConditionContext) -> TreeMetrics | None:
    """
    Port of the R per-condition metric calculations used by multiple regions.

    Returns None if required inputs are missing or the data are insufficient.
    """

    if not (ctx.condition_area_acres and ctx.condition_area_acres > 0):
        return None
    if ctx.trees is None or len(ctx.trees) == 0:
        return None

    trees = ctx.trees.copy()

    # Mature indices: dominant live trees (DIA >= 1, STATUSCD == 1, CCLCD in 1..3)
    mat_df = trees.loc[
        (trees.get("DIA") >= 1)
        & (trees.get("STATUSCD") == 1)
        & (trees.get("CCLCD").isin([1, 2, 3])),
        :,
    ]
    if len(mat_df) == 0:
        return None

    tpadom = float(len(mat_df) / ctx.condition_area_acres)

    # Ported R formula: sum(TPA_UNADJ * pi * (DIA/24) * 2)
    if "TPA_UNADJ" not in mat_df.columns or "DIA" not in mat_df.columns:
        return None
    badom = float((mat_df["TPA_UNADJ"] * _basal_term_r(mat_df["DIA"])).sum())

    if tpadom <= 0:
        return None
    qmd_dom = float(math.sqrt(badom / (tpadom * 0.005454)))

    # Diameter Diversity Index (DDI)
    dia = trees.get("DIA")
    if dia is None:
        return None
    # Ensure numeric dtype and handle pandas.NA gracefully
    dia = pd.to_numeric(dia, errors="coerce")
    dia_arr = dia.to_numpy()
    conds = [
        (dia_arr >= 2) & (dia_arr <= 9.8),
        (dia_arr >= 9.9) & (dia_arr <= 19.7),
        (dia_arr >= 19.8) & (dia_arr <= 39.4),
        (dia_arr >= 39.5),
    ]
    tpa_class = np.select(
        conds,
        ["class 0", "class 1", "class 2", "class 3"],
        default=None,
    )
    ddi_vec = pd.Series(tpa_class).dropna()
    if len(ddi_vec) == 0:
        return None
    p0 = (ddi_vec == "class 0").mean()
    p1 = (ddi_vec == "class 1").mean()
    p2 = (ddi_vec == "class 2").mean()
    p3 = (ddi_vec == "class 3").mean()
    ddi_score = float(1.0 - (p0**2 + p1**2 + p2**2 + p3**2))

    # Height vector weighted by floor(TPA_UNADJ) for living trees
    if "HT" not in trees.columns or "TPA_UNADJ" not in trees.columns or "STATUSCD" not in trees.columns:
        return None
    live = trees.loc[trees["STATUSCD"] == 1, ["HT", "TPA_UNADJ"]].dropna()
    if len(live) == 0:
        return None
    reps = np.floor(live["TPA_UNADJ"].to_numpy()).astype(int)
    reps = np.clip(reps, 0, None)
    if reps.sum() == 0:
        return None
    height_vec = np.repeat(live["HT"].to_numpy(), reps)
    if height_vec.size == 0:
        return None
    top_75 = float(np.quantile(height_vec, 0.75))
    ht_quart = float(height_vec[height_vec >= top_75].mean())
    ht_sd = float(height_vec.std(ddof=1)) if height_vec.size >= 2 else 0.0

    # Dead-tree basal area term (snags)
    dead = trees.loc[trees["STATUSCD"] == 2, :]
    snag_ba_tot = float((dead.get("TPA_UNADJ", 0) * _basal_term_r(dead.get("DIA", 0))).sum(skipna=True))

    # Southern: stand basal area for live trees DIA >= 5 using "forester's constant"
    live_basal = trees.loc[(trees["STATUSCD"] == 1) & (trees["DIA"] >= 5), ["DIA"]].copy()
    if len(live_basal) == 0:
        stand_basal_area_5in = 0.0
    else:
        basal = (live_basal["DIA"] ** 2) * 0.005454
        stand_basal_area_5in = float((basal * (1.0 / ctx.condition_area_acres)).sum())

    vals = [snag_ba_tot, ht_sd, ht_quart, ddi_score, qmd_dom, badom, tpadom]
    if any(pd.isna(v) for v in vals):
        return None

    return TreeMetrics(
        tpadom=tpadom,
        badom=badom,
        qmd_dom=qmd_dom,
        ddi_score=ddi_score,
        ht_quart=ht_quart,
        ht_sd=ht_sd,
        snag_ba_tot=snag_ba_tot,
        stand_basal_area_5in=stand_basal_area_5in,
    )


class RegionEvaluator:
    region: str

    def mog_vector(self, ctx: ConditionContext, metrics: TreeMetrics) -> List[float]:
        raise NotImplementedError


def _in_any(code: int, groups: Sequence[Iterable[int]]) -> bool:
    return any(code in g for g in groups)


class EasternEvaluator(RegionEvaluator):
    region = "eastern"

    # Old-growth rules: (forest_type_groups, stand_age>=, large_dia>=, large_density>=)
    _og_rules = [
        ([{805}], 140, 16, 10),
        ([{520, 801, 802, 809}], 140, 16, 10),
        ([{162, 163, 165, 167, 182, 184, 404, 405, 501, 502, 506, 507, 509, 510, 513, 515}], 100, 16, 20),
        ([{503, 504, 505, 511, 512, 516}], 160, 20, 5),
        ([{701, 702, 703, 704, 705, 706, 707, 708, 709}], 120, 18, 10),
        ([{104, 105, 401}], 140, 16, 10),
        ([{101, 102, 103}], 100, 12, 20),
        ([{121, 123, 124, 128, 129}], 140, 15, 10),
        ([{122, 125}], 140, 12, 10),
    ]

    # Mature rules: each is (forest_type_groups, contributions) where contributions are
    # tuples of (metric_name, comparator, threshold, weight)
    _mature_rules = [
        (
            [set([104, 105, 401]) | set(range(400, 410))],
            [
                ("qmd_dom", ">=", 14, 0.3),
                ("badom", ">=", 104.3, 0.22),
                ("snag_ba_tot", ">=", 14.5, 0.19),
                ("tpadom", "<=", 73.4, 0.16),
                ("ht_sd", ">=", 32, 0.13),
            ],
        ),
        (
            [
                set([520, 801, 802, 809, 805, 701, 702, 703, 704, 705, 706, 707, 708, 709])
                | set(range(900, 906))
                | set(range(600, 610))
                | set(range(500, 521))
                | set(range(960, 963))
            ],
            [
                ("qmd_dom", ">=", 9.9, 0.29),
                ("ht_quart", ">=", 43.3, 0.2),
                ("badom", ">=", 60.9, 0.18),
                ("tpadom", "<=", 97.6, 0.18),
                ("ht_sd", ">=", 32.9, 0.14),
            ],
        ),
        (
            [set([101, 102, 103]) | set(range(160, 169)) | set(range(380, 386)) | {390, 391}],
            [
                ("qmd_dom", ">=", 11.9, 0.3),
                ("ht_sd", ">=", 67.4, 0.22),
                ("ht_quart", ">=", 38, 0.21),
                ("tpadom", "<=", 83.2, 0.18),
                ("badom", ">=", 81.5, 0.09),  # R has a typo `>-`; ported as >=
            ],
        ),
        (
            [
                {
                    162,
                    173,  # present in the R maturity list
                    165,
                    167,
                    182,
                    184,
                    404,
                    405,
                    501,
                    502,
                    506,
                    507,
                    509,
                    510,
                    513,
                    515,
                    503,
                    504,
                    505,
                    511,
                    512,
                    516,
                }
            ],
            [
                ("qmd_dom", ">=", 12.7, 0.37),
                ("tpadom", "<=", 73.4, 0.22),
                ("ht_quart", ">=", 52.9, 0.21),
                ("ht_sd", ">=", 36.5, 0.2),
            ],
        ),
        (
            [set([122, 125, 121, 123, 124, 128, 129]) | set(range(120, 130))],
            [
                ("ddi_score", ">=", 22.2, 0.4),
                ("badom", ">=", 76.2, 0.36),
                ("ht_quart", ">=", 32, 0.24),
            ],
        ),
    ]

    def mog_vector(self, ctx: ConditionContext, metrics: TreeMetrics) -> List[float]:
        vec: List[float] = []

        # Old growth: specific rules, else default "all other forest types"
        matched_any = False
        for groups, age_thr, dia_thr, dens_thr in self._og_rules:
            if _in_any(ctx.forest_type, groups):
                matched_any = True
                status = self._evaluate_age_large_tree_og(
                    ctx=ctx, age_thr=age_thr, dia_thr=dia_thr, dens_thr=dens_thr
                )
                vec.append(float(status))

        if not matched_any:
            # R "all other forest types"
            status = self._evaluate_age_large_tree_og(ctx=ctx, age_thr=100, dia_thr=14, dens_thr=10)
            vec.append(float(status))

        if max(vec, default=0.0) >= 1.0:
            return vec

        # Mature assessment
        for groups, contributions in self._mature_rules:
            if _in_any(ctx.forest_type, groups):
                vec.append(float(_weighted_index(metrics, contributions)))
        return vec

    @staticmethod
    def _evaluate_age_large_tree_og(*, ctx: ConditionContext, age_thr: float, dia_thr: float, dens_thr: float) -> int:
        stand_age_ok = 1 if ctx.stand_age >= age_thr else 0
        live_large = ctx.trees.loc[(ctx.trees["DIA"] >= dia_thr) & (ctx.trees["STATUSCD"] == 1), :]
        dens = float(len(live_large) / ctx.condition_area_acres) if len(live_large) else 0.0
        dens_ok = 1 if dens >= dens_thr else 0
        return 1 if (stand_age_ok + dens_ok) == 2 else 0


class SouthernEvaluator(RegionEvaluator):
    region = "southern"

    # Crosswalk: FIA forest_type -> list of southern.type codes (Table 16)
    _crosswalk: List[tuple[Iterable[int], Sequence[int]]] = [
        ([104, 105, 123, 124], [2]),
        ([129], [31]),
        ([141], [26, 29]),
        ([142, 166, 407], [29]),
        ([161], [25]),
        ([162, 163, 404, 405, 409], [24, 25]),
        ([165, 167], [24]),
        ([400], [2, 24, 25, 26, 29]),
        ([401], [2]),
        ([403], [26]),
        ([406], [25]),
        ([500], [5, 13, 21, 22, 24, 27]),
        ([501], [22]),
        ([502, 515, 519], [21, 22]),
        ([504], [21, 27]),
        ([505], [21]),
        ([506, 511, 516], [5]),
        ([508], [13]),
        ([510], [21, 22, 24]),
        ([514], [22, 24]),
        ([517, 800, 801, 805], [1, 5]),
        ([520], [27]),
        ([600], [6, 10, 13, 22, 27, 28]),
        ([601, 602, 605, 706], [13]),
        ([607, 609], [14]),
        ([608, 809], [10]),
        ([700], [10, 28]),
        ([702, 703, 704], [28]),
        ([705], [13, 28]),
        ([708], [10, 13]),
        ([709], [28]),
        ([902], [31, 2]),
        ([962], [1, 5, 6, 10, 13, 21, 22, 27, 28]),
    ]

    # Old-growth thresholds per southern.type: age, large_dia, large_dens, basal_area, dead_dens
    _og_thresholds: Mapping[int, tuple[float, float, float, float, float]] = {
        1: (100, 14, 6, 40, 13),
        2: (140, 20, 6, 40, 6),
        5: (140, 30, 6, 40, 4),
        6: (120, 24, 6, 40, 4),
        10: (120, 20, 6, 40, 0),
        13: (100, 16, 6, 40, 0),
        14: (120, 8, 6, 40, 3),
        21: (130, 20, 6, 40, 26),
        22: (90, 8, 6, 10, 10),
        24: (100, 10, 6, 20, 6),
        25: (120, 19, 6, 40, 15),
        26: (80, 16, 6, 10, 0),
        27: (100, 20, 6, 40, 0),
        28: (100, 25, 6, 40, 6),
        29: (80, 9, 6, 10, 0),
        31: (120, 20, 6, 40, 14),
    }

    _mature_rules = [
        (
            [{105, 162, 402, 407, 401, 406, 409, 405}],
            [
                ("qmd_dom", ">=", 8.3, 0.42),
                ("tpadom", "<=", 111.6, 0.3),
                ("ht_quart", ">=", 39.2, 0.28),
            ],
        ),
        (
            [{141, 403}],
            [
                ("qmd_dom", ">=", 10.2, 0.31),
                ("ddi_score", ">=", 19, 0.23),
                ("tpadom", "<=", 54.7, 0.23),
                ("ht_sd", ">=", 24, 0.12),
                ("badom", ">=", 44.7, 0.12),
            ],
        ),
        (
            [{502, 510, 515, 514, 505, 504, 503, 501}],
            [
                ("qmd_dom", ">=", 9.5, 0.3),
                ("ddi_score", ">=", 22.8, 0.28),
                ("ht_quart", ">=", 44.1, 0.22),
                ("badom", ">=", 55, 0.2),
            ],
        ),
        (
            [{103, 104, 166, 142, 123, 165, 161, 164, 163, 167, 162}],
            [
                ("qmd_dom", ">=", 11.4, 0.38),
                ("tpadom", "<=", 60.4, 0.26),
                ("ht_quart", ">=", 65.8, 0.19),
                ("ht_sd", ">=", 38.6, 0.17),
            ],
        ),
        (
            [
                {
                    609,
                    520,
                    507,
                    516,
                    708,
                    608,
                    607,
                    962,
                    801,
                    703,
                    519,
                    602,
                    511,
                    802,
                    605,
                    706,
                    517,
                    809,
                    508,
                    506,
                    512,
                    905,
                    601,
                    704,
                    805,
                    702,
                    705,
                }
            ],
            [
                ("ddi_score", ">=", 30.1, 0.31),
                ("ht_quart", ">=", 43.8, 0.26),
                ("badom", ">=", 59.1, 0.22),
                ("ht_sd", ">=", 48, 0.21),
            ],
        ),
    ]

    def mog_vector(self, ctx: ConditionContext, metrics: TreeMetrics) -> List[float]:
        vec: List[float] = []

        southern_types = self._southern_types_for(ctx.forest_type)
        if not southern_types:
            return vec

        # In R, evaluation happens once per southern.type, appending each outcome.
        for st in southern_types:
            og = self._eval_southern_og(ctx=ctx, metrics=metrics, southern_type=st)
            vec.append(float(og))

            if og == 0:
                # Mature index uses FIA forest.type (not southern.type) in the R code.
                for groups, contributions in self._mature_rules:
                    if _in_any(ctx.forest_type, groups):
                        vec.append(float(_weighted_index(metrics, contributions)))
        return vec

    def _southern_types_for(self, forest_type: int) -> List[int]:
        out: List[int] = []
        for codes, types in self._crosswalk:
            if forest_type in set(codes):
                out.extend(list(types))
        # preserve duplicates? R can append duplicates; we dedupe to keep vector stable
        return sorted(set(out))

    def _eval_southern_og(self, *, ctx: ConditionContext, metrics: TreeMetrics, southern_type: int) -> int:
        thr = self._og_thresholds.get(southern_type)
        if thr is None:
            return 0
        age_thr, large_dia, large_dens, basal_thr, dead_dens_thr = thr

        stand_age_ok = 1 if ctx.stand_age >= age_thr else 0
        live = ctx.trees.loc[ctx.trees["STATUSCD"] == 1, :]
        dead = ctx.trees.loc[ctx.trees["STATUSCD"] == 2, :]
        dead_dens = float(len(dead) / ctx.condition_area_acres) if len(dead) else 0.0
        dead_ok = 1 if dead_dens >= dead_dens_thr else 0

        large = live.loc[live["DIA"] >= large_dia, :]
        large_d = float(len(large) / ctx.condition_area_acres) if len(large) else 0.0
        large_ok = 1 if large_d >= large_dens else 0

        basal_ok = 1 if metrics.stand_basal_area_5in >= basal_thr else 0

        return 1 if (stand_age_ok + dead_ok + large_ok + basal_ok) == 4 else 0


class RockyEvaluator(RegionEvaluator):
    region = "rocky"

    def mog_vector(self, ctx: ConditionContext, metrics: TreeMetrics) -> List[float]:
        vec: List[float] = []

        trees = ctx.trees.copy()
        raw_stand_age = ctx.stand_age
        trees["tree_age"] = np.where(
            ~trees.get("BHAGE").isna(),
            trees["BHAGE"],
            np.where(~trees.get("TOTAGE").isna(), trees["TOTAGE"], raw_stand_age),
        )

        n_trees_broken = float(
            len(trees.loc[(trees.get("CULLMSTOP", 0) > 0) | (trees.get("CULL", 0) > 0), :]) / ctx.condition_area_acres
        )
        n_dead_10 = float(len(trees.loc[(trees["DIA"] >= 10) & (trees["STATUSCD"] == 2), :]) / ctx.condition_area_acres)

        # Old growth, Table 8
        og_done = False
        og_done |= self._rocky_og(
            vec=vec,
            ctx=ctx,
            trees=trees,
            type_groups=[set([220, 221, 222, 224, 225, 226, 200, 201, 202, 203]) | set(range(120, 130)) | set(range(260, 272))],
            dia_thr=16,
            age_thr=200,
            dens_thr=10,
            broken_thr=1,
            dead_thr=2,
            n_trees_broken=n_trees_broken,
            n_dead=n_dead_10,
        )
        og_done |= self._rocky_og(
            vec=vec,
            ctx=ctx,
            trees=trees,
            type_groups=[set(range(900, 906))],
            dia_thr=14,
            age_thr=200,
            dens_thr=10,
            broken_thr=1,
            dead_thr=0,
            n_trees_broken=n_trees_broken,
            n_dead=n_dead_10,
        )
        og_done |= self._rocky_og(
            vec=vec,
            ctx=ctx,
            trees=trees,
            type_groups=[{280, 281}],
            dia_thr=10,
            age_thr=150,
            dens_thr=10,
            broken_thr=1,
            dead_thr=2,
            n_trees_broken=n_trees_broken,
            n_dead=n_dead_10,
        )
        og_done |= self._rocky_og(
            vec=vec,
            ctx=ctx,
            trees=trees,
            type_groups=[{180, 182, 184, 185}],
            dia_thr=12,
            age_thr=200,
            dens_thr=30,
            broken_thr=1,
            dead_thr=1,
            n_trees_broken=n_trees_broken,
            n_dead=n_dead_10,
        )
        og_done |= self._rocky_og(
            vec=vec,
            ctx=ctx,
            trees=trees,
            type_groups=[set(range(360, 370))],
            dia_thr=12,
            age_thr=200,
            dens_thr=10,
            broken_thr=0,
            dead_thr=0,
            n_trees_broken=n_trees_broken,
            n_dead=n_dead_10,
        )
        og_done |= self._rocky_og(
            vec=vec,
            ctx=ctx,
            trees=trees,
            type_groups=[set(range(970, 977))],
            dia_thr=4,
            age_thr=80,
            dens_thr=30,
            broken_thr=0,
            dead_thr=0,
            n_trees_broken=n_trees_broken,
            n_dead=n_dead_10,
        )
        og_done |= self._rocky_og(
            vec=vec,
            ctx=ctx,
            trees=trees,
            type_groups=[set(range(700, 723))],
            dia_thr=14,
            age_thr=100,
            dens_thr=20,
            broken_thr=0,
            dead_thr=0,
            n_trees_broken=n_trees_broken,
            n_dead=n_dead_10,
        )

        if max(vec, default=0.0) >= 1.0:
            return vec

        # Mature, Table 19 (subset implemented exactly as R block we extracted)
        ft = ctx.forest_type
        if ft in (set(range(900, 906)) | {703} | set(range(500, 521)) | {961, 962}):
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("ht_quart", ">=", 32.9, 0.31),
                            ("ddi_score", ">=", 18.6, 0.27),
                            ("badom", ">=", 55.1, 0.26),
                            ("ht_sd", ">=", 25.3, 0.15),
                        ],
                    )
                )
            )
        if ft in {201}:
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("ddi_score", ">=", 29.2, 0.3),
                            ("badom", ">=", 65.8, 0.21),
                            ("ht_quart", ">=", 40.6, 0.18),
                            ("qmd_dom", ">=", 9.3, 0.17),
                            ("snag_ba_tot", ">=", 21.3, 0.15),
                        ],
                    )
                )
            )
        if ft in set(range(970, 977)):
            score = _weighted_index(
                metrics,
                [
                    ("badom", ">=", 25.3, 0.3),
                    ("ddi_score", ">=", 8, 0.25),
                    ("ht_quart", ">=", 10.4, 0.24),
                    ("qmd_dom", ">=", 2.9, 0.21),
                ],
            )
            vec.append(0.5 if score >= 0.5 else 0.0)
        if ft in {280, 281}:
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("qmd_dom", ">=", 3.7, 0.56),
                            ("badom", ">=", 33.8, 0.38),
                            ("ht_sd", ">=", 17.5, 0.16),
                        ],
                    )
                )
            )
        if ft in (set(range(360, 370)) | {170, 171, 172}):
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("ddi_score", ">=", 24, 0.32),
                            ("qmd_dom", ">=", 6.5, 0.29),
                            ("ht_quart", ">=", 28.2, 0.24),
                            ("ht_sd", ">=", 21.6, 0.15),
                        ],
                    )
                )
            )
        if ft in {180, 182, 184, 185}:
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("ddi_score", ">=", 33.5, 0.55),
                            ("qmd_dom", ">=", 8.6, 0.45),
                        ],
                    )
                )
            )
        if ft in {220, 221, 222, 224, 225, 226}:
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("qmd_dom", ">=", 11.8, 0.33),
                            ("ddi_score", ">=", 31.6, 0.28),
                            ("ht_sd", ">=", 39, 0.21),
                            ("badom", ">=", 67.3, 0.18),
                        ],
                    )
                )
            )
        if ft in (set(range(260, 272)) | set(range(120, 130))):
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("ddi_score", ">=", 28.8, 0.31),
                            ("badom", ">=", 87.2, 0.27),
                            ("ht_quart", ">=", 43.5, 0.24),
                            ("ht_sd", ">=", 44.6, 0.18),
                        ],
                    )
                )
            )

        return vec

    @staticmethod
    def _rocky_og(
        *,
        vec: List[float],
        ctx: ConditionContext,
        trees: pd.DataFrame,
        type_groups: Sequence[Iterable[int]],
        dia_thr: float,
        age_thr: float,
        dens_thr: float,
        broken_thr: float,
        dead_thr: float,
        n_trees_broken: float,
        n_dead: float,
    ) -> bool:
        if not _in_any(ctx.forest_type, type_groups):
            return False
        large = trees.loc[(trees["DIA"] >= dia_thr) & (trees["STATUSCD"] == 1), :]
        large_old = large.loc[large["tree_age"] >= age_thr, :]
        dens = float(len(large_old) / ctx.condition_area_acres) if len(large_old) else 0.0
        ok_large = 1 if dens >= dens_thr else 0
        ok_broken = 1 if n_trees_broken >= broken_thr else 0
        ok_dead = 1 if n_dead >= dead_thr else 0
        status = 1 if (ok_large + ok_broken + ok_dead) == 3 else 0
        vec.append(float(status))
        return True


class _IntermountainEvaluator(RegionEvaluator):
    """
    Lightweight port of the intermountain rules.

    The original R code derives a string `intermountain.type` from several
    condition- and plot-level fields, then applies old-growth and mature
    thresholds by that label. Here we assume the caller has already
    computed and passed that label via `ConditionContext.condition_fortypcd`
    (as a string) or via an extra column on `trees` named
    ``INTERMOUNTAIN_TYPE``.
    """

    region = "intermountain"

    def mog_vector(self, ctx: ConditionContext, metrics: TreeMetrics) -> List[float]:
        # Derive intermountain.type using the R crosswalk logic.
        t = intermountain_type(
            statecd=ctx.plot_statecd or 0,
            fortypcd=ctx.condition_fortypcd or ctx.forest_type,
            physclcd=ctx.condition_physclcd,
            siteclcd=ctx.condition_siteclcd,
            adforcd=ctx.condition_adforcd,
            ecosubcd=ctx.ecosubcd,
            tree_species_codes=ctx.trees.get("SPCD", []),
        )
        if not t:
            return []

        vec: List[float] = []

        # Old-growth thresholds (Table 11) – ported as simple if/else using
        # the same labels as the R code. Only a subset is shown here; the
        # structure matches the others and can be extended as needed.
        if t in (
            "engelmann spruce - subalpine fir - warm - UT",
            "engelmann spruce - subalpine fir - warm - ID",
            "engelmann spruce - subalpine fir - cold",
            "engelmann spruce - subalpine fir - alpine",
        ):
            large = ctx.trees.loc[(ctx.trees["tree.age"] >= 220) & (ctx.trees["DIA"] >= (20 if "ID" not in t else 24)), :]
            dens = float(len(large) / ctx.condition_area_acres) if len(large) else 0.0
            vec.append(1.0 if dens >= 25 else 0.0)

        # Mature: use the same metric-threshold-weight pattern as other regions.
        # This is intentionally minimal; for most applications, northern / rocky /
        # southwest / coastal rules are more critical, and users can extend
        # intermountain rules locally following the same pattern.
        return vec


class _NorthernEvaluator(RegionEvaluator):
    """
    Northern Region (FIA Region 1) MOG rules from Green et al. (1992).

    Expects `northern_species_lookup` (SPCD → PLANTS) and `northern_veg_subplot`
    (P2VEG_SUBPLOT_SPP) on `ConditionContext` when habitat-based OG scoring is
    required. For Montana, set `northern_mt_east_of_divide` from a spatial join;
    if it is None, habitat OG is skipped for MT (maturity indices still run).
    """

    region = "northern"

    def mog_vector(self, ctx: ConditionContext, metrics: TreeMetrics) -> List[float]:
        vec: List[float] = []
        ft = int(ctx.forest_type)
        og = northern_og_forest_type(ft)
        st = statecd_to_abbrev(ctx.plot_statecd)
        sub = northern_subregion(
            st,
            mt_east_of_continental_divide=ctx.northern_mt_east_of_divide,
        )

        basal = northern_basal_area_per_acre(ctx.trees, ctx.condition_area_acres)
        plt_cn = ctx.trees["PLT_CN"].iloc[0] if "PLT_CN" in ctx.trees.columns else None
        condid = ctx.trees["CONDID"].iloc[0] if "CONDID" in ctx.trees.columns else None

        veg_code = northern_veg_code(
            ctx.trees,
            ctx.northern_species_lookup,
            ctx.northern_veg_subplot,
            plt_cn,
            condid,
        )
        letters = northern_habitat_letters(sub, veg_code)

        if sub == "northern Idaho zone":
            vec.extend(
                northern_idaho_og_vector(
                    letters,
                    og,
                    raw_stand_age=ctx.stand_age,
                    basal_area_per_acre=basal,
                    condition_area_acres=ctx.condition_area_acres,
                    trees=ctx.trees,
                )
            )
        elif sub == "western Montana zone":
            vec.extend(
                northern_west_mt_og_vector(
                    letters,
                    og,
                    raw_stand_age=ctx.stand_age,
                    basal_area_per_acre=basal,
                    condition_area_acres=ctx.condition_area_acres,
                    trees=ctx.trees,
                )
            )
        elif sub == "eastern Montana zone":
            vec.extend(
                northern_east_mt_og_vector(
                    letters,
                    og,
                    raw_stand_age=ctx.stand_age,
                    basal_area_per_acre=basal,
                    condition_area_acres=ctx.condition_area_acres,
                    trees=ctx.trees,
                )
            )

        # Mature (Table 19): Northern always appends these component scores in R.
        if og == "DF" or ft in set(range(200, 204)):
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("ddi_score", ">=", 32.6, 0.34),
                            ("badom", ">=", 82.5, 0.33),
                            ("qmd_dom", ">=", 10.3, 0.33),
                        ],
                    )
                )
            )

        if og in {"SAF", "GF", "WP"} or ft in set(range(260, 272)):
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("ddi_score", ">=", 24.0, 0.44),
                            ("ht_sd", ">=", 49.6, 0.3),
                            ("ht_quart", ">=", 39.2, 0.26),
                        ],
                    )
                )
            )

        if ft in set(range(910, 913)) | set(range(700, 723)) | set(range(900, 906)) | set(range(500, 521)) | set(
            range(970, 977)
        ):
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("ddi_score", ">=", 23.9, 0.31),
                            ("badom", ">=", 62.0, 0.28),
                            ("ht_quart", ">=", 38.4, 0.26),
                            ("ht_sd", ">=", 28.0, 0.15),
                        ],
                    )
                )
            )

        if ft in set(range(300, 306)):
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("ddi_score", ">=", 45.0, 0.38),
                            ("ht_sd", ">=", 74.4, 0.28),
                            ("ht_quart", ">=", 69.2, 0.21),
                            ("tpadom", "<=", 70.0, 0.13),
                        ],
                    )
                )
            )

        if og == "LP" or ft in (280, 281):
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("ht_quart", ">=", 25.0, 0.28),
                            ("ddi_score", ">=", 14.6, 0.26),
                            ("badom", ">=", 43.6, 0.26),
                            ("ht_sd", ">=", 24.0, 0.19),
                        ],
                    )
                )
            )

        # R uses `ddiscore >- 24`; interpreted as >= 24 (see MOG errata).
        if ft in set(range(360, 370)) | set(range(180, 186)):
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("ddi_score", ">=", 24.0, 0.3),
                            ("ht_quart", ">=", 28.6, 0.25),
                            ("qmd_dom", ">=", 7.0, 0.25),
                            ("ht_sd", ">=", 29.4, 0.2),
                        ],
                    )
                )
            )

        if og == "PP" or ft in set(range(220, 227)):
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("ddi_score", ">=", 31.5, 0.36),
                            ("qmd_dom", ">=", 13.0, 0.34),
                            ("ht_sd", ">=", 40.7, 0.3),
                        ],
                    )
                )
            )

        if og == "L":
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("qmd_dom", ">=", 15.8, 0.31),
                            ("ddi_score", ">=", 53.0, 0.31),
                            ("ht_sd", ">=", 80.9, 0.21),
                            ("tpadom", "<=", 69.0, 0.16),
                        ],
                    )
                )
            )

        return vec


class _SouthwestEvaluator(RegionEvaluator):
    """
    Southwest rules (``FUNCTION_mapMOG.R``): ERU from ``HABTYPCD1`` via
    :func:`fia_mog.southwest.core.southwest_eru` (also exposed on
    :mod:`fia_mog.crosswalk`), or override with ``ConditionContext.ecosubcd`` as a
    precomputed ERU string. Implementation: :mod:`fia_mog.southwest.evaluate`.
    """

    region = "southwest"

    def mog_vector(self, ctx: ConditionContext, metrics: TreeMetrics) -> List[float]:
        from .southwest.evaluate import southwest_mog_vector

        return southwest_mog_vector(ctx, metrics)


class _PacificSouthwestEvaluator(RegionEvaluator):
    """
    California / Region 5 MOG: Table 12 (``veg.type`` × site index) and Table 19 maturity.

    White fir uses NWFP polygon overlap on ``ConditionContext.pnw_inside_nwfp`` (R lines
    2497–2500). Other types do not require that flag.
    """

    region = "pacific southwest"

    def mog_vector(self, ctx: ConditionContext, metrics: TreeMetrics) -> List[float]:
        from .psw.evaluate import psw_mog_vector

        return psw_mog_vector(ctx, metrics)


class _PacificNorthwestEvaluator(RegionEvaluator):
    """
    Region 6 / Pacific Northwest: PAZ raster, NWFP polygon, Table 13/14 OG,
    then Table 19 maturity (see ``fia_mog.pnw``).
    """

    region = "pacific northwest"

    def mog_vector(self, ctx: ConditionContext, metrics: TreeMetrics) -> List[float]:
        from .pnw.evaluate import pacific_northwest_mog_vector

        # R always knows NWFP membership via ``st_filter``; without it we cannot
        # choose inside vs outside OG branches.
        if ctx.pnw_inside_nwfp is None:
            return []
        return pacific_northwest_mog_vector(ctx, metrics)


class MOGEngine:
    """
    Evaluate MOG vector for a single condition.

    This is the object you call from your plot/condition iteration code; it
    returns the full per-condition `MOG.vector` equivalent.
    """

    def __init__(self) -> None:
        self._evaluators: Mapping[str, RegionEvaluator] = {
            "eastern": EasternEvaluator(),
            "southern": SouthernEvaluator(),
            "rocky": RockyEvaluator(),
            # The remaining regions have specialized rules that depend on
            # additional landscape or vegetation inputs not yet wired into
            # ConditionContext. For now they are handled via dedicated
            # evaluators that assume the caller has precomputed the required
            # categorical labels.
            "intermountain": _IntermountainEvaluator(),
            "northern": _NorthernEvaluator(),
            "southwest": _SouthwestEvaluator(),
            "pacific southwest": _PacificSouthwestEvaluator(),
            "pacific northwest": _PacificNorthwestEvaluator(),
        }

    def mog_vector(self, ctx: ConditionContext) -> List[float]:
        metrics = compute_tree_metrics(ctx)
        if metrics is None:
            return []
        evaluator = self._evaluators.get(ctx.region)
        if evaluator is None:
            return []
        return evaluator.mog_vector(ctx, metrics)


def _weighted_index(metrics: TreeMetrics, contributions: Sequence[tuple[str, str, float, float]]) -> float:
    score = 0.0
    for attr, op, thr, w in contributions:
        val = getattr(metrics, attr)
        ok = False
        if op == ">=":
            ok = val >= thr
        elif op == "<=":
            ok = val <= thr
        else:
            raise ValueError(f"Unsupported comparator: {op}")
        score += w if ok else 0.0
    return float(score)

