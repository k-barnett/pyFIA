"""
Southwest MOG vector (``FUNCTION_mapMOG.R``): Table 9 old growth, Table 19 maturity.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, List

import numpy as np
import pandas as pd

from .core import southwest_eru
from ..engine import _basal_term_r, _weighted_index

if TYPE_CHECKING:
    from ..engine import ConditionContext, TreeMetrics


def southwest_mog_vector(ctx: "ConditionContext", metrics: "TreeMetrics") -> List[float]:
    """
    ERU from ``HABTYPCD1`` via :func:`fia_mog.southwest.core.southwest_eru`, or
    override with ``ConditionContext.ecosubcd`` as a precomputed ERU string.
    """

    eru = None
    if ctx.condition_habtypcd1 is not None and not pd.isna(ctx.condition_habtypcd1):
        eru = southwest_eru(ctx.condition_habtypcd1)
    elif isinstance(ctx.ecosubcd, str):
        eru = ctx.ecosubcd
    if not eru:
        return []

    vec: List[float] = []

    Zeide_b = 1.6064
    calc = ctx.trees.copy()
    calc["DIA"] = pd.to_numeric(calc.get("DIA"), errors="coerce")
    calc = calc.loc[~calc["DIA"].isna(), :]

    relative_sdi = 0.0
    if len(calc) > 0:
        dia_all = calc["DIA"].to_numpy()
        dr_all = (np.mean(dia_all**Zeide_b)) ** (1.0 / Zeide_b)
        n_all = dr_all ** (-Zeide_b) if dr_all > 0 else 0.0

        dia_large = calc.loc[calc["DIA"] >= 18, "DIA"].to_numpy()
        if len(dia_large) > 0:
            dr_large = (np.mean(dia_large**Zeide_b)) ** (1.0 / Zeide_b)
            n_large = dr_large ** (-Zeide_b) if dr_large > 0 else 0.0
        else:
            n_large = 0.0

        if not np.isfinite(n_large) or n_large < 0:
            n_large = 0.0
        if n_all > 0 and np.isfinite(n_all):
            relative_sdi = 100.0 * (n_large / n_all)
        else:
            relative_sdi = 0.0

    sw_qmd = 0.0
    sw_live = ctx.trees.copy()
    sw_live["DIA"] = pd.to_numeric(sw_live.get("DIA"), errors="coerce")
    sw_live = sw_live.loc[(sw_live.get("STATUSCD") == 1) & (sw_live["DIA"] >= 10), :]
    area = float(ctx.condition_area_acres)
    if len(sw_live) > 0 and area > 0 and "TPA_UNADJ" in sw_live.columns:
        sw_tpa = float(len(sw_live)) / area
        tpa_u = pd.to_numeric(sw_live["TPA_UNADJ"], errors="coerce").fillna(0.0)
        sw_ba = float((_basal_term_r(sw_live["DIA"]) * tpa_u).sum())
        if sw_tpa > 0:
            sw_qmd_val = math.sqrt(sw_ba / (sw_tpa * 0.005454))
            sw_qmd = float(sw_qmd_val) if np.isfinite(sw_qmd_val) else 0.0

    local_og = 0

    if eru in {
        "spruce-fir forest",
        "mixed conifer with aspen",
        "bristlecone pine",
        "pinyon juniper evergreen shrub",
        "pinyon juniper (persistent)",
        "pinyon juniper sagebrush",
        "pinyon juniper deciduous shrub",
        "gambel oak shrubland",
        "arizona walnut",
        "rio grande cottonwood/shrub",
        "narrowleaf cottonwood-spruce, narrowleaf cottonwood/shrub",
        "upper montane conifer/willow",
    }:
        local_og = 1 if sw_qmd >= 18 else 0
    elif eru == "mixed conifer - frequent fire":
        local_og = 1 if relative_sdi >= 56 else 0
    elif eru in {"ponderosa pine forest", "ponderosa pine/willow"}:
        local_og = 1 if relative_sdi >= 57 else 0
    elif eru == "ponderosa pine - evergreen oak":
        local_og = 1 if relative_sdi >= 26 else 0
    elif eru == "pinyon juniper grass":
        local_og = 1 if relative_sdi >= 29 else 0
    elif eru in {"juniper grass", "semi-desert grassland"}:
        local_og = 1 if relative_sdi >= 36 else 0
    elif eru in {"madrean pinyon-oak", "madrean encinal woodland"}:
        local_og = 1 if relative_sdi >= 20 else 0
    else:
        local_og = 0

    vec.append(float(local_og))

    if local_og == 0:
        ft = int(ctx.forest_type)

        if (
            eru
            in {
                "arizona walnut",
                "rio grande cottonwood/shrub",
                "gambel oak shrubland",
                "narrowleaf cottonwood-spruce, narrowleaf cottonwood/shrub",
                "upper montane conifer/willow",
                "other",
            }
            or ft in set(range(970, 977))
        ):
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("qmd_dom", ">=", 3.5, 0.34),
                            ("ddi_score", ">=", 7.7, 0.34),
                            ("ht_quart", ">=", 10.8, 0.2),
                            ("tpadom", "<=", 69.5, 0.12),
                        ],
                    )
                )
            )

        if eru == "juniper grass":
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("qmd_dom", ">=", 10.7, 0.3),
                            ("ht_quart", ">=", 11.2, 0.27),
                            ("ddi_score", ">=", 19, 0.27),
                            ("ht_sd", ">=", 4, 0.17),
                        ],
                    )
                )
            )

        if eru == "madrean encinal woodland":
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("qmd_dom", ">=", 8.8, 0.36),
                            ("ht_quart", ">=", 15.2, 0.3),
                            ("ddi_score", ">=", 16.8, 0.18),
                            ("tpadom", "<=", 56.4, 0.16),
                        ],
                    )
                )
            )

        if eru == "madrean pinyon-oak":
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("qmd_dom", ">=", 8.3, 0.32),
                            ("ht_quart", ">=", 14.4, 0.28),
                            ("ddi_score", ">=", 23.8, 0.23),
                            ("ht_sd", ">=", 10.4, 0.16),
                        ],
                    )
                )
            )

        if eru == "mixed conifer - frequent fire":
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("ddi_score", ">=", 21.4, 0.41),
                            ("qmd_dom", ">=", 13.3, 0.28),
                            ("ht_sd", ">=", 44.7, 0.21),
                        ],
                    )
                )
            )

        if eru in {"mixed conifer with aspen", "bristlecone pine"}:
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("ddi_score", ">=", 34.5, 0.39),
                            ("ht_sd", ">=", 41.2, 0.24),
                            ("ht_quart", ">=", 36.3, 0.22),
                            ("snag_ba_tot", "<=", 15, 0.15),
                        ],
                    )
                )
            )

        if eru in {"pinyon juniper grass", "pinyon juniper sagebrush", "semi-desert grassland"}:
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("ddi_score", ">=", 19.6, 0.29),
                            ("qmd_dom", ">=", 9.5, 0.26),
                            ("ht_quart", ">=", 12.8, 0.26),
                            ("ht_sd", ">=", 6.4, 0.19),
                        ],
                    )
                )
            )

        if eru in {
            "pinyon juniper (persistent)",
            "pinyon juniper deciduous shrub",
            "pinyon juniper evergreen shrub",
        } or ft in set(range(180, 186)):
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("ddi_score", ">=", 20.2, 0.46),
                            ("qmd_dom", ">=", 9.2, 0.34),
                            ("ht_quart", ">=", 13.3, 0.21),
                        ],
                    )
                )
            )

        if eru == "ponderosa pine forest":
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("ddi_score", ">=", 24.3, 0.45),
                            ("badom", ">=", 40, 0.28),
                            ("qmd_dom", ">=", 13.5, 0.27),
                        ],
                    )
                )
            )

        if eru in {"ponderosa pine - evergreen oak", "ponderosa pine/willow"}:
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("ddi_score", ">=", 32.4, 0.5),
                            ("qmd_dom", ">=", 9, 0.32),
                            ("ht_sd", ">=", 24.1, 0.18),
                        ],
                    )
                )
            )

        if eru == "spruce-fir forest" or ft in set(range(200, 204)):
            vec.append(
                float(
                    _weighted_index(
                        metrics,
                        [
                            ("ddi_score", ">=", 32.4, 0.24),
                            ("ht_sd", ">=", 51.8, 0.22),
                            ("qmd_dom", ">=", 11.4, 0.19),
                            ("ht_quart", ">=", 43.5, 0.19),
                            ("badom", ">=", 57.4, 0.17),
                        ],
                    )
                )
            )

    return vec
