"""
Outside-NWFP old-growth rules (R Table 14, ~3627–3941).
"""

from __future__ import annotations

import math
from typing import List

import pandas as pd

from .helpers import max_tree_and_stand_age


def _density_age_bin(
    trees: pd.DataFrame,
    acres: float,
    stand_age: float,
    *,
    dia_min: float,
    trees_per_acre_min: float,
    age_min: float,
) -> float:
    """Binary OG if large-tree density and max age both meet thresholds (R pattern)."""

    if acres <= 0 or trees is None or len(trees) == 0:
        return 0.0
    dia = pd.to_numeric(trees["DIA"], errors="coerce")
    n = int((dia >= dia_min).sum())
    dens_ok = 1.0 if (n / acres) >= trees_per_acre_min else 0.0
    mx = max_tree_and_stand_age(stand_age, trees)
    if math.isnan(mx):
        age_ok = 0.0
    else:
        age_ok = 1.0 if mx >= age_min else 0.0
    return 1.0 if (dens_ok + age_ok) == 2.0 else 0.0


def _site_int(site_class: float | None) -> int | None:
    if site_class is None or (isinstance(site_class, float) and math.isnan(site_class)):
        return None
    try:
        return int(round(float(site_class)))
    except (TypeError, ValueError):
        return None


def regional_geog_white_fir(ecosubcd: str | None, in_or_counties_layer: bool | None) -> str:
    """
    R lines 3634–3637: ``M242C`` prefix ⇒ central; intersect OR counties layer flips central → non-central.
    """

    is_central = False
    if ecosubcd is not None:
        s = str(ecosubcd).strip()
        if len(s) >= 5:
            is_central = s[:5] == "M242C"
    if is_central and in_or_counties_layer is True:
        return "non-central"
    return "central" if is_central else "non-central"


def outside_nwfp_og_vector(
    paz_group: str | None,
    trees: pd.DataFrame,
    condition_area_acres: float,
    stand_age: float,
    site_class: float | None,
    ecosubcd: str | None,
    in_or_counties_layer: bool | None,
) -> List[float]:
    """0/1 scores appended in R order (ponderosa may append two)."""

    vec: List[float] = []
    if not paz_group or trees is None or len(trees) == 0:
        return vec

    acres = float(condition_area_acres)
    sc = _site_int(site_class)

    if paz_group == "white fir - grand fir":
        geo = regional_geog_white_fir(ecosubcd, in_or_counties_layer)
        if sc is None:
            return vec
        if geo == "central":
            if sc in (1, 2):
                vec.append(
                    _density_age_bin(
                        trees,
                        acres,
                        stand_age,
                        dia_min=21.0,
                        trees_per_acre_min=15.0,
                        age_min=150.0,
                    )
                )
            else:
                vec.append(
                    _density_age_bin(
                        trees,
                        acres,
                        stand_age,
                        dia_min=21.0,
                        trees_per_acre_min=10.0,
                        age_min=150.0,
                    )
                )
        else:  # non-central
            if sc in (1, 2):
                vec.append(
                    _density_age_bin(
                        trees,
                        acres,
                        stand_age,
                        dia_min=21.0,
                        trees_per_acre_min=20.0,
                        age_min=150.0,
                    )
                )
            else:
                vec.append(
                    _density_age_bin(
                        trees,
                        acres,
                        stand_age,
                        dia_min=21.0,
                        trees_per_acre_min=10.0,
                        age_min=150.0,
                    )
                )
        return vec

    if paz_group == "douglas fir":
        vec.append(
            _density_age_bin(
                trees,
                acres,
                stand_age,
                dia_min=21.0,
                trees_per_acre_min=8.0,
                age_min=150.0,
            )
        )
        return vec

    if paz_group == "lodgepole pine":
        vec.append(
            _density_age_bin(
                trees,
                acres,
                stand_age,
                dia_min=12.0,
                trees_per_acre_min=60.0,
                age_min=120.0,
            )
        )
        return vec

    if paz_group == "silver fir":
        if sc is None:
            return vec
        if sc == 5:
            vec.append(
                _density_age_bin(
                    trees,
                    acres,
                    stand_age,
                    dia_min=22.0,
                    trees_per_acre_min=9.0,
                    age_min=260.0,
                )
            )
        elif sc == 6:
            vec.append(
                _density_age_bin(
                    trees,
                    acres,
                    stand_age,
                    dia_min=22.0,
                    trees_per_acre_min=1.0,
                    age_min=360.0,
                )
            )
        elif sc in (2, 3):
            vec.append(
                _density_age_bin(
                    trees,
                    acres,
                    stand_age,
                    dia_min=26.0,
                    trees_per_acre_min=6.0,
                    age_min=180.0,
                )
            )
        elif sc == 4:
            vec.append(
                _density_age_bin(
                    trees,
                    acres,
                    stand_age,
                    dia_min=25.0,
                    trees_per_acre_min=7.0,
                    age_min=200.0,
                )
            )
        return vec

    if paz_group == "ponderosa pine":
        if sc is None:
            return vec
        if sc in (1, 2, 3, 4):
            vec.append(
                _density_age_bin(
                    trees,
                    acres,
                    stand_age,
                    dia_min=21.0,
                    trees_per_acre_min=13.0,
                    age_min=150.0,
                )
            )
            vec.append(
                _density_age_bin(
                    trees,
                    acres,
                    stand_age,
                    dia_min=31.0,
                    trees_per_acre_min=3.0,
                    age_min=200.0,
                )
            )
        elif sc > 4:
            vec.append(
                _density_age_bin(
                    trees,
                    acres,
                    stand_age,
                    dia_min=21.0,
                    trees_per_acre_min=10.0,
                    age_min=150.0,
                )
            )
            vec.append(
                _density_age_bin(
                    trees,
                    acres,
                    stand_age,
                    dia_min=31.0,
                    trees_per_acre_min=2.0,
                    age_min=200.0,
                )
            )
        return vec

    if paz_group == "subalpine fir":
        if sc is None:
            return vec
        if sc in (1, 2):
            vec.append(
                _density_age_bin(
                    trees,
                    acres,
                    stand_age,
                    dia_min=21.0,
                    trees_per_acre_min=10.0,
                    age_min=150.0,
                )
            )
        elif sc > 4:
            vec.append(
                _density_age_bin(
                    trees,
                    acres,
                    stand_age,
                    dia_min=13.0,
                    trees_per_acre_min=10.0,
                    age_min=150.0,
                )
            )
        return vec

    if paz_group == "western hemlock":
        if sc is None:
            return vec
        if sc == 1:
            vec.append(
                _density_age_bin(
                    trees,
                    acres,
                    stand_age,
                    dia_min=42.0,
                    trees_per_acre_min=8.0,
                    age_min=200.0,
                )
            )
        elif sc == 2:
            vec.append(
                _density_age_bin(
                    trees,
                    acres,
                    stand_age,
                    dia_min=35.0,
                    trees_per_acre_min=8.0,
                    age_min=200.0,
                )
            )
        elif sc == 3:
            vec.append(
                _density_age_bin(
                    trees,
                    acres,
                    stand_age,
                    dia_min=31.0,
                    trees_per_acre_min=8.0,
                    age_min=200.0,
                )
            )
        elif sc in (4, 5):
            vec.append(
                _density_age_bin(
                    trees,
                    acres,
                    stand_age,
                    dia_min=21.0,
                    trees_per_acre_min=8.0,
                    age_min=200.0,
                )
            )
        return vec

    return vec
