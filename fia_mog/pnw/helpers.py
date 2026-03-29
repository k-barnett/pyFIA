"""Shared numeric helpers for PNW OGSI / Table 14 logic (ported from R)."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd


def acres_to_hectares(acres: float) -> float:
    return float(acres) / 2.471


def tree_age_years(trees: pd.DataFrame) -> pd.Series:
    """R ``tree.age``: ``BHAGE`` if present else ``TOTAGE``."""

    if "BHAGE" in trees.columns:
        bh = pd.to_numeric(trees["BHAGE"], errors="coerce")
    else:
        bh = pd.Series(np.nan, index=trees.index)
    if "TOTAGE" in trees.columns:
        tot = pd.to_numeric(trees["TOTAGE"], errors="coerce")
    else:
        tot = pd.Series(np.nan, index=trees.index)
    return bh.where(~bh.isna(), tot)


def p_live_stems(trees: pd.DataFrame) -> float:
    """R ``nrow(live)/nrow(all)``."""

    if trees is None or len(trees) == 0:
        return 0.0
    n = len(trees)
    live = trees.loc[trees.get("STATUSCD") == 1]
    return float(len(live) / n)


def mean_live_diameter(trees: pd.DataFrame) -> float:
    live = trees.loc[trees.get("STATUSCD") == 1]
    if len(live) == 0:
        return float("nan")
    dia = pd.to_numeric(live["DIA"], errors="coerce")
    return float(dia.mean())


def density_per_ha_stems(df: pd.DataFrame, cond_ha: float) -> float:
    if cond_ha <= 0:
        return 0.0
    return float(len(df) / cond_ha)


def diam_class_densities_per_ha(trees: pd.DataFrame, cond_ha: float) -> tuple[float, float, float, float]:
    """R inside-NWFP diameter classes (all stems), densities per hectare."""

    dia = pd.to_numeric(trees["DIA"], errors="coerce")
    c1 = ((dia >= 2) & (dia < 9.9)).sum()
    c2 = ((dia >= 9.9) & (dia < 19.7)).sum()
    c3 = ((dia >= 19.7) & (dia < 39.4)).sum()
    c4 = (dia >= 39.4).sum()
    inv = 1.0 / cond_ha if cond_ha > 0 else 0.0
    return (float(c1 * inv), float(c2 * inv), float(c3 * inv), float(c4 * inv))


def diameter_diversity_score(c1: float, c2: float, c3: float, c4: float) -> float:
    """
    R ``diamDiv.score`` (lines 3065–3069 and analogs).

    Note: the class-3 branch uses **transformed** class-1 in ``ifelse(class.1.trees < 40, ...)``.
    """

    c1t = 0.005 * c1 if c1 < 200 else 1.0
    c2t = (1.0 / 75.0) * c2 if c2 < 75 else 1.0
    c3t = 0.025 * c3 if c1t < 40 else 1.0
    c4t = (1.0 / 30.0) * c4 if c4 < 30 else 1.0
    return 10.0 * (c1t + 2.0 * c2t + 3.0 * c3t + 4.0 * c4t)


def woody_debris_cover_pct(woody: pd.DataFrame | None) -> float:
    """Sum ``COVER_PCT`` for rows with ``TRANSDIA >= 9.8`` (R inside-NWFP)."""

    if woody is None or woody.empty:
        return 0.0
    td = pd.to_numeric(woody.get("TRANSDIA", 0), errors="coerce")
    sub = woody.loc[td >= 9.8]
    if sub.empty:
        return 0.0
    cov = pd.to_numeric(sub.get("COVER_PCT", 0), errors="coerce").fillna(0.0)
    return float(cov.sum())


def pw3(x: float, b1: float, b2: float, s0: float, k0: float, s1: float, k1: float, s2: float, k2: float) -> float:
    """Three-segment piecewise linear function (R ``if / if / if`` chain)."""

    if x < b1:
        return s0 + k0 * x
    if x < b2:
        return s1 + k1 * x
    return s2 + k2 * x


def clip100(x: float) -> float:
    return float(min(100.0, x))


def max_tree_and_stand_age(stand_age: float, trees: pd.DataFrame) -> float:
    ages = [stand_age]
    ta = tree_age_years(trees)
    ages.extend(float(v) for v in ta.dropna().tolist() if not math.isnan(float(v)))
    return max(ages) if ages else float("nan")
