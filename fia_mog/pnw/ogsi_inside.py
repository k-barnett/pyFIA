"""
Inside-NWFP OGSI200 logic (R ``FUNCTION_mapMOG.R`` ~3010–3624, Table 13).
"""

from __future__ import annotations

from typing import List

import pandas as pd

from .helpers import (
    acres_to_hectares,
    clip100,
    diam_class_densities_per_ha,
    diameter_diversity_score,
    mean_live_diameter,
    p_live_stems,
    pw3,
    woody_debris_cover_pct,
)


def inside_nwfp_og_vector(
    paz_group: str | None,
    trees: pd.DataFrame,
    condition_area_acres: float,
    woody_debris: pd.DataFrame | None,
) -> List[float]:
    """Return 0/1 scores to append to ``MOG.vector`` (one element when rules run)."""

    vec: List[float] = []
    if not paz_group or trees is None or len(trees) == 0:
        return vec
    if p_live_stems(trees) < 0.1:
        return vec

    ha = acres_to_hectares(condition_area_acres)
    mean_live = mean_live_diameter(trees)
    if pd.isna(mean_live):
        return vec

    live = trees.loc[trees.get("STATUSCD") == 1]
    dead = trees.loc[trees.get("STATUSCD") == 2]
    live_dia = pd.to_numeric(live["DIA"], errors="coerce")
    dead_dia = pd.to_numeric(dead["DIA"], errors="coerce")
    debris_cov = woody_debris_cover_pct(woody_debris)
    c1, c2, c3, c4 = diam_class_densities_per_ha(trees, ha)
    diam_div = diameter_diversity_score(c1, c2, c3, c4)

    def live_ld(thr: float) -> float:
        return float((live_dia >= thr).sum() / ha)

    def dead_ld(thr: float) -> float:
        return float((dead_dia >= thr).sum() / ha)

    # --- white fir - grand fir (threshold 48)
    if paz_group == "white fir - grand fir":
        if mean_live >= 0.5 * 29.5:
            lld = live_ld(29.5)
            tree = clip100(pw3(lld, 14.814, 41.973, 0.0, 3.3751, 36.3636, 0.9205, 54.7831, 0.4816))
            dld = dead_ld(19.7)
            if dld < 9.876:
                snag = 0.0 + 5.0627 * lld
            elif dld >= 9.876 and lld < 19.8088:
                snag = 25.1429 + 2.5169 * lld
            else:
                snag = 66.9783 + 0.4049 * lld
            snag = clip100(snag)
            deb = clip100(pw3(debris_cov, 1.6277, 3.4429, 0.0, 30.7181, 27.5823, 13.7725, 63.8187, 3.2476))
            ogsi = (tree + snag + deb + diam_div) / 4.0
            vec.append(1.0 if ogsi >= 48.0 else 0.0)
        return vec

    # --- juniper (70.95, tree only)
    if paz_group == "juniper":
        if mean_live >= 0.5 * 19.7:
            lld = live_ld(19.7)
            tree = clip100(pw3(lld, 0.1, 14.8708, 0.0, 500.0, 49.8307, 1.6925, 66.5006, 0.5715))
            vec.append(1.0 if tree >= 70.95 else 0.0)
        return vec

    # --- mountain hemlock (41.34)
    if paz_group == "mountain hemlock":
        if mean_live >= 0.5 * 29.5:
            lld = live_ld(29.5)
            tree = clip100(pw3(lld, 9.876, 34.5663, 0.0, 5.0627, 40.0001, 1.0125, 56.6049, 0.5321))
            dld = dead_ld(19.7)
            if dld < 12.345:
                snag = 0.0 + 4.0502 * lld
            elif dld >= 12.345 and lld < 24.7468:
                snag = 25.1145 + 2.0158 * lld
            else:
                snag = 62.6278 + 0.4999 * lld
            snag = clip100(snag)
            deb = clip100(pw3(debris_cov, 2.1817, 4.3632, 0.0, 22.9179, 24.9982, 11.4597, 59.4239, 3.5698))
            ogsi = (tree + snag + deb + diam_div) / 4.0
            vec.append(1.0 if ogsi >= 41.34 else 0.0)
        return vec

    # --- oak woodland (62.46, tree only)
    if paz_group == "oak woodland":
        if mean_live >= 0.5 * 19.7:
            lld = live_ld(19.7)
            tree = clip100(pw3(lld, 2.469, 23.0339, 0.0, 20.2511, 46.9985, 1.2156, 66.9116, 0.3511))
            vec.append(1.0 if tree >= 62.46 else 0.0)
        return vec

    # --- ponderosa pine inside (67.95, tree only)
    if paz_group == "ponderosa pine":
        if mean_live >= 0.5 * 29.5:
            lld = live_ld(29.5)
            tree = clip100(pw3(lld, 2.469, 12.345, 0.0, 20.2511, 43.75, 2.5313, 61.1419, 1.1225))
            vec.append(1.0 if tree >= 67.95 else 0.0)
        return vec

    # --- port orford cedar (45.01)
    if paz_group == "port orford cedar":
        if mean_live >= 0.5 * 29.5:
            lld = live_ld(29.5)
            tree = clip100(pw3(lld, 12.345, 32.097, 0.0, 4.0502, 34.375, 1.2656, 58.5609, 0.5121))
            dld = dead_ld(19.7)
            if dld < 14.8708:
                snag = 0.0 + 3.3622 * lld
            elif dld >= 14.8708 and lld < 31.5365:
                snag = 27.6925 + 1.5 * lld
            else:
                snag = 30.0436 + 1.4255 * lld
            snag = clip100(snag)
            deb = clip100(pw3(debris_cov, 0.9715, 2.6162, 0.0, 51.4641, 35.2318, 15.2005, 61.9655, 4.9821))
            ogsi = (tree + snag + deb + diam_div) / 4.0
            vec.append(1.0 if ogsi >= 45.01 else 0.0)
        return vec

    # --- shasta red fir (52.81)
    if paz_group == "shasta red fir":
        if mean_live >= 0.5 * 29.5:
            lld = live_ld(29.5)
            tree = clip100(pw3(lld, 19.752, 46.911, 0.0, 2.5313, 31.8181, 0.9205, 38.4615, 0.7788))
            dld = dead_ld(19.7)
            if dld < 7.407:
                snag = 0.0 + 6.7503 * lld
            elif dld >= 7.407 and lld < 14.814:
                snag = 25.0 + 3.3751 * lld
            else:
                snag = 64.8879 + 0.6826 * lld
            snag = clip100(snag)
            deb = clip100(pw3(debris_cov, 0.6477, 2.2328, 0.0, 77.1962, 39.7847, 15.7716, 66.1432, 3.9666))
            ogsi = (tree + snag + deb + diam_div) / 4.0
            vec.append(1.0 if ogsi >= 52.81 else 0.0)
        return vec

    # --- silver fir (43.39)
    if paz_group == "silver fir":
        if mean_live >= 0.5 * 29.5:
            lld = live_ld(29.5)
            tree = clip100(pw3(lld, 24.2677, 54.318, 0.0, 2.0603, 29.8106, 0.8319, 46.3446, 0.5275))
            dld = dead_ld(19.7)
            if dld < 21.42:
                snag = 0.0 + 2.3342 * lld
            elif dld >= 21.42 and lld < 37.1711:
                snag = 16.0024 + 1.5871 * lld
            else:
                snag = 55.5568 + 0.523 * lld
            snag = clip100(snag)
            deb = clip100(pw3(debris_cov, 4.3974, 7.7551, 0.0, 11.3703, 17.2593, 7.4454, 60.9019, 1.8179))
            ogsi = (tree + snag + deb + diam_div) / 4.0
            vec.append(1.0 if ogsi >= 43.39 else 0.0)
        return vec

    # --- sitka spruce (59.96), large live DIA >= 39.4
    if paz_group == "sitka spruce":
        if mean_live >= 0.5 * 39.4:
            lld = live_ld(39.4)
            tree = clip100(pw3(lld, 9.2587, 28.3935, 0.0, 5.4002, 37.9032, 1.3065, 44.086, 1.0887))
            dld = dead_ld(19.7)
            if dld < 9.876:
                snag = 0.0 + 5.0627 * lld
            elif dld >= 9.876 and lld < 16.0485:
                snag = 10.0 + 4.0502 * lld
            else:
                snag = 51.4492 + 1.4674 * lld
            snag = clip100(snag)
            deb = clip100(pw3(debris_cov, 4.6871, 7.1243, 0.0, 10.6674, 1.9207, 10.2576, 53.8812, 2.9643))
            ogsi = (tree + snag + deb + diam_div) / 4.0
            vec.append(1.0 if ogsi >= 59.96 else 0.0)
        return vec

    # --- subalpine fir (45.08), large live >= 19.7
    if paz_group == "subalpine fir":
        if mean_live >= 0.5 * 19.7:
            lld = live_ld(19.7)
            tree = clip100(pw3(lld, 14.814, 59.3242, 0.0, 3.3751, 41.6794, 0.5616, 55.1433, 0.3347))
            dld = dead_ld(19.7)
            if dld < 2.6195:
                snag = 0.0 + 19.0876 * lld
            elif dld >= 2.6195 and lld < 14.8708:
                snag = 44.6546 + 2.0406 * lld
            else:
                snag = 68.6455 + 0.4273 * lld
            snag = clip100(snag)
            deb = clip100(pw3(debris_cov, 1.2612, 3.0509, 0.0, 39.6447, 32.383, 13.9684, 66.5385, 2.7733))
            ogsi = (tree + snag + deb + diam_div) / 4.0
            vec.append(1.0 if ogsi >= 45.08 else 0.0)
        return vec

    # --- tanoak (48.22), snag large >= 39.4
    if paz_group == "tanoak":
        if mean_live >= 0.5 * 39.4:
            lld = live_ld(39.4)
            tree = clip100(pw3(lld, 12.345, 24.69, 0.0, 4.0502, 25.0, 2.0251, 46.7832, 1.1428))
            dld = dead_ld(39.4)
            if dld < 4.938:
                snag = 0.0 + 10.1255 * lld
            elif dld >= 4.938 and lld < 14.9374:
                snag = 37.6542 + 2.5001 * lld
            else:
                snag = 64.8584 + 0.6789 * lld
            snag = clip100(snag)
            deb = clip100(pw3(debris_cov, 1.3635, 3.4982, 0.0, 36.6703, 34.032, 11.7109, 61.7148, 3.7976))
            ogsi = (tree + snag + deb + diam_div) / 4.0
            vec.append(1.0 if ogsi >= 48.22 else 0.0)
        return vec

    # --- western hemlock (44.63), live & snag large >= 39.4
    if paz_group == "western hemlock":
        if mean_live >= 0.5 * 39.4:
            lld = live_ld(39.4)
            tree = clip100(pw3(lld, 9.876, 27.254, 0.0, 5.0627, 35.7923, 1.4386, 55.2341, 0.7252))
            dld = dead_ld(39.4)
            if dld < 6.7116:
                snag = 0.0 + 7.4497 * lld
            elif dld >= 6.7116 and lld < 12.345:
                snag = 20.2151 + 4.4378 * lld
            else:
                snag = 62.5286 + 1.0102 * lld
            snag = clip100(snag)
            deb = clip100(pw3(debris_cov, 3.9884, 7.1244, 0.0, 12.5363, 18.2047, 7.9719, 59.9428, 2.1134))
            ogsi = (tree + snag + deb + diam_div) / 4.0
            vec.append(1.0 if ogsi >= 44.63 else 0.0)
        return vec

    # --- douglas fir (50.81)
    if paz_group == "douglas fir":
        if mean_live >= 0.5 * 29.5:
            lld = live_ld(29.5)
            tree = clip100(pw3(lld, 4.938, 24.69, 0.0, 10.1255, 43.75, 1.2656, 59.5297, 0.6265))
            dld = dead_ld(19.7)
            if dld < 2.469:
                snag = 0.0 + 20.2511 * lld
            elif dld >= 2.469 and lld < 7.5099:
                snag = 37.7551 + 4.9594 * lld
            else:
                snag = 70.3665 + 0.6169 * lld
            snag = clip100(snag)
            deb = clip100(pw3(debris_cov, 0.784, 1.909, 0.0, 63.7755, 32.5777, 22.2222, 65.8906, 4.7717))
            ogsi = (tree + snag + deb + diam_div) / 4.0
            vec.append(1.0 if ogsi >= 50.81 else 0.0)
        return vec

    # --- lodgepole pine (61.46, tree only, large >= 9.8)
    if paz_group == "lodgepole pine":
        if mean_live >= 0.5 * 9.8:
            lld = live_ld(9.8)
            tree = clip100(pw3(lld, 55.984, 247.752, 0.0, 0.8931, 42.7016, 0.1303, 51.7325, 0.0939))
            vec.append(1.0 if tree >= 61.46 else 0.0)
        return vec

    # --- jeffrey pine - knobcone pine (61.77, tree only)
    if paz_group == "jeffrey pine - knobcone pine":
        if mean_live >= 0.5 * 29.5:
            lld = live_ld(29.5)
            tree = clip100(pw3(lld, 4.938, 17.283, 0.0, 10.1255, 40.0, 2.0251, 57.9434, 0.9868))
            vec.append(1.0 if tree >= 61.77 else 0.0)
        return vec

    # --- redwood (50.94)
    if paz_group == "redwood":
        if mean_live >= 0.5 * 39.4:
            lld = live_ld(39.4)
            tree = clip100(pw3(lld, 16.6657, 37.6522, 0.0, 3.0001, 30.147, 1.1912, 47.2727, 0.7364))
            dld = dead_ld(39.4)
            if dld < 2.469:
                snag = 0.0 + 20.2511 * lld
            elif dld >= 2.469 and lld < 5.5319:
                snag = 29.8478 + 8.162 * lld
            else:
                snag = 63.2309 + 2.1274 * lld
            snag = clip100(snag)
            deb = clip100(pw3(debris_cov, 2.0109, 4.5045, 0.0, 24.8644, 29.8397, 10.0254, 65.7197, 2.0601))
            ogsi = (tree + snag + deb + diam_div) / 4.0
            vec.append(1.0 if ogsi >= 50.94 else 0.0)
        return vec

    return vec
