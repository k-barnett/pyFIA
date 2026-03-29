"""
PNW mature weighted indices (R ``FUNCTION_mapMOG.R`` ~3944–4125, Table 19).
"""

from __future__ import annotations

from typing import Any, List


def _rng(a: int, b: int) -> frozenset[int]:
    return frozenset(range(a, b + 1))


_HARDWOODS_FT = _rng(700, 722) | _rng(900, 905) | _rng(960, 976) | _rng(920, 935)
_SOFTWOODS_FT = _rng(360, 369) | _rng(180, 185)


def pnw_mature_vector(
    forest_type: int,
    paz_group: str | None,
    inside_nwfp: bool,
    metrics: Any,
) -> List[float]:
    """
    Each matching R ``if`` appends one weighted maturity score (sum of 0–1 contributions).

    ``metrics`` must provide: ``ddi_score``, ``ht_quart``, ``badom``, ``snag_ba_tot``,
    ``qmd_dom``, ``ht_sd``, ``tpadom``.
    """

    vec: List[float] = []
    ft = int(forest_type)
    ddi = float(metrics.ddi_score)
    htq = float(metrics.ht_quart)
    bad = float(metrics.badom)
    snag = float(metrics.snag_ba_tot)
    qmd = float(metrics.qmd_dom)
    htsd = float(metrics.ht_sd)
    tpa = float(metrics.tpadom)

    if ft in _HARDWOODS_FT:
        t_ddi = 0.31 if ddi >= 47.5 else 0.0
        t_htq = 0.31 if htq >= 68.5 else 0.0
        t_bad = 0.21 if bad >= 150.0 else 0.0
        t_snag = 0.16 if snag >= 24.9 else 0.0
        vec.append(float(t_ddi + t_htq + t_bad + t_snag))

    if ft in _SOFTWOODS_FT:
        t_qmd = 0.54 if qmd >= 14.2 else 0.0
        t_bad = 0.46 if bad >= 30.9 else 0.0
        vec.append(float(t_qmd + t_bad))

    if paz_group == "douglas fir" and not inside_nwfp:
        t_qmd = 0.42 if qmd >= 11.1 else 0.0
        t_ddi = 0.38 if ddi >= 30.2 else 0.0
        t_bad = 0.21 if bad >= 60.1 else 0.0
        vec.append(float(t_qmd + t_ddi + t_bad))

    if ft in range(200, 204):
        t_qmd = 0.42 if qmd >= 11.1 else 0.0
        t_ddi = 0.38 if ddi >= 30.2 else 0.0
        t_bad = 0.21 if bad >= 60.1 else 0.0
        vec.append(float(t_qmd + t_ddi + t_bad))

    if paz_group == "douglas fir" and inside_nwfp:
        t_qmd = 0.45 if qmd >= 12.7 else 0.0
        t_ddi = 0.33 if ddi >= 32.6 else 0.0
        t_bad = 0.23 if bad >= 42.3 else 0.0
        vec.append(float(t_qmd + t_ddi + t_bad))

    if paz_group == "mountain hemlock" or ft in _rng(260, 271):
        t_qmd = 0.29 if qmd >= 13.1 else 0.0
        t_bad = 0.2 if bad >= 126.6 else 0.0
        t_htsd = 0.2 if htsd >= 58.5 else 0.0
        t_htq = 0.19 if htq >= 42.7 else 0.0
        t_tpa = 0.12 if tpa <= 77.4 else 0.0
        vec.append(float(t_qmd + t_bad + t_htsd + t_htq + t_tpa))

    if paz_group in {"ponderosa pine", "jeffrey pine - knobcone pine", "lodgepole pine"} or ft in _rng(220, 226):
        t_qmd = 0.34 if qmd >= 7.7 else 0.0
        t_ddi = 0.28 if ddi >= 15.3 else 0.0
        t_htsd = 0.22 if htsd >= 31.2 else 0.0
        t_tpa = 0.16 if tpa <= 31.7 else 0.0
        vec.append(float(t_qmd + t_ddi + t_htsd + t_tpa))

    if paz_group in {"port orford cedar", "redwood"}:
        t_ddi = 0.62 if ddi >= 44.4 else 0.0
        t_qmd = 0.38 if qmd >= 13.0 else 0.0
        vec.append(float(t_ddi + t_qmd))

    if paz_group in {"shasta red fir", "silver fir"}:
        t_qmd = 0.29 if qmd >= 17.1 else 0.0
        t_htsd = 0.2 if htsd >= 72.2 else 0.0  # R typo ``t.HTsd``; use live height SD
        t_bad = 0.19 if bad >= 161.6 else 0.0
        t_snag = 0.18 if snag >= 39.7 else 0.0
        t_tpa = 0.14 if tpa <= 53.1 else 0.0
        vec.append(float(t_qmd + t_htsd + t_bad + t_snag + t_tpa))

    if paz_group == "sitka spruce":
        t_qmd = 0.3 if qmd >= 24.3 else 0.0
        t_htsd = 0.22 if htsd >= 63.5 else 0.0
        t_bad = 0.2 if bad >= 184.6 else 0.0
        t_snag = 0.13 if snag >= 54.5 else 0.0
        t_tpa = 0.15 if tpa <= 37.6 else 0.0
        vec.append(float(t_qmd + t_htsd + t_bad + t_snag + t_tpa))

    if paz_group == "subalpine fir":
        t_ddi = 0.4 if ddi >= 27.8 else 0.0
        t_htq = 0.32 if htq >= 39.8 else 0.0
        t_htsd = 0.29 if htsd >= 41.3 else 0.0
        vec.append(float(t_ddi + t_htq + t_htsd))

    if ft in (265, 266):
        t_ddi = 0.42 if ddi >= 33.2 else 0.0
        t_qmd = 0.35 if qmd >= 8.8 else 0.0
        t_htsd = 0.23 if htsd >= 42.9 else 0.0
        vec.append(float(t_ddi + t_qmd + t_htsd))

    if paz_group == "tanoak":
        t_qmd = 0.29 if qmd >= 15.3 else 0.0
        t_ddi = 0.24 if ddi >= 56.0 else 0.0
        t_htq = 0.16 if htq >= 51.7 else 0.0
        t_tpa = 0.16 if tpa <= 55.9 else 0.0
        t_htsd = 0.15 if htsd >= 64.0 else 0.0
        vec.append(float(t_qmd + t_ddi + t_htq + t_tpa + t_htsd))

    if paz_group == "western hemlock":
        t_qmd = 0.33 if qmd >= 19.9 else 0.0
        t_bad = 0.2 if bad >= 156.2 else 0.0
        t_htsd = 0.17 if htsd >= 25.9 else 0.0
        t_snag = 0.17 if snag >= 63.2 else 0.0
        t_tpa = 0.14 if tpa <= 42.0 else 0.0
        vec.append(float(t_qmd + t_bad + t_htsd + t_snag + t_tpa))

    if paz_group == "white fir - grand fir" or ft in (267, 261):
        t_qmd = 0.33 if qmd >= 12.3 else 0.0
        t_ddi = 0.31 if ddi >= 40.1 else 0.0
        t_htsd = 0.2 if htsd >= 46.8 else 0.0
        t_snag = 0.16 if snag >= 8.6 else 0.0
        vec.append(float(t_qmd + t_ddi + t_htsd + t_snag))

    return vec
