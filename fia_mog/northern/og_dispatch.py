"""
Northern Region old-growth 0/1 rule blocks (Green et al. 1992).

Mirrors the R inner loops for Idaho, western Montana, and eastern Montana zones.
Tree densities use nrow(DIA >= x) / condition.area (no TPA), matching the R script.
"""

from __future__ import annotations

import pandas as pd

_IDAHO_CEDAR_SPCD = frozenset({41, 42, 43, 67, 68, 81, 241, 242, 978})
_IDAHO_SPEC_SPCD = frozenset({201, 202, 17, 73, 263, 114, 119, 122})
_IDAHO_ALL_EX = _IDAHO_CEDAR_SPCD | _IDAHO_SPEC_SPCD


def _n_large_per_acre(trees: pd.DataFrame, dia_min: float, area: float) -> float:
    if area <= 0:
        return 0.0
    dia = pd.to_numeric(trees["DIA"], errors="coerce")
    return float((dia >= dia_min).sum() / area)


def _idaho_group7_density(trees: pd.DataFrame, area: float) -> float:
    if area <= 0:
        return 0.0
    dia = pd.to_numeric(trees["DIA"], errors="coerce")
    sp = pd.to_numeric(trees["SPCD"], errors="coerce")
    n_ced = int(((sp.isin(_IDAHO_CEDAR_SPCD)) & (dia >= 25)).sum())
    n_spec = int(((sp.isin(_IDAHO_SPEC_SPCD)) & (dia >= 21)).sum())
    n_oth = int(((~sp.isin(_IDAHO_ALL_EX)) & (dia >= 17)).sum())
    return float((n_ced + n_spec + n_oth) / area)


def northern_idaho_og_vector(
    letters: list[str],
    og_type: str | None,
    *,
    raw_stand_age: float,
    basal_area_per_acre: float,
    condition_area_acres: float,
    trees: pd.DataFrame,
) -> list[float]:
    vec: list[float] = []
    if not og_type or condition_area_acres <= 0:
        return vec
    area = condition_area_acres
    age = float(raw_stand_age)
    basal = float(basal_area_per_acre)

    for v in letters:
        # group 1
        if v in {"A", "B"} and og_type in {"PP", "DF", "L"}:
            t_age = 1.0 if age >= 150 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 21, area) >= 8 else 0.0
            t_bas = 1.0 if basal >= 40 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        # group 2
        if v in {"B", "C", "D", "E", "G", "H", "I", "J", "K"} and og_type == "LP":
            t_age = 1.0 if age >= 120 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 13, area) >= 10 else 0.0
            t_bas = 1.0 if basal >= 60 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        # group 3 (Y + yew)
        if v in {"C", "C1", "G1"} and og_type == "Y":
            t_age = 1.0 if age >= 150 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 21, area) >= 3 else 0.0
            t_bas = 1.0 if basal >= 80 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        # group 4A
        if v in {"C", "C1", "D", "E"} and og_type in {"DF", "GF", "SAF", "WP", "PP"}:
            t_age = 1.0 if age >= 150 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 21, area) >= 10 else 0.0
            t_bas = 1.0 if basal >= 80 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        # group 4B
        if v in {"F", "G", "G1", "H", "I"} and og_type in {"DF", "GF", "WH", "WP", "PP"}:
            t_age = 1.0 if age >= 150 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 21, area) >= 10 else 0.0
            if v in {"F", "G", "G1"}:
                t_bas = 1.0 if basal >= 120 else 0.0
            else:
                t_bas = 1.0 if basal >= 80 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        # group 5
        if v in {"F", "G", "G1", "H", "I"} and og_type in {"SAF", "MAF"}:
            t_age = 1.0 if age >= 150 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 17, area) >= 10 else 0.0
            t_bas = 1.0 if basal >= 80 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        # group 6 (WBP)
        if v in {"I", "J", "K"} and og_type == "WBP":
            t_age = 1.0 if age >= 150 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 13, area) >= 5 else 0.0
            if v in {"I", "J"}:
                t_bas = 1.0 if basal >= 60 else 0.0
            else:
                t_bas = 1.0 if basal >= 40 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        # group 7 (C / cedars)
        if v in {"F", "G", "G1"} and og_type == "C":
            t_age = 1.0 if age >= 150 else 0.0
            t_dens = 1.0 if _idaho_group7_density(trees, area) >= 10 else 0.0
            t_bas = 1.0 if basal >= 120 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        # group 8
        if v == "J" and og_type in {"DF", "L", "SAF", "MAF", "WP"}:
            t_age = 1.0 if age >= 150 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 17, area) >= 10 else 0.0
            t_bas = 1.0 if basal >= 60 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        # group 9
        if v == "K" and og_type in {"SAF", "MAF"}:
            t_age = 1.0 if age >= 150 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 13, area) >= 5 else 0.0
            t_bas = 1.0 if basal >= 40 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

    return vec


def northern_west_mt_og_vector(
    letters: list[str],
    og_type: str | None,
    *,
    raw_stand_age: float,
    basal_area_per_acre: float,
    condition_area_acres: float,
    trees: pd.DataFrame,
) -> list[float]:
    vec: list[float] = []
    if not og_type or condition_area_acres <= 0:
        return vec
    area = condition_area_acres
    age = float(raw_stand_age)
    basal = float(basal_area_per_acre)

    for v in letters:
        # group 1
        if v in {"A", "B"} and og_type in {"PP", "DF", "L", "GF", "LP"}:
            t_age = 1.0 if age >= 170 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 21, area) >= 8 else 0.0
            t_bas = 1.0 if basal >= 60 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        # group 2
        if v == "C" and og_type in {"DF", "L", "PP", "SAF", "GF"}:
            t_age = 1.0 if age >= 170 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 21, area) >= 8 else 0.0
            t_bas = 1.0 if basal >= 80 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        # group 3 (LP)
        if v in {"C", "D", "E", "F", "G", "H"} and og_type == "LP":
            t_age = 1.0 if age >= 140 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 13, area) >= 10 else 0.0
            if v == "E":
                t_bas = 1.0 if basal >= 60 else 0.0
            elif v in {"C", "H"}:
                t_bas = 1.0 if basal >= 70 else 0.0
            else:
                t_bas = 1.0 if basal >= 80 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        # group 4
        if v in {"D", "E", "F"} and og_type in {"SAF", "DF", "GF", "C", "L", "MAF", "PP", "WP", "WH", "WSL"}:
            t_age = 1.0 if age >= 180 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 21, area) >= 10 else 0.0
            t_bas = 1.0 if basal >= 80 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        # group 5
        if v in {"G", "H"} and og_type in {"SAF", "DF", "GF", "L", "MAF", "PP", "WP", "WSL"}:
            t_age = 1.0 if age >= 180 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 17, area) >= 10 else 0.0
            if v == "H" and og_type == "SAF":
                t_bas = 1.0 if basal >= 70 else 0.0
            else:
                t_bas = 1.0 if basal >= 80 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        # group 6
        if v == "I" and og_type in {"SAF", "WSL", "DF", "L"}:
            t_age = 1.0 if age >= 180 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 13, area) >= 10 else 0.0
            t_bas = 1.0 if basal >= 60 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        # group 7 (LP + I)
        if v == "I" and og_type == "LP":
            t_age = 1.0 if age >= 140 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 9, area) >= 30 else 0.0
            t_bas = 1.0 if basal >= 70 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        # group 8
        if v == "J" and og_type in {"SAF", "WSL"}:
            t_age = 1.0 if age >= 180 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 13, area) >= 20 else 0.0
            t_bas = 1.0 if basal >= 80 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

    return vec


def northern_east_mt_og_vector(
    letters: list[str],
    og_type: str | None,
    *,
    raw_stand_age: float,
    basal_area_per_acre: float,
    condition_area_acres: float,
    trees: pd.DataFrame,
) -> list[float]:
    vec: list[float] = []
    if not og_type or condition_area_acres <= 0:
        return vec
    area = condition_area_acres
    age = float(raw_stand_age)
    basal = float(basal_area_per_acre)

    for v in letters:
        if v == "A" and og_type == "DF":
            t_age = 1.0 if age >= 200 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 17, area) >= 4 else 0.0
            t_bas = 1.0 if basal >= 60 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        if v in {"B", "C", "D", "E", "F", "H"} and og_type == "DF":
            t_age = 1.0 if age >= 200 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 19, area) >= 5 else 0.0
            t_bas = 1.0 if basal >= 60 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        if v == "G" and og_type == "DF":
            t_age = 1.0 if age >= 180 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 17, area) >= 10 else 0.0
            t_bas = 1.0 if basal >= 80 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        # R uses "K" in set but no habitat letter K exists in eastern MT lists
        if v in {"A", "B", "C", "K"} and og_type == "PP":
            t_age = 1.0 if age >= 180 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 17, area) >= 4 else 0.0
            t_bas = 1.0 if basal >= 40 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        if v in {"A", "B"} and og_type == "PF":
            t_age = 1.0 if age >= 120 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 9, area) >= 6 else 0.0
            t_bas = 1.0 if basal >= 50 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        if v in {"A", "B", "C", "D", "E", "F", "G", "H", "I"} and og_type == "LP":
            t_age = 1.0 if age >= 150 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 10, area) >= 12 else 0.0
            t_bas = 1.0 if basal >= 50 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        if v == "C" and og_type == "SAF":
            t_age = 1.0 if age >= 160 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 17, area) >= 12 else 0.0
            t_bas = 1.0 if basal >= 80 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        if v in {"D", "E"} and og_type == "SAF":
            t_age = 1.0 if age >= 160 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 17, area) >= 7 else 0.0
            t_bas = 1.0 if basal >= 80 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        if v in {"F", "G", "H", "I"} and og_type == "SAF":
            t_age = 1.0 if age >= 160 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 13, area) >= 10 else 0.0
            t_bas = 1.0 if basal >= 60 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        if v == "J" and og_type == "SAF":
            t_age = 1.0 if age >= 135 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 13, area) >= 8 else 0.0
            t_bas = 1.0 if basal >= 40 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        if v in {"D", "E", "F", "G", "H", "I"} and og_type == "WBP":
            t_age = 1.0 if age >= 150 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 13, area) >= 11 else 0.0
            t_bas = 1.0 if basal >= 60 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

        if v == "J" and og_type == "WBP":
            t_age = 1.0 if age >= 135 else 0.0
            t_dens = 1.0 if _n_large_per_acre(trees, 13, area) >= 7 else 0.0
            t_bas = 1.0 if basal >= 40 else 0.0
            vec.append(1.0 if (t_age + t_dens + t_bas) == 3 else 0.0)

    return vec
