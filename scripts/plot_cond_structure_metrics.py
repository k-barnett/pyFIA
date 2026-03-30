"""
Build per-condition (PLT_CN, CONDID) structure metrics: TPA/BAA strata from
``tpa(..., tree_list=True)`` plus QMD and height summaries from ``TREE``.

Mirrors common rFIA-style dplyr workflows; point ``FIADB_DIR`` at a folder of
FIA CSVs (e.g. from ``get_fia``), then run::

    python scripts/plot_cond_structure_metrics.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Project root (parent of ``scripts/``)
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from fia import clip_fia, read_fia, tpa

# ---------------------------------------------------------------------------
# User: set this to the directory containing FIA table CSVs (e.g. MT_PLOT.csv).
# ---------------------------------------------------------------------------
FIADB_DIR = Path("/Users/kevinbarnett/Documents/FIADB").expanduser()

TABLES = [
    "COND",
    "COND_DWM_CALC",
    "INVASIVE_SUBPLOT_SPP",
    "PLOT",
    "POP_ESTN_UNIT",
    "POP_EVAL",
    "POP_EVAL_GRP",
    "POP_EVAL_TYP",
    "POP_PLOT_STRATUM_ASSGN",
    "POP_STRATUM",
    "SUBPLOT",
    "TREE",
    "TREE_GRM_COMPONENT",
    "TREE_GRM_MIDPT",
    "TREE_GRM_BEGIN",
    "SUBP_COND_CHNG_MTRX",
    "SEEDLING",
    "SURVEY",
    "SUBP_COND",
    "P2VEG_SUBP_STRUCTURE",
    "PLOTGEOM",
]


def _qmd(dia: pd.Series, tpa_unadj: pd.Series, mask: pd.Series) -> float:
    dia = pd.to_numeric(dia, errors="coerce")
    w = pd.to_numeric(tpa_unadj, errors="coerce")
    m = mask.fillna(False) & dia.notna() & w.notna() & (w > 0)
    if not m.any():
        return float("nan")
    num = ((dia[m] ** 2) * w[m]).sum()
    den = w[m].sum()
    if den == 0:
        return float("nan")
    return float(np.sqrt(num / den))


def qmd_height_by_condition(tree: pd.DataFrame) -> pd.DataFrame:
    """QMD (live, live >5\") and height stats by PLT_CN, CONDID (from raw TREE)."""

    def _summarise_condition(g: pd.DataFrame) -> pd.Series:
        ht = pd.to_numeric(g["HT"], errors="coerce")
        dia = pd.to_numeric(g["DIA"], errors="coerce")
        tpa_u = pd.to_numeric(g["TPA_UNADJ"], errors="coerce")
        live = g["STATUSCD"].eq(1)
        dead = g["STATUSCD"].eq(2)
        live_gt5 = live & dia.gt(5)

        return pd.Series(
            {
                "QMD_live": _qmd(dia, tpa_u, live),
                "QMD_live_gt5": _qmd(dia, tpa_u, live_gt5),
                "MAXHT_all": ht.max(skipna=True),
                "MAXHT_live": ht.where(live).max(skipna=True),
                "MAXHT_dead": ht.where(dead).max(skipna=True),
                "AVGHT_all": ht.mean(skipna=True),
                "AVGHT_live": ht.where(live).mean(skipna=True),
                "AVGHT_dead": ht.where(dead).mean(skipna=True),
            }
        )

    return (
        tree.groupby(["PLT_CN", "CONDID"], sort=False)
        .apply(_summarise_condition)
        .reset_index()
    )


def tpa_baa_by_condition(
    db_clipped,
    tree: pd.DataFrame,
    *,
    land_type: str = "forest",
    tree_type: str = "all",
) -> pd.DataFrame:
    """
    TPA/BAA sums by PLT_CN, CONDID with live/dead and diameter cuts (>5\", >17\"),
    using ``tpa(..., tree_list=True)`` joined back to TREE for STATUSCD/DIA.
    """
    tl = tpa(
        db_clipped,
        land_type=land_type,
        tree_type=tree_type,
        tree_list=True,
    )

    tree_keys = tree[["PLT_CN", "CONDID", "SUBP", "TREE", "STATUSCD", "DIA"]].copy()
    m = tl.merge(
        tree_keys,
        on=["PLT_CN", "CONDID", "SUBP", "TREE"],
        how="inner",
    )

    d = pd.to_numeric(m["DIA"], errors="coerce")
    live = m["STATUSCD"].eq(1)
    dead = m["STATUSCD"].eq(2)
    gt5 = d.gt(5)
    gt17 = d.gt(17)

    m = m.assign(
        tp_gt5=m["TPA"].where(gt5),
        tp_live=m["TPA"].where(live),
        tp_live_gt5=m["TPA"].where(live & gt5),
        tp_dead=m["TPA"].where(dead),
        tp_dead_gt5=m["TPA"].where(dead & gt5),
        tp_gt17=m["TPA"].where(gt17),
        tp_live_gt17=m["TPA"].where(live & gt17),
        tp_dead_gt17=m["TPA"].where(dead & gt17),
        ba_gt5=m["BAA"].where(gt5),
        ba_live=m["BAA"].where(live),
        ba_live_gt5=m["BAA"].where(live & gt5),
        ba_dead=m["BAA"].where(dead),
        ba_dead_gt5=m["BAA"].where(dead & gt5),
        ba_gt17=m["BAA"].where(gt17),
        ba_live_gt17=m["BAA"].where(live & gt17),
        ba_dead_gt17=m["BAA"].where(dead & gt17),
    )

    return (
        m.groupby(["PLT_CN", "CONDID"], sort=False)
        .agg(
            TPA_all=("TPA", "sum"),
            TPA_all_gt5=("tp_gt5", "sum"),
            TPA_live=("tp_live", "sum"),
            TPA_live_gt5=("tp_live_gt5", "sum"),
            TPA_dead=("tp_dead", "sum"),
            TPA_dead_gt5=("tp_dead_gt5", "sum"),
            TPA_all_gt17=("tp_gt17", "sum"),
            TPA_live_gt17=("tp_live_gt17", "sum"),
            TPA_dead_gt17=("tp_dead_gt17", "sum"),
            BAA_all=("BAA", "sum"),
            BAA_all_gt5=("ba_gt5", "sum"),
            BAA_live=("ba_live", "sum"),
            BAA_live_gt5=("ba_live_gt5", "sum"),
            BAA_dead=("ba_dead", "sum"),
            BAA_dead_gt5=("ba_dead_gt5", "sum"),
            BAA_all_gt17=("ba_gt17", "sum"),
            BAA_live_gt17=("ba_live_gt17", "sum"),
            BAA_dead_gt17=("ba_dead_gt17", "sum"),
        )
        .reset_index()
    )


def build_condition_structure_table(
    db_clipped,
    *,
    land_type: str = "forest",
    tree_type: str = "all",
) -> pd.DataFrame:
    """
    Full outer merge of TPA/BAA block and QMD/height block on PLT_CN, CONDID.

    Re-runs ``tpa`` with the given land/tree domain; QMD/height always uses the
    raw ``TREE`` table in ``db_clipped`` (all rows in that table).
    """
    tree = db_clipped["TREE"]
    tpa_block = tpa_baa_by_condition(
        db_clipped, tree, land_type=land_type, tree_type=tree_type
    )
    qmd_block = qmd_height_by_condition(tree)
    return tpa_block.merge(qmd_block, on=["PLT_CN", "CONDID"], how="outer")


def main() -> None:
    pd.set_option("display.float_format", "{:.4f}".format)

    if not FIADB_DIR.is_dir():
        print(f"Set FIADB_DIR to a folder of FIA CSVs; not found: {FIADB_DIR}", file=sys.stderr)
        sys.exit(1)

    db = read_fia(str(FIADB_DIR), tables=TABLES)
    db_current = clip_fia(db, most_recent=True)

    tree = db_current["TREE"]
    cond = db_current["COND"]
    plot_df = db_current["PLOT"]
    plot_geom = db_current["PLOTGEOM"]

    cond_metrics = build_condition_structure_table(
        db_current,
        land_type="forest",
        tree_type="all",
    )

    print(f"Rows (PLT_CN, CONDID): {len(cond_metrics)}")
    print(cond_metrics.head())


if __name__ == "__main__":
    main()
