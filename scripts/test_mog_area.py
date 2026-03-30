"""
Run from the ``fia_py`` repo root, e.g. ``python scripts/test_mog_area.py``,
or install the project in editable mode (``pip install -e .``).

CSV exports (default directory ``diagnostics_export/`` under the repo root):

- ``cond_mog_az.csv`` — full condition-level table (includes ``SW_*`` columns for AZ).
- ``sw_eru_summary_az.csv`` — southwest Table 9 / SDI vs QMD summary by ERU (all AZ).
- ``cond_mog_adforcd307.csv`` and ``sw_eru_summary_adforcd307.csv`` — Kaibab subset, if present.

Use ``--out-dir /path/to/folder`` to choose another directory.

Import note: use ``fia.data_io`` / ``fia.clip`` instead of ``import fia`` so this
script does not pull in ``fia.__init__`` (which binds ``old_growth_area`` early).
In notebooks, restart the kernel after editing ``fia_mog`` or run
``import importlib, fia_mog.estimators as e; importlib.reload(e)``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Repo root must win over any other ``fia`` / ``fia_mog`` on sys.path.
_rp = Path(__file__).resolve().parent.parent
sys.path[:] = [str(_rp)] + [p for p in sys.path if p != str(_rp)]

from fia.clip import clip_fia
from fia.data_io import get_fia, read_fia
from fia_mog.estimators import old_growth_area
from fia_mog.southwest.diagnostics import summarize_southwest_og_by_eru

# ``get_fia`` requires an explicit ``tables=`` list in the Python port.
_TABLES_FOR_MOG = [
    "COND",
    "PLOT",
    "TREE",
    "POP_ESTN_UNIT",
    "POP_EVAL",
    "POP_EVAL_GRP",
    "POP_EVAL_TYP",
    "POP_PLOT_STRATUM_ASSGN",
    "POP_STRATUM",
    "PLOTGEOM",
    "P2VEG_SUBPLOT_SPP",
    "DWM_COARSE_WOODY_DEBRIS",
]


def main() -> None:
    parser = argparse.ArgumentParser(description="MOG old-growth area test (AZ) + CSV diagnostics.")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Directory for CSV exports (default: <repo>/diagnostics_export)",
    )
    args = parser.parse_args()
    out_dir = args.out_dir if args.out_dir is not None else (_rp / "diagnostics_export")
    out_dir = out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    get_fia(
        states="AZ",
        dir="./tests/fiadb/AZ",
        tables=_TABLES_FOR_MOG,
        load=False,
    )
    db = read_fia("./tests/fiadb/AZ", tables=_TABLES_FOR_MOG)
    db_current = clip_fia(db)

    res = old_growth_area(
        db_current,
        states=["AZ"],
        grp_by=["ADFORCD"],
        eval_typ="CURR",
        method="TI",
        totals=True,
        variance=True,
    )

    print(res.df)

    cm = res.cond_mog
    cond_path = out_dir / "cond_mog_az.csv"
    cm.to_csv(cond_path, index=False)
    print(f"Wrote condition-level diagnostics: {cond_path}")

    if "SW_ERU" in cm.columns:
        sw_mask = cm["SW_ERU"].notna() & (cm["SW_ERU"].astype(str).str.len() > 0)
        sw_cm = cm.loc[sw_mask]
        if len(sw_cm):
            sum_az = summarize_southwest_og_by_eru(sw_cm)
            sum_path = out_dir / "sw_eru_summary_az.csv"
            sum_az.to_csv(sum_path, index=False)
            print(f"Wrote southwest ERU summary (AZ): {sum_path}")

    if "ADFORCD" in cm.columns and {"MOG_SCORE", "OG_FLAG"}.issubset(cm.columns):
        k = cm[cm["ADFORCD"] == 307]
        if len(k):
            n_mature_only = int(((k["MOG_SCORE"] >= 1.0 - 1e-9) & (k["OG_FLAG"] < 0.5)).sum())
            print(
                f"Kaibab (307): conditions={len(k)}, "
                f"mean OG_FLAG={k['OG_FLAG'].mean():.4f}, "
                f"rows with MOG_SCORE≈1 but OG_FLAG=0 (mature-not-OG)={n_mature_only}"
            )
            k_path = out_dir / "cond_mog_adforcd307.csv"
            k.to_csv(k_path, index=False)
            print(f"Wrote Kaibab condition-level diagnostics: {k_path}")
            if "SW_ERU" in k.columns:
                print("\nKaibab (307) — Table 9 / SDI vs QMD by ERU:")
                sk = summarize_southwest_og_by_eru(k)
                print(sk.to_string(index=False))
                sk_path = out_dir / "sw_eru_summary_adforcd307.csv"
                sk.to_csv(sk_path, index=False)
                print(f"Wrote Kaibab ERU summary: {sk_path}")

    print(res.cond_mog.head())
    print(res.forest_area_df)


if __name__ == "__main__":
    main()
