"""
Montana MOG diagnostics: old-growth area + **NR_*** vegetation crosswalk gaps.

Run from the ``fia_py`` repo root, e.g. ``python scripts/test_mog_area_mt.py``,
or install the project in editable mode (``pip install -e .``).

Requires a Montana FIADB extract under **``tests/fiadb/MT/``** (gitignored; see README).

CSV exports (default **``diagnostics_export/``**):

- **``cond_mog_mt.csv``** — condition-level ``mog_condition_scores`` (includes ``NR_*``).
- **``mt_northern_veg_gap_summary.csv``** — condition/plot counts by gap type.
- **``mt_veg_not_in_crosswalk_by_code.csv``** — top ``NR_VEG_CODE`` where the code is built
  but **no habitat letters** remain after exact + prefix fallback (``NR_FLAG_VEG_NOT_IN_CROSSWALK``).

Use ``--fiadb-dir`` if your MT extract lives elsewhere; ``--out-dir`` for exports.

Import note: use ``fia.data_io`` / ``fia.clip`` instead of ``import fia`` so this
script does not pull in ``fia.__init__`` early (see ``scripts/test_mog_area.py``).

**Prove you are running this repo's code**

- Run ``python scripts/test_mog_area_mt.py --verify-mog-code`` (no FIADB): must print
  ``MOG_RADICAL_VERIFY_OK`` and absolute paths under your ``fia_py`` checkout.
- A normal run prints ``MOG_MT_ACTIVE_CODE_STAMP`` and writes ``mog_code_fingerprint.txt``
  next to the CSV exports (unless ``--quiet``).
"""

from __future__ import annotations

import argparse
import inspect
import sys
from pathlib import Path

import pandas as pd

_rp = Path(__file__).resolve().parent.parent
sys.path[:] = [str(_rp)] + [p for p in sys.path if p != str(_rp)]

from fia.clip import clip_fia
from fia.data_io import get_fia, read_fia
from fia_mog.estimators import (
    MOG_ESTIMATORS_REVISION,
    mog_condition_scores,
    old_growth_area,
    _mog_forest_type_from_cond_row,
)

# Must equal ``fia_mog.estimators.MOG_ESTIMATORS_REVISION`` (``--verify-mog-code`` checks this).
MOG_MT_SCRIPT_BUILD_STAMP = MOG_ESTIMATORS_REVISION
from fia_mog.northern import diagnostics as _nr_diag
from fia_mog.northern.core import _MT_LON_HEURISTIC_EASTERN_IF_GTE, northern_subregion
from fia_mog.northern.diagnostics import (
    filter_montana_mog_rows,
    nr_flag_true,
    summarize_mt_northern_veg_gaps,
)

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


def _print_active_code_banner(*, quiet: bool) -> None:
    if quiet:
        return
    import fia_mog
    import fia_mog.estimators as _est
    import fia_mog.engine as _eng

    print()
    print("=" * 72)
    print("MOG_MT_ACTIVE_CODE_STAMP:", MOG_MT_SCRIPT_BUILD_STAMP)
    print("python executable:", sys.executable)
    print("repo root (sys.path[0]):", _rp)
    print("fia_mog/__init__.py:  ", Path(fia_mog.__file__).resolve())
    print("fia_mog/estimators.py:", Path(_est.__file__).resolve())
    print("fia_mog/engine.py:    ", Path(_eng.__file__).resolve())
    print("=" * 72)
    print()


def run_verify_mog_code() -> None:
    """
    Radical self-test: no FIADB required. Fails loudly if you are not importing this repo's ``fia_mog``.
    """

    import fia_mog
    import fia_mog.estimators as est
    import fia_mog.engine as eng

    print("MOG_RADICAL_VERIFY_START", MOG_MT_SCRIPT_BUILD_STAMP)

    if not hasattr(est, "MOG_ESTIMATORS_REVISION"):
        raise SystemExit(
            "OLD fia_mog.estimators: no MOG_ESTIMATORS_REVISION — wrong install or PYTHONPATH"
        )
    if est.MOG_ESTIMATORS_REVISION != MOG_ESTIMATORS_REVISION:
        raise SystemExit(
            f"REVISION_MISMATCH imported={est.MOG_ESTIMATORS_REVISION!r} script={MOG_ESTIMATORS_REVISION!r}"
        )

    # FORTYPCD-only row (ponderosa) must resolve — proves FLDTYPCD-only filter is not the only path.
    row = pd.Series({"FLDTYPCD": pd.NA, "FORTYPCD": 221})
    if _mog_forest_type_from_cond_row(row) != 221:
        raise SystemExit("RADICAL_FAIL: _mog_forest_type_from_cond_row(FORTYPCD=221) != 221")

    # Northern: unmapped FLDTYPCD must not hide mappable FORTYPCD (common MT / FIADB pattern).
    row_mixed = pd.Series({"FLDTYPCD": 999, "FORTYPCD": 221})
    if _mog_forest_type_from_cond_row(row_mixed, region="northern") != 221:
        raise SystemExit(
            "RADICAL_FAIL: northern row FLDTYPCD=999 FORTYPCD=221 should use 221 for OG type"
        )

    # Longitude heuristic when divide flag is unknown.
    if _MT_LON_HEURISTIC_EASTERN_IF_GTE != -110.5:
        raise SystemExit(f"Unexpected divide lon constant: {_MT_LON_HEURISTIC_EASTERN_IF_GTE}")
    east = northern_subregion("MT", mt_east_of_continental_divide=None, plot_lon=-107.0)
    west = northern_subregion("MT", mt_east_of_continental_divide=None, plot_lon=-114.0)
    if east != "eastern Montana zone":
        raise SystemExit(f"RADICAL_FAIL: expected eastern zone for lon=-107, got {east!r}")
    if west != "western Montana zone":
        raise SystemExit(f"RADICAL_FAIL: expected western zone for lon=-114, got {west!r}")

    bundle_src = inspect.getsource(_nr_diag.compute_northern_habitat_og_bundle)
    for needle in (
        "infer_northern_og_type_from_species",
        "fallback_northern_habitat_letters",
        "plot_lon",
        "_extend_habitat_og_for_zone",
        "_NORTHERN_OG_ZONES",
    ):
        if needle not in bundle_src:
            raise SystemExit(f"RADICAL_FAIL: compute_northern_habitat_og_bundle source missing {needle!r}")

    print("MOG_RADICAL_VERIFY_OK — imported from:")
    print(" ", Path(fia_mog.__file__).resolve())
    print(" ", Path(est.__file__).resolve())
    print(" ", Path(eng.__file__).resolve())
    print(" ", Path(_nr_diag.__file__).resolve())
    print("If OG acres still look unchanged, the code is loading; thresholds / data may explain it.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="MOG old-growth area (MT) + northern vegetation crosswalk gap report."
    )
    parser.add_argument(
        "--verify-mog-code",
        action="store_true",
        help="Run import-path + assertion self-test (no FIADB), then exit 0.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress the MOG_MT_ACTIVE_CODE_STAMP banner (normal run only).",
    )
    parser.add_argument(
        "--fiadb-dir",
        type=Path,
        default=None,
        help="Directory with MT FIADB CSVs (default: <repo>/tests/fiadb/MT)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Directory for CSV exports (default: <repo>/diagnostics_export)",
    )
    args = parser.parse_args()

    if args.verify_mog_code:
        run_verify_mog_code()
        return

    _print_active_code_banner(quiet=args.quiet)

    fiadb_dir = args.fiadb_dir if args.fiadb_dir is not None else (_rp / "tests" / "fiadb" / "MT")
    fiadb_dir = fiadb_dir.resolve()
    if not fiadb_dir.is_dir():
        raise SystemExit(
            f"Montana FIADB directory not found: {fiadb_dir}\n"
            "Place an MT extract there or pass --fiadb-dir."
        )

    out_dir = args.out_dir if args.out_dir is not None else (_rp / "diagnostics_export")
    out_dir = out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    get_fia(
        states="MT",
        dir=str(fiadb_dir),
        tables=_TABLES_FOR_MOG,
        load=False,
    )
    db = read_fia(str(fiadb_dir), tables=_TABLES_FOR_MOG)
    db_current = clip_fia(db)

    cm = mog_condition_scores(
        db_current,
        states=["MT"],
        eval_typ="CURR",
        land_type="forest",
    )
    if cm.empty:
        print("No Montana forest conditions returned by mog_condition_scores.")
        return

    cond_path = out_dir / "cond_mog_mt.csv"
    cm.to_csv(cond_path, index=False)
    print(f"Wrote condition-level table: {cond_path}")

    fp_path = out_dir / "mog_code_fingerprint.txt"
    import fia_mog.estimators as _est_fp

    fp_path.write_text(
        f"MOG_MT_SCRIPT_BUILD_STAMP={MOG_MT_SCRIPT_BUILD_STAMP}\n"
        f"MOG_ESTIMATORS_REVISION={MOG_ESTIMATORS_REVISION}\n"
        f"fia_mog.estimators={Path(_est_fp.__file__).resolve()}\n",
        encoding="utf-8",
    )
    print(f"Wrote code fingerprint (proves which estimators.py ran): {fp_path}")

    if "STATECD" in cm.columns and "NR_FLAG_MT_SUBREGION_MISSING" in cm.columns:
        summary = summarize_mt_northern_veg_gaps(cm)
        sum_path = out_dir / "mt_northern_veg_gap_summary.csv"
        summary.to_csv(sum_path, index=False)
        print(f"Wrote Montana northern veg gap summary: {sum_path}")
        print("\n--- Montana northern vegetation / crosswalk gaps ---")
        print(summary.to_string(index=False))

        mt = filter_montana_mog_rows(cm)
        cx_path = out_dir / "mt_veg_not_in_crosswalk_by_code.csv"
        if "NR_VEG_CODE" in mt.columns and "NR_FLAG_VEG_NOT_IN_CROSSWALK" in mt.columns:
            cross = nr_flag_true(mt["NR_FLAG_VEG_NOT_IN_CROSSWALK"])
            ser = mt.loc[cross, "NR_VEG_CODE"].astype(str).replace("", "(empty)")
            if len(ser) == 0:
                vc = pd.DataFrame(columns=["NR_VEG_CODE", "n_conditions"])
            else:
                vc = (
                    ser.value_counts()
                    .rename_axis("NR_VEG_CODE")
                    .reset_index(name="n_conditions")
                )
        else:
            vc = pd.DataFrame(columns=["NR_VEG_CODE", "n_conditions"])
        vc.to_csv(cx_path, index=False)
        print(f"\nWrote veg-code counts (not in letter crosswalk): {cx_path}")

        if "NR_FLAG_NO_HABITAT_OG_BLOCKS" in mt.columns:
            no_blk = nr_flag_true(mt["NR_FLAG_NO_HABITAT_OG_BLOCKS"])
            br_cols = [
                c
                for c in (
                    "NR_SUBREGION",
                    "NR_OG_FOREST_TYPE",
                    "NR_N_HABITAT_LETTERS",
                    "FLDTYPCD",
                    "FORTYPCD",
                    "NR_VEG_CODE",
                )
                if c in mt.columns
            ]
            out_br = br_cols + ["n_conditions"]
            if br_cols and no_blk.any():
                br = (
                    mt.loc[no_blk, br_cols]
                    .fillna("")
                    .astype(str)
                    .groupby(br_cols, dropna=False)
                    .size()
                    .reset_index(name="n_conditions")
                    .sort_values("n_conditions", ascending=False)
                )
            else:
                br = pd.DataFrame(columns=out_br)
            br_path = out_dir / "mt_nr_no_blocks_by_sub_og_letters_type.csv"
            br.to_csv(br_path, index=False)
            print(f"Wrote NR_FLAG_NO_HABITAT_OG_BLOCKS breakdown: {br_path}")

        if "PLT_CN" in mt.columns:
            any_gap = nr_flag_true(mt["NR_FLAG_MT_SUBREGION_MISSING"])
            for col in (
                "NR_FLAG_VEG_CODE_MISSING",
                "NR_FLAG_VEG_NOT_IN_CROSSWALK",
                "NR_FLAG_NO_HABITAT_OG_BLOCKS",
            ):
                if col in mt.columns:
                    any_gap = any_gap | nr_flag_true(mt[col])
            n_p = int(mt.loc[any_gap, "PLT_CN"].nunique())
            n_c = int(any_gap.sum())
            n_pt = int(mt["PLT_CN"].nunique())
            print(
                f"\nPlots with ≥1 forest condition having any NR_* gap: {n_p} / {n_pt} "
                f"(conditions flagged: {n_c} / {len(mt)})."
            )
    else:
        print("NR_* columns missing from mog_condition_scores output; upgrade fia_mog.")

    res = old_growth_area(
        db_current,
        states=["MT"],
        grp_by=["ADFORCD"],
        eval_typ="CURR",
        land_type="forest",
        method="TI",
        totals=True,
        variance=True,
    )
    print("\n--- old_growth_area by ADFORCD (MT) ---")
    print(res.df)
    print(res.cond_mog.head())
    print(res.forest_area_df)


if __name__ == "__main__":
    main()
