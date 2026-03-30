"""
Southwest (Region 3) Table 9 / Table 19 diagnostics on ``mog_condition_scores`` output.
"""

from __future__ import annotations

import pandas as pd

# Table 9 diagnostics (no maturity); maturity is appended in ``mog_condition_scores``.
SOUTHWEST_OG_DIAGNOSTIC_KEYS: tuple[str, ...] = (
    "SW_ERU",
    "SW_TABLE9_RULE",
    "SW_TABLE9_FAMILY",
    "SW_REL_SDI",
    "SW_SDI_OFE",
    "SW_SDI_RMRS",
    "SW_QMD",
    "SW_OG_TABLE9",
    "SW_PASS_QMD18",
    "SW_PASS_SDI_GE_57",
    "SW_PASS_SDI_GE_56",
    "SW_PASS_SDI_GE_26",
    "SW_PASS_SDI_GE_29",
    "SW_PASS_SDI_GE_36",
    "SW_PASS_SDI_GE_20",
    "SW_APPLIES_QMD_TABLE9",
    "SW_APPLIES_SDI_TABLE9",
)

SW_MATURITY_SCORE_KEY = "SW_MATURITY_SCORE"

# Columns appended by :func:`fia_mog.estimators.mog_condition_scores` for southwest rows.
MOG_CONDITION_SW_COLUMNS: tuple[str, ...] = SOUTHWEST_OG_DIAGNOSTIC_KEYS + (SW_MATURITY_SCORE_KEY,)

_STAT_COLS: tuple[str, ...] = (
    "SW_OG_TABLE9",
    "SW_PASS_QMD18",
    "SW_PASS_SDI_GE_57",
    "SW_PASS_SDI_GE_56",
    "SW_PASS_SDI_GE_26",
    "SW_PASS_SDI_GE_29",
    "SW_PASS_SDI_GE_36",
    "SW_PASS_SDI_GE_20",
    "SW_APPLIES_QMD_TABLE9",
    "SW_APPLIES_SDI_TABLE9",
)


def summarize_southwest_og_by_eru(cond_df: pd.DataFrame) -> pd.DataFrame:
    """
    For each ``SW_ERU``, report condition counts, mean relative SDI / QMD, and for
    each 0/1 diagnostic column the **count** and **proportion** of conditions with 1.

    ``cond_df`` should be the output of :func:`fia_mog.estimators.mog_condition_scores`
    (optionally filtered to Arizona / New Mexico or a single administrative forest).
    Rows with empty ``SW_ERU`` are dropped.
    """

    if "SW_ERU" not in cond_df.columns:
        raise ValueError(
            "cond_df must include southwest columns from mog_condition_scores (missing SW_ERU)."
        )
    if "PLT_CN" not in cond_df.columns:
        raise ValueError("cond_df must include PLT_CN.")

    mask = cond_df["SW_ERU"].notna() & (cond_df["SW_ERU"].astype(str).str.len() > 0)
    sub = cond_df.loc[mask].copy()
    if sub.empty:
        return pd.DataFrame()

    rows: list[dict[str, object]] = []
    stat_cols = [c for c in _STAT_COLS if c in sub.columns]

    for eru, part in sub.groupby("SW_ERU", sort=True):
        n = int(len(part))
        row: dict[str, object] = {
            "SW_ERU": eru,
            "n_conditions": n,
            "mean_SW_REL_SDI": float(part["SW_REL_SDI"].mean()) if "SW_REL_SDI" in part.columns else float("nan"),
            "mean_SW_QMD": float(part["SW_QMD"].mean()) if "SW_QMD" in part.columns else float("nan"),
        }
        if SW_MATURITY_SCORE_KEY in part.columns:
            m = pd.to_numeric(part[SW_MATURITY_SCORE_KEY], errors="coerce")
            row["mean_SW_MATURITY_SCORE"] = float(m.mean()) if m.notna().any() else float("nan")

        for c in stat_cols:
            s = pd.to_numeric(part[c], errors="coerce").fillna(0).astype(int)
            short = c.replace("SW_", "").lower()
            row[f"n_{short}"] = int(s.sum())
            row[f"prop_{short}"] = float(s.mean()) if n else 0.0

        # How often each Table 9 rule id appears within this ERU (should be 1.0 for pure ERU groups).
        if "SW_TABLE9_RULE" in part.columns:
            vc = part["SW_TABLE9_RULE"].astype(str).value_counts(normalize=True)
            row["dominant_table9_rule"] = vc.index[0] if len(vc) else ""
            row["p_dominant_table9_rule"] = float(vc.iloc[0]) if len(vc) else 0.0

        rows.append(row)

    out = pd.DataFrame(rows)
    if not out.empty and "n_conditions" in out.columns:
        out = out.sort_values("n_conditions", ascending=False).reset_index(drop=True)
    return out


__all__ = [
    "MOG_CONDITION_SW_COLUMNS",
    "SOUTHWEST_OG_DIAGNOSTIC_KEYS",
    "SW_MATURITY_SCORE_KEY",
    "summarize_southwest_og_by_eru",
]
