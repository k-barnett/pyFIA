"""
Southwest MOG vector: **Table 9** old growth, **Table 19** maturity.

Old-growth minimum criteria for Southwestern Region (Region 3) follow
*Mature and Old-Growth Forests: Definition, Identification, and Initial Inventory*
(FS-1215a, April 2024), **Table 9**—ecological response units and minimum
**% SDI from trees ≥ 18 in** compared to **COND.SDI_RMRS**, or **minimum QMD**
for **live** trees **≥ 10 in** where SDI is n/a
(`mature-and-old-growth-forests-tech.pdf`). **Percent SDI** matches the OFE / rFIA-style
recipe: ``SDI_OFE = Σ(TPA_UNADJ × (DIA/10)^1.6)`` on **live** trees with **DIA ≥ 18** (inches),
then ``100 × SDI_OFE / SDI_RMRS`` (``SDI_RMRS`` from the condition table). Table 9 **QMD** uses conventional basal
area ``π (DIA/24)²`` ft² per tree (DBH in inches) and ``sqrt((BA_acres / ΣTPA) / 0.005454)``
with ``BA_acres = Σ(TPA_UNADJ × BA_tree)``, matching typical FIA / OFE-style helpers—not
``FUNCTION_mapMOG.R``’s linear ``π (DIA/24) × 2`` term.

Habitat → ERU mapping matches **Table 10** in that document (ported from R
``FUNCTION_mapMOG.R`` / ``southwest_eru``).
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, List

import numpy as np
import pandas as pd

from .core import southwest_eru
from .diagnostics import SOUTHWEST_OG_DIAGNOSTIC_KEYS
from ..engine import _basal_term_r, _weighted_index

if TYPE_CHECKING:
    from ..engine import ConditionContext, TreeMetrics

# Table 9 ERUs that use live-tree QMD (≥10 in) ≥ 18 in, not relative SDI.
_QMD_TABLE9_ERUS: frozenset[str] = frozenset(
    {
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
    }
)

def _basal_area_sq_ft_standard(dia_inches: pd.Series | np.ndarray) -> np.ndarray:
    """
    Tree basal area (ft²) from DBH in inches: ``π × (DIA/24)²``.

    Same as ``pi*(DIA/12/2)^2`` in R; contrasts with :func:`fia_mog.engine._basal_term_r`
    (mapMOG / maturity metrics).
    """

    d = np.asarray(pd.to_numeric(dia_inches, errors="coerce"), dtype=float)
    return np.pi * (d / 24.0) ** 2


def _resolve_southwest_eru(ctx: "ConditionContext") -> str | None:
    if ctx.condition_habtypcd1 is not None and not pd.isna(ctx.condition_habtypcd1):
        return southwest_eru(ctx.condition_habtypcd1)
    if isinstance(ctx.ecosubcd, str):
        return ctx.ecosubcd
    return None


# Exponent in ``(DIA/10)^exp`` for condition-level SDI_OFE (live, DIA ≥ 18 in).
# Matches ``rFIA-ofe-region3.R`` ``variable_sdi_calc`` (not mapMOG Zeide 1.6064).
_SOUTHWEST_SDI_OFE_EXPONENT = 1.6


def _southwest_sdi_ofe_live_18(trees: pd.DataFrame) -> float:
    """
    ``Σ TPA_UNADJ × (DIA/10)^1.6`` for **live** (``STATUSCD == 1``) trees with ``DIA >= 18``.
    """

    t = trees.copy()
    t["DIA"] = pd.to_numeric(t.get("DIA"), errors="coerce")
    mask = t["DIA"] >= 18
    if "STATUSCD" in t.columns:
        sc = pd.to_numeric(t["STATUSCD"], errors="coerce")
        mask = mask & (sc == 1)
    sub = t.loc[mask]
    if sub.empty or "TPA_UNADJ" not in sub.columns:
        return 0.0
    dia = sub["DIA"].to_numpy(dtype=float)
    tpa = pd.to_numeric(sub["TPA_UNADJ"], errors="coerce").fillna(0.0).to_numpy(dtype=float)
    terms = (dia / 10.0) ** _SOUTHWEST_SDI_OFE_EXPONENT
    total = float((terms * tpa).sum())
    return total if math.isfinite(total) else 0.0


def _southwest_relative_sdi_and_qmd(ctx: "ConditionContext") -> tuple[float, float, float]:
    """
    Table 9 **percent SDI** (``100 × SDI_OFE / SDI_RMRS``), Southwest QMD (live ≥10 in),
    and **SDI_OFE** (same numerator).

    If ``condition_sdi_rmrs`` is missing or not positive, **percent SDI** is ``nan`` (no OG pass
    on SDI rules).
    """

    sdi_ofe = _southwest_sdi_ofe_live_18(ctx.trees)
    denom = ctx.condition_sdi_rmrs
    if denom is None or not math.isfinite(denom) or denom <= 0:
        relative_sdi = float("nan")
    else:
        pct = 100.0 * sdi_ofe / denom
        relative_sdi = float(pct) if math.isfinite(pct) else float("nan")

    sw_qmd = 0.0
    sw_live = ctx.trees.copy()
    sw_live["DIA"] = pd.to_numeric(sw_live.get("DIA"), errors="coerce")
    sw_live = sw_live.loc[(sw_live.get("STATUSCD") == 1) & (sw_live["DIA"] >= 10), :]
    if len(sw_live) > 0 and "TPA_UNADJ" in sw_live.columns:
        tpa_u = pd.to_numeric(sw_live["TPA_UNADJ"], errors="coerce").fillna(0.0)
        # Align with OFE-style R: BA_tree = π*(DIA/24)² ft²; BA_ac = Σ(TPA_UNADJ×BA_tree);
        # QMD = sqrt((BA_ac / ΣTPA_UNADJ) / 0.005454). No rounding (R `round(...,1)` is for display).
        sw_tpa = float(tpa_u.sum())
        ba_tree = _basal_area_sq_ft_standard(sw_live["DIA"].to_numpy())
        sw_ba = float((ba_tree * tpa_u.to_numpy()).sum())
        if sw_tpa > 0:
            sw_qmd_val = math.sqrt((sw_ba / sw_tpa) / 0.005454)
            sw_qmd = float(sw_qmd_val) if np.isfinite(sw_qmd_val) else 0.0

    return relative_sdi, float(sw_qmd), sdi_ofe


def _southwest_table9_family_rule_og(
    relative_sdi: float, sw_qmd: float, eru: str | None
) -> tuple[str, str, int]:
    """
    Return ``(family, rule_id, local_og)`` matching :func:`southwest_mog_vector` Table 9 logic.

    ``family`` is ``qmd``, ``sdi``, ``none`` (ERU not in Table 9), or ``missing_eru``.
    """

    if not eru:
        return ("missing_eru", "no_eru", 0)
    if eru in _QMD_TABLE9_ERUS:
        og = 1 if sw_qmd >= 18 else 0
        return ("qmd", "qmd_ge_18", og)
    if eru == "mixed conifer - frequent fire":
        return ("sdi", "sdi_ge_56_mcff", 1 if relative_sdi >= 56 else 0)
    if eru in {"ponderosa pine forest", "ponderosa pine/willow"}:
        return ("sdi", "sdi_ge_57_pp", 1 if relative_sdi >= 57 else 0)
    if eru == "ponderosa pine - evergreen oak":
        return ("sdi", "sdi_ge_26_ppo", 1 if relative_sdi >= 26 else 0)
    if eru == "pinyon juniper grass":
        return ("sdi", "sdi_ge_29", 1 if relative_sdi >= 29 else 0)
    if eru in {"juniper grass", "semi-desert grassland"}:
        return ("sdi", "sdi_ge_36", 1 if relative_sdi >= 36 else 0)
    if eru in {"madrean pinyon-oak", "madrean encinal woodland"}:
        return ("sdi", "sdi_ge_20", 1 if relative_sdi >= 20 else 0)
    return ("none", "no_table9", 0)


def southwest_og_diagnostic_row(ctx: "ConditionContext") -> dict[str, object]:
    """
    Per-condition Table 9 diagnostics for Region 3 (Southwest).

    Returns 0/1 flags for each published SDI threshold (independent of ERU) plus
    which Table 9 rule applies for this ERU. Use with :func:`summarize_southwest_og_by_eru`.
    """

    eru = _resolve_southwest_eru(ctx)
    rsdi, sqmd, sdi_ofe = _southwest_relative_sdi_and_qmd(ctx)
    family, rule, og = _southwest_table9_family_rule_og(rsdi, sqmd, eru)
    rmrs_out: float | object = ctx.condition_sdi_rmrs if ctx.condition_sdi_rmrs is not None else pd.NA
    out = {
        "SW_ERU": eru if eru else "",
        "SW_TABLE9_RULE": rule,
        "SW_TABLE9_FAMILY": family,
        "SW_REL_SDI": rsdi,
        "SW_SDI_OFE": float(sdi_ofe),
        "SW_SDI_RMRS": rmrs_out,
        "SW_QMD": sqmd,
        "SW_OG_TABLE9": int(og),
        "SW_PASS_QMD18": 1 if sqmd >= 18 else 0,
        "SW_PASS_SDI_GE_57": 1 if rsdi >= 57 else 0,
        "SW_PASS_SDI_GE_56": 1 if rsdi >= 56 else 0,
        "SW_PASS_SDI_GE_26": 1 if rsdi >= 26 else 0,
        "SW_PASS_SDI_GE_29": 1 if rsdi >= 29 else 0,
        "SW_PASS_SDI_GE_36": 1 if rsdi >= 36 else 0,
        "SW_PASS_SDI_GE_20": 1 if rsdi >= 20 else 0,
        "SW_APPLIES_QMD_TABLE9": 1 if family == "qmd" else 0,
        "SW_APPLIES_SDI_TABLE9": 1 if family == "sdi" else 0,
    }
    if set(out) != set(SOUTHWEST_OG_DIAGNOSTIC_KEYS):
        raise RuntimeError("southwest_og_diagnostic_row keys drifted from SOUTHWEST_OG_DIAGNOSTIC_KEYS")
    return out


def southwest_mog_vector(ctx: "ConditionContext", metrics: "TreeMetrics") -> List[float]:
    """
    ERU from ``HABTYPCD1`` via :func:`fia_mog.southwest.core.southwest_eru`, or
    override with ``ConditionContext.ecosubcd`` as a precomputed ERU string.
    """

    eru = _resolve_southwest_eru(ctx)
    if not eru:
        return []

    vec: List[float] = []

    relative_sdi, sw_qmd, _ = _southwest_relative_sdi_and_qmd(ctx)
    _, _, local_og = _southwest_table9_family_rule_og(relative_sdi, sw_qmd, eru)

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
