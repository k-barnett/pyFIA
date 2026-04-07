"""
Northern (Region 1) vegetation / habitat crosswalk diagnostics on ``mog_condition_scores`` output.

Used to count Montana (and other northern) conditions where **Green et al.** habitat OG rules
cannot run because the continental divide, species→PLANTS lookup, P2VEG understory, or
**veg_code → habitat letter** crosswalk is missing or unmatched.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import pandas as pd

from fia_mog.crosswalk import statecd_to_abbrev

from .core import (
    fallback_northern_habitat_letters,
    infer_northern_og_type_from_species,
    northern_basal_area_per_acre,
    northern_habitat_letters,
    northern_og_forest_type,
    northern_subregion,
    northern_veg_code,
    refine_habitat_letters_environmental,
)
from .og_dispatch import (
    northern_east_mt_og_vector,
    northern_idaho_og_vector,
    northern_west_mt_og_vector,
)

# Columns appended by :func:`fia_mog.estimators.mog_condition_scores` for northern rows.
NORTHERN_VEG_DIAGNOSTIC_KEYS: tuple[str, ...] = (
    "NR_SUBREGION",
    "NR_OG_FOREST_TYPE",
    "NR_VEG_CODE",
    "NR_HABITAT_LETTERS",
    "NR_N_HABITAT_LETTERS",
    "NR_N_OG_BLOCKS",
    "NR_FLAG_MT_SUBREGION_MISSING",
    "NR_FLAG_VEG_CODE_MISSING",
    "NR_FLAG_VEG_NOT_IN_CROSSWALK",
    "NR_FLAG_NO_HABITAT_OG_BLOCKS",
    "NR_OG_TYPE_INFERRED",
    "NR_HABITAT_LETTERS_FALLBACK",
)

NR_MATURITY_SCORE_KEY = "NR_MATURITY_SCORE"

MOG_CONDITION_NR_COLUMNS: tuple[str, ...] = NORTHERN_VEG_DIAGNOSTIC_KEYS + (NR_MATURITY_SCORE_KEY,)


def nr_flag_true(series: pd.Series) -> pd.Series:
    """True where a numeric NR_* flag column equals 1; NA and non-numeric treated as 0."""

    return pd.to_numeric(series, errors="coerce").fillna(0) == 1


def _northern_state_abbrev_from_ctx(ctx: object) -> str | None:
    """
    State for ``northern_subregion``: ``plot_statecd`` on context, else first TREE ``STATECD``.
    """

    s = statecd_to_abbrev(getattr(ctx, "plot_statecd", None))
    if s is not None:
        return s
    trees = getattr(ctx, "trees", None)
    if trees is not None and len(trees) > 0 and "STATECD" in trees.columns:
        return statecd_to_abbrev(trees["STATECD"].iloc[0])
    return None


def _optional_plot_float(val: object) -> float | None:
    if val is None:
        return None
    try:
        if pd.isna(val):
            return None
    except TypeError:
        return None
    try:
        x = float(val)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


_NORTHERN_OG_ZONES: frozenset[str] = frozenset(
    {"northern Idaho zone", "western Montana zone", "eastern Montana zone"}
)

# When the plot's zone yields no rule blocks, try other Green et al. tables in this order.
_CROSS_ZONE_OG_ORDER: dict[str, tuple[str, ...]] = {
    "western Montana zone": ("eastern Montana zone", "northern Idaho zone"),
    "eastern Montana zone": ("western Montana zone", "northern Idaho zone"),
    "northern Idaho zone": ("western Montana zone", "eastern Montana zone"),
}


def _extend_habitat_og_for_zone(
    vec: list[float],
    zone: str,
    letters_list: list[str],
    og_type: str | None,
    *,
    ctx: object,
    basal: float,
) -> None:
    kw = dict(
        raw_stand_age=ctx.stand_age,
        basal_area_per_acre=basal,
        condition_area_acres=ctx.condition_area_acres,
        trees=ctx.trees,
    )
    if zone == "northern Idaho zone":
        vec.extend(northern_idaho_og_vector(letters_list, og_type, **kw))
    elif zone == "western Montana zone":
        vec.extend(northern_west_mt_og_vector(letters_list, og_type, **kw))
    elif zone == "eastern Montana zone":
        vec.extend(northern_east_mt_og_vector(letters_list, og_type, **kw))


def filter_montana_mog_rows(cond_df: pd.DataFrame) -> pd.DataFrame:
    """
    Subset to Montana (``STATECD`` 30), robust to string/float encodings (e.g. CSV round-trips).

    If every ``STATECD`` is missing, returns ``cond_df`` unchanged (caller already restricted
    the FIADB extract to Montana). If ``STATECD`` is absent, returns an empty frame.
    """

    if cond_df.empty:
        return cond_df
    if "STATECD" not in cond_df.columns:
        return cond_df.iloc[0:0].copy()
    sc = pd.to_numeric(cond_df["STATECD"], errors="coerce")
    if sc.isna().all():
        return cond_df.copy()
    return cond_df.loc[sc == 30].copy()


@dataclass(frozen=True)
class NorthernHabitatOGBundle:
    """Habitat-layer OG inputs shared by :class:`fia_mog.engine._NorthernEvaluator` and diagnostics."""

    subregion: str | None
    og_forest_type: str | None
    og_type_inferred: bool
    veg_code: str | None
    letters: tuple[str, ...]
    letters_from_exact_crosswalk: bool
    habitat_letters_fallback: bool
    og_scores: tuple[float, ...]


def compute_northern_habitat_og_bundle(ctx: object) -> NorthernHabitatOGBundle:
    """
    Replicate habitat OG vector construction (no Table 19 maturity).

    ``ctx`` is a :class:`fia_mog.engine.ConditionContext`; typed as ``object`` to avoid a
    circular import at module load.

    **Eastern Montana zone:** Green et al. only pairs ``PP`` with habitat letters A/B/C/K,
    while the east veg table also assigns D–J (shared with fir/spruce rows). If FIA maps
    ponderosa to ``PP`` but letters are D–J, the east dispatch emits **no** blocks. We then
    try ``PF`` (east ponderosa shading) and, if still empty, **western Montana** ``PP``/DF/LP
    rule blocks for the same letters (same numeric thresholds as the west-zone R branch).

    **Species-based OG override:** if habitat letters exist but the primary zone still emits
    no blocks, we re-run species composition inference and dispatch again when it disagrees
    with the type-code OG (e.g. ponderosa type with fir-dominated trees).

    **Cross-zone habitat OG dispatch:** if blocks are still empty, we try the other northern
    zone tables for the same ``og`` and letters (NR_SUBREGION is unchanged; this only fills
    the score vector when a neighboring table has a matching letter/type pair).
    """

    ft = int(ctx.forest_type)
    og = northern_og_forest_type(ft)
    og_type_inferred = False
    if og is None:
        og = infer_northern_og_type_from_species(ctx.trees, ctx.northern_species_lookup)
        og_type_inferred = og is not None

    st = _northern_state_abbrev_from_ctx(ctx)
    sub = northern_subregion(
        st,
        mt_east_of_continental_divide=ctx.northern_mt_east_of_divide,
        plot_lon=_optional_plot_float(getattr(ctx, "plot_lon", None)),
    )

    basal = northern_basal_area_per_acre(ctx.trees, ctx.condition_area_acres)
    plt_cn = ctx.trees["PLT_CN"].iloc[0] if "PLT_CN" in ctx.trees.columns else None
    condid = ctx.trees["CONDID"].iloc[0] if "CONDID" in ctx.trees.columns else None

    veg_code = northern_veg_code(
        ctx.trees,
        ctx.northern_species_lookup,
        ctx.northern_veg_subplot,
        plt_cn,
        condid,
    )
    letters_list = northern_habitat_letters(sub, veg_code)
    letters_from_exact_crosswalk = len(letters_list) > 0
    habitat_letters_fallback = False

    if not letters_list and sub:
        letters_list = fallback_northern_habitat_letters(
            sub,
            veg_code,
            ctx.trees,
            ctx.northern_species_lookup,
        )
        habitat_letters_fallback = len(letters_list) > 0

    if habitat_letters_fallback and len(letters_list) > 2:
        letters_list = refine_habitat_letters_environmental(
            sub,
            letters_list,
            siteclcd=_optional_plot_float(getattr(ctx, "condition_siteclcd", None)),
            elev_ft=_optional_plot_float(getattr(ctx, "plot_elev_ft", None)),
        )

    letters = tuple(letters_list)

    vec: list[float] = []
    if sub == "northern Idaho zone":
        vec.extend(
            northern_idaho_og_vector(
                letters_list,
                og,
                raw_stand_age=ctx.stand_age,
                basal_area_per_acre=basal,
                condition_area_acres=ctx.condition_area_acres,
                trees=ctx.trees,
            )
        )
    elif sub == "western Montana zone":
        vec.extend(
            northern_west_mt_og_vector(
                letters_list,
                og,
                raw_stand_age=ctx.stand_age,
                basal_area_per_acre=basal,
                condition_area_acres=ctx.condition_area_acres,
                trees=ctx.trees,
            )
        )
    elif sub == "eastern Montana zone":
        vec.extend(
            northern_east_mt_og_vector(
                letters_list,
                og,
                raw_stand_age=ctx.stand_age,
                basal_area_per_acre=basal,
                condition_area_acres=ctx.condition_area_acres,
                trees=ctx.trees,
            )
        )

    og_out = og
    area_ac = float(ctx.condition_area_acres or 0.0)

    if (
        not vec
        and letters_list
        and sub == "eastern Montana zone"
        and og == "PP"
        and area_ac > 0.0
    ):
        pf_vec = northern_east_mt_og_vector(
            letters_list,
            "PF",
            raw_stand_age=ctx.stand_age,
            basal_area_per_acre=basal,
            condition_area_acres=ctx.condition_area_acres,
            trees=ctx.trees,
        )
        if pf_vec:
            vec.extend(pf_vec)
            og_out = "PF"

    if (
        not vec
        and letters_list
        and sub == "eastern Montana zone"
        and og
        and area_ac > 0.0
    ):
        vec.extend(
            northern_west_mt_og_vector(
                letters_list,
                og,
                raw_stand_age=ctx.stand_age,
                basal_area_per_acre=basal,
                condition_area_acres=ctx.condition_area_acres,
                trees=ctx.trees,
            )
        )

    # Species-based OG override: type code can disagree with stand composition.
    if not vec and letters_list and sub in _NORTHERN_OG_ZONES and area_ac > 0.0:
        og_retry = infer_northern_og_type_from_species(ctx.trees, ctx.northern_species_lookup)
        if og_retry and og_retry != og_out:
            before = len(vec)
            _extend_habitat_og_for_zone(
                vec, sub, letters_list, og_retry, ctx=ctx, basal=basal
            )
            if len(vec) > before:
                og_out = og_retry
                og_type_inferred = True

    # Cross-zone habitat OG dispatch (NR_SUBREGION stays the plot's true zone).
    if not vec and letters_list and og_out and area_ac > 0.0 and sub in _NORTHERN_OG_ZONES:
        for alt in _CROSS_ZONE_OG_ORDER.get(sub, ()):
            before = len(vec)
            _extend_habitat_og_for_zone(
                vec, alt, letters_list, og_out, ctx=ctx, basal=basal
            )
            if len(vec) > before:
                break

    return NorthernHabitatOGBundle(
        subregion=sub,
        og_forest_type=og_out,
        og_type_inferred=og_type_inferred,
        veg_code=veg_code,
        letters=letters,
        letters_from_exact_crosswalk=letters_from_exact_crosswalk,
        habitat_letters_fallback=habitat_letters_fallback,
        og_scores=tuple(vec),
    )


def northern_veg_diagnostic_row(ctx: object, mog_vec: list[float]) -> dict[str, object]:
    """
    Per-condition northern vegetation / habitat diagnostics.

    ``NR_FLAG_VEG_NOT_IN_CROSSWALK``: subregion assigned, ``veg_code`` built, but **no**
    habitat letters after exact table + prefix fallback (true crosswalk / composition gap).
    For an exact-table miss with successful fallback, see ``NR_HABITAT_LETTERS_FALLBACK``
    and ``NR_N_HABITAT_LETTERS`` while this flag stays **0**.

    ``NR_FLAG_NO_HABITAT_OG_BLOCKS``: **zero** habitat OG 0/1 rule blocks were emitted
    after OG-type and letter resolution (unmapped composition types such as aspen
    still yield no blocks).

    ``NR_OG_TYPE_INFERRED`` / ``NR_HABITAT_LETTERS_FALLBACK``: QA flags when
    ``FLDTYPCD`` or the printed veg crosswalk did not fully specify habitat OG inputs.
    """

    b = compute_northern_habitat_og_bundle(ctx)
    st = _northern_state_abbrev_from_ctx(ctx)
    is_mt = st == "MT"

    letters_joined = ",".join(b.letters) if b.letters else ""
    n_let = len(b.letters)
    n_og = len(b.og_scores)

    flag_mt_sub = 1 if (is_mt and b.subregion is None) else 0
    flag_veg_miss = 1 if (b.subregion is not None and b.veg_code is None) else 0
    # No habitat letters after exact match + prefix fallback (actionable crosswalk gap).
    flag_cross = 1 if (b.subregion is not None and b.veg_code is not None and n_let == 0) else 0
    flag_no_blocks = 1 if n_og == 0 else 0

    n_og_len = len(b.og_scores)
    mat_max: float | None
    if len(mog_vec) > n_og_len:
        mat_max = max(float(x) for x in mog_vec[n_og_len:])
    else:
        mat_max = None

    out: dict[str, object] = {
        "NR_SUBREGION": b.subregion if b.subregion else "",
        "NR_OG_FOREST_TYPE": b.og_forest_type if b.og_forest_type else "",
        "NR_VEG_CODE": b.veg_code if b.veg_code else "",
        "NR_HABITAT_LETTERS": letters_joined,
        "NR_N_HABITAT_LETTERS": int(n_let),
        "NR_N_OG_BLOCKS": int(n_og),
        "NR_FLAG_MT_SUBREGION_MISSING": int(flag_mt_sub),
        "NR_FLAG_VEG_CODE_MISSING": int(flag_veg_miss),
        "NR_FLAG_VEG_NOT_IN_CROSSWALK": int(flag_cross),
        "NR_FLAG_NO_HABITAT_OG_BLOCKS": int(flag_no_blocks),
        "NR_OG_TYPE_INFERRED": int(b.og_type_inferred),
        "NR_HABITAT_LETTERS_FALLBACK": int(b.habitat_letters_fallback),
        NR_MATURITY_SCORE_KEY: float(mat_max) if mat_max is not None else pd.NA,
    }
    if set(out) != set(MOG_CONDITION_NR_COLUMNS):
        raise RuntimeError("northern_veg_diagnostic_row keys drifted from MOG_CONDITION_NR_COLUMNS")
    return out


def summarize_mt_northern_veg_gaps(cond_df: pd.DataFrame) -> pd.DataFrame:
    """
    Summarize Montana forest conditions where habitat OG crosswalk / inputs are incomplete.

    Expects ``mog_condition_scores`` output with **STATECD** (Montana = 30).
    """

    if cond_df.empty:
        return pd.DataFrame()
    need = {"PLT_CN", "STATECD", "NR_FLAG_MT_SUBREGION_MISSING"}
    if not need.issubset(cond_df.columns):
        raise ValueError(
            "cond_df must include STATECD and NR_* columns from mog_condition_scores (Montana run)."
        )

    mt = filter_montana_mog_rows(cond_df)
    if mt.empty:
        return pd.DataFrame()

    rows: list[dict[str, object]] = []

    def _add_row(label: str, mask: pd.Series) -> None:
        sub = mt.loc[mask]
        n_c = int(len(sub))
        n_p = int(sub["PLT_CN"].nunique()) if n_c and "PLT_CN" in sub.columns else 0
        rows.append({"issue": label, "n_conditions": n_c, "n_plots": n_p})

    all_mt = mt["PLT_CN"].notna()
    _add_row("all MT forest conditions (mog rows)", all_mt)

    m_sub = nr_flag_true(mt["NR_FLAG_MT_SUBREGION_MISSING"])
    _add_row("MT: continental divide / subregion missing (NR_FLAG_MT_SUBREGION_MISSING)", m_sub)

    if "NR_FLAG_VEG_CODE_MISSING" in mt.columns:
        _add_row("MT: veg_code missing (species/P2VEG/lookup)", nr_flag_true(mt["NR_FLAG_VEG_CODE_MISSING"]))

    if "NR_FLAG_VEG_NOT_IN_CROSSWALK" in mt.columns:
        _add_row(
            "MT: veg_code present but no habitat letters after exact + prefix fallback (NR_FLAG_VEG_NOT_IN_CROSSWALK)",
            nr_flag_true(mt["NR_FLAG_VEG_NOT_IN_CROSSWALK"]),
        )

    if "NR_FLAG_NO_HABITAT_OG_BLOCKS" in mt.columns:
        _add_row(
            "MT: zero habitat OG rule blocks (NR_FLAG_NO_HABITAT_OG_BLOCKS)",
            nr_flag_true(mt["NR_FLAG_NO_HABITAT_OG_BLOCKS"]),
        )

    combined = m_sub.copy()
    for c in (
        "NR_FLAG_VEG_CODE_MISSING",
        "NR_FLAG_VEG_NOT_IN_CROSSWALK",
        "NR_FLAG_NO_HABITAT_OG_BLOCKS",
    ):
        if c in mt.columns:
            combined = combined | nr_flag_true(mt[c])

    _add_row("MT: any of the above habitat-input / crosswalk gaps", combined)

    return pd.DataFrame(rows)


__all__ = [
    "MOG_CONDITION_NR_COLUMNS",
    "NORTHERN_VEG_DIAGNOSTIC_KEYS",
    "NR_MATURITY_SCORE_KEY",
    "NorthernHabitatOGBundle",
    "compute_northern_habitat_og_bundle",
    "filter_montana_mog_rows",
    "nr_flag_true",
    "northern_veg_diagnostic_row",
    "summarize_mt_northern_veg_gaps",
]
