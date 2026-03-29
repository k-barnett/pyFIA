"""
Orchestrate PNW MOG vector: inside-NWFP OGSI, outside Table 14, then mature (R gating).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, List

from .mature import pnw_mature_vector
from .ogsi_inside import inside_nwfp_og_vector
from .ogsi_outside import outside_nwfp_og_vector
from .paz import paz_value_to_group

if TYPE_CHECKING:
    from ..engine import ConditionContext, TreeMetrics


def pacific_northwest_mog_vector(ctx: "ConditionContext", metrics: "TreeMetrics") -> List[float]:
    """
    Full ``region %in% "pacific northwest"`` branch.

    Requires :class:`~fia_mog.engine.ConditionContext` fields populated for PNW
    (PAZ value, inside-NWFP flag, woody debris, site class, optional ECOSUBCD /
    OR-county overlay). See ``mog_condition_scores`` for automatic filling when
    auxiliary rasters/shapefiles and tables are available.
    """

    paz_raw = ctx.pnw_paz_raster_value
    paz_group = paz_value_to_group(paz_raw) if paz_raw is not None else None

    og: List[float] = []
    inside = ctx.pnw_inside_nwfp

    woody = ctx.pnw_woody_debris
    if woody is not None and woody.empty:
        woody = None

    if inside is True:
        og.extend(
            inside_nwfp_og_vector(
                paz_group,
                ctx.trees,
                ctx.condition_area_acres,
                woody,
            )
        )
    elif inside is False:
        site_cls = ctx.pnw_site_class_max
        og.extend(
            outside_nwfp_og_vector(
                paz_group,
                ctx.trees,
                ctx.condition_area_acres,
                ctx.stand_age,
                site_cls,
                ctx.ecosubcd,
                ctx.pnw_plot_in_or_counties_layer,
            )
        )

    last_og = og[-1] if og else 0.0
    out = list(og)
    if last_og == 0.0:
        out.extend(
            pnw_mature_vector(
                ctx.forest_type,
                paz_group,
                bool(inside) if inside is not None else False,
                metrics,
            )
        )

    return out
