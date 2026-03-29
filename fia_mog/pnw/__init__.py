"""
Pacific Northwest (R ``region %in% "pacific northwest"``) MOG helpers.

Spatial layers (under ``mog_auxillary/``): ``utility_R5_PAZ.tif``,
``utility_NWFPboundary/NWFPboundary.shp``, optional ``utility_ORcounties/counties.shp``.

The full MOG vector is assembled by :class:`fia_mog.engine.MOGEngine` via
``pnw.evaluate.pacific_northwest_mog_vector`` (not imported here to keep
``import fia_mog.crosswalk`` lightweight).

Raster and vector sampling use ``rasterio`` and ``geopandas`` (declared as core dependencies).
"""

from __future__ import annotations

from .paz import paz_value_to_group
from .spatial import (
    oregon_plot_in_counties_layer,
    plot_inside_nwfp,
    sample_r5_paz_value,
)

__all__ = [
    "oregon_plot_in_counties_layer",
    "paz_value_to_group",
    "plot_inside_nwfp",
    "sample_r5_paz_value",
]
