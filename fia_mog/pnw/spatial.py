"""
Sample ``utility_R5_PAZ.tif`` and intersect NWFP / Oregon county layers (R ``terra`` / ``sf``).

Uses ``rasterio`` and ``geopandas`` (core ``fia_py`` dependencies) when layers are present.
"""

from __future__ import annotations

from pathlib import Path


def sample_r5_paz_value(
    lon: float,
    lat: float,
    raster_path: Path | str,
    *,
    plot_crs: int = 4326,
) -> float | None:
    """
    Extract a single PAZ cell at ``(lon, lat)`` (WGS84 by default), matching R
    ``project(vect(local.plot), paz)`` then ``terra::extract``.

    Returns ``None`` if ``rasterio`` is missing, the file is absent, or the
    value is nodata.
    """

    path = Path(raster_path)
    if not path.is_file():
        return None
    try:
        import numpy as np
        import rasterio
        from rasterio.warp import transform as warp_transform
    except ImportError:
        return None

    with rasterio.open(path) as ds:
        xs, ys = warp_transform(plot_crs, ds.crs, [float(lon)], [float(lat)])
        try:
            samples = list(ds.sample([(float(xs[0]), float(ys[0]))]))
        except ValueError:
            return None
        if not samples:
            return None
        val = float(samples[0][0])
        nodata = ds.nodata
        if nodata is not None and (
            val == nodata or (isinstance(nodata, float) and np.isnan(nodata) and np.isnan(val))
        ):
            return None
        if np.isnan(val):
            return None
    return val


def plot_inside_nwfp(
    lon: float,
    lat: float,
    nwfp_shp: Path | str,
    *,
    plot_crs: int = 4326,
) -> bool | None:
    """
    R ``nrow(st_filter(local.plot, nwfp)) > 0``.

    Returns ``None`` if ``geopandas`` is unavailable or the shapefile is missing.
    """

    path = Path(nwfp_shp)
    if not path.is_file():
        return None
    try:
        import geopandas as gpd  # type: ignore[import-not-found]
        from shapely.geometry import Point  # type: ignore[import-untyped]
    except ImportError:
        return None

    nwfp = gpd.read_file(path)
    pt = gpd.GeoDataFrame(geometry=[Point(float(lon), float(lat))], crs=plot_crs)
    pt = pt.to_crs(nwfp.crs)
    joined = gpd.sjoin(pt, nwfp, predicate="intersects", how="inner")
    return len(joined) > 0


def oregon_plot_in_counties_layer(
    lon: float,
    lat: float,
    counties_shp: Path | str,
    *,
    plot_crs: int = 4326,
) -> bool | None:
    """
    R ``nrow(st_filter(local.plot, counties)) > 0`` for Hood River / Wasco override.

    Returns ``None`` if geospatial stack or file is unavailable.
    """

    path = Path(counties_shp)
    if not path.is_file():
        return None
    try:
        import geopandas as gpd  # type: ignore[import-not-found]
        from shapely.geometry import Point  # type: ignore[import-untyped]
    except ImportError:
        return None

    counties = gpd.read_file(path)
    pt = gpd.GeoDataFrame(geometry=[Point(float(lon), float(lat))], crs=plot_crs)
    pt = pt.to_crs(counties.crs)
    joined = gpd.sjoin(pt, counties, predicate="intersects", how="inner")
    return len(joined) > 0
