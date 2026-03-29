"""
Loaders for files referenced by ``FUNCTION_mapMOG.R`` (``source.path``).

These match the R script’s use of ``read.csv`` / ``st_read`` on the auxiliary folder.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .paths import default_mog_auxiliary_dir, master_tree_species_csv, montana_continental_divide_shp


def load_usfs_master_tree_species_list(auxiliary_dir: Path | None = None) -> pd.DataFrame:
    """
    National Master Tree Species List (R ``USFStrees``).

    CSV columns are normalized to dot form like R ``read.csv`` (e.g. ``FIA Code`` → ``FIA.Code``)
    so merges can use ``by.x = 'SPCD', by.y = 'FIA.Code'`` as in mapMOG.
    """

    path = master_tree_species_csv(auxiliary_dir)
    df = pd.read_csv(path, dtype_backend="numpy_nullable")
    df.columns = [str(c).strip().replace(" ", ".") for c in df.columns]
    return df


def species_lookup_from_master_list(df: pd.DataFrame | None) -> pd.DataFrame | None:
    """Return a two-column frame suitable for :func:`fia_mog.crosswalk.northern_veg_code`."""

    if df is None or df.empty:
        return None
    if "FIA.Code" not in df.columns:
        return None
    plants_col = "PLANTS.Code" if "PLANTS.Code" in df.columns else None
    if plants_col is None:
        for c in df.columns:
            if str(c).upper().replace(".", "_") in {"PLANTS_CODE", "PLANTS_SYMBOL"}:
                plants_col = c
                break
    if plants_col is None:
        return None
    out = df[["FIA.Code", plants_col]].copy()
    out["FIA.Code"] = pd.to_numeric(out["FIA.Code"], errors="coerce")
    return out.rename(columns={"FIA.Code": "SPCD", plants_col: "PLANTS.Code"}).dropna(subset=["SPCD"])


def montana_plot_in_cont_divide_east(
    lon: float,
    lat: float,
    auxiliary_dir: Path | None = None,
    *,
    crs_plot: int | str = 4326,
) -> bool | None:
    """
    Whether an FIA plot point intersects the ``utility_MTcontDivide`` polygon layer.

    R uses ``st_filter(local.plot, ContDivideEast)`` and then ``nrow(...) > 1`` for
    “eastern Montana zone” vs ``< 1`` for western. For ordinary point plots the
    filtered row count is 0 or 1; treating **``>= 1``** as “in / east-side test
    polygon” matches the intended split when a single plot row is returned.

    Returns ``None`` if the shapefile is missing or ``geopandas`` cannot be imported.
    """

    shp = montana_continental_divide_shp(auxiliary_dir)
    if not shp.is_file():
        return None
    try:
        import geopandas as gpd  # type: ignore[import-not-found]
        from shapely.geometry import Point  # type: ignore[import-untyped]
    except ImportError:
        return None

    divide = gpd.read_file(shp)
    pt = gpd.GeoDataFrame(geometry=[Point(float(lon), float(lat))], crs=crs_plot)
    pt = pt.to_crs(divide.crs)
    joined = gpd.sjoin(pt, divide, predicate="intersects", how="inner")
    n = len(joined)
    if n >= 1:
        return True
    return False


def montana_northern_subregion_flags(
    lon: float,
    lat: float,
    auxiliary_dir: Path | None = None,
) -> tuple[bool | None, bool | None]:
    """
    Return ``(is_eastern_montana_zone, is_western_montana_zone)`` for MT plots.

    If the divide test is unavailable, returns ``(None, None)`` (caller should skip
    habitat OG or supply flags another way).
    """

    inside = montana_plot_in_cont_divide_east(lon, lat, auxiliary_dir)
    if inside is None:
        return None, None
    if inside:
        return True, False
    return False, True


def resolve_mog_auxiliary_dir(explicit: Path | str | None) -> Path:
    """Use an explicit path, else default under ``fia_py``."""

    if explicit is not None:
        return Path(explicit).expanduser().resolve()
    return default_mog_auxiliary_dir()
