"""
Default locations for MOG auxiliary data shipped with ``fia_py`` (R ``source.path``).

Layout mirrors ``FUNCTION_mapMOG.R`` expectations under ``mog_auxillary/``:
``utility.v9-5_*_Natl_MasterTreeSpeciesList.csv``, ``utility_MTcontDivide/``, etc.
"""

from __future__ import annotations

from pathlib import Path

_FIA_MOG_ROOT = Path(__file__).resolve().parent
_FIA_PY_ROOT = _FIA_MOG_ROOT.parent


def default_mog_auxiliary_dir() -> Path:
    """Sibling of ``fia_mog``: ``<fia_py>/mog_auxillary``."""

    return _FIA_PY_ROOT / "mog_auxillary"


def master_tree_species_csv(auxiliary_dir: Path | None = None) -> Path:
    root = auxiliary_dir if auxiliary_dir is not None else default_mog_auxiliary_dir()
    return root / "utility.v9-5_2024-10_Natl_MasterTreeSpeciesList.csv"


def montana_continental_divide_shp(auxiliary_dir: Path | None = None) -> Path:
    root = auxiliary_dir if auxiliary_dir is not None else default_mog_auxiliary_dir()
    return root / "utility_MTcontDivide" / "utility_MTcontDivide.shp"


def pacific_northwest_paz_raster(auxiliary_dir: Path | None = None) -> Path:
    """R ``utility_R5_PAZ.tif`` (plant association zones for Region 5 / PNW)."""

    root = auxiliary_dir if auxiliary_dir is not None else default_mog_auxiliary_dir()
    return root / "utility_R5_PAZ.tif"


def pacific_northwest_nwfp_boundary_shp(auxiliary_dir: Path | None = None) -> Path:
    """R ``utility_NWFPboundary/NWFPboundary.shp``."""

    root = auxiliary_dir if auxiliary_dir is not None else default_mog_auxiliary_dir()
    return root / "utility_NWFPboundary" / "NWFPboundary.shp"


def oregon_counties_shp(auxiliary_dir: Path | None = None) -> Path:
    """R ``utility_ORcounties/counties.shp`` (Hood River / Wasco override for white fir)."""

    root = auxiliary_dir if auxiliary_dir is not None else default_mog_auxiliary_dir()
    return root / "utility_ORcounties" / "counties.shp"
