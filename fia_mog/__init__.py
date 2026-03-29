"""
Stand-alone MOG (mature / old-growth) logic aligned with ``FUNCTION_mapMOG.R``.

Auxiliary tables and shapefiles live in ``fia_py/mog_auxillary`` (R ``source.path``).
Use :func:`fia_mog.auxiliary.load_usfs_master_tree_species_list` for the same
``USFStrees`` merge as mapMOG (``FIA.Code`` → ``PLANTS.Code``).
"""

from __future__ import annotations

from .auxiliary import (
    load_usfs_master_tree_species_list,
    montana_plot_in_cont_divide_east,
    resolve_mog_auxiliary_dir,
    species_lookup_from_master_list,
)
from .crosswalk import classify_region, statecd_to_abbrev
from .engine import ConditionContext, MOGEngine, compute_tree_metrics
from .estimators import MOGAreaResult, mog_condition_scores, old_growth_area
from .paths import default_mog_auxiliary_dir, master_tree_species_csv

__all__ = [
    "MOGAreaResult",
    "MOGEngine",
    "ConditionContext",
    "compute_tree_metrics",
    "classify_region",
    "statecd_to_abbrev",
    "default_mog_auxiliary_dir",
    "master_tree_species_csv",
    "resolve_mog_auxiliary_dir",
    "load_usfs_master_tree_species_list",
    "species_lookup_from_master_list",
    "montana_plot_in_cont_divide_east",
    "mog_condition_scores",
    "old_growth_area",
]
