"""
Run from the ``fia_py`` repo root, e.g. ``python scripts/test_mog_area.py``,
or install the project in editable mode (``pip install -e .``).
"""

from __future__ import annotations

import sys
from pathlib import Path

_repo_root = Path(__file__).resolve().parent.parent
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

from fia import clip_fia, read_fia
from fia_mog import old_growth_area

db = read_fia("./tests/fiadb/AZ")
db_current = clip_fia(db)

res = old_growth_area(
    db_current,
    states=["AZ"],        # optional
    grp_by=["ADFORCD"],          # or e.g. ["STATECD"]
    eval_typ="CURR",
    method="TI",
    totals=True,
    variance=True,
)

# Design-based area totals/vars live here:
print(res.df)

# Condition-level diagnostics (MOG_SCORE, OG_FLAG, etc.):
print(res.cond_mog.head())

# Forest area df to compare in R
print(res.forest_area_df)