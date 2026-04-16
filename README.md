## FIA (Python)

`pyFIA` is an in-progress Python port of the R package [`rFIA`](https://github.com/doserjef/rFIA), providing a user-friendly interface to the USFS Forest Inventory and Analysis (FIA) Database.

The goal is to mirror the core public API of `rFIA` while adopting Pythonic conventions and the scientific Python stack (e.g., `pandas`).

### Recent updates

- **`area(..., cond_list=True)`** now returns **`AREA_BASIS`** (from `COND.PROP_BASIS`, normalized to `SUBP` / `MACR`) together with **`PLT_CN`**, **`CONDID`**, **`YEAR`**, optional **`grp_by`** columns, and **`PROP_FOREST`**, so condition lists plug directly into **`custom_pse`** with correct non-response adjustment.
- **`custom_pse`** for **`method="TI"`** and condition-level **`AREA_BASIS`** inputs (no **`TREE_BASIS`**) uses the same stratified pipeline as **`area()`** (`_rfia_sum_to_plot` → `_rfia_sum_to_eu` → **`combine_mr`**), so totals such as **`PROP_FOREST_TOTAL`** match **`AREA_TOTAL`** from **`area()`** for the same grouping. Tree-level inputs and non-TI methods still use the legacy plot-mean × **`AREA_USED`** shortcut inside each estimation unit.
- **`get_fia(..., load=False)`** removes downloaded **`.zip`** archives after extracting CSVs into **`dir`**, leaving a clean directory of tables.
- **Old-growth / MOG** logic is implemented in the stand-alone **`fia_mog`** package (aligned with **`mog_auxillary/FUNCTION_mapMOG.R`**). **`mog_condition_scores`**, **`old_growth_area`**, and **`ConditionContext` / `MOGEngine`** remain importable from **`fia`** for compatibility; prefer **`fia_mog`** for new code.

### Install

**From a clone or sdist / wheel** (project root):

```bash
pip install .
# editable while developing:
pip install -e .
```

**Optional dev dependencies** (tests, building wheels):

```bash
pip install -e ".[dev]"
```

**After you publish to PyPI**, install with:

```bash
pip install fia_py
```

*(The distribution name on PyPI is **`fia_py`**; you still `import fia` and `import fia_mog` in code.)*

Then in Python:

```python
from fia import (
    FiaDatabase,
    get_fia,
    read_fia,
    clip_fia,
    tpa,
    area,
    custom_pse,
    # Old-growth / MOG (optional)
    ConditionContext,
    MOGEngine,
    compute_tree_metrics,
    MOGAreaResult,
    mog_condition_scores,
    old_growth_area,
)
```

### Current functionality

#### Core data access

- **`get_fia(states, dir=None, common=True, tables=None, load=True)`**
  - Downloads FIA CSVs from the FIA Datamart and optionally loads them into memory.
  - **`states`** may be a single string (e.g. **`"AZ"`**) or a sequence — a string is treated as one state, not as individual characters.
  - When **`dir`** is set, the folder is **created if it does not exist** (including parents).
  - **`tables` is required** in the Python port (full-state ZIP-only workflow is not implemented yet). **`tables`** may be one table name as a string or a list of names.
  - Merges per-state tables into a single logical table per name (e.g. `PLOT`, `TREE`).
  - Returns a **`FiaDatabase`** (a mapping from table name → `pandas.DataFrame`) when **`load=True`**; returns **`None`** when **`load=False`** after writing CSVs to **`dir`**.
  - When **`load=False`**, **`.zip`** files are deleted after extraction so **`dir`** contains only CSVs.
  - If the expected CSV for a requested state/table (e.g. **`AZ_PLOT.csv`**) already exists in **`dir`**, that table is **not** re-downloaded; **`get_fia`** prints a skip message. With **`load=True`**, it still loads from the existing file.

- **`read_fia(dir, tables=None, common=True, states=None)`**
  - Reads into a **`FiaDatabase`** from either a **directory of CSVs** (Datamart-style names such as **`AZ_PLOT.csv`**) **or** a **SQLite FIADB file** path (e.g. **`.sqlite`** / **`.db`**). If **`dir`** points to a file, it is opened with **`sqlite3`**; table and view names are matched case-insensitively and stored under **uppercase** keys like **`PLOT`**, **`TREE`**. Column names in each loaded table are uppercased so they match CSV-style FIADB field names.
  - Supports selecting a specific list of tables (**`tables`**) or a common subset when **`common=True`** (similar to **`readFIA(..., common=TRUE)`** in R).
  - Optional **`states`**: for **CSV** folders, only loads prefixed files for those codes (e.g. **`states=["MT"]`** reads **`MT_PLOT.csv`** but skips **`AZ_PLOT.csv`**); unprefixed national **`PLOT.csv`**-style files load only if **`"ENTIRE"`** is included; **`REF_SPECIES`**-style files require **`"REF"`**. For **SQLite**, filters **`PLOT`** by **`STATECD`** and child tables by **`PLT_CN`** / **`STATECD`** ( **`PLOT`** must be in **`tables`** or use **`common=True`** ).

- **Unclipped FIADB and estimators**
  - **`tpa`**, **`area`**, **`cond_height_percentiles`**, **`cond_mean_crown_ratio`**, and design helpers normalize **`PLOT` / `COND` / `TREE`** columns (uppercase names; if **`PLOT`** has **`CN`** but not **`PLT_CN`**, **`PLT_CN`** is set from **`CN`**) so you can run them **without** calling **`clip_fia`** first. Clipping is still recommended when you want the **most recent evaluation** subset only.

#### Database clipping

- **`clip_fia(db, most_recent=True)`**
  - In-memory analogue of **`clipFIA`** (temporal only, no spatial mask yet).
  - When POP design tables are present, restricts **`db`** to plots and design records belonging to the **most recent inventory per state**.
  - Returns a new **`FiaDatabase`** suitable for “current” snapshots like **`az_current`** in the R examples.

#### Design-based engine (TI)

Implemented in **`fia.design`**:

- **`get_design_info(db, eval_types=("VOL",), most_recent=True)`** — port of **`getDesignInfo()`**, joining **`POP_EVAL`**, **`POP_EVAL_TYP`**, **`POP_ESTN_UNIT`**, **`POP_STRATUM`**, **`POP_PLOT_STRATUM_ASSGN`**.
- **`handle_pops(db, eval_types, method="TI", most_recent=True)`** — TI-focused port of **`handlePops()`**; recodes **`ESTN_METHOD`** and returns the population / design frame used by estimators.
- **`ratio_var(...)`**, **`combine_mr(...)`** — ports of **`ratioVar()`** and **`combineMR()`** from rFIA’s utilities.

#### Estimators

##### **`tpa`** — tree abundance and basal area

```python
tpa(
    db,
    grp_by=None,
    land_type="forest",
    tree_type="live",
    by_species=False,
    by_size_class=None,
    by_plot=False,
    by_cond=False,
    tree_list=False,
)
```

- **`land_type`**: `"forest"`, `"timber"`, `"all"` (domain via **`COND_STATUS_CD`**, **`SITECLCD`**, **`RESERVCD`**).
- **`tree_type`**: `"live"`, `"dead"`, `"gs"`, `"all"`.
- **`by_species=True`**: groups by species (**`SPECIES`** table when available; else **`SPCD`**).
- **`by_size_class`**: ``None`` (default) skips diameter classes; otherwise a **non-empty list/tuple of inch breakpoints** (positive, strictly increasing). With *k* cutoffs you get **k + 1** half-open bins: ``[0, c0)``, ``[c0, c1)``, …, ``[c_{k-2}, c_{k-1})``, and a **final open-ended bin** ``[c_{k-1}, ∞)`` (trees at or above the largest cutoff go here). Column **`sizeClass`** uses pandas **`Interval`** labels. Example: ``(2, 4, 6, 8)`` → five bins: 0–2, 2–4, 4–6, 6–8, and 8–∞.
- **`by_plot=True`**: per-plot **`TPA`**, **`BAA`**, **`PROP_FOREST`**.
- **`by_cond=True`**: per **condition** (**`PLT_CN`**, **`CONDID`**, **`YEAR`**, optional **`grp_by`**) **`TPA`**, **`BAA`**, and condition-level **`PROP_FOREST`**. Do not set together with **`by_plot=True`**.
- **`tree_list=True`**: per-tree rows with **`PLT_CN`**, **`CONDID`**, **`SUBP`**, **`TREE`**, **`YEAR`**, **`TPA`**, **`BAA`**, **`PROP_FOREST`** (and related columns as implemented).

With POP design tables, uses a **TI design-based pipeline** (stratum → estimation unit → totals and variances): **`TPA`**, **`BAA`**, **`TREE_TOTAL`**, **`BA_TOTAL`**, **`AREA_TOTAL`**, variance / SE columns, and **`nPlots_*`**, **`N`**. Without POP tables, falls back to a simple plot-mean estimator.

##### **`cond_height_percentiles`** — height percentiles by plot–condition

```python
cond_height_percentiles(
    db,
    land_type="forest",
    tree_type="live",
    percentile_method="raw",
    weight_by="tpa",
    percentiles=(5, 25, 50, 75, 95),
)
```

- One row per **`(PLT_CN, CONDID, YEAR)`** with **`n_trees`** and columns **`HT_P{p}`** for each requested percentile (heights in feet per **`TREE.HT`**).
- **`tree_type`**: **`"all"`**, **`"live"`**, or **`"dead"`** (same domain flags as **`tpa`** via **`tDI`**).
- **`percentile_method`**: **`"raw"`** (unweighted over in-domain trees) or **`"weighted"`** (cumulative-weight interpolation as in the internal **`_weighted_percentile`** helper).
- **`weight_by`** (weighted only): **`"tpa"`** or **`"baa"`** using the same tree-level **`TPA`** / **`BAA`** contributions as **`tpa`** (i.e. **`TPA_UNADJ × tDI`** and **`basal_area × TPA_UNADJ × tDI`**).

##### **`cond_mean_crown_ratio`** — mean crown ratio by plot–condition

```python
cond_mean_crown_ratio(
    db,
    land_type="forest",
    tree_type="live",
    weight_by="tpa",
    weighted=True,
)
```

- One row per **`(PLT_CN, CONDID, YEAR)`** with **`n_trees`** and **`mean_crown_ratio`** ( **`TREE.CR`** as a **0–1 proportion**: FIADB’s usual **0–100** crown ratio is divided by 100).
- **`land_type`** / **`tree_type`**: same domains as **`cond_height_percentiles`** (**`"all"`**, **`"live"`**, **`"dead"`** for trees).
- **`weighted=True`** (default): **`numpy.average(CR, weights=…)`** with **`weight_by="tpa"`** or **`"baa"`** (same tree-level weights as **`tpa`**). **`weighted=False`**: simple mean of **`CR`** over in-domain trees, still scaled to proportion.

##### **`area`** — land area and percentage of area

```python
area(
    db,
    grp_by=None,
    land_type="forest",
    cond_list=False,
    by_plot=False,
)
```

- **`land_type`**: `"forest"`, `"timber"`, `"all"`.
- **`cond_list=True`**: one row per condition with **`PLT_CN`**, **`CONDID`**, **`YEAR`**, **`AREA_BASIS`**, **`PROP_FOREST`** ( **`CONDPROP_UNADJ × landD`** ), plus any **`grp_by`** columns present on the joined plot–condition frame. Use with **`custom_pse`** and the same **`EVAL_TYP`** as in **`handle_pops`** (e.g. **`CURR`**).
- **`by_plot=True`**: per-plot **`PROP_FOREST`** by **`PLT_CN`**, **`YEAR`**, and **`grp_by`**.

With POP tables, the TI path returns **`PERC_AREA`**, **`AREA_TOTAL`** (acres), variance / SE columns, **`nPlots_AREA`**, **`N`**, etc. **`PERC_AREA`** is a ratio of numerator / denominator EU estimates (denominator groups by **`YEAR`** only when no extra domain grouping is applied, matching the rFIA-style **`area()`** split). Without POP tables, uses a simple plot-mean fallback.

##### **`custom_pse`** — custom population summaries

```python
custom_pse(
    db,
    x,
    x_vars,
    x_grp_by=None,
    x_transform=None,
    y=None,
    y_vars=None,
    y_grp_by=None,
    y_transform=None,
    method="TI",
    lambda_=0.5,
    totals=True,
    variance=True,
)
```

- **`x`** / optional **`y`** must include **`PLT_CN`**, **`EVAL_TYP`**, and exactly one of **`TREE_BASIS`** or **`AREA_BASIS`**.
- Condition-level **`AREA_BASIS`** rows should include **`CONDID`** (e.g. from **`area(..., cond_list=True)`**).
- **TI + `AREA_BASIS` (no `TREE_BASIS`)**: uses the **same stratified estimators as `area()`** so **`{x_var}_TOTAL`** matches **`area()`**’s acre-style totals for the same variable and groups. With a denominator **`y`**, ratio variance uses **`ratio_var`** with the summed covariance term (**`fa_cv`**) like **`area()`**, not a zero covariance.
- **`TREE_BASIS`**, or non-TI **`method`**: still aggregated with **`_sum_to_plot`**, then an **unstratified** EU shortcut (mean of plot values × **`AREA_USED`** per estimation unit). Those paths are useful for trees but **do not** reproduce full FIA area expansion for arbitrary condition variables.
- **`totals` / `variance`**: control whether **`_TOTAL`** / **`_VAR`** / **`_SE`** columns are retained in the returned frame (see implementation for exact pruning rules).

##### Old-growth / MOG (mature and old-growth structure)

Implementation package: **`fia_mog`** (sibling of **`fia`**). Auxiliary tables and shapefiles ship under **`mog_auxillary/`** (R **`source.path`**), including the National Master Tree Species list (**`USFStrees`** merge: **`FIA.Code` → `PLANTS.Code`**), **`utility_MTcontDivide`** for Montana eastern vs western subregions, and (for WA/OR) **`utility_R5_PAZ.tif`**, **`utility_NWFPboundary/`**, and optional **`utility_ORcounties/`**. **`geopandas`** and **`rasterio`** are core dependencies so **`mog_condition_scores`** can use those layers automatically when present. California (Pacific Southwest / Table 12) does **not** use the PAZ raster in R; it uses the same NWFP polygon only for **white fir** old-growth splits, which the Python pipeline mirrors via **`pnw_inside_nwfp`** on the condition context.

- **`mog_condition_scores(db, ..., mog_auxiliary_dir=None, use_montana_divide_shapefile=True)`** — condition-level MOG scores; defaults to **`fia_py/mog_auxillary`** for species CSV and MT divide when present.
- **`old_growth_area(...)`** — same kwargs; returns **`MOGAreaResult`** with **`df`** (**`OG_PROP`**), **`forest_area_df`** (**`FOREST_PROP`**), **`cond_mog`**.
- **`fia.mog`** re-exports **`fia_mog.engine`** (**`ConditionContext`**, **`MOGEngine`**, …). MOG scoring and area helpers are on the top-level **`fia`** import (**`mog_condition_scores`**, **`old_growth_area`**, **`MOGAreaResult`**) from **`fia_mog.estimators`**.

**Note (Southwest Region 3, Table 9 SDI%):** Table 9 SDI thresholds are interpreted as **percent of COND.SDI_RMRS** using the same numerator as OFE / ``rFIA-ofe-region3.R``: **SDI_OFE** = Σ``TPA_UNADJ × (DIA/10)^1.6`` on **live** trees with **DIA ≥ 18** in, then **SW_REL_SDI** = ``100 × SDI_OFE / SDI_RMRS``. **`mog_condition_scores`** merges **``SDI_RMRS``** from **COND** (when present). If **``SDI_RMRS``** is missing or not positive, **SW_REL_SDI** is undefined (**NaN**) and SDI-based Table 9 rules do not pass. Diagnostics also emit **``SW_SDI_OFE``** and **``SW_SDI_RMRS``**. Example thresholds: **Ponderosa Pine -- Evergreen Oak** **≥ 26%**; *Ponderosa Pine Forest* and *Ponderosa Pine/Willow* **≥ 57%**; *Mixed Conifer -- Frequent Fire* **≥ 56%**. **`SW_PASS_SDI_GE_26`** (and similar) flag ``SW_REL_SDI`` against fixed cutoffs for summaries.

**Note (Southwest Table 9 QMD, live trees ≥ 10 in):** **`FUNCTION_mapMOG.R`** uses a **non-standard** basal term ``π × (DIA/24) × 2`` (linear in DBH) and, in the QMD ratio, **sample stem density** ``nrow / condition.area`` while the numerator uses **expanded** ``TPA_UNADJ``. **`fia_mog.southwest`** follows **conventional** FIA-style QMD: basal area per tree ``π × (DIA/24)²`` ft² (DBH in inches), basal area per acre ``Σ(TPA_UNADJ × BA_tree)``, then ``QMD = sqrt((BA_ac / ΣTPA_UNADJ) / 0.005454)``—the same structure as typical OFE/MOG R helpers that use ``sum(TPA_UNADJ)`` in the denominator. This **differs from bundled mapMOG** but matches standard BA–TPA–QMD algebra; OG flags use **full precision** (no ``round(..., 1)``).

**Note (`OG_PROP` / `OG_FLAG`):** The MOG vector is **old-growth first, then Table 19 maturity** (same pattern as R). Maturity components are weighted scores that can reach **1.0** when every criterion is met. **`mog_condition_scores`** sets **`OG_FLAG`** from **regional OG rules only** (Table 9 / 11 / 12 / 13–14 as implemented in **`fia_mog.engine`**) via **`MOGEngine.old_growth_flag`**, so a stand that is **mature but not OG** does **not** get **`OG_PROP`** credit. **`MOG_SCORE`** remains the **max** of the full vector (maturity can still drive it to **1.0** for diagnostics).

**If `OG_PROP` totals do not change after updating this logic:** Python may be using a **cached** `fia_mog.estimators` (e.g. Jupyter imported `fia` before editing the repo, or `import fia` bound `old_growth_area` from an older install). **Restart the kernel** or run `import importlib, fia_mog.estimators as e; importlib.reload(e)` (and reload `fia` if you use `from fia import old_growth_area`). Prefer **`from fia_mog.estimators import old_growth_area`** or **`from fia.data_io import …`** / **`from fia.clip import clip_fia`** so **`fia.__init__`** does not pin an old MOG binding; **`scripts/test_mog_area.py`** follows that pattern.

**Note (removed compatibility modules, March 2026):** The stub modules **`fia.mog_crosswalk`** and **`fia.mog_estimators`** were deleted as redundant with **`fia_mog.crosswalk`** and **`fia_mog.estimators`**. **`MOGAreaResult`**, **`mog_condition_scores`**, and **`old_growth_area`** remain on the top-level **`fia`** import; crosswalk helpers (**`classify_region`**, …) live only under **`fia_mog.crosswalk`**. If old code used **`from fia.mog_crosswalk import …`** or **`from fia.mog_estimators import …`**, use those **`fia_mog`** modules (or **`from fia import …`** for the three estimators above). To revisit: restore two small files that **`from fia_mog.crosswalk import *`** and **`from fia_mog.estimators import MOGAreaResult, mog_condition_scores, old_growth_area`** respectively, matching the previous layout.

### Tests

Tests use **`pytest`**.

| Location | What it covers | Network / data |
|----------|----------------|------------------|
| **`tests/test_custom_pse_area_match.py`** | **`custom_pse`** **`PROP_FOREST_TOTAL`** vs **`area`** **`AREA_TOTAL`** by **`ADFORCD`** on bundled **AZ** FIADB | **No network**; requires **`tests/fiadb/AZ/`** CSVs |
| **`tests/test_rfia_montana.py`** | **`tpa`**: tree-list columns + **FIADB-API** check (**`snum=1004`**, ``BA_TOTAL`` within **1%** relative error of API); **`area(cond_list=True)`**; **`custom_pse`**; **`area`** vs **FIADB-API** forest acres (**`snum=2`**, integer-rounded match) | **FIADB**: place a Montana extract under **`tests/fiadb/MT/`** locally (that tree is **`.gitignore`d**—files exceed GitHub’s size limit). **FIADB-API** needs network for API-backed tests |
| **`tests/test_area.py`**, **`tests/test_mog_area.py`**, **`tests/test_mog.py`**, **`tests/validation_figure.py`** | Example / exploratory scripts (not necessarily collected by **`pytest`**) | Local **`tests/fiadb/...`** as referenced in each script |

Run everything:

```bash
pytest
```

Run only tests that avoid the network (example: AZ regression only):

```bash
pytest tests/test_custom_pse_area_match.py
```

Montana integration tests are marked **`@pytest.mark.integration`**; you can select or deselect them with **`pytest -m integration`** or **`-m "not integration"`** (markers are registered in **`pyproject.toml`**).

**Large MOG assets not in git:** **`mog_auxillary/utility_R5_PAZ.tif`** (PAZ raster for Pacific Northwest) is also **`.gitignore`d** for the same reason. Keep it next to the other **`mog_auxillary/`** layers on disk, or use **Git LFS** / release downloads if you need it in a shared checkout.

### Dependencies

Core runtime dependencies are declared in **`pyproject.toml`**:

- **`pandas[pyarrow]>=1.5`** (PyArrow-backed dtypes for large tables)
- **`requests>=2.31`**
- **`tqdm>=4.66`**
- **`geopandas>=0.14`**, **`rasterio>=1.3`** (MOG: Montana divide; PNW PAZ + NWFP + OR counties; Pacific Southwest NWFP for California white fir)

```bash
pip install -e .
```

installs these automatically, including spatial stack for **`mog_condition_scores`** when auxiliary shapefiles and rasters are available.

### Status and roadmap

**Implemented**

- Core data I/O (**`get_fia`**, **`read_fia`**, **`FiaDatabase`**).
- Temporal clipping (**`clip_fia`**, most recent per state).
- TI design helpers (**`get_design_info`**, **`handle_pops`**, **`ratio_var`**, **`combine_mr`**).
- Estimators: **`tpa`**, **`area`**, **`custom_pse`** (stratified TI path for condition **`AREA_BASIS`**; tree and legacy paths for other cases).
- Condition list output with **`AREA_BASIS`** for **`custom_pse`** workflows.
- MOG / old-growth (**`fia_mog`**, **`mog_auxillary`**, re-exported on **`fia`**).
- Tests: bundled **AZ** FIADB regression; **Montana** + FIADB-API integration.

**Planned / future**

- Additional estimators (**`biomass`**, **`carbon`**, **`volume`**, **`diversity`**, …).
- Spatial clipping (**`mask`** support akin to **`clipFIA(mask=...)`**).
- Space–time summaries and plotting helpers.
- Remote / lazy-loading database support similar to **`Remote.FIA.Database`**.

Expect breaking changes until a stable **`1.0.0`** release.
