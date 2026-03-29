## FIA (Python)

`fia_py` is an in-progress Python port of the R package [`rFIA`](https://github.com/doserjef/rFIA), providing a user-friendly interface to the USFS Forest Inventory and Analysis (FIA) Database.

The goal is to mirror the core public API of `rFIA` while adopting Pythonic conventions and the scientific Python stack (e.g., `pandas`).

### Recent updates

- **`area(..., cond_list=True)`** now returns **`AREA_BASIS`** (from `COND.PROP_BASIS`, normalized to `SUBP` / `MACR`) together with **`PLT_CN`**, **`CONDID`**, **`YEAR`**, optional **`grp_by`** columns, and **`PROP_FOREST`**, so condition lists plug directly into **`custom_pse`** with correct non-response adjustment.
- **`custom_pse`** for **`method="TI"`** and condition-level **`AREA_BASIS`** inputs (no **`TREE_BASIS`**) uses the same stratified pipeline as **`area()`** (`_rfia_sum_to_plot` → `_rfia_sum_to_eu` → **`combine_mr`**), so totals such as **`PROP_FOREST_TOTAL`** match **`AREA_TOTAL`** from **`area()`** for the same grouping. Tree-level inputs and non-TI methods still use the legacy plot-mean × **`AREA_USED`** shortcut inside each estimation unit.
- **`get_fia(..., load=False)`** removes downloaded **`.zip`** archives after extracting CSVs into **`dir`**, leaving a clean directory of tables.
- **Old-growth / MOG** logic is implemented in the stand-alone **`fia_mog`** package (aligned with **`mog_auxillary/FUNCTION_mapMOG.R`**). **`mog_condition_scores`**, **`old_growth_area`**, and **`ConditionContext` / `MOGEngine`** remain importable from **`fia`** for compatibility; prefer **`fia_mog`** for new code.

### Install (local, editable)

From the project root:

```bash
pip install -e .
```

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
  - **`tables` is required** in the Python port (full-state ZIP-only workflow is not implemented yet).
  - Merges per-state tables into a single logical table per name (e.g. `PLOT`, `TREE`).
  - Returns a **`FiaDatabase`** (a mapping from table name → `pandas.DataFrame`) when **`load=True`**; returns **`None`** when **`load=False`** after writing CSVs to **`dir`**.
  - When **`load=False`**, **`.zip`** files are deleted after extraction so **`dir`** contains only CSVs.

- **`read_fia(dir, tables=None, common=True)`**
  - Reads FIA CSVs from a local directory into a **`FiaDatabase`**.
  - Supports selecting a specific list of tables (**`tables`**) or a common subset when **`common=True`** (similar to **`readFIA(..., common=TRUE)`** in R).

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
    by_size_class=False,
    by_plot=False,
    tree_list=False,
)
```

- **`land_type`**: `"forest"`, `"timber"`, `"all"` (domain via **`COND_STATUS_CD`**, **`SITECLCD`**, **`RESERVCD`**).
- **`tree_type`**: `"live"`, `"dead"`, `"gs"`, `"all"`.
- **`by_species=True`**: groups by species (**`SPECIES`** table when available; else **`SPCD`**).
- **`by_size_class=True`**: 2-inch **`sizeClass`** groups (odd integers: 1, 3, 5, …).
- **`by_plot=True`**: per-plot **`TPA`**, **`BAA`**, **`PROP_FOREST`**.
- **`tree_list=True`**: per-tree rows with **`PLT_CN`**, **`CONDID`**, **`SUBP`**, **`TREE`**, **`YEAR`**, **`TPA`**, **`BAA`**, **`PROP_FOREST`** (and related columns as implemented).

With POP design tables, uses a **TI design-based pipeline** (stratum → estimation unit → totals and variances): **`TPA`**, **`BAA`**, **`TREE_TOTAL`**, **`BA_TOTAL`**, **`AREA_TOTAL`**, variance / SE columns, and **`nPlots_*`**, **`N`**. Without POP tables, falls back to a simple plot-mean estimator.

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
