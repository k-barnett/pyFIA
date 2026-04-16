from __future__ import annotations

import shutil
import sqlite3
import tempfile
from contextlib import closing
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

import pandas as pd
import requests
from tqdm import tqdm


FIA_BASE_URL = "https://apps.fs.usda.gov/fia/datamart/CSV"

# Default table subset for ``read_fia(..., common=True)`` (CSV or SQLite).
READ_FIA_COMMON_TABLES: Set[str] = {
    "COND",
    "COND_DWM_CALC",
    "INVASIVE_SUBPLOT_SPP",
    "PLOT",
    "POP_ESTN_UNIT",
    "POP_EVAL",
    "POP_EVAL_GRP",
    "POP_EVAL_TYP",
    "POP_PLOT_STRATUM_ASSGN",
    "POP_STRATUM",
    "SUBPLOT",
    "TREE",
    "TREE_GRM_COMPONENT",
    "TREE_GRM_MIDPT",
    "TREE_GRM_BEGIN",
    "SUBP_COND_CHNG_MTRX",
    "SEEDLING",
    "SURVEY",
    "SUBP_COND",
    "P2VEG_SUBP_STRUCTURE",
    "PLOTGEOM",
}

# List of known table names, ported from rFIA::getFIA
KNOWN_TABLES: List[str] = [
    "BOUNDARY",
    "COND_DWM_CALC",
    "COND",
    "COUNTY",
    "DWM_COARSE_WOODY_DEBRIS",
    "DWM_DUFF_LITTER_FUEL",
    "DWM_FINE_WOODY_DEBRIS",
    "DWM_MICROPLOT_FUEL",
    "DWM_RESIDUAL_PILE",
    "DWM_TRANSECT_SEGMENT",
    "DWM_VISIT",
    "GRND_CVR",
    "INVASIVE_SUBPLOT_SPP",
    "LICHEN_LAB",
    "LICHEN_PLOT_SUMMARY",
    "LICHEN_VISIT",
    "OZONE_BIOSITE_SUMMARY",
    "OZONE_PLOT_SUMMARY",
    "OZONE_PLOT",
    "OZONE_SPECIES_SUMMARY",
    "OZONE_VALIDATION",
    "OZONE_VISIT",
    "P2VEG_SUBP_STRUCTURE",
    "P2VEG_SUBPLOT_SPP",
    "PLOT_REGEN",
    "PLOT",
    "PLOTGEOM",
    "PLOTSNAP",
    "POP_ESTN_UNIT",
    "POP_EVAL_ATTRIBUTE",
    "POP_EVAL_GRP",
    "POP_EVAL_TYP",
    "POP_EVAL",
    "POP_PLOT_STRATUM_ASSGN",
    "POP_STRATUM",
    "REF_SPECIES",
    "SEEDLING_REGEN",
    "SEEDLING",
    "SITETREE",
    "SOILS_EROSION",
    "SOILS_LAB",
    "SOILS_SAMPLE_LOC",
    "SOILS_VISIT",
    "SUBP_COND_CHNG_MTRX",
    "SUBP_COND",
    "SUBPLOT_REGEN",
    "SUBPLOT",
    "SURVEY",
    "TREE_GRM_BEGIN",
    "TREE_GRM_COMPONENT",
    "TREE_GRM_ESTN",
    "TREE_GRM_MIDPT",
    "TREE_REGIONAL_BIOMASS",
    "TREE_WOODLAND_STEMS",
    "TREE",
    "VEG_PLOT_SPECIES",
    "VEG_QUADRAT",
    "VEG_SUBPLOT_SPP",
    "VEG_SUBPLOT",
    "VEG_VISIT",
    "CITATION",
    "DIFFERENCE_TEST_PER_ACRE",
    "DIFFERENCE_TEST_TOTALS",
    "FIADB_VERSION",
    "FOREST_TYPE",
    "FOREST_TYPE_GROUP",
    "GRM_TYPE",
    "HABTYP_DESCRIPTION",
    "HABTYP_PUBLICATION",
    "INVASIVE_SPECIES",
    "LICHEN_SPECIES",
    "LICHEN_SPP_COMMENTS",
    "NVCS_HEIRARCHY_STRCT",
    "NVCS_LEVEL_1_CODES",
    "NVCS_LEVEL_2_CODES",
    "NVCS_LEVEL_3_CODES",
    "NVCS_LEVEL_4_CODES",
    "NVCS_LEVEL_5_CODES",
    "NVCS_LEVEL_6_CODES",
    "NVCS_LEVEL_7_CODES",
    "NVCS_LEVEL_8_CODES",
    "OWNGRPCD",
    "PLANT_DICTIONARY",
    "POP_ATTRIBUTE",
    "POP_EVAL_TYP_DESCR",
    "RESEARCH_STATION",
    "SPECIES",
    "SPECIES_GROUP",
    "STATE_ELEV",
    "UNIT",
    "FVS_PLOTINIT_PLOT",
    "FVS_GROUPADDFILESANDKEYWORDS",
    "FVS_STANDINIT_COND",
    "FVS_TREEINIT_PLOT",
    "FVS_STANDINIT_PLOT",
    "FVS_TREEINIT_COND",
]

ALL_STATES: List[str] = [
    "AL",
    "AK",
    "AZ",
    "AR",
    "CA",
    "CO",
    "CT",
    "DE",
    "FL",
    "GA",
    "HI",
    "ID",
    "IL",
    "IN",
    "IA",
    "KS",
    "KY",
    "LA",
    "ME",
    "MD",
    "MA",
    "MI",
    "MN",
    "MS",
    "MO",
    "MT",
    "NE",
    "NV",
    "NH",
    "NJ",
    "NM",
    "NY",
    "NC",
    "ND",
    "OH",
    "OK",
    "OR",
    "PA",
    "RI",
    "SC",
    "SD",
    "TN",
    "TX",
    "UT",
    "VT",
    "VA",
    "WA",
    "WV",
    "WI",
    "WY",
    "AS",
    "FM",
    "GU",
    "MP",
    "PW",
    "PR",
    "VI",
    "ENTIRE",
    "REF",
]

# USPS / territory abbrev → FIA ``STATECD`` (plot-level filter for SQLite extracts).
# Inverted from ``fia_mog.crosswalk.STATECD_TO_ST_ABBREV`` plus common extras.
ST_ABBREV_TO_STATECD: Dict[str, int] = {
    "AL": 1,
    "AK": 2,
    "AZ": 4,
    "AR": 5,
    "CA": 6,
    "CO": 8,
    "CT": 9,
    "DE": 10,
    "DC": 11,
    "FL": 12,
    "GA": 13,
    "HI": 15,
    "ID": 16,
    "IL": 17,
    "IN": 18,
    "IA": 19,
    "KS": 20,
    "KY": 21,
    "LA": 22,
    "ME": 23,
    "MD": 24,
    "MA": 25,
    "MI": 26,
    "MN": 27,
    "MS": 28,
    "MO": 29,
    "MT": 30,
    "NE": 31,
    "NV": 32,
    "NH": 33,
    "NJ": 34,
    "NM": 35,
    "NY": 36,
    "NC": 37,
    "ND": 38,
    "OH": 39,
    "OK": 40,
    "OR": 41,
    "PA": 42,
    "RI": 44,
    "SC": 45,
    "SD": 46,
    "TN": 47,
    "TX": 48,
    "UT": 49,
    "VT": 50,
    "VA": 51,
    "WA": 53,
    "WV": 54,
    "WI": 55,
    "WY": 56,
    "AS": 60,
    "FM": 64,
    "GU": 66,
    "MP": 69,
    "PW": 70,
    "PR": 72,
    "VI": 78,
}


def _normalize_read_fia_states(
    states: Optional[str | Sequence[str]],
) -> Optional[Tuple[str, ...]]:
    """
    Validate and normalize ``states`` for :func:`read_fia`.

    Returns ``None`` when no filter should be applied.
    """
    if states is None:
        return None
    if isinstance(states, str):
        seq: Sequence[str] = [states]
    else:
        seq = list(states)
    out = tuple(s.strip().upper() for s in seq if str(s).strip())
    if not out:
        return None
    bad = [s for s in out if s not in ALL_STATES]
    if bad:
        raise ValueError(
            "read_fia: unrecognized state(s): "
            + ", ".join(bad)
            + ". Use two-letter USPS-style codes (e.g. 'MT', 'AZ') or 'REF' / 'ENTIRE' where applicable."
        )
    return out


def _csv_stem_state_prefix(stem: str) -> Optional[str]:
    """Return the Datamart-style state prefix (e.g. ``MT``) or ``None`` if unprefixed."""
    parts = stem.split("_", 1)
    if len(parts) == 2 and len(parts[0]) in (2, 3):
        return parts[0].upper()
    return None


def _sqlite_table_columns_upper(conn: sqlite3.Connection, table_real: str) -> Set[str]:
    cur = conn.execute(f"PRAGMA table_info({_sqlite_quote_ident(table_real)})")
    return {str(row[1]).upper() for row in cur.fetchall()}


def _read_sqlite_chunked_pltcn(
    conn: sqlite3.Connection,
    quoted_table: str,
    plt_cns: Sequence[int],
    chunk_size: int = 500,
) -> pd.DataFrame:
    """Load rows with ``PLT_CN`` in ``plt_cns`` (SQLite parameter limit–safe)."""
    ids = [int(x) for x in plt_cns]
    if not ids:
        probe = pd.read_sql_query(
            f"SELECT * FROM {quoted_table} WHERE 0",
            conn,
            dtype_backend="pyarrow",
        )
        return probe.iloc[0:0]
    frames: List[pd.DataFrame] = []
    for i in range(0, len(ids), chunk_size):
        chunk = ids[i : i + chunk_size]
        ph = ",".join("?" * len(chunk))
        frames.append(
            pd.read_sql_query(
                f"SELECT * FROM {quoted_table} WHERE PLT_CN IN ({ph})",
                conn,
                params=chunk,
                dtype_backend="pyarrow",
            )
        )
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


@dataclass
class FiaDatabase(Mapping[str, pd.DataFrame]):
    """
    Minimal in-memory representation of an FIA database.

    This mirrors the R `FIA.Database` class conceptually: a mapping
    from uppercase table names (e.g., 'PLOT', 'TREE') to pandas
    DataFrames.
    """

    tables: Dict[str, pd.DataFrame]

    def __getitem__(self, key: str) -> pd.DataFrame:
        return self.tables[key.upper()]

    def __iter__(self):
        return iter(self.tables)

    def __len__(self) -> int:
        return len(self.tables)

    def keys(self) -> Iterable[str]:  # type: ignore[override]
        return self.tables.keys()


def normalize_fiadb_dataframe_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Uppercase all column names so SQLite FIADB and CSV exports share the same
    identifiers (e.g. ``plt_cn`` → ``PLT_CN``).
    """
    out = df.copy()
    out.columns = pd.Index([str(c).upper() for c in out.columns])
    return out


def ensure_plot_plt_cn(plot_df: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure ``PLOT`` has a ``PLT_CN`` column for joins with COND/TREE.

    Official FIADB uses ``PLT_CN``; some SQLite extracts only expose ``CN`` as
    the plot control number (same value). :func:`clip_fia` applies the same rule.
    """
    plot_df = plot_df.copy()
    if "PLT_CN" not in plot_df.columns and "CN" in plot_df.columns:
        plot_df["PLT_CN"] = plot_df["CN"]
    return plot_df


def _ensure_dir(path: Path) -> None:
    """
    Ensure ``path`` exists as a directory, creating it (and parents) if needed.

    Raises
    ------
    NotADirectoryError
        If ``path`` exists but is not a directory.
    """
    if path.exists():
        if not path.is_dir():
            raise NotADirectoryError(f"Path exists but is not a directory: {path}")
        return
    path.mkdir(parents=True, exist_ok=True)


def _download_file(url: str, dest: Path, timeout: int = 3600) -> None:
    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        total = int(r.headers.get("Content-Length", 0))
        with open(dest, "wb") as f, tqdm(
            total=total, unit="B", unit_scale=True, desc=dest.name
        ) as pbar:
            for chunk in r.iter_content(chunk_size=8192):
                if not chunk:
                    continue
                f.write(chunk)
                pbar.update(len(chunk))


def get_fia(
    states: str | Sequence[str],
    dir: Optional[str] = None,
    common: bool = True,
    tables: Optional[str | Sequence[str]] = None,
    load: bool = True,
) -> Optional[FiaDatabase]:
    """
    Download FIA data from the FIA Datamart, loosely mirroring rFIA::getFIA.

    Parameters
    ----------
    states:
        State / territory codes. Pass a single string (e.g. ``"AZ"``) or a
        sequence (e.g. ``["AZ", "NM"]``). A bare string is **not** split into
        characters — this matches ``rFIA::getFIA(states = "AZ")``.
    dir:
        Destination folder for downloaded archives and extracted CSVs. If this
        path does not exist, it is created (including parent directories). If
        ``dir`` is omitted and ``load=True``, a temporary directory is used.
    tables:
        Table names to download. Pass one name as a string or several as a
        sequence (e.g. ``\"PLOT\"`` or ``[\"PLOT\", \"TREE\"]``).

    Differences from the R implementation:
    - Only the "download specific tables" pathway is currently implemented
      (i.e., `tables` must be provided).
    - Parallelism / nCores is not yet implemented; downloads are sequential.
    - Error messages are simplified but follow the same logic.
    - If the expected extracted CSV for a state/table pair already exists in
      ``dir``, that download is skipped and a short message is printed.
    """

    if dir is None and not load:
        raise ValueError('Must specify `dir` when `load=False`.')

    if isinstance(states, str):
        states = [states]
    else:
        states = list(states)
    states = [s.upper() for s in states]

    if "REF" in states and len(set(states)) > 1:
        raise ValueError('Download reference tables and state subsets separately.')

    if any(s not in ALL_STATES for s in states):
        missing = [s for s in states if s not in ALL_STATES]
        raise ValueError(
            "Data unavailable for: "
            + ", ".join(missing)
            + ". Use state/territory abbreviations like 'AL', 'CT', 'PR', etc."
        )

    if tables is None:
        raise NotImplementedError(
            "`get_fia` in Python currently requires `tables` to be specified. "
            "The full-state ZIP workflow has not been ported yet."
        )

    if isinstance(tables, str):
        tables = [tables]
    else:
        tables = list(tables)
    tables = [t.upper() for t in tables]
    if any(t not in KNOWN_TABLES for t in tables):
        missing = [t for t in tables if t not in KNOWN_TABLES]
        raise ValueError(
            "Tables unavailable or unrecognized: "
            + ", ".join(missing)
            + ". Check the FIA Datamart for valid table names."
        )

    if dir is None:
        dest_dir = Path(tempfile.mkdtemp(prefix="rfia_"))
    else:
        dest_dir = Path(dir).expanduser().absolute()
        _ensure_dir(dest_dir)

    result_tables: Dict[str, pd.DataFrame] = {}

    # ENTIRE naming convention uses just TABLE.zip; otherwise STATE_TABLE.zip.
    state_prefixes: List[str]
    if "ENTIRE" in states:
        state_prefixes = [""]
    else:
        state_prefixes = [f"{s}_" for s in states]

    for s_prefix in state_prefixes:
        for t in tables:
            zip_name = f"{s_prefix}{t}.zip"
            url = f"{FIA_BASE_URL}/{zip_name}"
            zip_path = dest_dir / zip_name
            # Datamart CSV basename matches the zip stem (e.g. AZ_PLOT.zip → AZ_PLOT.csv).
            expected_csv = dest_dir / f"{s_prefix}{t}.csv"

            if expected_csv.is_file():
                print(
                    f"get_fia: skip download (already present): {expected_csv}",
                    flush=True,
                )
                if load:
                    df = pd.read_csv(expected_csv, dtype_backend="pyarrow")
                    result_tables[expected_csv.stem] = df
                continue

            _download_file(url, zip_path)

            # Extract the CSV inside
            with tempfile.TemporaryDirectory() as tmpdir:
                tmpdir_path = Path(tmpdir)
                shutil.unpack_archive(str(zip_path), str(tmpdir_path))
                # Assume exactly one CSV inside, matching the rFIA convention
                csv_files = list(tmpdir_path.glob("*.csv"))
                if not csv_files:
                    raise RuntimeError(f"No CSV found in archive {zip_name}")
                csv_path = csv_files[0]

                # Optionally copy the CSV to dest_dir if user requested dir
                if dir is not None:
                    target_path = dest_dir / csv_path.name
                    shutil.move(str(csv_path), target_path)
                    csv_path = target_path

                if load:
                    df = pd.read_csv(csv_path, dtype_backend="pyarrow")
                    key = csv_path.stem  # e.g. 'AL_PLOT'
                    result_tables[key] = df

            # When load=False, we do not need to retain the ZIP archives on disk;
            # remove them so get_fia leaves only the extracted CSVs in `dir`.
            if not load and zip_path.exists():
                zip_path.unlink()

    if not load:
        return None

    # Merge across states into a single table per logical name, similar to rFIA
    merged: Dict[str, pd.DataFrame] = {}
    for key, df in result_tables.items():
        # Drop state prefix if present: 'AL_PLOT' -> 'PLOT'
        name_parts = key.split("_", 1)
        logical_name = name_parts[1] if len(name_parts) == 2 else name_parts[0]
        logical_name = logical_name.upper()
        if logical_name in merged:
            merged[logical_name] = pd.concat([merged[logical_name], df], ignore_index=True)
        else:
            merged[logical_name] = df

    return FiaDatabase(merged)


def _sqlite_quote_ident(name: str) -> str:
    """Double-quote a SQLite identifier (table or view name)."""
    return '"' + name.replace('"', '""') + '"'


def _read_fia_sqlite(
    path: Path,
    tables: Optional[Sequence[str]],
    common: bool,
    states: Optional[str | Sequence[str]] = None,
) -> FiaDatabase:
    """Load selected FIADB tables from a SQLite file into a :class:`FiaDatabase`."""
    want_states_t = _normalize_read_fia_states(states)
    geo_states = [s for s in (want_states_t or ()) if s not in ("REF", "ENTIRE")]
    use_geo = bool(geo_states)
    statecds: List[int] = []
    if use_geo:
        missing = [s for s in geo_states if s not in ST_ABBREV_TO_STATECD]
        if missing:
            raise ValueError(
                "read_fia SQLite state filter: no FIA STATECD mapping for: "
                + ", ".join(missing)
            )
        statecds = [ST_ABBREV_TO_STATECD[s] for s in geo_states]

    with closing(sqlite3.connect(str(path))) as conn:
        cur = conn.cursor()
        cur.execute(
            "SELECT name FROM sqlite_master WHERE type IN ('table', 'view') "
            "AND name NOT LIKE 'sqlite_%' ORDER BY name"
        )
        raw_names = [row[0] for row in cur.fetchall()]

        by_upper: Dict[str, str] = {}
        for n in raw_names:
            u = n.upper()
            if u not in by_upper:
                by_upper[u] = n

        if tables is not None:
            wanted = {t.upper() for t in tables}
        elif common:
            wanted = set(READ_FIA_COMMON_TABLES)
        else:
            wanted = set(by_upper.keys()).intersection({t.upper() for t in KNOWN_TABLES})

        plot_filtered: Optional[pd.DataFrame] = None
        plot_cn_ids: Optional[List[int]] = None
        if use_geo:
            if "PLOT" not in wanted:
                raise ValueError(
                    "read_fia(..., states=...) with a SQLite FIADB requires the PLOT table "
                    "so rows can be restricted by STATECD. Add 'PLOT' to `tables` or use common=True."
                )
            real_plot = by_upper["PLOT"]
            qplot = _sqlite_quote_ident(real_plot)
            cols_p = _sqlite_table_columns_upper(conn, real_plot)
            if "STATECD" not in cols_p:
                raise ValueError(
                    "read_fia(..., states=...): PLOT has no STATECD column in this database; "
                    "use Datamart CSV extracts or a standard FIADB SQLite file."
                )
            ph_cd = ",".join("?" * len(statecds))
            plot_filtered = pd.read_sql_query(
                f"SELECT * FROM {qplot} WHERE STATECD IN ({ph_cd})",
                conn,
                params=statecds,
                dtype_backend="pyarrow",
            )
            plot_filtered = normalize_fiadb_dataframe_columns(plot_filtered)
            if "CN" in plot_filtered.columns:
                ser = pd.to_numeric(plot_filtered["CN"], errors="coerce").dropna()
            elif "PLT_CN" in plot_filtered.columns:
                ser = pd.to_numeric(plot_filtered["PLT_CN"], errors="coerce").dropna()
            else:
                raise ValueError(
                    "read_fia(..., states=...): PLOT has neither CN nor PLT_CN; cannot link TREE/COND."
                )
            plot_cn_ids = ser.astype("int64").tolist()

        tables_dict: Dict[str, pd.DataFrame] = {}
        for logical in sorted(wanted):
            if logical not in by_upper:
                continue
            real = by_upper[logical]
            q = _sqlite_quote_ident(real)
            cols = _sqlite_table_columns_upper(conn, real)

            if use_geo and logical == "PLOT":
                df = plot_filtered  # type: ignore[assignment]
            elif use_geo and "PLT_CN" in cols and logical != "PLOT":
                df = _read_sqlite_chunked_pltcn(conn, q, plot_cn_ids or [])
            elif use_geo and "STATECD" in cols:
                ph = ",".join("?" * len(statecds))
                df = pd.read_sql_query(
                    f"SELECT * FROM {q} WHERE STATECD IN ({ph})",
                    conn,
                    params=statecds,
                    dtype_backend="pyarrow",
                )
            else:
                df = pd.read_sql_query(f"SELECT * FROM {q}", conn, dtype_backend="pyarrow")

            df = normalize_fiadb_dataframe_columns(df)
            if logical in tables_dict:
                tables_dict[logical] = pd.concat(
                    [tables_dict[logical], df], ignore_index=True
                )
            else:
                tables_dict[logical] = df

    return FiaDatabase(tables_dict)


def read_fia(
    dir: str,
    tables: Optional[Sequence[str]] = None,
    common: bool = True,
    states: Optional[str | Sequence[str]] = None,
) -> FiaDatabase:
    """
    Read FIA data into a :class:`FiaDatabase` from CSV files or a SQLite FIADB.

    This is a partial port of `rFIA::readFIA` for the in-memory case.

    Parameters
    ----------
    dir:
        Path to a **directory** of FIA CSV files (Datamart extract), **or** a path
        to a **SQLite** FIADB file (``.sqlite``, ``.db``, etc.). If ``dir`` is a
        file, it is opened with :mod:`sqlite3`; if that fails, a :exc:`ValueError`
        is raised.
    tables:
        Optional sequence of table base names (e.g. ['PLOT', 'TREE']). If
        omitted and `common=True`, a default "common" subset will be loaded.
    common:
        If True and `tables` is None, load a common subset of frequently used
        tables (PLOT, TREE, COND, etc.), mirroring the R package behavior.
    states:
        Optional state / territory codes (e.g. ``"MT"`` or ``["AZ", "NM"]``).

        **CSV directory:** only reads files whose Datamart prefix matches one of
        the codes (e.g. ``MT_PLOT.csv`` when ``states=["MT"]``). Unprefixed
        national files such as ``PLOT.csv`` are included only if ``"ENTIRE"``
        appears in ``states``. Reference-style prefixes such as ``REF_SPECIES``
        are included only if ``"REF"`` is listed.

        **SQLite FIADB:** when you pass real geographic codes (not ``REF`` /
        ``ENTIRE`` alone), loads **PLOT** rows with ``STATECD`` in the requested
        states, then restricts other tables that have ``PLT_CN`` to those plot
        control numbers (and applies ``STATECD`` filtering on tables that have
        that column but not ``PLT_CN``). You must include **PLOT** in ``tables``
        (or use ``common=True``). ``REF`` / ``ENTIRE`` do not apply row filtering
        on SQLite; pass ``states=None`` to load the full file.
    """

    path = Path(dir).expanduser().absolute()
    if not path.exists():
        raise FileNotFoundError(f"Path does not exist: {path}")

    if path.is_file():
        try:
            return _read_fia_sqlite(path, tables=tables, common=common, states=states)
        except sqlite3.DatabaseError as e:
            raise ValueError(
                f"Could not read SQLite FIADB at {path}: {e}"
            ) from e

    if not path.is_dir():
        raise NotADirectoryError(f"Expected a directory or SQLite file: {path}")

    want_states_t = _normalize_read_fia_states(states)
    want_state_set: Optional[Set[str]] = set(want_states_t) if want_states_t else None

    csv_files = sorted(p for p in path.iterdir() if p.suffix.lower() == ".csv")
    if not csv_files:
        raise FileNotFoundError(f"Directory {path} contains no .csv files.")

    all_tables = {p.stem for p in csv_files}

    # Derive base names (drop state prefixes like 'AL_', 'REF_')
    base_names: List[str] = []
    for name in all_tables:
        parts = name.split("_", 1)
        if len(parts) == 2 and len(parts[0]) in (2, 3):
            base_names.append(parts[1])
        else:
            base_names.append(name)

    base_names_set = {b.upper() for b in base_names}

    if tables is not None:
        wanted = {t.upper() for t in tables}
    elif common:
        wanted = set(READ_FIA_COMMON_TABLES)
    else:
        wanted = base_names_set.intersection({t.upper() for t in KNOWN_TABLES})

    selected_files: List[Path] = []
    for p in csv_files:
        stem = p.stem
        parts = stem.split("_", 1)
        if len(parts) == 2 and len(parts[0]) in (2, 3):
            logical = parts[1].upper()
        else:
            logical = stem.upper()
        if logical not in wanted:
            continue
        if want_state_set is not None:
            pref = _csv_stem_state_prefix(stem)
            if pref is None:
                if "ENTIRE" not in want_state_set:
                    continue
            elif pref not in want_state_set:
                continue
        selected_files.append(p)

    tables_dict: Dict[str, pd.DataFrame] = {}
    for p in selected_files:
        df = pd.read_csv(p, dtype_backend="pyarrow")
        stem = p.stem
        parts = stem.split("_", 1)
        if len(parts) == 2 and len(parts[0]) in (2, 3):
            logical = parts[1].upper()
        else:
            logical = stem.upper()

        if logical in tables_dict:
            tables_dict[logical] = pd.concat([tables_dict[logical], df], ignore_index=True)
        else:
            tables_dict[logical] = df

    return FiaDatabase(tables_dict)

