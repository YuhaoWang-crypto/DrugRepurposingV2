from __future__ import annotations

import gzip
import io
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

import numpy as np

def int_like_fraction(df: pd.DataFrame, max_cells: int = 20000) -> float:
    """Estimate fraction of values that are (almost) integers.

    Used to guess whether a matrix looks like raw counts vs normalized expression.
    Returns 0.0 if no numeric values can be sampled.
    """
    if df is None or df.empty:
        return 0.0
    # sample a small block for speed
    sub = df.iloc[: min(200, df.shape[0]), : min(20, df.shape[1])]
    # convert to numeric; ignore non-numeric
    vals = pd.to_numeric(sub.stack(), errors="coerce").dropna().values
    if vals.size == 0:
        return 0.0
    if vals.size > max_cells:
        # deterministic subsample
        rng = np.random.default_rng(0)
        vals = rng.choice(vals, size=max_cells, replace=False)
    frac = float(np.mean(np.isclose(vals, np.round(vals), atol=1e-6)))
    return frac


from .ftp import download_geo_suppl_file, get_geo_matrix_files, get_geo_suppl_files, rank_bulk_matrix_candidates


def _open_text(path: Path):
    """Open plain text or .gz text with best-effort decoding."""
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8', errors='ignore')
    return open(path, 'rt', encoding='utf-8', errors='ignore')


def read_geo_series_matrix(path: Path) -> pd.DataFrame:
    """Parse a GEO *series matrix* (GSE*_series_matrix.txt(.gz)).

    The table is delimited by:
      !series_matrix_table_begin
      ...tab-delimited table...
      !series_matrix_table_end

    Returns a dataframe with rows=features (ID_REF) and columns=GSM*.
    """
    begin_seen = False
    end_seen = False
    table_lines: List[str] = []

    # Some GEO series_matrix files are malformed (e.g. extra header lines,
    # inconsistent row widths). We:
    #   1) Prefer the official markers.
    #   2) Fall back to detecting a "GSM" header line if markers are missing.
    #   3) Parse with a tolerant CSV reader.
    fallback_mode = False

    with _open_text(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            low = line.lower()

            if not begin_seen:
                if low.startswith("!series_matrix_table_begin"):
                    begin_seen = True
                    fallback_mode = False
                    continue

                # Fallback: detect a likely header line (no marker present).
                # Typical header: "ID_REF\tGSM...".
                if ("gsm" in low) and ("\t" in line or "," in line) and (not low.startswith("!")):
                    if low.startswith("id_ref") or low.startswith("id"):
                        begin_seen = True
                        fallback_mode = True
                        table_lines.append(line)
                continue

            # begin seen
            if low.startswith("!series_matrix_table_end"):
                end_seen = True
                break
            if fallback_mode and low.startswith("!"):
                # In fallback mode, stop if metadata resumes.
                end_seen = True
                break
            table_lines.append(line)

    # If we got here with no table content, return empty (caller will treat as unsuitable).
    if not begin_seen or not table_lines:
        return pd.DataFrame()

    # Choose delimiter based on the header line.
    header = table_lines[0]
    sep = "\t" if header.count("\t") >= header.count(",") else ","

    text = "\n".join(table_lines)
    try:
        df = pd.read_csv(
            io.StringIO(text),
            sep=sep,
            engine="python",
            on_bad_lines="skip",
        )
    except TypeError:
        # Older pandas: no on_bad_lines
        df = pd.read_csv(io.StringIO(text), sep=sep, engine="python")
    except Exception:
        # Last resort: return empty and let the caller try other candidates.
        return pd.DataFrame()

    if df.shape[1] < 2:
        return pd.DataFrame()

    first = df.columns[0]
    df = df.set_index(first)
    df = df.dropna(axis=1, how="all")
    return df


def read_matrix_generic(path: Path) -> pd.DataFrame:
    """Read a generic expression matrix in common tabular formats.

    Supported:
      - .csv/.tsv/.txt with optional gzip
      - GEO series matrix (.series_matrix...)

    Output: rows=gene/probe IDs, columns=samples
    """
    name_low = path.name.lower()
    if 'series_matrix' in name_low:
        return read_geo_series_matrix(path)

    # Excel matrices occasionally appear in /suppl
    if path.suffix.lower() in ('.xlsx', '.xls'):
        # pandas uses openpyxl for .xlsx; include it in requirements.txt
        df = pd.read_excel(path, sheet_name=0)
        # If the first column looks like an ID column, use it as index
        if df.shape[1] > 1:
            df = df.set_index(df.columns[0])
        df = df.dropna(axis=1, how='all')
        return df

    # Heuristic: some series-matrix-like files might start with '!'
    try:
        with _open_text(path) as fh:
            head = fh.readline()
        if head.lstrip().startswith('!'):
            # might be series matrix despite filename
            return read_geo_series_matrix(path)
    except Exception:
        pass

    sep = '\t'
    if path.suffix.lower() == '.csv' or name_low.endswith('.csv.gz'):
        sep = ','

    df = pd.read_csv(path, sep=sep)
    # choose first column as index if it looks like gene id
    if df.shape[1] > 1:
        df = df.set_index(df.columns[0])
    # drop columns that are all NA
    df = df.dropna(axis=1, how='all')
    return df


def download_and_load_matrix(
    gse: str,
    raw_gse_dir: Path,
    timeout: int = 240,
    min_genes: int = 1000,
    min_samples: int = 4,
) -> Tuple[str, Path, pd.DataFrame]:
    """Download and load an expression matrix for a GSE.

    Strategy:
      1) Try GEO /suppl directory (processed matrices are often here)
      2) Try GEO /matrix directory (GSE*_series_matrix.txt.gz lives here)
      3) Try multiple candidate files and keep the first that looks like an actual
         expression matrix (enough genes/samples).

    Returns (mode, path, df_expr).
    """
    raw_gse_dir.mkdir(parents=True, exist_ok=True)

    files: List[str] = []
    errs: List[str] = []

    try:
        files.extend(get_geo_suppl_files(gse, timeout=timeout))
    except Exception as e:
        errs.append(f'suppl: {repr(e)}')

    try:
        files.extend(get_geo_matrix_files(gse, timeout=timeout))
    except Exception as e:
        errs.append(f'matrix: {repr(e)}')

    # de-dup while keeping order
    seen = set()
    files2: List[str] = []
    for f in files:
        if f not in seen:
            seen.add(f)
            files2.append(f)
    files = files2

    if not files:
        msg = 'No GEO files listed in suppl/matrix directories.'
        if errs:
            msg += ' Errors: ' + '; '.join(errs)
        raise RuntimeError(msg)

    candidates = rank_bulk_matrix_candidates(files)
    if not candidates:
        raise RuntimeError('No candidate matrix-like files after ranking.')

    tried: List[str] = []
    for rel_path in candidates[:25]:
        try:
            local = download_geo_suppl_file(rel_path, raw_gse_dir, timeout=timeout)
            df = read_matrix_generic(local)
            # basic sanity checks
            if df.shape[0] >= min_genes and df.shape[1] >= min_samples:
                return 'bulk', local, df
            tried.append(f"{Path(rel_path).name} -> shape={df.shape}")
        except Exception as e:
            tried.append(f"{Path(rel_path).name} -> {repr(e)}")
            continue

    raise RuntimeError(
        'No suitable bulk matrix after trying candidates. Tried: ' + ' | '.join(tried[:12])
    )


def map_columns_to_gsm_by_text(df_expr: pd.DataFrame, df_meta: pd.DataFrame) -> Dict[str, str]:
    """Map expression matrix columns to GSM IDs using sample metadata.

    - If column already contains GSM\d+, extract that.
    - Else, try exact matching to sample titles.

    Returns {old_col -> GSMxxxx}.
    """
    import re

    def _norm(s: str) -> str:
        # Lowercase, strip, collapse whitespace/punctuation for looser matching
        s = str(s).strip().lower()
        s = re.sub(r"\s+", " ", s)
        s = re.sub(r"[^a-z0-9 ]+", "", s)
        return s.strip()

    mapping: Dict[str, str] = {}

    # Build a unique mapping from various metadata text fields to GSM IDs
    text_fields = [c for c in ["title", "source_name_ch1", "characteristics_ch1"] if c in df_meta.columns]
    norm_to_gsms: Dict[str, List[str]] = {}
    raw_to_gsms: Dict[str, List[str]] = {}
    if text_fields:
        for gsm, row in df_meta[text_fields].astype(str).iterrows():
            for field in text_fields:
                val = str(row.get(field, "")).strip()
                if not val:
                    continue
                raw_to_gsms.setdefault(val, []).append(str(gsm))
                nv = _norm(val)
                if nv:
                    norm_to_gsms.setdefault(nv, []).append(str(gsm))

    def _unique_lookup(d: Dict[str, List[str]], key: str) -> str | None:
        gs = d.get(key)
        if not gs:
            return None
        gs = list(dict.fromkeys(gs))  # preserve order, de-dup
        return gs[0] if len(gs) == 1 else None

    for col in df_expr.columns:
        c = str(col)

        # 1) Column already contains GSM id
        if "GSM" in c:
            m = re.search(r"(GSM\d+)", c)
            if m:
                mapping[col] = m.group(1)
                continue

        c_strip = c.strip()

        # 2) Exact raw match to one of the text fields (title/source_name/characteristics)
        gsm = _unique_lookup(raw_to_gsms, c_strip)
        if gsm:
            mapping[col] = gsm
            continue

        # 3) Normalized match
        cn = _norm(c_strip)
        if cn:
            gsm = _unique_lookup(norm_to_gsms, cn)
            if gsm:
                mapping[col] = gsm
                continue

            # 4) Substring match (only if unique and reasonably informative)
            if len(cn) >= 4:
                hits = [k for k in norm_to_gsms.keys() if cn in k]
                if len(hits) == 1:
                    gsm = _unique_lookup(norm_to_gsms, hits[0])
                    if gsm:
                        mapping[col] = gsm

    # keep only mappings that are unique target GSMs
    inv = {}
    out = {}
    for k, v in mapping.items():
        inv.setdefault(v, []).append(k)
    for v, ks in inv.items():
        if len(ks) == 1:
            out[ks[0]] = v
    return out


def align_expr_and_meta(
    df_expr: pd.DataFrame,
    df_meta: pd.DataFrame,
    gse: str | None = None,
    config: dict | None = None,
    **_kwargs,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Align expression columns with sample metadata.

    This function tries multiple strategies because GEO matrices often use
    non-GSM sample names.

    Returns (df_expr_aligned, meta_aligned).

    Important: meta_aligned will include a column 'expr_col' that corresponds to
    the aligned df_expr_aligned columns. Downstream DE should use that instead of
    assuming meta index == expression columns.
    
    Notes:
      - `gse`/`config` are accepted for backward/forward compatibility with
        other pipeline versions (they may be used for logging or per-GSE
        heuristics). This implementation does not require them.
    """

    # Ensure meta index is gsm when possible
    if 'gsm' in df_meta.columns and df_meta.index.name != 'gsm':
        df_meta = df_meta.set_index('gsm')

    # Strategy 1: direct intersection
    common = [c for c in df_expr.columns if c in df_meta.index]

    # Strategy 2: rename columns via text mapping
    if not common:
        mapping = map_columns_to_gsm_by_text(df_expr, df_meta)
        if mapping:
            df_expr2 = df_expr.rename(columns=mapping)
            common = [c for c in df_expr2.columns if c in df_meta.index]
            df_expr = df_expr2

    # Strategy 3: fallback by position (same number of samples)
    if not common:
        if df_expr.shape[1] == df_meta.shape[0]:
            meta2 = df_meta.copy()
            df_expr2 = df_expr.copy()
            df_expr2.columns = meta2.index.tolist()
            meta2 = meta2.copy()
            meta2['expr_col'] = df_expr2.columns.tolist()
            return df_expr2, meta2
        raise RuntimeError('Could not align expression columns with sample metadata.')

    # Subset and reorder
    df_expr2 = df_expr[common].copy()
    meta2 = df_meta.loc[common].copy()
    meta2['expr_col'] = df_expr2.columns.tolist()
    return df_expr2, meta2
