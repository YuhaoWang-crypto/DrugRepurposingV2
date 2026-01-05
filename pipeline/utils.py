from __future__ import annotations

import gzip
import io
import json
import os
import re
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Union

import pandas as pd
import requests


def safe_mkdir(p: Union[str, Path]) -> Path:
    p = Path(p)
    p.mkdir(parents=True, exist_ok=True)
    return p


def write_json(path: Union[str, Path], obj: Any) -> None:
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False))


def read_text_maybe_gz(path: Union[str, Path]) -> str:
    path = Path(path)
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as f:
            return f.read()
    return path.read_text(encoding="utf-8", errors="ignore")


def requests_get(
    url: str,
    timeout: int = 120,
    retries: int = 4,
    backoff: float = 1.7,
    headers: Optional[Dict[str, str]] = None,
    binary: bool = False,
) -> Union[str, bytes]:
    """GET with retries. Returns text (default) or bytes (binary=True)."""
    last = None
    for i in range(retries):
        try:
            r = requests.get(url, timeout=timeout, headers=headers)
            if r.status_code >= 500:
                raise RuntimeError(f"ServerError {r.status_code}: {r.text[:200]}")
            r.raise_for_status()
            return r.content if binary else r.text
        except Exception as e:
            last = e
            time.sleep(backoff ** i)
    raise last


def parse_ftp_dir_listing(html: str) -> List[str]:
    """Parse a simple NCBI ftp web directory listing page and return linked names."""
    # NCBI's directory listing is plain HTML with <a href="filename">filename</a>
    # We only need the href target names.
    hrefs = re.findall(r'href="([^"]+)"', html, flags=re.IGNORECASE)
    out = []
    for h in hrefs:
        h = h.strip()
        if h in ("../", "./"):
            continue
        # The listing sometimes includes absolute or query URLs; keep only basename-ish.
        h = h.split("?")[0]
        # filter anchors and oddities
        if not h or h.startswith("#"):
            continue
        out.append(h)
    # keep order but deduplicate
    seen = set()
    dedup = []
    for x in out:
        if x not in seen:
            seen.add(x)
            dedup.append(x)
    return dedup


def read_table_auto(path: Union[str, Path], nrows: Optional[int] = None) -> pd.DataFrame:
    """Read txt/tsv/csv(.gz) into DataFrame with delimiter sniffing."""
    path = Path(path)
    # pandas can read gz directly
    # delimiter sniff: try tab then comma
    for sep in ("\t", ","):
        try:
            df = pd.read_csv(path, sep=sep, nrows=nrows, low_memory=False)
            if df.shape[1] > 1:
                return df
        except Exception:
            pass
    # fallback: whitespace
    return pd.read_csv(path, sep=r"\s+", nrows=nrows, engine="python", low_memory=False)


def normalize_gene_symbol(s: str) -> str:
    s = (s or "").strip()
    if not s:
        return s
    return s.upper()


def normalize_gse(s: str) -> str:
    s = (s or "").strip().upper()
    if not s.startswith("GSE"):
        return s
    return s
