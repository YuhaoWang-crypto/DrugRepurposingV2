from __future__ import annotations

import gzip
import json
import re
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import pandas as pd
import requests


# -------------------------
# FS helpers
# -------------------------
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


# -------------------------
# HTTP helpers
# -------------------------
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


def requests_post_json(
    url: str,
    payload: Dict[str, Any],
    timeout: int = 120,
    retries: int = 4,
    backoff: float = 1.7,
    headers: Optional[Dict[str, str]] = None,
) -> Dict[str, Any]:
    """POST json with retries; returns parsed json."""
    hdrs = {"content-type": "application/json"}
    if headers:
        hdrs.update(headers)

    last = None
    for i in range(retries):
        try:
            r = requests.post(url, json=payload, timeout=timeout, headers=hdrs)
            if r.status_code >= 500:
                raise RuntimeError(f"ServerError {r.status_code}: {r.text[:200]}")
            r.raise_for_status()
            return r.json()
        except Exception as e:
            last = e
            time.sleep(backoff ** i)
    raise last


def download_file(
    url: str,
    out_path_or_dir: Union[str, Path],
    fname: Optional[str] = None,
    timeout: int = 240,
    retries: int = 4,
    backoff: float = 1.7,
    overwrite: bool = False,
) -> Path:
    """
    Download url to a local file.

    Compatible with multiple call styles:
      - download_file(url, "/tmp/a.txt.gz")
      - download_file(url, "/tmp/")  (auto filename from url)
      - download_file(url, "/tmp", fname="a.txt.gz")
    """
    out = Path(out_path_or_dir)

    # decide if out is a dir-like or file-like target
    out_str = str(out_path_or_dir)
    dir_like = False

    if fname is not None:
        dir_like = True
    elif out_str.endswith("/"):
        dir_like = True
    elif out.exists() and out.is_dir():
        dir_like = True
    else:
        # heuristic: no suffix => likely a directory (e.g. "raw/GSE216834")
        # but be careful: some dirs may contain dots; still ok for our use
        if out.suffix == "":
            dir_like = True

    if dir_like:
        if fname is None:
            fname = url.split("/")[-1].split("?")[0]
        out = out / fname

    out.parent.mkdir(parents=True, exist_ok=True)

    if out.exists() and out.stat().st_size > 0 and (not overwrite):
        return out

    data = requests_get(url, timeout=timeout, retries=retries, backoff=backoff, binary=True)
    out.write_bytes(data)
    return out


# -------------------------
# GEO FTP dir listing parser
# -------------------------
def parse_ftp_dir_listing(html: str) -> List[str]:
    """Parse a simple NCBI ftp web directory listing page and return linked names."""
    hrefs = re.findall(r'href="([^"]+)"', html, flags=re.IGNORECASE)
    out = []
    for h in hrefs:
        h = (h or "").strip()
        if not h or h in ("../", "./") or h.startswith("#"):
            continue
        h = h.split("?")[0]
        out.append(h)

    # de-dup preserve order
    seen = set()
    dedup = []
    for x in out:
        if x not in seen:
            seen.add(x)
            dedup.append(x)
    return dedup


# -------------------------
# Light data helpers
# -------------------------
def read_table_auto(path: Union[str, Path], nrows: Optional[int] = None) -> pd.DataFrame:
    """Read txt/tsv/csv(.gz) into DataFrame with delimiter sniffing."""
    path = Path(path)
    for sep in ("\t", ","):
        try:
            df = pd.read_csv(path, sep=sep, nrows=nrows, low_memory=False)
            if df.shape[1] > 1:
                return df
        except Exception:
            pass
    return pd.read_csv(path, sep=r"\s+", nrows=nrows, engine="python", low_memory=False)


def normalize_gene_symbol(s: str) -> str:
    s = (s or "").strip()
    return s.upper() if s else s


def normalize_gse(s: str) -> str:
    s = (s or "").strip().upper()
    return s
