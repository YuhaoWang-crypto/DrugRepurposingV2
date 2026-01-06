\
from __future__ import annotations

import re
import time
from typing import Any, Dict, List, Optional

import pandas as pd
import requests

from .utils import requests_get


EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _eutils_get(url: str, params: dict, timeout: int = 60, retries: int = 6, backoff: float = 1.8) -> requests.Response:
    """GET wrapper with retry/backoff.

    In practice, NCBI eUtils can intermittently return transient errors (429/5xx)
    or momentarily drop connections. We retry those cases to make the Streamlit
    app resilient on Streamlit Cloud.
    """

    headers = {
        # A lightweight UA helps some hosting environments avoid being treated
        # as an anonymous/bot-like client.
        "User-Agent": "geo-drug-repurposing-explorer/1.0 (+https://streamlit.io)",
        "Accept": "application/json",
    }

    last_err: Exception | None = None
    for i in range(max(1, int(retries))):
        try:
            r = requests.get(url, params=params, timeout=timeout, headers=headers)
            # Retry common transient statuses.
            if r.status_code in (429, 500, 502, 503, 504):
                raise requests.HTTPError(
                    f"HTTP {r.status_code} for {url}", response=r
                )
            r.raise_for_status()
            return r
        except Exception as e:
            last_err = e
            # Exponential backoff with a small cap.
            time.sleep(min(30.0, backoff ** i))

    assert last_err is not None
    raise last_err


def geo_esearch(term: str, retmax: int = 40) -> List[str]:
    """
    Search GEO DataSets (GDS) database via eutils, return list of internal IDs.
    """
    params = {
        "db": "gds",
        "term": term,
        "retmode": "json",
        "retmax": int(retmax),
    }
    r = _eutils_get(f"{EUTILS}/esearch.fcgi", params=params, timeout=60)
    js = r.json()
    ids = js.get("esearchresult", {}).get("idlist", []) or []
    return [str(x) for x in ids]


def geo_esummary(ids: List[str]) -> List[Dict[str, Any]]:
    if not ids:
        return []

    # Be conservative: chunk IDs to avoid overly long URLs and intermittent 414/4xx.
    chunk_size = 200
    out: List[Dict[str, Any]] = []
    for i in range(0, len(ids), chunk_size):
        chunk = ids[i : i + chunk_size]
        params = {
            "db": "gds",
            "id": ",".join(chunk),
            "retmode": "json",
        }
        r = _eutils_get(f"{EUTILS}/esummary.fcgi", params=params, timeout=60)
        js = r.json()
        result = js.get("result", {})
        for k, v in result.items():
            if k == "uids":
                continue
            if isinstance(v, dict):
                out.append(v)
    return out


def extract_gse_accession(accession: str) -> Optional[str]:
    """
    accession can be like "GSE216834" or "GDSxxxx"; we only keep GSE.
    """
    if not accession:
        return None
    m = re.search(r"(GSE\d+)", accession.upper())
    if m:
        return m.group(1)
    return None


def geo_search_candidates(queries: List[str], retmax_each: int = 40) -> pd.DataFrame:
    """
    Run multiple queries, return a de-duplicated table of GSE hits.
    """
    rows: List[Dict[str, Any]] = []
    errors: List[Dict[str, Any]] = []
    for q in queries:
        try:
            ids = geo_esearch(q, retmax=retmax_each)
            summ = geo_esummary(ids)
        except Exception as e:
            # Don't crash the whole UI on transient GEO/eUtils hiccups.
            errors.append({"query": q, "error": repr(e)})
            continue
        for v in summ:
            acc = extract_gse_accession(str(v.get("accession", "")))
            if not acc:
                continue
            rows.append({
                "accession": acc,
                "title": v.get("title", ""),
                "summary": v.get("summary", ""),
                "gdsType": v.get("gdstype", ""),
                "n_samples": v.get("n_samples", ""),
                "PDAT": v.get("pdat", ""),
                "taxon": v.get("taxon", ""),
                "entryType": v.get("entrytype", ""),
                "query": q,
            })
    if not rows:
        df_empty = pd.DataFrame()
        df_empty.attrs["errors"] = errors
        return df_empty
    df = pd.DataFrame(rows)
    df = df.drop_duplicates(subset=["accession"]).sort_values("PDAT", ascending=False)
    df.attrs["errors"] = errors
    return df
