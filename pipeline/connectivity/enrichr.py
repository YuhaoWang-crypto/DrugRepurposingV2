from __future__ import annotations

import time
from typing import Dict, List

import requests


def _clean_gene_symbols(genes: List[str]) -> List[str]:
    """Sanitize gene symbols for external APIs."""
    out: List[str] = []
    for g in genes or []:
        if g is None:
            continue
        s = str(g).strip()
        if not s or s.lower() == "nan":
            continue
        # Common junk tokens
        if s in {"NA", "N/A", "None", "NULL"}:
            continue
        out.append(s)
    # Keep order but de-duplicate
    seen = set()
    uniq: List[str] = []
    for s in out:
        if s not in seen:
            uniq.append(s)
            seen.add(s)
    return uniq


def enrichr_add_list(base_url: str, genes: List[str], description: str = "signature") -> int:
    """Create an Enrichr list and return userListId.

    Enrichr can be finicky about content-type; multipart/form-data (requests `files=`)
    is the most compatible.
    """
    genes = _clean_gene_symbols(genes)
    if not genes:
        raise ValueError("Empty gene list for Enrichr")

    list_str = "\n".join(genes)

    # Try both historical and newer paths.
    endpoints = [
        f"{base_url.rstrip('/')}/addList",
        f"{base_url.rstrip('/')}/api/addList",
    ]

    last_err = None
    for url in endpoints:
        try:
            payload = {
                "list": (None, list_str),
                "description": (None, description),
            }
            r = requests.post(url, files=payload, timeout=120, headers={"User-Agent": "Mozilla/5.0"})
            if r.status_code >= 500:
                raise RuntimeError(f"Enrichr server error {r.status_code}: {(r.text or '')[:200]}")
            r.raise_for_status()
            js = r.json()
            # Typical response contains `userListId`
            if "userListId" not in js:
                raise RuntimeError(f"Enrichr addList response missing userListId: {str(js)[:200]}")
            return int(js["userListId"])
        except Exception as e:
            last_err = e
            time.sleep(0.5)
            continue

    raise RuntimeError(f"Enrichr addList failed on all endpoints. Last error: {repr(last_err)}")


def enrichr_enrich(base_url: str, user_list_id: int, library: str) -> Dict:
    url = f"{base_url.rstrip('/')}/enrich"
    r = requests.post(url, data={"userListId": user_list_id, "backgroundType": library}, timeout=120)
    r.raise_for_status()
    return r.json()


def enrichr_run_dual(
    base_url: str,
    up: List[str],
    dn: List[str],
    libraries: List[str],
    desc_prefix: str | None = None,
) -> Dict[str, Dict]:
    """Run Enrichr for UP and DN lists across libraries.

    Note: some callers may want to tag submissions with a dataset prefix (e.g. GSE).
    The Enrichr API supports a short 'description' field for uploaded gene sets.
    """
    results: Dict[str, Dict] = {}
    up_desc = f"{desc_prefix} UP" if desc_prefix else "UP"
    dn_desc = f"{desc_prefix} DN" if desc_prefix else "DN"
    up_id = enrichr_add_list(base_url, up, description=up_desc)
    dn_id = enrichr_add_list(base_url, dn, description=dn_desc)
    for lib in libraries:
        results[f"UP::{lib}"] = enrichr_enrich(base_url, up_id, lib)
        results[f"DN::{lib}"] = enrichr_enrich(base_url, dn_id, lib)
    return results
