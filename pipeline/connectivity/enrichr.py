"""Enrichr connectivity utilities.

This module provides a small Enrichr client and a helper `enrichr_run_dual`
used by the main pipeline.

The pipeline expects **two** pandas DataFrames:

1) df_long: long-form enrichment results for each queried library and each
   direction (UP / DN query lists).
2) df_rev:  an aggregated *reversal* table with columns:
      - term (drug / perturbation signature)
      - enrichr_reversal_score (higher = stronger reversal support)
      - n_support (how many library-direction hits contributed)

Reversal scoring logic
----------------------
Enrichr libraries used here are typically of the form:
  - *_up   : genes UP-regulated by a perturbation (drug)
  - *_down : genes DOWN-regulated by a perturbation

Given a disease signature split into UP and DN genes (disease vs control),
we treat a perturbation as *reversal-supporting* when:
  - disease UP genes enrich in drug DOWN gene sets
  - disease DN genes enrich in drug UP gene sets

We use Enrichr's "combined score" as the per-hit contribution.
"""

from __future__ import annotations

import re
from typing import Any, Dict, Iterable, List, Optional, Tuple

import pandas as pd
import requests


def _clean_gene_list(genes: Iterable[str], max_n: int = 5000) -> List[str]:
    """Return a de-duplicated, sanitized gene list for Enrichr."""

    out: List[str] = []
    seen = set()
    for g in genes:
        if g is None:
            continue
        s = str(g).strip()
        if not s:
            continue
        # Enrichr expects one gene per line. Drop anything with newlines.
        s = s.replace("\n", " ").replace("\r", " ").strip()
        if not s:
            continue
        # Common messy identifiers: remove version suffixes (ENSG...\.12)
        s = re.sub(r"\.(\d+)$", "", s)
        # Upper-case symbols (safe for Enrichr; Ensembl IDs are already upper)
        s = s.upper()
        if s in seen:
            continue
        seen.add(s)
        out.append(s)
        if len(out) >= max_n:
            break
    return out


def _post_add_list(base_url: str, genes: List[str], description: str, timeout: int = 120) -> int:
    """POST gene list to Enrichr and return userListId.

    Enrichr historically supported both:
      - /addList
      - /api/addList
    depending on deployment.
    """

    if not genes:
        raise RuntimeError("Empty gene list provided to Enrichr addList")

    list_str = "\n".join(genes)

    # Try the modern canonical endpoint first.
    urls = [base_url.rstrip("/") + "/addList", base_url.rstrip("/") + "/api/addList"]
    last_err: Optional[Exception] = None
    for url in urls:
        try:
            r = requests.post(
                url,
                files={
                    # NOTE: Enrichr expects multipart/form-data
                    "list": (None, list_str),
                    "description": (None, description),
                },
                timeout=timeout,
            )
            r.raise_for_status()
            js = r.json()
            if "userListId" not in js:
                raise RuntimeError(f"Enrichr addList missing userListId: keys={list(js.keys())}")
            return int(js["userListId"])
        except Exception as e:  # noqa: BLE001
            last_err = e
            continue
    raise RuntimeError(f"Enrichr addList failed on all endpoints. Last error: {last_err}")


def _get_enrich(
    base_url: str,
    user_list_id: int,
    library: str,
    timeout: int = 120,
) -> List[List[Any]]:
    """Call Enrichr /enrich and return the raw results list."""

    url = base_url.rstrip("/") + "/enrich"
    params = {"userListId": str(user_list_id), "backgroundType": library}
    r = requests.get(url, params=params, timeout=timeout)
    r.raise_for_status()
    js = r.json()
    if library not in js:
        # Some deployments return {"results": [...]} — handle that too.
        if "results" in js and isinstance(js["results"], list):
            return js["results"]
        raise RuntimeError(f"Enrichr enrich response missing library key '{library}'. keys={list(js.keys())}")
    return js[library]


def _is_up_library(lib: str) -> bool:
    s = lib.lower()
    return s.endswith("_up") or s.endswith("-up") or "_up_" in s or "pert_up" in s


def _is_down_library(lib: str) -> bool:
    s = lib.lower()
    return s.endswith("_down") or s.endswith("-down") or "_down_" in s or "pert_down" in s


def enrichr_run_dual(
    base_url: str,
    up: List[str],
    dn: List[str],
    libraries: List[str],
    *,
    desc_prefix: str = "SIG",
    top_k: int = 100,
    timeout: int = 120,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Run Enrichr on UP and DN lists across multiple libraries.

    Returns:
      (df_long, df_rev)

    df_long columns:
      - library, direction, rank, term, pval, zscore, combined_score, adj_pval,
        overlap_genes

    df_rev columns:
      - term, enrichr_reversal_score, n_support
    """

    up_clean = _clean_gene_list(up)
    dn_clean = _clean_gene_list(dn)
    if not up_clean and not dn_clean:
        return (
            pd.DataFrame(
                columns=[
                    "library",
                    "direction",
                    "rank",
                    "term",
                    "pval",
                    "zscore",
                    "combined_score",
                    "adj_pval",
                    "overlap_genes",
                ]
            ),
            pd.DataFrame(columns=["term", "enrichr_reversal_score", "n_support"]),
        )

    # Register lists.
    up_id: Optional[int] = None
    dn_id: Optional[int] = None
    if up_clean:
        up_id = _post_add_list(base_url, up_clean, f"{desc_prefix}__UP", timeout=timeout)
    if dn_clean:
        dn_id = _post_add_list(base_url, dn_clean, f"{desc_prefix}__DN", timeout=timeout)

    rows: List[Dict[str, Any]] = []
    for lib in libraries:
        if up_id is not None:
            res = _get_enrich(base_url, up_id, lib, timeout=timeout)
            for rec in res[:top_k]:
                # Enrichr record format (as documented):
                # [rank, term, pval, zscore, combined_score, overlapping_genes,
                #  adjusted_pval, old_pval, old_adjusted_pval]
                rows.append(
                    {
                        "library": lib,
                        "direction": "UP",
                        "rank": rec[0],
                        "term": rec[1],
                        "pval": rec[2],
                        "zscore": rec[3],
                        "combined_score": rec[4],
                        "overlap_genes": rec[5],
                        "adj_pval": rec[6] if len(rec) > 6 else None,
                    }
                )
        if dn_id is not None:
            res = _get_enrich(base_url, dn_id, lib, timeout=timeout)
            for rec in res[:top_k]:
                rows.append(
                    {
                        "library": lib,
                        "direction": "DN",
                        "rank": rec[0],
                        "term": rec[1],
                        "pval": rec[2],
                        "zscore": rec[3],
                        "combined_score": rec[4],
                        "overlap_genes": rec[5],
                        "adj_pval": rec[6] if len(rec) > 6 else None,
                    }
                )

    df_long = pd.DataFrame(rows)
    if df_long.empty:
        return df_long, pd.DataFrame(columns=["term", "enrichr_reversal_score", "n_support"])

    # Compute reversal score.
    def support_score(row: pd.Series) -> float:
        lib = str(row["library"])
        direction = str(row["direction"])  # UP / DN query list
        cs = float(row.get("combined_score") or 0.0)

        # disease UP -> drug DOWN libraries
        if direction == "UP" and _is_down_library(lib):
            return cs
        # disease DN -> drug UP libraries
        if direction == "DN" and _is_up_library(lib):
            return cs
        return 0.0

    df_long["reversal_support"] = df_long.apply(support_score, axis=1)
    df_sup = df_long[df_long["reversal_support"] > 0].copy()
    if df_sup.empty:
        return df_long, pd.DataFrame(columns=["term", "enrichr_reversal_score", "n_support"])

    df_rev = (
        df_sup.groupby("term", as_index=False)
        .agg(enrichr_reversal_score=("reversal_support", "sum"), n_support=("reversal_support", "size"))
        .sort_values(["enrichr_reversal_score", "n_support"], ascending=False)
        .reset_index(drop=True)
    )
    return df_long, df_rev
