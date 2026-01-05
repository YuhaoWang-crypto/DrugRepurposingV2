from __future__ import annotations

import re
from pathlib import Path
from typing import List

from .utils import parse_ftp_dir_listing, requests_get

NCBI_FTP = "https://ftp.ncbi.nlm.nih.gov"


def geo_series_dir(gse: str) -> str:
    """Return the GEO FTP directory path for a GSE.

    Correct grouping is:
      GSE216834 -> /geo/series/GSE216nnn/GSE216834
      GSE146111 -> /geo/series/GSE146nnn/GSE146111
      GSE2430   -> /geo/series/GSE2nnn/GSE2430
      GSE999    -> /geo/series/GSEnnn/GSE999
    """
    gse_num = int(re.findall(r"\d+", gse)[0])
    prefix = gse_num // 1000
    prefix_str = str(prefix) if prefix > 0 else ""
    return f"/geo/series/GSE{prefix_str}nnn/{gse}"


def gse_soft_url(gse: str) -> str:
    base = geo_series_dir(gse)
    return NCBI_FTP + base + f"/soft/{gse}_family.soft.gz"


def get_geo_suppl_files(gse: str, timeout: int = 120) -> List[str]:
    base = geo_series_dir(gse)
    url = NCBI_FTP + base + "/suppl/"
    html = requests_get(url, timeout=timeout)
    names = parse_ftp_dir_listing(html)
    names = [n for n in names if not n.endswith("/")]
    return [base + "/suppl/" + n for n in names]


def get_geo_matrix_files(gse: str, timeout: int = 120) -> List[str]:
    base = geo_series_dir(gse)
    url = NCBI_FTP + base + "/matrix/"
    html = requests_get(url, timeout=timeout)
    names = parse_ftp_dir_listing(html)
    names = [n for n in names if not n.endswith("/")]
    return [base + "/matrix/" + n for n in names]


def download_geo_suppl_file(file_rel: str, out_dir: Path, timeout: int = 240) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    url = NCBI_FTP + file_rel
    fname = file_rel.split("/")[-1]
    out_path = out_dir / fname
    if out_path.exists() and out_path.stat().st_size > 0:
        return out_path
    data = requests_get(url, timeout=timeout, binary=True)
    out_path.write_bytes(data)
    return out_path


def _candidate_score(name: str) -> int:
    n = name.lower()
    score = 0

    if "series_matrix" in n:
        score += 120
    if "txi" in n:
        score += 100
    if "counts" in n or "count" in n:
        score += 80
    if "fpkm" in n or "tpm" in n or "cpm" in n:
        score += 60
    if "expr" in n or "expression" in n:
        score += 40

    neg = ["readme", "meta", "metadata", "annotation", "annot", "design", "sample", "pheno", "phenotype", "sra"]
    if any(k in n for k in neg):
        score -= 80

    good_ext = (".txt", ".tsv", ".csv", ".gz")
    if n.endswith(good_ext) or any(n.endswith(e + ".gz") for e in (".txt", ".tsv", ".csv")):
        score += 10

    if n.endswith((".tar", ".tar.gz", ".zip", ".7z")):
        score -= 100

    return score


def rank_bulk_matrix_candidates(files: List[str], top_k: int = 50) -> List[str]:
    scored = []
    for rel in files:
        fname = rel.split("/")[-1]
        s = _candidate_score(fname)
        scored.append((s, rel))
    scored.sort(key=lambda x: x[0], reverse=True)
    ranked = [rel for s, rel in scored if s > 0]
    return ranked[:top_k]


def pick_bulk_matrix_file(files: List[str]) -> str:
    ranked = rank_bulk_matrix_candidates(files, top_k=1)
    if not ranked:
        raise RuntimeError("No candidate bulk matrix file found in GEO directory.")
    return ranked[0]
