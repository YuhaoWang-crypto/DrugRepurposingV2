import json
import os
from pathlib import Path

import pandas as pd
import streamlit as st

from pipeline.config import default_config
from pipeline.geo import geo_search_candidates
from pipeline.soft import download_soft, parse_soft
from pipeline.label import label_conditions
from pipeline.pipeline import run_one_gse


st.set_page_config(page_title="GEO Drug Repurposing Explorer", layout="wide")


def deep_update(base, upd):
    """Recursively merge `upd` into `base`.

    Notes:
      - dicts are merged key-by-key
      - lists/scalars are REPLACED (not extended)
    """
    if isinstance(base, dict) and isinstance(upd, dict):
        out = dict(base)
        for k, v in upd.items():
            if k in out:
                out[k] = deep_update(out[k], v)
            else:
                out[k] = v
        return out
    # replace lists/scalars
    return upd


def _init_state():
    if "config" not in st.session_state:
        st.session_state["config"] = default_config()

    # Keep a separate editor buffer so it is *not* overwritten on reruns.
    if "cfg_editor" not in st.session_state:
        st.session_state["cfg_editor"] = json.dumps(st.session_state["config"], ensure_ascii=False, indent=2)

    if "candidates" not in st.session_state:
        st.session_state["candidates"] = None
    if "geo_queries" not in st.session_state:
        st.session_state["geo_queries"] = []
    if "validation" not in st.session_state:
        st.session_state["validation"] = None
    if "run_results" not in st.session_state:
        st.session_state["run_results"] = []
    if "logs" not in st.session_state:
        st.session_state["logs"] = []


_init_state()


def add_log(msg: str):
    st.session_state["logs"].append(msg)


def _quote_term_if_needed(t: str) -> str:
    t = (t or "").strip()
    if not t:
        return ""
    # If the user already included quotes, keep.
    if '"' in t:
        return t
    # Quote phrases with whitespace
    if any(ch.isspace() for ch in t):
        return f'"{t}"'
    return t


def build_geo_queries(genes, phenotype, extra_terms, exclude_terms=None, max_queries: int = 3):
    """Build GEO search queries.

    Notes:
    - `exclude_terms` is taken from config JSON (geo_exclude_terms).
    - Keep the query list small; GEO's esearch may rate-limit.
    """
    genes = [g.strip() for g in genes if g.strip()]
    phenotype = (phenotype or "").strip()
    extra_terms = (extra_terms or "").strip()

    if genes:
        base = " OR ".join(genes)
        base = f"({base})"
    else:
        base = ""

    if phenotype:
        base = f"({base}) AND ({phenotype})" if base else f"({phenotype})"

    if extra_terms:
        base = f"({base}) AND ({extra_terms})" if base else f"({extra_terms})"

    # Encourage expression datasets
    base = f"({base}) AND (RNA-seq OR microarray OR transcriptome OR expression OR 'single cell' OR 'single nucleus' OR scRNAseq OR snRNAseq)" if base else "(RNA-seq OR microarray OR transcriptome OR expression OR 'single cell' OR 'single nucleus' OR scRNAseq OR snRNAseq)"

    exclude_terms = exclude_terms or []
    exclude_terms = [_quote_term_if_needed(t) for t in exclude_terms]
    exclude_terms = [t for t in exclude_terms if t]
    if exclude_terms:
        ex = " OR ".join(exclude_terms)
        base = f"({base}) NOT ({ex})"

    # Query variants: sometimes adding tissue hints helps
    queries = [base]
    hints = ["human", "mouse", "patient", "biopsy", "cell line"]
    for h in hints:
        if len(queries) >= max_queries:
            break
        queries.append(f"({base}) AND ({h})")

    # Deduplicate while preserving order
    seen = set()
    out = []
    for q in queries:
        if q not in seen:
            out.append(q)
            seen.add(q)
    return out[:max_queries]


def summarize_characteristics_kv(df_meta: pd.DataFrame, top_n: int = 30) -> pd.DataFrame:
    """Summarize which `characteristics_kv` keys exist and how frequent they are.

    Returns a dataframe with columns: key, n_samples, examples.
    """
    if df_meta is None or df_meta.empty or "characteristics_kv" not in df_meta.columns:
        return pd.DataFrame(columns=["key", "n_samples", "examples"])

    from collections import Counter

    cnt: Counter = Counter()
    examples: Dict[str, List[str]] = {}
    for d in df_meta["characteristics_kv"].tolist():
        if not isinstance(d, dict):
            continue
        for k, v in d.items():
            if not k:
                continue
            cnt[k] += 1
            if k not in examples:
                examples[k] = []
            if v and (len(examples[k]) < 3) and (str(v) not in examples[k]):
                examples[k].append(str(v)[:120])

    rows = []
    for k, n in cnt.most_common(top_n):
        rows.append({
            "key": k,
            "n_samples": n,
            "examples": " | ".join(examples.get(k, [])),
        })
    return pd.DataFrame(rows)


def build_characteristics_view(df_meta: pd.DataFrame, keys: List[str]) -> pd.DataFrame:
    """Build a user-friendly preview table for sample annotations."""
    if df_meta is None or df_meta.empty:
        return pd.DataFrame()
    keys = [k for k in keys if k]
    base_cols = [c for c in ["gsm", "condition", "title", "source", "description"] if c in df_meta.columns]
    out = df_meta[base_cols].copy()

    # Add selected characteristic keys as columns
    if "characteristics_kv" in df_meta.columns:
        for k in keys:
            out[k] = df_meta["characteristics_kv"].apply(lambda d: (d or {}).get(k, "") if isinstance(d, dict) else "")

    # Add the raw characteristics block as a readable string
    if "characteristics" in df_meta.columns:
        out["characteristics"] = df_meta["characteristics"].apply(
            lambda xs: "\n".join(xs) if isinstance(xs, list) else ("" if xs is None else str(xs))
        )
    return out


# ---------------- Sidebar ----------------
st.sidebar.header("Inputs")

# These are *UI* inputs; they are not written into config unless you choose to.
genes_in = st.sidebar.text_input("Target gene(s) / protein(s) (comma separated)", value="CLCN2")
phenotype_in = st.sidebar.text_input("Phenotype keywords (optional)", value="lysosome OR endosome")
extra_geo_terms = st.sidebar.text_input("Extra GEO terms (optional)", value="")

st.sidebar.markdown("---")

st.sidebar.subheader("Config (editable JSON)")

# Streamlit quirk: do not mutate st.session_state[...] keys bound to widgets
# directly in the main body after widget creation. Use callbacks.
def _reset_config_callback():
    cfg0 = default_config()
    st.session_state["config"] = cfg0
    # Don't mutate widget-backed state (cfg_editor) directly inside callbacks.
    # Instead, set a pending buffer that we apply BEFORE rendering the widget.
    st.session_state["_cfg_editor_pending"] = json.dumps(cfg0, ensure_ascii=False, indent=2)
    st.session_state.pop("_cfg_error", None)
    add_log("[UI] Reset config JSON to defaults")


def _apply_config_callback():
    try:
        cfg_new = json.loads(st.session_state.get("cfg_editor", "{}"))
        # Merge on top of defaults so missing keys are filled (keeps the pipeline robust)
        cfg_merged = deep_update(default_config(), cfg_new)
        st.session_state["config"] = cfg_merged
        st.session_state["_cfg_editor_pending"] = json.dumps(cfg_merged, ensure_ascii=False, indent=2)
        st.session_state.pop("_cfg_error", None)
        add_log("[UI] Applied config JSON")
    except Exception as e:
        st.session_state["_cfg_error"] = str(e)


colA, colB = st.sidebar.columns(2)
colA.button("Apply config JSON", use_container_width=True, on_click=_apply_config_callback)
colB.button("Reset defaults", use_container_width=True, on_click=_reset_config_callback)

if st.session_state.get("_cfg_error"):
    st.sidebar.error(f"Invalid JSON: {st.session_state['_cfg_error']}")

# Apply any pending editor update BEFORE the widget is created.
if "_cfg_editor_pending" in st.session_state:
    st.session_state["cfg_editor"] = st.session_state.pop("_cfg_editor_pending")

# Editor
st.sidebar.text_area("Config JSON", key="cfg_editor", height=380)

# Make it obvious when the editor buffer differs from applied config.
try:
    applied = json.dumps(st.session_state["config"], ensure_ascii=False, sort_keys=True)
    edited = json.dumps(json.loads(st.session_state["cfg_editor"]), ensure_ascii=False, sort_keys=True)
    if applied != edited:
        st.sidebar.warning("Config JSON has unapplied changes. Click **Apply config JSON** to use them.")
except Exception:
    # editor may not be valid JSON yet
    pass

st.sidebar.markdown("---")

st.sidebar.subheader("Optional PPI / module")
ppi_file = st.sidebar.file_uploader("Upload PPI edges CSV (two columns a,b; gene symbols)", type=["csv"])
if ppi_file is not None:
    # Save to a stable temp path (per-session) to avoid overwriting other runs.
    tmp_dir = Path("./_streamlit_uploads")
    tmp_dir.mkdir(parents=True, exist_ok=True)
    ppi_path = tmp_dir / ("ppi_edges.csv")
    ppi_path.write_bytes(ppi_file.getvalue())
    st.session_state["config"].setdefault("ppi_edges_csv", str(ppi_path))
    st.sidebar.success(f"Uploaded: {ppi_path}")

module_genes = st.sidebar.text_area(
    "Module gene symbols (one per line; used only if you enable PPI edges)",
    value="\n".join(st.session_state["config"].get("module_gene_symbols", [])),
    height=160,
)
if st.sidebar.button("Apply module genes"):
    genes = [g.strip().upper() for g in module_genes.splitlines() if g.strip()]
    st.session_state["config"]["module_gene_symbols"] = genes
    st.sidebar.success(f"Saved {len(genes)} module genes.")


# ---------------- Main UI ----------------
st.title("GEO Drug Repurposing Explorer")

# Tabs
geo_tab, val_tab, run_tab, res_tab, logs_tab = st.tabs(
    ["GEO Search", "Validate & Select GSE", "Run Pipeline", "Results", "Logs"]
)


# -------- Tab 1: GEO Search --------
with geo_tab:
    st.subheader("Build query and search GEO")
    cfg = st.session_state["config"]

    genes = [g.strip() for g in genes_in.split(",") if g.strip()]

    # If the user provided explicit GEO queries in the Config JSON, honor them.
    # Otherwise build queries from the UI inputs.
    cfg_queries = cfg.get("geo_query_list") or []
    cfg_queries = [str(q).strip() for q in cfg_queries if str(q).strip()]
    if cfg_queries:
        queries = cfg_queries[:3]
    else:
        queries = build_geo_queries(
            genes,
            phenotype_in,
            extra_geo_terms,
            exclude_terms=cfg.get("geo_exclude_terms", []),
            max_queries=3,
        )

    @st.cache_data(show_spinner=False, ttl=3600)
    def cached_geo_search(queries_tuple: tuple, retmax_each: int):
        # Cached wrapper to reduce NCBI request volume (helps avoid 429/rate-limit).
        return geo_search_candidates(list(queries_tuple), retmax_each=retmax_each)

    if st.button("Build & Search GEO"):
        add_log("[UI] GEO search started")
        add_log("Queries: " + " | ".join(queries))

        try:
            with st.spinner("Searching GEO..."):
                df = cached_geo_search(tuple(queries), retmax_each=int(cfg.get("geo_retmax_each", 40)))

            # Surface partial failures (e.g. transient NCBI 429/5xx) without crashing the app.
            errs = []
            try:
                errs = (df.attrs or {}).get("errors", []) if df is not None else []
            except Exception:
                errs = []

            if errs:
                st.warning(
                    f"Some GEO queries failed ({len(errs)}). Showing partial results. "
                    "If this persists, try lowering geo_retmax_each or simplifying keywords."
                )
                # Show only a few to keep UI readable.
                st.code("\n".join([f"{e.get('query','')}: {e.get('error','')}" for e in errs[:5]]))

        except Exception as e:
            add_log(f"[ERROR] GEO search failed: {repr(e)}")
            st.error(
                "GEO search failed (network/rate-limit or query issue). "
                "Please try again, or simplify query terms.\n\n" + repr(e)
            )
            df = pd.DataFrame()

        # Normalise candidate schema: some GEO utilities return 'accession' rather than 'gse'.
        if df is not None and not df.empty:
            df = df.copy()
            if "gse" not in df.columns:
                if "accession" in df.columns:
                    df["gse"] = df["accession"]
                elif "Accession" in df.columns:
                    df["gse"] = df["Accession"]

        st.session_state["candidates"] = df
        st.session_state["geo_queries"] = queries

        # New GEO candidates invalidate downstream state (validation / run results)
        # because option sets change across searches.
        st.session_state.pop("validation", None)
        st.session_state.pop("run_results", None)
        st.session_state.pop("ranked", None)

    df_cand = st.session_state.get("candidates")
    if df_cand is not None and not df_cand.empty:
        st.write(f"Candidates: {df_cand.shape[0]}")
        st.dataframe(df_cand, use_container_width=True)

        csv = df_cand.to_csv(index=False).encode("utf-8")
        st.download_button("Download candidates CSV", csv, file_name="geo_candidates.csv")
    else:
        st.info("No GEO candidates yet. Click Build & Search GEO.")


# -------- Tab 2: Validate --------
with val_tab:
    st.subheader("Download SOFT and validate case/control feasibility")
    cfg = st.session_state["config"]

    df_cand = st.session_state.get("candidates")
    if df_cand is None or df_cand.empty:
        st.info("Run GEO search first.")
    else:
        # Defensive: some candidate frames may still carry 'accession' instead of 'gse'.
        if "gse" not in df_cand.columns:
            if "accession" in df_cand.columns:
                df_cand = df_cand.copy()
                df_cand["gse"] = df_cand["accession"].astype(str)
            elif "Accession" in df_cand.columns:
                df_cand = df_cand.copy()
                df_cand["gse"] = df_cand["Accession"].astype(str)

        gse_options = (
            sorted(df_cand["gse"].astype(str).unique().tolist())
            if "gse" in df_cand.columns
            else []
        )
        selected = st.multiselect("Select GSE to validate", gse_options, default=gse_options[:10])
        if st.button("Validate selected GSE"):
            # Cache parsed SOFT metadata for later inspection (characteristics preview)
            meta_cache: Dict[str, pd.DataFrame] = st.session_state.setdefault("meta_cache", {})
            rows = []
            for gse in selected:
                try:
                    soft_path = download_soft(gse, Path(cfg.get("project_dir", "./repurpose_pipeline")) / "raw" / gse)
                    meta = parse_soft(soft_path)
                    meta = label_conditions(meta, gse, cfg)
                    # Save labeled metadata so users can inspect characteristics and refine config.
                    meta_cache[gse] = meta
                    c = meta["condition"].value_counts().to_dict()
                    rows.append(
                        {
                            "gse": gse,
                            "pass_strict": c.get("unknown", 0) / len(meta) <= float(cfg.get("max_unknown_frac_strict", 0.1)),
                            "pass_relaxed": c.get("unknown", 0) / len(meta) <= float(cfg.get("max_unknown_frac_relaxed", 0.35)),
                            "case": int(c.get("case", 0)),
                            "control": int(c.get("control", 0)),
                            "ambiguous": int(c.get("ambiguous", 0)),
                            "unknown": int(c.get("unknown", 0)),
                            "unknown_frac": c.get("unknown", 0) / len(meta),
                            "total": len(meta),
                        }
                    )
                    add_log(f"[Validate] {gse} counts {c}")
                except Exception as e:
                    rows.append({"gse": gse, "pass_strict": False, "pass_relaxed": False, "case": 0, "control": 0, "ambiguous": 0, "unknown": 0, "unknown_frac": 1.0, "total": 0, "error": repr(e)})
                    add_log(f"[Validate][ERROR] {gse}: {repr(e)}")

            df_val = pd.DataFrame(rows)
            st.session_state["validation"] = df_val

        df_val = st.session_state.get("validation")
        if df_val is not None and not df_val.empty:
            st.dataframe(df_val, use_container_width=True)
            st.download_button("Download validation CSV", df_val.to_csv(index=False).encode("utf-8"), file_name="validation.csv")

            # ---- Sample metadata preview (characteristics) ----
            meta_cache = st.session_state.get("meta_cache", {})
            cached_gses = [g for g in df_val["gse"].astype(str).tolist() if g in meta_cache]
            if cached_gses:
                st.markdown("#### Preview sample annotations (characteristics)")
                gse_preview = st.selectbox(
                    "Choose a GSE to preview",
                    options=cached_gses,
                    index=0,
                    key="gse_preview_select",
                )
                df_meta_preview = meta_cache.get(gse_preview)
                if isinstance(df_meta_preview, pd.DataFrame) and not df_meta_preview.empty:
                    # Quick counts
                    if "condition" in df_meta_preview.columns:
                        vc = df_meta_preview["condition"].value_counts(dropna=False).rename_axis("condition").reset_index(name="n")
                        st.caption("Condition label counts")
                        st.dataframe(vc, use_container_width=True, hide_index=True)

                    # Key frequency summary
                    df_keys = summarize_characteristics_kv(df_meta_preview, top_n=30)
                    with st.expander("Most frequent characteristics keys", expanded=False):
                        st.dataframe(df_keys, use_container_width=True, hide_index=True)

                    # Build a view with selected keys
                    all_keys = df_keys["key"].astype(str).tolist() if not df_keys.empty else []
                    prefer = cfg.get("condition", {}).get("prefer_keys", [])
                    default_keys = [k for k in prefer if k in set(all_keys)] or all_keys[:8]

                    selected_keys = st.multiselect(
                        "Show these characteristics keys as columns",
                        options=all_keys,
                        default=default_keys,
                        key=f"char_keys_{gse_preview}",
                    )
                    df_view = build_characteristics_view(df_meta_preview, selected_keys)
                    st.dataframe(df_view, use_container_width=True, height=360)

                    st.download_button(
                        f"Download {gse_preview} sample annotations CSV",
                        df_view.to_csv(index=False).encode("utf-8"),
                        file_name=f"{gse_preview}_sample_annotations.csv",
                    )
                else:
                    st.info("No cached sample metadata for this GSE. Click 'Validate selected GSE' again.")

            # For running the pipeline, limit options to validated GSEs.
            gse_run_options = sorted(df_val["gse"].tolist())

            # Streamlit will raise if default contains values not present in options
            # (this can happen if you re-run GEO Search and keep an old validation table).
            gse_run_default = [
                x for x in df_val[df_val["pass_relaxed"] == True]["gse"].tolist()
                if x in set(gse_run_options)
            ]
            if not gse_run_default and gse_run_options:
                gse_run_default = gse_run_options[:3]

            gse_run = st.multiselect(
                "Select GSE to run pipeline",
                gse_run_options,
                default=gse_run_default,
            )

            st.session_state["gse_to_run"] = gse_run


# -------- Tab 3: Run Pipeline --------
with run_tab:
    st.subheader("Run pipeline")
    cfg = st.session_state["config"]

    gse_run = st.session_state.get("gse_to_run", [])
    if not gse_run:
        st.info("Select GSEs in the Validate tab first.")
    else:
        st.write("Will run:", ", ".join(gse_run))

        if st.button("Run selected GSE"):
            st.session_state["run_results"] = []
            add_log("[UI] Run pipeline started")

            # Use a *copy* of the config for this run, so we can safely set transient params.
            cfg_run = json.loads(json.dumps(cfg))

            # If user didn't change project_dir, keep everything under ./repurpose_pipeline
            project_dir = Path(cfg_run.get("project_dir", "./repurpose_pipeline"))
            project_dir.mkdir(parents=True, exist_ok=True)

            for i, gse in enumerate(gse_run, 1):
                add_log(f"=== RUN {gse} ({i}/{len(gse_run)}) ===")
                try:
                    out = run_one_gse(gse, cfg_run, logger=add_log)
                    st.session_state["run_results"].append(out)
                except Exception as e:
                    st.session_state["run_results"].append({"gse": gse, "error": repr(e)})
                    add_log(f"[ERROR] {gse}: {repr(e)}")


# -------- Tab 4: Results --------
with res_tab:
    st.subheader("Ranked compounds")
    results = st.session_state.get("run_results", [])
    if not results:
        st.info("No ranked results yet. Run pipeline first.")
    else:
        df_sum = pd.DataFrame([{k: v for k, v in r.items() if k in ("gse", "compound_ranked_csv", "error", "l1000cds2_top_csv", "l1000fwd_error", "enrichr_error")} for r in results])
        st.write("Run summary:")
        st.dataframe(df_sum, use_container_width=True)

        # show per-GSE compound tables
        for r in results:
            gse = r.get("gse")
            csv_path = r.get("compound_ranked_csv")
            if csv_path and os.path.exists(csv_path):
                st.markdown(f"#### {gse}")
                df = pd.read_csv(csv_path)
                st.dataframe(df.head(50), use_container_width=True)
                st.download_button(
                    f"Download {gse} compounds CSV",
                    df.to_csv(index=False).encode("utf-8"),
                    file_name=f"{gse}_compound_ranked.csv",
                )


# -------- Tab 5: Logs --------
with logs_tab:
    st.subheader("Logs")
    st.text("\n".join(st.session_state["logs"][-400:]))
