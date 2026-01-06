from __future__ import annotations

from typing import Any, Dict


def default_config() -> Dict[str, Any]:
    """Default configuration shown in Streamlit as editable JSON.

    This should be *generic* (no disease-specific branding), so the public app can
    be reused across targets/indications. You can always paste a disease-specific
    template into the Config JSON editor.
    """
    return {
        # Project dirs
        "project_dir": "./repurpose_pipeline",

        # GEO search
        "geo_retmax_each": 40,
        # Terms to exclude from GEO query (set [] to disable)
        "geo_exclude_terms": [],

        # If you already know the GSEs, you can paste them here; otherwise leave empty
        "gse_list": [],

        # Query templates (auto-built in UI; you can override here)
        "geo_query_list": [],

        # Condition inference
        "condition": {
            "prefer_keys": [
                "diagnosis",
                "disease",
                "disease state",
                "condition",
                "group",
                "treatment",
                "genotype",
                "phenotype",
                "infection",
                "time",
            ],
            # These are *heuristics*. For best results, set disease-specific terms.
            "case_terms": [
                "case",
                "patient",
                "disease",
                "tumor",
                "cancer",
                "infected",
                "treated",
                "stimulated",
                "mutation",
                "mutant",
                "knockout",
                "ko",
                "overexpression",
                "oe",
            ],
            "control_terms": [
                "control",
                "healthy",
                "normal",
                "baseline",
                "untreated",
                "vehicle",
                "wild type",
                "wild-type",
                "wt",
            ],
            # ambiguous | case | control
            "ambiguous_policy": "ambiguous",
        },

        # Per-GSE overrides (optional). Regex strings are used for title/column fallbacks.
        "per_gse_overrides": {},

        # Validation thresholds
        "min_reps_per_group": 2,
        "max_unknown_frac_strict": 0.1,
        "max_unknown_frac_relaxed": 0.35,

        # Bulk DE parameters
        "de_min_samples_total": 4,
        "de_eps": 1e-6,

        # Signature parameters
        "signature_top_n": 150,
        "sig_min_abs_log2fc": 0.25,

        # Connectivity toggles
        "run_l1000cds2": True,
        "run_l1000fwd": True,
        "run_enrichr": True,

        # L1000CDS2
        "l1000_mode": "reverse",
        "l1000_urls": [
            "https://maayanlab.cloud/L1000CDS2/query",
            "https://amp.pharm.mssm.edu/L1000CDS2/query",
        ],
        "l1000_timeout": 120,
        "l1000_limit": 50,

        # L1000FWD
        "l1000fwd_base_urls": [
            "https://maayanlab.cloud/l1000fwd",
            "https://maayanlab.cloud/L1000FWD",
        ],
        "l1000fwd_timeout": 120,
        "l1000fwd_limit": 50,

        # Enrichr
        "enrichr_base_url": "https://maayanlab.cloud/Enrichr",
        "enrichr_libraries": [
            "LINCS_L1000_Chem_Pert_up",
            "LINCS_L1000_Chem_Pert_down",
            "Drug_Perturbations_from_GEO_up",
            "Drug_Perturbations_from_GEO_down",
        ],

        # Ranking weights
        "rank_weights": {"cds2": 1.0, "l1000fwd": 1.0, "enrichr": 0.7},

        # Optional PPI / module parameters
        "module_gene_symbols": [],
        "module_go_terms": [],
        "quickgo_taxon": "9606",
        "ppi_edges_csv": "",
        "ppi_max_dist": 2,
    }
