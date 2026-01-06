# CLCN GEO → Connectivity (L1000CDS2 / L1000FWD / Enrichr) → (optional) PPI → Ranked Drug List

This is a Streamlit wrapper around a Python pipeline that:

1) builds GEO queries from your **gene(s)** + optional phenotype keywords  
2) searches GEO (GDS) to get candidate **GSE**  
3) downloads **GSE family SOFT** and infers **case/control** labels (with editable rules)  
4) downloads a **bulk** expression matrix from GEO FTP (heuristic selection)  
5) computes a simple DE table (log2FC + Welch t-test + BH-FDR)  
6) builds UP/DN gene signatures  
7) runs connectivity / enrichment:
   - L1000CDS2 (Characteristic Direction search engine)
   - L1000FWD (signature similarity search)
   - Enrichr drug/perturbation libraries
8) merges & re-ranks compounds across GSE, and (optionally) runs PPI module checks.

> Notes
- This implementation focuses on **bulk matrices** (TSV/CSV/TXT). 10x/scRNA support is not fully automated here.
- For reproducibility and speed, everything is cached under `project_dir/raw` and `project_dir/proc`.

## Run locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

## What you need to prepare (optional)

- A PPI edges CSV with two columns: `a,b` (gene symbols). Upload in the UI.
- GO terms for your module(s), default includes lysosome/autophagy/ion transport.

## File structure

- `app.py`: Streamlit UI
- `pipeline/`: reusable pipeline modules
