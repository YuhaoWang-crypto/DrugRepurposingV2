from __future__ import annotations

import numpy as np
import pandas as pd


def compute_de_bulk_ttest(df_expr: pd.DataFrame, meta: pd.DataFrame, eps: float = 1e-6) -> pd.DataFrame:
    """Compute differential expression with a simple t-test.

    Important: GEO /supplementary matrices often have columns that don't
    exactly match GSM IDs. We therefore prefer a pre-aligned meta column
    `expr_col` (added by align_expr_and_meta) to map conditions -> columns.
    """
    if 'condition' not in meta.columns:
        raise RuntimeError('meta is missing required column: condition')

    # Prefer explicit mapping from meta->expression columns
    if 'expr_col' in meta.columns:
        case_cols = meta.loc[meta['condition'] == 'case', 'expr_col'].astype(str).tolist()
        ctrl_cols = meta.loc[meta['condition'] == 'control', 'expr_col'].astype(str).tolist()
        # Keep only columns that actually exist (safety)
        case_cols = [c for c in case_cols if c in df_expr.columns]
        ctrl_cols = [c for c in ctrl_cols if c in df_expr.columns]
    else:
        cond = meta['condition'].astype(str)
        case_cols = [c for c in df_expr.columns if cond.get(c, 'unknown') == 'case']
        ctrl_cols = [c for c in df_expr.columns if cond.get(c, 'unknown') == 'control']

    if len(case_cols) < 2 or len(ctrl_cols) < 2:
        raise RuntimeError(f"Not enough samples for DE: case={len(case_cols)} control={len(ctrl_cols)}")

    X = df_expr[case_cols].to_numpy(dtype=float) + eps
    Y = df_expr[ctrl_cols].to_numpy(dtype=float) + eps

    # log2 fold change
    lfc = np.log2(X.mean(axis=1)) - np.log2(Y.mean(axis=1))

    # Welch t-test (vectorized)
    nx, ny = X.shape[1], Y.shape[1]
    vx = X.var(axis=1, ddof=1)
    vy = Y.var(axis=1, ddof=1)
    denom = np.sqrt(vx / nx + vy / ny)
    denom[denom == 0] = np.nan
    tstat = (X.mean(axis=1) - Y.mean(axis=1)) / denom

    # approximate two-sided p-values using normal approximation
    # (good enough for ranking; avoids scipy dependency)
    from math import erf, sqrt

    def norm_sf(z):
        # survival function for standard normal
        return 0.5 * (1 - erf(z / sqrt(2)))

    pvals = np.array([2 * norm_sf(abs(z)) if np.isfinite(z) else 1.0 for z in tstat])

    df = pd.DataFrame({
        'gene_raw': df_expr.index.astype(str),
        'log2fc': lfc,
        'pval': pvals,
    })
    df = df.sort_values('pval', ascending=True).reset_index(drop=True)
    return df
