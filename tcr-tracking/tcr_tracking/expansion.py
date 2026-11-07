"""Per-clonotype expansion / contraction tests between pre- and on-treatment repertoires.

For clonotype i with a_i templates in the pre sample (total A) and b_i in the on sample (total B),
the one-sided Fisher exact test asks whether its share of the on-treatment repertoire is larger
(expansion) or smaller (contraction) than its share at baseline. Sequencing depth enters only through
the sample totals A and B, so the comparison is depth-normalized. The p-value is the hypergeometric tail
    P(X >= b_i), X ~ Hypergeom(population A+B, successes B, draws a_i+b_i)
which is identical to scipy.stats.fisher_exact on [[b_i, B-b_i], [a_i, A-a_i]] (alternative="greater").

This is the same family of test the immunoSEQ Analyzer applies (DeWitt et al. 2015, as used by Riaz et al.
2017). Differences: Benjamini-Hochberg FDR here, versus the Analyzer's positive FDR, and an explicit
minimum-count filter (see ``min_total``).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import hypergeom


def bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (q-values)."""
    p = np.asarray(p, dtype=float)
    n = p.size
    if n == 0:
        return p.copy()
    order = np.argsort(p, kind="mergesort")
    ranked = p[order] * n / np.arange(1, n + 1)
    q_sorted = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.minimum(q_sorted, 1.0)
    return q


def expansion_pvalues(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """One-sided Fisher p-values (expansion, contraction) for every clonotype."""
    A, B = int(a.sum()), int(b.sum())
    M, N = A + B, a + b
    p_expand = hypergeom.sf(b - 1, M, B, N)      # P(X >= b)
    p_contract = hypergeom.cdf(b, M, B, N)       # P(X <= b)
    return p_expand, p_contract


def compute_expansion(a: np.ndarray, b: np.ndarray, min_total: int = 2, alpha: float = 0.05,
                min_log2fc: float = 0.0) -> dict[str, np.ndarray]:
    """Run the tests on aligned count vectors. Returns arrays (NaN where a clonotype was not tested)."""
    a, b = np.asarray(a, dtype=np.int64), np.asarray(b, dtype=np.int64)
    A, B = a.sum(), b.sum()
    p_exp, p_con = expansion_pvalues(a, b)
    total = a + b
    t_exp = (b >= 1) & (total >= min_total)      # only clonotypes seen on treatment can expand
    t_con = (a >= 1) & (total >= min_total)      # only clonotypes seen at baseline can contract
    q_exp = np.full(a.shape, np.nan)
    q_con = np.full(a.shape, np.nan)
    q_exp[t_exp] = bh_fdr(p_exp[t_exp])
    q_con[t_con] = bh_fdr(p_con[t_con])
    log2fc = np.log2(((b + 0.5) / B) / ((a + 0.5) / A))      # descriptive only; not used by the test
    with np.errstate(invalid="ignore"):
        expanded = t_exp & (q_exp < alpha) & (log2fc >= min_log2fc)
        contracted = t_con & (q_con < alpha) & (log2fc <= -min_log2fc)
    return {
        "p_expand": np.where(t_exp, p_exp, np.nan), "q_expand": q_exp,
        "p_contract": np.where(t_con, p_con, np.nan), "q_contract": q_con,
        "log2fc": log2fc, "expanded": expanded, "contracted": contracted,
    }


def expansion_test(pre: pd.Series, on: pd.Series, min_total: int = 2, alpha: float = 0.05,
                   min_log2fc: float = 0.0) -> pd.DataFrame:
    """Test every clonotype in the union of two repertoires (Series of template counts, same clonotype index)."""
    idx = pre.index.union(on.index)
    a = pre.reindex(idx, fill_value=0).to_numpy(dtype=np.int64)
    b = on.reindex(idx, fill_value=0).to_numpy(dtype=np.int64)
    res = compute_expansion(a, b, min_total, alpha, min_log2fc)
    df = pd.DataFrame({"pre": a, "on": b, "freq_pre": a / a.sum(), "freq_on": b / b.sum()}, index=idx)
    for k, v in res.items():
        df[k] = v
    df["novel"] = df["pre"] == 0           # not detected at baseline
    return df


def summarize_expansion(df: pd.DataFrame) -> dict[str, float]:
    A, B = int(df["pre"].sum()), int(df["on"].sum())
    e, c = df[df["expanded"]], df[df["contracted"]]
    return {
        "n_tested_expand": int(df["p_expand"].notna().sum()),
        "n_expanded": int(len(e)),
        "n_expanded_preexisting": int((~e["novel"]).sum()),
        "n_expanded_novel": int(e["novel"].sum()),
        "frac_on_templates_in_expanded": float(e["on"].sum() / B),
        "top_expanded_freq_on": float(e["freq_on"].max()) if len(e) else 0.0,
        "n_tested_contract": int(df["p_contract"].notna().sum()),
        "n_contracted": int(len(c)),
        "frac_pre_templates_in_contracted": float(c["pre"].sum() / A),
    }
