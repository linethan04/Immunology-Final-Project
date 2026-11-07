"""Depth-robust diversity of a single repertoire and group comparison of pre -> on changes.

Observed richness (number of unique clonotypes) grows with sequencing depth, so a raw pre/on richness ratio
partly measures how deeply each sample was sequenced. ``rarefied_diversity`` downsamples a repertoire to a
fixed number of templates (without replacement) and averages the metrics over repeated draws, so that two
samples can be compared at equal depth.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from .metrics import d_fraction, evenness, shannon
from .tracking import rarefy


def diversity_summary(counts: np.ndarray) -> dict[str, float]:
    s = pd.Series(np.asarray(counts)[np.asarray(counts) > 0].astype(np.int64))
    p = s.to_numpy(dtype=float) / s.sum()
    ev = evenness(s)
    return {
        "templates": int(s.sum()),
        "richness": int(len(s)),
        "shannon_bits": shannon(s),
        "evenness": ev,
        "clonality": 1.0 - ev,
        "inv_simpson": float(1.0 / (p * p).sum()),
        "top_clone_share": float(p.max()),
        "D90": d_fraction(s, 0.9),
    }


def rarefied_diversity(counts: np.ndarray, depth: int, n_iter: int = 20, seed: int = 1) -> dict[str, float]:
    """Mean over ``n_iter`` random downsamplings to ``depth`` templates (and ``<metric>_sd``)."""
    rng = np.random.default_rng(seed)
    draws = pd.DataFrame([diversity_summary(rarefy(np.asarray(counts, dtype=np.int64), depth, rng))
                          for _ in range(n_iter)]).drop(columns="templates")
    out = draws.mean().to_dict()
    out.update({f"{k}_sd": float(v) for k, v in draws.std(ddof=1).items()})
    return out


def mwu_p(a, b) -> float:
    """Two-sided Mann-Whitney p. Exact when there are no ties (as R's wilcox.test does for n < 50),
    otherwise the normal approximation. scipy's default only goes exact when both groups have <= 8 values."""
    pooled = np.concatenate([np.asarray(a, float), np.asarray(b, float)])
    method = "exact" if len(np.unique(pooled)) == len(pooled) else "asymptotic"
    return float(mannwhitneyu(a, b, alternative="two-sided", method=method).pvalue)


def compare_groups(df: pd.DataFrame, value_col: str) -> list[dict]:
    """Benefit vs no benefit on a per-patient change, separately for each cohort (two-sided Mann-Whitney).

    ``df`` needs columns Cohort and benefit ("benefit" / "no_benefit") plus ``value_col``.
    """
    rows = []
    for cohort, g in df.groupby("Cohort"):
        a = g.loc[g["benefit"] == "benefit", value_col].dropna()
        b = g.loc[g["benefit"] == "no_benefit", value_col].dropna()
        if len(a) >= 2 and len(b) >= 2:
            rows.append({"cohort": cohort, "metric": value_col, "n_benefit": len(a), "n_no_benefit": len(b),
                         "median_benefit": float(a.median()), "median_no_benefit": float(b.median()),
                         "p_mannwhitney": mwu_p(a, b)})
    return rows
