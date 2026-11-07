"""Single-repertoire diversity metrics, following Riaz et al. (2017) / Sims et al. (2016).

All functions take a pandas Series of template counts (one entry per clonotype).
Entropy is in bits (log2); evenness = H / log2(N); clonality = 1 - evenness.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def richness(counts: pd.Series) -> int:
    """Number of unique clonotypes."""
    return int((counts > 0).sum())


def shannon(counts: pd.Series) -> float:
    c = counts[counts > 0].to_numpy(dtype=float)
    p = c / c.sum()
    return float(-(p * np.log2(p)).sum())


def evenness(counts: pd.Series, single_clone_value: float = 0.0) -> float:
    """H / log2(N). With N == 1 the ratio is 0/0; the authors' per-VJ files report 0 there."""
    n = richness(counts)
    if n <= 1:
        return single_clone_value
    return shannon(counts) / np.log2(n)


def clonality(counts: pd.Series) -> float:
    return 1.0 - evenness(counts)


def d_fraction(counts: pd.Series, frac: float = 0.9, include_crossing_clone: bool = False) -> float:
    """Fraction of unique clonotypes that account for ``frac`` of all templates (D90 for 0.9).

    Clonotypes are ranked by abundance. By default the count is the number of top clonotypes whose
    cumulative share is still below ``frac``; this reproduces the values in the authors' TCRmetrics.csv.
    With ``include_crossing_clone=True`` the clonotype that carries the cumulative share to or past
    ``frac`` is counted too (the strict "minimum number needed to reach 90%"), which is one clone higher.
    """
    c = np.sort(counts[counts > 0].to_numpy(dtype=float))[::-1]
    cum = np.cumsum(c) / c.sum()
    k = int(np.searchsorted(cum, frac - 1e-12))  # clonotypes with cumulative share < frac
    if include_crossing_clone:
        k += 1
    return k / len(c)


def summarize(counts: pd.Series) -> dict[str, float]:
    return {
        "templates": int(counts.sum()),
        "richness": richness(counts),
        "shannon_bits": shannon(counts),
        "evenness": evenness(counts),
        "clonality": clonality(counts),
        "D90": d_fraction(counts, 0.9),
    }
