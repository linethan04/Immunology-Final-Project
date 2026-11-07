"""Track clonotypes between the pre- and on-treatment repertoires of one patient.

Replacement metrics answer: how much of the on-treatment repertoire is made of clonotypes that
were already detectable at baseline? Because detection depends on sequencing depth, every metric
is available in two forms:

  raw       computed on the full tables as sequenced
  rarefied  both samples are randomly downsampled (without replacement, at the template level) to
            a common depth and the metric is averaged over repeated draws
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def align(pre: pd.Series, on: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    """Template counts for the union of clonotypes (0 where a clonotype is absent), same order in both."""
    idx = pre.index.union(on.index)
    a = pre.reindex(idx, fill_value=0).to_numpy(dtype=np.int64)
    b = on.reindex(idx, fill_value=0).to_numpy(dtype=np.int64)
    return a, b


def morisita_horn(a: np.ndarray, b: np.ndarray) -> float:
    """Morisita-Horn overlap on abundances (0 = no overlap, 1 = identical composition)."""
    x, y = a.astype(float), b.astype(float)
    X, Y = x.sum(), y.sum()
    if X == 0 or Y == 0:
        return float("nan")
    dx, dy = (x * x).sum() / X**2, (y * y).sum() / Y**2
    return float(2.0 * (x * y).sum() / ((dx + dy) * X * Y))


def replacement_metrics(a: np.ndarray, b: np.ndarray) -> dict[str, float]:
    in_pre, in_on = a > 0, b > 0
    shared = in_pre & in_on
    n_pre, n_on, n_shared = int(in_pre.sum()), int(in_on.sum()), int(shared.sum())
    A, B = a.sum(), b.sum()
    return {
        "pre_templates": int(A),
        "on_templates": int(B),
        "pre_clonotypes": n_pre,
        "on_clonotypes": n_on,
        "shared_clonotypes": n_shared,
        # fraction of the on-treatment repertoire that comes from baseline clones
        "frac_on_templates_from_baseline": float(b[in_pre].sum() / B) if B else float("nan"),
        "frac_on_clonotypes_from_baseline": n_shared / n_on if n_on else float("nan"),
        # converse: how much of the baseline repertoire is still detectable on treatment
        "frac_pre_templates_persisting": float(a[in_on].sum() / A) if A else float("nan"),
        "frac_pre_clonotypes_persisting": n_shared / n_pre if n_pre else float("nan"),
        "jaccard": n_shared / int((in_pre | in_on).sum()),
        "morisita_horn": morisita_horn(a, b),
    }


def rarefy(x: np.ndarray, depth: int, rng: np.random.Generator) -> np.ndarray:
    """Downsample a count vector to ``depth`` templates without replacement."""
    if depth >= x.sum():
        return x.copy()
    return rng.multivariate_hypergeometric(x, int(depth))


RAREFIED_KEYS = (
    "frac_on_templates_from_baseline", "frac_on_clonotypes_from_baseline",
    "frac_pre_templates_persisting", "frac_pre_clonotypes_persisting",
    "jaccard", "morisita_horn", "pre_clonotypes", "on_clonotypes", "shared_clonotypes",
)


def rarefied_metrics(a: np.ndarray, b: np.ndarray, depth: int, n_iter: int = 20, seed: int = 1) -> dict[str, float]:
    """Mean and SD (over ``n_iter`` random downsamplings to ``depth``) of the replacement metrics."""
    rng = np.random.default_rng(seed)
    draws = [replacement_metrics(rarefy(a, depth, rng), rarefy(b, depth, rng)) for _ in range(n_iter)]
    df = pd.DataFrame(draws)[list(RAREFIED_KEYS)]
    out = {f"rar_{k}": float(v) for k, v in df.mean().items()}
    out.update({f"rar_sd_{k}": float(v) for k, v in df.std(ddof=1).items()})
    out["rar_depth"] = int(depth)
    return out
