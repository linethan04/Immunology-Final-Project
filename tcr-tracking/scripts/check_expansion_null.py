"""Calibration check for the expansion test on real clone-count structure.

Each sample is split into two random halves by binomial thinning of every clonotype's count (p = 0.5)
and the two halves are compared with the same test. With sampling noise as the only difference, the
test should call (almost) nothing at q < alpha.

This is a best-case null: it contains no PCR, library-prep or biopsy-to-biopsy variation, so passing it
does not show the test is calibrated on real pre/on pairs (see README, assumption A12).

Usage: python scripts/check_expansion_null.py [--level aa_vj]
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tcr_tracking.expansion import compute_expansion  # noqa: E402
from tcr_tracking.load import find_pairs, load_sample, productive_filter, repertoire  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", default=os.environ.get("TCR_DATA_DIR"))
    ap.add_argument("--level", default="aa_vj", choices=["aa", "aa_vj"])
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--min-total", type=int, default=2)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    if not args.data_dir:
        ap.error("--data-dir is required (or set TCR_DATA_DIR)")

    rng = np.random.default_rng(args.seed)
    rows = []
    for patient, files in find_pairs(args.data_dir).items():
        for tp in ("pre", "on"):
            x = repertoire(productive_filter(load_sample(files[tp])), args.level).to_numpy(dtype=np.int64)
            h1 = rng.binomial(x, 0.5)
            h2 = x - h1
            r = compute_expansion(h1, h2, args.min_total, args.alpha)
            tested = ~np.isnan(r["p_expand"])
            rows.append({"patient": patient, "sample": tp, "templates": int(x.sum()), "n_tested": int(tested.sum()),
                         "expanded_calls": int(r["expanded"].sum()), "contracted_calls": int(r["contracted"].sum()),
                         "frac_p_expand_below_0.05": float((r["p_expand"][tested] < 0.05).mean())})
    out = pd.DataFrame(rows)
    pd.set_option("display.width", 200)
    print(f"Null split-half check (level={args.level}, alpha={args.alpha}, {len(out)} samples)")
    print(f"total calls: expanded={out.expanded_calls.sum()}, contracted={out.contracted_calls.sum()} "
          f"(expected ~0); samples with any call: {int(((out.expanded_calls + out.contracted_calls) > 0).sum())}")
    print(f"fraction of tested clonotypes with raw p < 0.05: mean={out['frac_p_expand_below_0.05'].mean():.4f} "
          f"max={out['frac_p_expand_below_0.05'].max():.4f} (should not exceed ~0.05)")
    Path("outputs").mkdir(exist_ok=True)
    out.to_csv("outputs/expansion_null_check.csv", index=False)


if __name__ == "__main__":
    main()
