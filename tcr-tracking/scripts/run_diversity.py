"""Pre vs on-treatment diversity per patient: raw, and rarefied to a common depth.

Usage:
    python scripts/run_diversity.py                 # clonotype = CDR3 amino acid (as in the authors' diversity metrics)
    python scripts/run_diversity.py --depth 3500    # fixed depth for the cross-patient version (default 3,500)
    python scripts/run_diversity.py --level aa_vj

Three versions of every metric are produced:
    raw        full tables as sequenced
    rar_pair   both samples downsampled to the smaller of the pair (equal depth within a patient)
    rar_fixed  both samples downsampled to --depth templates (equal depth across patients; pairs below it are NaN)

Outputs: outputs/diversity_summary_<level>.csv (per patient) and outputs/diversity_group_tests_<level>.csv.
The group tests are exploratory (5-9 patients per group). Paths: --data-dir / --metrics-csv or TCR_DATA_DIR / TCR_METRICS_CSV.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tcr_tracking.diversity import compare_groups, diversity_summary, rarefied_diversity  # noqa: E402
from tcr_tracking.expansion import bh_fdr  # noqa: E402
from tcr_tracking.load import find_pairs, load_sample, productive_filter, repertoire  # noqa: E402

RATIO_METRICS = ("richness", "inv_simpson")                       # reported as on / pre
DIFF_METRICS = ("shannon_bits", "evenness", "clonality", "top_clone_share", "D90")   # reported as on - pre
TEST_COLS = ("richness_fc", "evenness_delta", "shannon_bits_delta", "inv_simpson_fc", "top_clone_share_delta")


def add_changes(rec: dict, version: str, pre: dict, on: dict) -> None:
    for m in RATIO_METRICS:
        rec[f"{version}_{m}_pre"], rec[f"{version}_{m}_on"] = pre[m], on[m]
        rec[f"{version}_{m}_fc"] = on[m] / pre[m]
    for m in DIFF_METRICS:
        rec[f"{version}_{m}_pre"], rec[f"{version}_{m}_on"] = pre[m], on[m]
        rec[f"{version}_{m}_delta"] = on[m] - pre[m]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", default=os.environ.get("TCR_DATA_DIR"))
    ap.add_argument("--metrics-csv", default=os.environ.get("TCR_METRICS_CSV"))
    ap.add_argument("--level", default="aa", choices=["aa", "aa_vj"])
    ap.add_argument("--depth", type=int, default=3500)
    ap.add_argument("--n-iter", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    if not args.data_dir or not args.metrics_csv:
        ap.error("--data-dir and --metrics-csv are required (or set TCR_DATA_DIR / TCR_METRICS_CSV)")

    ref = pd.read_csv(args.metrics_csv, sep=None, engine="python")
    ref["patient"] = ref["Patient"].str.lower().str.replace("pt", "", regex=False).astype(int)
    ref = ref.set_index("patient")[["Cohort", "myBOR"]]

    rows = []
    for patient, files in find_pairs(args.data_dir).items():
        pre = repertoire(productive_filter(load_sample(files["pre"])), args.level).to_numpy(dtype=np.int64)
        on = repertoire(productive_filter(load_sample(files["on"])), args.level).to_numpy(dtype=np.int64)
        rec: dict = {"patient": patient, "pre_templates": int(pre.sum()), "on_templates": int(on.sum())}
        rec["log2_depth_ratio_on_over_pre"] = float(np.log2(on.sum() / pre.sum()))
        add_changes(rec, "raw", diversity_summary(pre), diversity_summary(on))
        d_pair = int(min(pre.sum(), on.sum()))
        add_changes(rec, "rar_pair", rarefied_diversity(pre, d_pair, args.n_iter, args.seed),
                    rarefied_diversity(on, d_pair, args.n_iter, args.seed + 1))
        if min(pre.sum(), on.sum()) >= args.depth:
            add_changes(rec, "rar_fixed", rarefied_diversity(pre, args.depth, args.n_iter, args.seed),
                        rarefied_diversity(on, args.depth, args.n_iter, args.seed + 1))
        rows.append(rec)
        print(f"Pt{patient}: richness pre/on raw {rec['raw_richness_pre']:>6}/{rec['raw_richness_on']:>6}  "
              f"rarefied-pair {rec['rar_pair_richness_pre']:>8.0f}/{rec['rar_pair_richness_on']:>8.0f}", flush=True)

    out = pd.DataFrame(rows).set_index("patient").join(ref)
    out["benefit"] = out["myBOR"].map({"PRCR": "benefit", "SD": "benefit", "PD": "no_benefit"})
    Path("outputs").mkdir(exist_ok=True)
    out.to_csv(f"outputs/diversity_summary_{args.level}.csv")

    pd.set_option("display.width", 220)
    print(f"\nlevel={args.level}; patients={len(out)}; fixed-depth version has {out['rar_fixed_richness_fc'].notna().sum()}")

    print("\nHow much does the richness change just reflect depth? Spearman(richness fold change, log2 depth ratio on/pre):")
    for v in ("raw", "rar_pair", "rar_fixed"):
        s = out[[f"{v}_richness_fc", "log2_depth_ratio_on_over_pre"]].dropna()
        r, p = spearmanr(np.log2(s[f"{v}_richness_fc"]), s["log2_depth_ratio_on_over_pre"])
        print(f"  {v:10s} rho = {r:+.2f}  (n={len(s)})")

    print("\nOverall pre -> on shift (Wilcoxon signed-rank, all patients with data):")
    shift = []
    for v in ("raw", "rar_pair", "rar_fixed"):
        for m in ("richness", "evenness", "shannon_bits", "inv_simpson", "top_clone_share"):
            s = out[[f"{v}_{m}_pre", f"{v}_{m}_on"]].dropna()
            shift.append({"version": v, "metric": m, "n": len(s), "median_pre": s.iloc[:, 0].median(),
                          "median_on": s.iloc[:, 1].median(), "p_wilcoxon": wilcoxon(s.iloc[:, 0], s.iloc[:, 1]).pvalue})
    print(pd.DataFrame(shift).round(4).to_string(index=False))

    tests = []
    for v in ("raw", "rar_pair", "rar_fixed"):
        rows_v = []
        for col in TEST_COLS:
            for r in compare_groups(out, f"{v}_{col}"):
                r["version"], r["metric"] = v, col
                rows_v.append(r)
        if rows_v:
            q = bh_fdr(np.array([r["p_mannwhitney"] for r in rows_v]))
            for r, qi in zip(rows_v, q):
                r["q_bh_within_version"] = float(qi)
            tests += rows_v
    tests = pd.DataFrame(tests)[["version", "metric", "cohort", "n_benefit", "n_no_benefit", "median_benefit",
                                 "median_no_benefit", "p_mannwhitney", "q_bh_within_version"]]
    tests.to_csv(f"outputs/diversity_group_tests_{args.level}.csv", index=False)
    print("\nBenefit vs no benefit on the pre->on change, per cohort (Mann-Whitney; exploratory):")
    print(tests.round(3).to_string(index=False))
    print(f"\nWrote outputs/diversity_summary_{args.level}.csv and outputs/diversity_group_tests_{args.level}.csv")


if __name__ == "__main__":
    main()
