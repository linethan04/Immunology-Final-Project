"""Per-patient clonotype expansion tests (pre vs on treatment), with a rarefied-depth sensitivity run.

Usage:
    python scripts/run_expansion.py                       # clonotype = CDR3 aa + V + J
    python scripts/run_expansion.py --level aa
    python scripts/run_expansion.py --depth 0             # skip the rarefied sensitivity run
    python scripts/run_expansion.py --min-log2fc 1        # also require a 2-fold increase to call a clone expanded

Outputs (in outputs/):
    expansion_summary_<level>.csv   one row per patient (full data + rarefied sensitivity + cohort/response)
    expansion_clones_<level>.csv    only the clonotypes called expanded or contracted (q < alpha)
Paths come from --data-dir / --metrics-csv or TCR_DATA_DIR / TCR_METRICS_CSV.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tcr_tracking.expansion import expansion_test, summarize_expansion, compute_expansion  # noqa: E402
from tcr_tracking.load import find_pairs, load_sample, productive_filter, repertoire  # noqa: E402
from tcr_tracking.tracking import rarefy  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", default=os.environ.get("TCR_DATA_DIR"))
    ap.add_argument("--metrics-csv", default=os.environ.get("TCR_METRICS_CSV"))
    ap.add_argument("--level", default="aa_vj", choices=["aa", "aa_vj"])
    ap.add_argument("--alpha", type=float, default=0.05, help="FDR threshold")
    ap.add_argument("--min-total", type=int, default=2, help="test only clonotypes with pre+on templates >= this")
    ap.add_argument("--min-log2fc", type=float, default=0.0, help="extra fold-change requirement (log2) for a call")
    ap.add_argument("--depth", type=int, default=3500, help="sensitivity: downsample both samples to this many templates (0 = skip)")
    ap.add_argument("--n-iter", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    if not args.data_dir or not args.metrics_csv:
        ap.error("--data-dir and --metrics-csv are required (or set TCR_DATA_DIR / TCR_METRICS_CSV)")

    ref = pd.read_csv(args.metrics_csv, sep=None, engine="python")
    ref["patient"] = ref["Patient"].str.lower().str.replace("pt", "", regex=False).astype(int)
    ref = ref.set_index("patient")[["Cohort", "myBOR"]]

    summary_rows, sig_tables = [], []
    for patient, files in find_pairs(args.data_dir).items():
        pre = repertoire(productive_filter(load_sample(files["pre"])), args.level)
        on = repertoire(productive_filter(load_sample(files["on"])), args.level)
        df = expansion_test(pre, on, args.min_total, args.alpha, args.min_log2fc)
        rec = {"patient": patient, "pre_templates": int(df["pre"].sum()), "on_templates": int(df["on"].sum()),
               **summarize_expansion(df)}

        # sensitivity: repeat on randomly downsampled data at a common depth
        if args.depth and rec["pre_templates"] >= args.depth and rec["on_templates"] >= args.depth:
            a, b = df["pre"].to_numpy(), df["on"].to_numpy()
            rng = np.random.default_rng(args.seed)
            counts = [int(compute_expansion(rarefy(a, args.depth, rng), rarefy(b, args.depth, rng),
                                      args.min_total, args.alpha, args.min_log2fc)["expanded"].sum())
                      for _ in range(args.n_iter)]
            rec.update({"rar_depth": args.depth, "rar_n_expanded_mean": float(np.mean(counts)),
                        "rar_n_expanded_sd": float(np.std(counts, ddof=1))})
        summary_rows.append(rec)

        sig = df[df["expanded"] | df["contracted"]].copy()
        sig["call"] = np.where(sig["expanded"], "expanded", "contracted")
        sig = sig.reset_index()
        sig.insert(0, "patient", patient)
        sig_tables.append(sig.drop(columns=["expanded", "contracted"]))
        print(f"Pt{patient}: expanded={rec['n_expanded']:>4} (novel {rec['n_expanded_novel']:>3}) "
              f"contracted={rec['n_contracted']:>4} | on-templates in expanded clones={rec['frac_on_templates_in_expanded']:.3f}"
              + (f" | rarefied n_expanded={rec['rar_n_expanded_mean']:.1f}" if "rar_n_expanded_mean" in rec else ""), flush=True)

    out = pd.DataFrame(summary_rows).set_index("patient").join(ref)
    out["benefit"] = out["myBOR"].map({"PRCR": "benefit", "SD": "benefit", "PD": "no_benefit"})
    Path("outputs").mkdir(exist_ok=True)
    out.to_csv(f"outputs/expansion_summary_{args.level}.csv")
    pd.concat(sig_tables, ignore_index=True).to_csv(f"outputs/expansion_clones_{args.level}.csv", index=False)

    pd.set_option("display.width", 220)
    cols = ["n_expanded", "rar_n_expanded_mean", "frac_on_templates_in_expanded", "n_contracted"]
    g = out.groupby(["Cohort", "benefit"])
    print(f"\nlevel={args.level} alpha={args.alpha} min_total={args.min_total} min_log2fc={args.min_log2fc}")
    print("Medians by cohort and benefit group (descriptive only; n is small):")
    print(pd.concat([g.size().rename("n"), g[cols].median().round(3)], axis=1).to_string())
    print("\nWrote outputs/expansion_summary_%s.csv and outputs/expansion_clones_%s.csv" % (args.level, args.level))


if __name__ == "__main__":
    main()
