"""Per-patient clonotype tracking and baseline-replacement metrics (raw and depth-matched).

Usage:
    python scripts/run_tracking.py                       # clonotype = CDR3 aa + V + J, per-pair common depth
    python scripts/run_tracking.py --level aa            # clonotype = CDR3 amino acid only
    python scripts/run_tracking.py --depth 5000          # fixed common depth; pairs with a sample below it get NaN

Paths come from --data-dir / --metrics-csv or the TCR_DATA_DIR / TCR_METRICS_CSV environment variables.
Response labels are taken from the authors' TCRmetrics.csv (Cohort, myBOR); "benefit" = CR/PR or SD.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tcr_tracking.load import find_pairs, load_sample, productive_filter, repertoire  # noqa: E402
from tcr_tracking.tracking import align, rarefied_metrics, replacement_metrics  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", default=os.environ.get("TCR_DATA_DIR"))
    ap.add_argument("--metrics-csv", default=os.environ.get("TCR_METRICS_CSV"))
    ap.add_argument("--level", default="aa_vj", choices=["aa", "aa_vj"])
    ap.add_argument("--depth", type=int, default=None, help="fixed rarefaction depth in templates (default: min of each pair)")
    ap.add_argument("--n-iter", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    if not args.data_dir or not args.metrics_csv:
        ap.error("--data-dir and --metrics-csv are required (or set TCR_DATA_DIR / TCR_METRICS_CSV)")
    out_path = args.out or f"outputs/tracking_{args.level}.csv"

    ref = pd.read_csv(args.metrics_csv, sep=None, engine="python")
    ref["patient"] = ref["Patient"].str.lower().str.replace("pt", "", regex=False).astype(int)
    ref = ref.set_index("patient")[["Cohort", "myBOR"]]

    rows = []
    for patient, files in find_pairs(args.data_dir).items():
        pre = repertoire(productive_filter(load_sample(files["pre"])), args.level)
        on = repertoire(productive_filter(load_sample(files["on"])), args.level)
        a, b = align(pre, on)
        rec = {"patient": patient, **{f"raw_{k}": v for k, v in replacement_metrics(a, b).items()}}
        rec["depth_ratio_on_over_pre"] = rec["raw_on_templates"] / rec["raw_pre_templates"]
        depth = args.depth if args.depth else int(min(a.sum(), b.sum()))
        if a.sum() >= depth and b.sum() >= depth:
            rec.update(rarefied_metrics(a, b, depth, n_iter=args.n_iter, seed=args.seed))
        else:
            rec["rar_depth"] = depth  # a sample is shallower than the requested depth: rarefied values left empty
        rows.append(rec)
        print(f"Pt{patient}: pre={rec['raw_pre_templates']:>7,} on={rec['raw_on_templates']:>7,} templates; "
              f"on-from-baseline raw={rec['raw_frac_on_templates_from_baseline']:.3f} "
              f"rarefied={rec.get('rar_frac_on_templates_from_baseline', float('nan')):.3f}", flush=True)

    out = pd.DataFrame(rows).set_index("patient").join(ref)
    out["benefit"] = out["myBOR"].map({"PRCR": "benefit", "SD": "benefit", "PD": "no_benefit"})
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path)

    pd.set_option("display.width", 200)
    print(f"\nlevel={args.level}  patients={len(out)}  (rarefied: {out['rar_frac_on_templates_from_baseline'].notna().sum()})")
    cols = ["raw_frac_on_templates_from_baseline", "rar_frac_on_templates_from_baseline",
            "rar_frac_on_clonotypes_from_baseline", "rar_morisita_horn"]
    print("\nMedians by cohort and benefit group:")
    g = out.groupby(["Cohort", "benefit"])
    print(pd.concat([g.size().rename("n"), g[cols].median().round(3)], axis=1).to_string())
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
