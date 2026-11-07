"""Compare locally computed repertoire metrics with the authors' TCRmetrics.csv.

Usage:
    python scripts/validate_metrics.py --data-dir ~/data/bms038_raw \
        --metrics-csv /path/to/bms038_analysis/data/tcr/TCRmetrics.csv

Both arguments can instead be set with the environment variables TCR_DATA_DIR and TCR_METRICS_CSV.
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
from tcr_tracking.metrics import summarize  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", default=os.environ.get("TCR_DATA_DIR"))
    ap.add_argument("--metrics-csv", default=os.environ.get("TCR_METRICS_CSV"))
    ap.add_argument("--out", default="outputs/validation_aa.csv")
    ap.add_argument("--level", default="aa", choices=["aa", "aa_vj"])
    args = ap.parse_args()
    if not args.data_dir or not args.metrics_csv:
        ap.error("--data-dir and --metrics-csv are required (or set TCR_DATA_DIR / TCR_METRICS_CSV)")

    # TCRmetrics.csv uses old-Mac (CR) line endings, so let the python engine sniff the format.
    ref = pd.read_csv(args.metrics_csv, sep=None, engine="python")
    ref["patient"] = ref["Patient"].str.lower().str.replace("pt", "", regex=False).astype(int)
    ref = ref.set_index("patient")

    rows = []
    for patient, files in find_pairs(args.data_dir).items():
        rec = {"patient": patient}
        for tp, label in (("pre", "Pre"), ("on", "On")):
            counts = repertoire(productive_filter(load_sample(files[tp])), args.level)
            m = summarize(counts)
            rec[f"templates_{label}"] = m["templates"]
            rec[f"richness_{label}"] = m["richness"]
            rec[f"clonality_{label}"] = m["clonality"]
            rec[f"D90_{label}"] = m["D90"]
            if patient in ref.index:
                rec[f"ref_richness_{label}"] = ref.loc[patient, f"Num_CDR3_{label}"]
                rec[f"ref_clonality_{label}"] = ref.loc[patient, f"CLcdr3_{label}"]
                rec[f"ref_D90_{label}"] = ref.loc[patient, f"D90_{label}"]
        rows.append(rec)
    out = pd.DataFrame(rows).set_index("patient")

    for label in ("Pre", "On"):
        out[f"richness_reldiff_{label}"] = (out[f"richness_{label}"] - out[f"ref_richness_{label}"]) / out[f"ref_richness_{label}"]
        out[f"clonality_absdiff_{label}"] = out[f"clonality_{label}"] - out[f"ref_clonality_{label}"]
        out[f"D90_absdiff_{label}"] = out[f"D90_{label}"] - out[f"ref_D90_{label}"]

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out)

    pd.set_option("display.width", 200)
    print(f"level={args.level}; patients compared: {len(out)}\n")
    for metric in ("richness", "clonality", "D90"):
        for label in ("Pre", "On"):
            a, b = out[f"{metric}_{label}"], out[f"ref_{metric}_{label}"]
            ok = a.notna() & b.notna()
            r = np.corrcoef(a[ok], b[ok])[0, 1]
            mad = (a[ok] - b[ok]).abs().max()
            print(f"{metric:10s} {label:3s}  n={int(ok.sum()):2d}  pearson r={r:.5f}  max |diff|={mad:.4g}")
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
