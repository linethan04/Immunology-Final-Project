"""Reproduce the group statistics of Riaz et al. (2017) Fig. 5B from the authors' TCRmetrics.csv.

Benefit (CR/PR/SD) vs no benefit (PD), two-sided Mann-Whitney, separately for each cohort:
  * richness = fold change in number of unique CDR3s (all patients)
  * evenness = change in CDR3 evenness, excluding the two patients the authors' createTCRfigures.r removes
    ("two outliers per Grubbs": pt103 in the Ipi-progressed cohort, pt28 in the Ipi-naive cohort)

Published values are copied from the paper text (Results, "On-Therapy Changes in Intratumoral T Cell Repertoire
Diversity"). Usage: python scripts/reproduce_fig5b_stats.py --metrics-csv /path/to/TCRmetrics.csv
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tcr_tracking.diversity import mwu_p  # noqa: E402

PUBLISHED = {
    ("richness fold change", "NIV3-PROG"): 0.016, ("richness fold change", "NIV3-NAIVE"): 0.489,
    ("evenness change", "NIV3-NAIVE"): 0.036, ("evenness change", "NIV3-PROG"): 0.594,
}
EXCLUDE = {"NIV3-PROG": "pt103", "NIV3-NAIVE": "pt28"}      # from the authors' createTCRfigures.r (evenness only)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metrics-csv", default=os.environ.get("TCR_METRICS_CSV"))
    args = ap.parse_args()
    if not args.metrics_csv:
        ap.error("--metrics-csv is required (or set TCR_METRICS_CSV)")
    t = pd.read_csv(args.metrics_csv, sep=None, engine="python")
    t["pid"] = t["Patient"].str.lower()
    t["benefit"] = np.where(t["myBOR"] == "PD", "no_benefit", "benefit")

    rows = []
    for (label, cohort), pub in PUBLISHED.items():
        col = "Num_CDR3_FoldChange" if label.startswith("richness") else "evenness.cdr3_On-Pre"
        s = t[t["Cohort"] == cohort]
        if label.startswith("evenness"):
            s = s[s["pid"] != EXCLUDE[cohort]]
        a, b = s.loc[s["benefit"] == "benefit", col], s.loc[s["benefit"] == "no_benefit", col]
        p = mwu_p(a, b)
        rows.append({"analysis": label, "cohort": cohort, "n_benefit": len(a), "n_no_benefit": len(b),
                     "p_reproduced": round(p, 3), "p_published": pub, "abs_diff": round(abs(p - pub), 3)})
    print(pd.DataFrame(rows).to_string(index=False))

    # sensitivity to the outlier exclusion (not part of the published analysis)
    print("\nEvenness change WITHOUT the exclusions:")
    for cohort in EXCLUDE:
        s = t[t["Cohort"] == cohort]
        a, b = s.loc[s["benefit"] == "benefit", "evenness.cdr3_On-Pre"], s.loc[s["benefit"] == "no_benefit", "evenness.cdr3_On-Pre"]
        print(f"  {cohort}: p = {mwu_p(a, b):.3f}")


if __name__ == "__main__":
    main()
