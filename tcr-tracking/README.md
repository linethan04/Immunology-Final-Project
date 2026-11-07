# tcr-tracking

The authors' data and code are **not** included here. Need to be cloned separately and point this
project at them:

    git clone https://github.com/riazn/bms038_analysis      # reference only (no LICENSE file in that repo)

## Setup

    python3 -m venv .venv && source .venv/bin/activate
    pip install -r requirements.txt

Unpack the raw immunoSEQ tables (inspect first with `tar -tzf`):

    mkdir -p ~/data/bms038_raw && chmod 700 ~/data/bms038_raw
    tar -xzf bms038_analysis/data/tcr/Immunoseq.tgz -C ~/data/bms038_raw --strip-components=7 --exclude='._*'
    ls ~/data/bms038_raw/Pt*_*.tsv | wc -l      # expect 58 (29 patients x pre/on)

Set the paths (or pass `--data-dir` / `--metrics-csv` on the command line):

    export TCR_DATA_DIR=~/data/bms038_raw
    export TCR_METRICS_CSV=/path/to/bms038_analysis/data/tcr/TCRmetrics.csv

## Step 1: validate against the authors' published metrics

    python scripts/validate_metrics.py

Compares richness (`Num_CDR3`), clonality (`CLcdr3`) and `D90` with `TCRmetrics.csv` for the 29
patients that have raw tables (the table lists 34; patients 84, 89, 92, 93, 94 have no raw data in
the repo). Expected result on the repo data: richness identical; clonality within 0.0005 (the
reference table is rounded to 3 decimals); D90 within about 0.0002.

## Conventions (re-implemented from the authors' Python 2.7 scripts)

- Counts: `count (templates)` column.
- Productive filter ("only"): amino-acid CDR3 non-empty, starts with `C`, contains no `*`.
- Entropy in bits; evenness = H / log2(N); clonality = 1 - evenness.
- D90: number of top clonotypes with cumulative template share still below 90%, divided by the
  number of unique clonotypes (`include_crossing_clone=True` gives the strict minimum, one higher).
- Clonotype levels: `aa` (authors' aaCDR3) and `aa_vj` (CDR3 amino acid + V gene + J gene).

## Layout

    tcr_tracking/load.py       read tables, filters, clonotype collapsing, pre/on pairing
    tcr_tracking/metrics.py    richness, Shannon, evenness, clonality, D90
    scripts/validate_metrics.py
    scripts/make_all_figures.py is the one that makes the final figures used in powerpoint, and those are found under output -> figures folder 

