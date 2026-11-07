"""Load Adaptive immunoSEQ tables and apply the filtering used by Riaz et al. (2017).

The authors' JS_splitter_Imsqv4.5.py (Python 2.7) is re-implemented here:
  * counts come from the ``count (templates)`` column
  * V and J names resolve gene -> gene ties -> family -> family ties
  * the "only" productive filter keeps a CDR3 amino-acid string if it is non-empty,
    starts with "C" and contains no "*"
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

SAMPLE_RE = re.compile(r"^Pt(\d+)_(pre|on)\.tsv$")

USECOLS = [
    "nucleotide", "aminoAcid", "count (templates)", "sequenceStatus",
    "vGeneName", "vGeneNameTies", "vFamilyName", "vFamilyTies",
    "jGeneName", "jGeneNameTies", "jFamilyName", "jFamilyTies",
]


def find_pairs(data_dir: str | Path) -> dict[int, dict[str, Path]]:
    """Return {patient_number: {"pre": path, "on": path}} for patients with both timepoints.

    Names that do not match ``Pt<N>_pre.tsv`` / ``Pt<N>_on.tsv`` (for example the macOS
    ``._*`` files) are ignored.
    """
    found: dict[int, dict[str, Path]] = {}
    for f in sorted(Path(data_dir).iterdir()):
        m = SAMPLE_RE.match(f.name)
        if m:
            found.setdefault(int(m.group(1)), {})[m.group(2)] = f
    return {p: d for p, d in sorted(found.items()) if {"pre", "on"} <= d.keys()}


def _first_tie(s: pd.Series) -> pd.Series:
    return s.fillna("").astype(str).str.split(",").str[0].str.strip().str.strip('"')


def resolve_gene(gene: pd.Series, gene_ties: pd.Series, family: pd.Series, family_ties: pd.Series) -> pd.Series:
    """Most specific TCR gene designation available, as in the authors' splitter."""
    out = gene.fillna("").astype(str)
    ok = out.str.startswith("TCR")
    for cand in (_first_tie(gene_ties), family.fillna("").astype(str), _first_tie(family_ties)):
        use = ~ok & cand.str.startswith("TCR")
        out = out.where(~use, cand)
        ok = ok | use
    return out.where(ok, "unresolved")


def load_sample(path: str | Path) -> pd.DataFrame:
    """Read one immunoSEQ table (all rows, unfiltered) with tidy columns.

    Columns: nt, aa, v, j, templates, status
    """
    df = pd.read_csv(path, sep="\t", usecols=USECOLS, low_memory=False)
    return pd.DataFrame({
        "nt": df["nucleotide"],
        "aa": df["aminoAcid"].fillna("").astype(str),
        "v": resolve_gene(df["vGeneName"], df["vGeneNameTies"], df["vFamilyName"], df["vFamilyTies"]),
        "j": resolve_gene(df["jGeneName"], df["jGeneNameTies"], df["jFamilyName"], df["jFamilyTies"]),
        "templates": df["count (templates)"].astype(int),
        "status": df["sequenceStatus"],
    })


def productive_filter(df: pd.DataFrame) -> pd.DataFrame:
    """Authors' "only" filter on the amino-acid CDR3: non-empty, starts with C, no stop codon."""
    aa = df["aa"]
    keep = (aa != "") & aa.str.startswith("C") & ~aa.str.contains("*", regex=False)
    return df[keep]


def repertoire(df: pd.DataFrame, level: str = "aa") -> pd.Series:
    """Collapse to clonotypes and sum templates.

    level="aa"    -> clonotype = CDR3 amino-acid sequence (authors' aaCDR3 resolution)
    level="aa_vj" -> clonotype = (CDR3 amino acid, V gene, J gene) (authors' totCDR3 resolution)
    Returns a Series of template counts indexed by clonotype, sorted descending.
    """
    if level == "aa":
        keys = ["aa"]
    elif level == "aa_vj":
        keys = ["aa", "v", "j"]
    else:
        raise ValueError(f"unknown level: {level!r}")
    s = df.groupby(keys)["templates"].sum()
    s = s[s > 0]
    return s.sort_values(ascending=False)
