# Bulk RNA-seq in nivolumab-treated melanoma: a reanalysis of Riaz et al. (2017)

Bulk RNA sequencing part of our group project integrating RNA, WES and TCR data from
Riaz et al., *Cell* 2017 ([DOI 10.1016/j.cell.2017.09.028](https://doi.org/10.1016/j.cell.2017.09.028)).

**Analysis by:** Ethan (RNA lead). Group: Krithika (TCR), Ethan (RNA), Kwaku (WES).

## Question

What changes in a melanoma's gene expression when nivolumab works? We compare tumor biopsies taken
before treatment and about 4 weeks into treatment, across response groups:
PR/CR (partial or complete response), SD (stable disease) and PD (progressive disease).

## Findings

| Measure, median | PR/CR | SD | PD | PR/CR − PD (95% CI) | Wilcoxon p |
|---|---|---|---|---|---|
| Baseline T-cell-inflamed score (n = 10 / 16 / 23) | 0.19 | −0.47 | −0.34 | +0.53 (−0.73 to +1.33) | 0.20 |
| Share of genes changing 2-fold (n = 9 / 15 / 18) | 8.6% | 4.7% | 4.5% | +4.2 points (0.4 to 16.0) | 0.048 |
| Change in tumor score (n = 9 / 15 / 18) | −1.07 | −0.17 | −0.01 | −1.05 (−1.77 to −0.23) | 0.0002 |
| Change in CD8 score (n = 9 / 15 / 18) | +0.47 | +0.31 | 0.00 | +0.47 (0.00 to +1.20) | 0.07 |
| On-treatment CD8 score (n = 12 / 18 / 24) | 1.35 | 0.20 | −0.31 | +1.66 (+0.75 to +2.15) | 0.0004 |

Genes changing On vs Pre within each group (DESeq2, FDR < 0.1): **3,582 PR/CR, 28 SD, 0 PD**.
A direct test finds 1,721 genes whose change differs between PR/CR and PD.

- Baseline immune scores do not separate responders: the interval crosses zero.
- Responders lose melanoma and proliferation programs. The tumor score tracks Kwaku's mutation
  contraction (Spearman rho = −0.69, n = 25).
- Immune pathways rise in every group, including SD. On-treatment CD8 separates PR/CR from PD, and
  the RNA CD8 score tracks the TCR-seq T-cell fraction (rho = 0.88, n = 30).
- These are exploratory results from small groups. Stable disease is not part of the tests.
  Confidence intervals are percentile bootstraps (10,000 resamples, seed 2017), as in the WES analysis.

## How to run

Requires R (tested with 4.5.2). From this folder (`RNA_Ethan/`), with the shared `../data/` and
`../WES_Kwaku/` folders in place:

```bash
Rscript rna/setup.R
Rscript rna/run_rna.R
python3 slides/build_slides.py
```

Or open `RNA_analysis.Rmd` in RStudio and click Knit. It runs every step and writes
`RNA_analysis.html` (included here). The first run takes about 20 minutes; later runs reuse the
DESeq2 fits cached in `rna/results/cache/` (about 1 GB, not committed).

## Folder layout

| Path | Contents |
|---|---|
| `rna/setup.R` | One-time package install |
| `rna/run_rna.R` | Runs all steps in order |
| `rna/01`–`05_*.R` | QC and normalization, DESeq2, pathways, signature scores, links to WES/TCR/outcome |
| `rna/rna_bootstrap_ci.R` | Bootstrap CIs, Wilcoxon tests, tumor-dilution check |
| `rna/paper_comparison.R` | Reruns the authors' designs and compares gene counts with the paper |
| `rna/results/` | Output tables (`.csv`) and slide figures `01`–`05` (`.png`, `.pdf`) |
| `rna/results/supplementary_figures/` | QC, volcano plots, heatmaps, trajectories, AUC and survival (`S01`–`S13`) |
| `slides/RNA_11slides_Ethan.pptx` | Slides with speaker notes, built on the WES deck's theme |
| `slides/presentation_script.md` | The presenter script |
| `slides/build_slides.py` | Builds both slide files from `rna/results/` |

Key tables: `rna_patient_features.csv` (one row per patient, for the group's integration),
`de_summary.csv`, `gsea_hallmark.csv`, `bootstrap_median_ci.csv`, `cross_modality_correlations.csv`,
and `dilution_check.csv`.

## How figures match the WES analysis

| WES figure | RNA figure | Idea |
|---|---|---|
| `01` baseline mutations by response | `01_baseline_immune_by_response` | Does the starting state predict response? |
| `02` mutation changes per patient | `02_gene_changes_per_patient` | What changed between biopsies, patient by patient? |
| `03` contraction by response | Stat panels on slides 6–7 | Summary of that change by group |
| `04` neoantigen ranking | `03_tumor_immune_change_per_patient` | Ranked per-patient result, with a reference line |
| No equivalent | `04_pathways_by_response` | Which biological programs moved? |
| No equivalent | `05_rna_vs_wes_tcr` | Do RNA, WES and TCR agree? |

## Comparison with the paper

| Finding | Paper | Here |
|---|---|---|
| Immune genes up, cell-cycle and melanin programs down in responders | Yes | Yes |
| Baseline immune activity only marginally predictive | Yes | Yes |
| "Hot" baseline in responders who progressed on ipilimumab | Yes | Yes |
| Pre PR/CR vs PD, genes at q < 0.2 | 189 | 397 with authors' design and gene filter; 460 with ours |
| On vs Pre, genes at q < 0.2 | 475 | 157 with authors' design; 227 with ours |
| Change on treatment, benefit vs no benefit, q < 0.2 | 2,670 | 7 (design from the paper's Methods); 3,479 for PR/CR vs PD |

We could not match the paper's gene counts. Differences include DESeq2 version (1.12.4 vs 1.50.2)
and IHW, which was not used here for Table S6A. The repository has no script for Table S6D.

## Choices and differences from the authors

- Patient IDs come from the sample name, because the table's `PatientID` mislabels Pt7_Pre as Pt72.
- `.2` libraries in the count table are not in the sample table, so they are left out of the analysis.
  We use them only as a noise check, assuming they are repeat libraries of the same biopsy (median
  difference 3.1% of genes, up to 15%). The repository does not document them.
- Pt23 has no response call, so it is left out. 25 of the 26 shared patients are analyzed.
- The paper removed Pt3 from transcriptome analysis, though the authors' public scripts keep it.
  We keep it and report tests without it in `wilcoxon_PRCR_vs_PD.csv`.
- We use marker-gene scores rather than the authors' CIBERSORT deconvolution.
- We use MSigDB Hallmark sets with fgsea rather than the authors' clusterProfiler and IPA.

## Limitations

- Bulk RNA is a mixture. Tumor loss raises the share of every other cell type: adjusting for the tumor
  score cuts the PR/CR CD8 rise from +0.67 to +0.21. The on-treatment CD8 level holds (+1.30 to +1.21).
- Response is defined by tumor shrinkage, so the tumor score partly restates the outcome.
- This analysis starts from the authors' count table. Reads were not realigned.
- Enrichment p-values are optimistic because genes in a pathway are correlated.
- Figures use all 42 patients with paired RNA, not only the 26 shared across all three data types.
