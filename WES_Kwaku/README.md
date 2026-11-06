# WES in nivolumab-treated melanoma: a reanalysis of Riaz et al. (2017)

Whole exome sequencing (WES) part of our group project integrating RNA, WES and TCR data from
Riaz et al., *Cell* 2017 ([DOI 10.1016/j.cell.2017.09.028](https://doi.org/10.1016/j.cell.2017.09.028)).

**Analysis by:** Kwaku Darkwah (WES lead). Group: Krithika (TCR), Ethan (RNA), Kwaku (WES).

## Question

When nivolumab works, do melanomas lose their mutations? We compare tumor biopsies taken before
treatment and about 4 weeks into treatment, across response groups:
PR/CR (partial or complete response), SD (stable disease) and PD (progressive disease).

## Findings

| Measure | PR/CR | SD | PD | PR/CR − PD (95% CI) |
|---|---|---|---|---|
| Baseline nonsynonymous mutations, median (n = 15 / 23 / 30) | 217 | 214 | 162 | +55 (−95 to +393) |
| Share of mutations contracting on treatment, median (n = 9 / 13 / 19) | 0.91 | 0.26 | 0.09 | +0.83 (0.48 to 0.88) |
| Neoantigens observed ÷ expected on treatment, median (n = 9 / 13 / 19) | 0.43 | 1.01 | 1.03 | −0.60 (−0.70 to −0.40) |

- Baseline mutation burden alone does not separate responders: the interval crosses zero.
- In responders, most mutations contract (are lost, or now found in fewer cancer cells), and
  fewer neoantigens remain than expected from silent mutations. Wilcoxon PR/CR vs PD:
  p ≈ 0.0001 (contraction) and p ≈ 0.0008 (neoantigens).
- These are exploratory results from small groups. Stable disease is not part of the tests.
  Confidence intervals are percentile bootstraps (10,000 resamples, seed 2017).

## How to run

Requires R (tested with 4.6.1) and git. From the repository root:

```bash
Rscript wes/setup.R
Rscript wes/run_wes.R
Rscript wes/wes_bootstrap_ci.R
```

Or run everything as one report: open `WES_analysis.Rmd` in RStudio and click Knit. It does the
setup, analysis and confidence intervals, and writes `WES_analysis.html` (included here).

`setup.R` downloads the authors' repository [riazn/bms038_analysis](https://github.com/riazn/bms038_analysis)
at commit `137111c38a32f8626ee6a436e14ac667bd31b1e7` into `external/`, extracts their `nrstats`
helper package, and installs any missing CRAN packages (ggplot2, survival, RColorBrewer, ks).
The data are not copied into this repository.

## Repository layout

| Path | Contents |
|---|---|
| `wes/setup.R` | One-time download of the authors' code and data, plus package install |
| `wes/run_wes.R` | Runs the authors' exome scripts, then writes the figures and per-patient tables |
| `wes/wes_bootstrap_ci.R` | Bootstrap confidence intervals and the Wilcoxon tests quoted on the slides |
| `wes/results/` | Output tables (`.csv`), slide figures `01`–`04` (`.png`, `.pdf`) and the regenerated original figures |
| `slides/WES_10slides_Kwaku.pptx` | Presentation slides with speaker notes |

Key tables: `paired_mutation_features.csv` (per-patient contraction), `neoantigen_patient_features.csv`
(observed/expected neoantigens), `bootstrap_median_ci.csv`, and `all_modality_paired_patient_ids.csv`
(the 26 patients with WES, RNA and TCR at both timepoints, for the group's integration).

## Changes from the authors' code

The upstream scripts are run unchanged on disk. `run_wes.R` adjusts them in memory only:

- Points `PARENTDIR` at the local copy and removes `rm(list = ls())` so results can be collected.
- Replaces the two `sqldf` `COUNT(*) ... GROUP BY` queries with a plain-R equivalent.
  `sqldf`'s compiled dependencies are blocked by Windows Smart App Control on some machines.
- Sets string sorting to the C locale. Upstream code depends on `merge()` row order, which
  otherwise varies by system locale and breaks a column assignment in `baseline_clonality.r`.
- Uses a current ggplot2 replacement for the old `nrstats` boxplot helper (keeping its Welch and
  Wilcoxon tests), sources only the `nrstats` survival functions, and raises PNG resolution.

## Limitations

- This reproduces the analysis from the authors' processed data. It does not re-call mutations
  or re-predict neoantigens from raw reads.
- "Contraction" uses the authors' cancer cell fraction categories at a 0.95 threshold. It includes
  mutations that went from clonal to subclonal, not only those that disappeared.
- More immune cells in an on-treatment biopsy can make mutations look lost. Cancer cell fractions
  adjust for tumor purity, but low purity and read depth remain caveats.
- The figures use all 41 patients with paired WES, not only the 26 shared across all three data types.
