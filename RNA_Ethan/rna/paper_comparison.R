# Rerun the authors' RNA comparisons as closely as this environment allows and compare gene counts.
# Authors' scripts (rnaseq/TableS6.A.R, TableS6.C.R) use their gene filter, DESeq2 1.12.4 and, for S6A, IHW.
# Table S6D has no script in their repository, so that row uses the design from the paper's Methods.
source("rna/common.R")
suppressPackageStartupMessages({library(genefilter); library(matrixStats)})

mat <- read.delim(file.path(data_dir, "CountData.BMS038.txt"))
tab <- read.csv(file.path(data_dir, "SampleTableCorrected.9.19.16.csv"), row.names = 1)
inter <- intersect(colnames(mat), rownames(tab)); tab <- tab[inter, ]
mat <- as.matrix(mat[, inter])
base <- DESeqDataSetFromMatrix(mat, tab, ~ 1)

# Authors' filter, condensed from filterthegenes() in TableS6.A.R
authors_filter <- function(x) {
  x <- estimateSizeFactors(x[, !is.na(x$Response)])
  low <- rownames(x)[!genefilter(counts(x), filterfun(kOverA(k = 4, A = 2, na.rm = TRUE)))]
  m <- estimateSizeFactors(x[!rownames(x) %in% low, ]); n <- counts(m, normalized = TRUE)
  dl <- n[, (ncol(m) / 2 + 1):ncol(m)] - n[, 1:(ncol(m) / 2)]
  v1 <- rownames(n)[log10(rowVars(dl)) < 0]
  m <- estimateSizeFactors(x[!rownames(x) %in% c(low, v1), ]); n <- counts(m, normalized = TRUE)
  v2 <- rownames(n)[log10(rowVars(n)) < 0 | log10(rowMads(n)) < 0]
  unique(c(low, v1, v2))
}

# Table S6A: Pre, PR/CR vs PD
pre <- base[, base$PreOn == "Pre" & !is.na(base$Response) & !colnames(base) %in% c("Pt48_Pre", "Pt67_Pre", "Pt52_Pre", "Pt27_Pre")]
pre <- pre[!rownames(pre) %in% authors_filter(pre), ]
pre$Response <- factor(pre$Response, c("PD", "SD", "PRCR")); design(pre) <- ~ Response
fit_S6A <- cached("fit_S6A", DESeq(pre, quiet = TRUE))
n_S6A <- sum(results(fit_S6A, name = "Response_PRCR_vs_PD")$padj < 0.2, na.rm = TRUE)

# Table S6C: On vs Pre with a patient term, all samples, authors' PatientID column as supplied
onpre <- base; onpre$PreOn <- factor(onpre$PreOn, c("Pre", "On")); onpre$PatientID <- factor(onpre$PatientID)
design(onpre) <- ~ PatientID + PreOn
fit_S6C <- cached("fit_S6C", DESeq(onpre, quiet = TRUE))
n_S6C <- sum(results(fit_S6C, contrast = c("PreOn", "On", "Pre"))$padj < 0.2, na.rm = TRUE)

ours <- read_csv(file.path(out, "de_summary.csv"), show_col_types = FALSE)
get <- function(n) ours$fdr_0.2[ours$contrast == n]
comparison <- tibble(
  comparison = c("Pre: PR/CR vs PD", "On vs Pre, all patients", "On vs Pre change: benefit vs no benefit", "On vs Pre change: PR/CR vs PD"),
  paper_q0.2 = c(189, 475, 2670, NA),
  authors_design_this_R = c(n_S6A, n_S6C, NA, NA),
  our_pipeline_fdr0.2 = c(get("A_pre_PRCR_vs_PD"), get("B_on_vs_pre_paired"), get("D_on_vs_pre_benefit_vs_none"), get("C_on_vs_pre_PRCR_vs_PD")))
write_csv(comparison, file.path(out, "paper_comparison.csv"))
print(comparison)
