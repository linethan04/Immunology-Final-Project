# Ethan sandbox: bulk RNA-seq, Riaz et al. 2017 (BMS038, nivolumab)
# Run from the project root (open Immunology-Final-Project.Rproj).
# 1) DESeq2 pre-treatment PRCR vs PD (paper Table S6A)
# 2) DESeq2 on- vs pre-treatment, paired by patient (paper Table S6C)
# 3) Immune signature scores per sample -> per-patient table for integration

suppressPackageStartupMessages({
  library(DESeq2)
  library(org.Hs.eg.db)
  library(tidyverse)
})
out <- "output/rna"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# ---- Load ----
counts <- read.delim("data/CountData.BMS038.txt", check.names = FALSE)  # rows = Entrez IDs

samples <- read.csv("data/SampleTableCorrected.9.19.16.csv") |>
  transmute(sample   = Sample,
            patient  = sub("_.*", "", Sample),  # PatientID column mislabels Pt7_Pre as Pt72
            time     = factor(PreOn, c("Pre", "On")),
            response = factor(Response, c("PD", "SD", "PRCR")),
            benefit  = Response.1,              # LB / NB
            cohort   = Cohort) |>
  filter(!is.na(response), sample %in% colnames(counts))

# The 26 patients with WES + TCR-seq + RNA-seq at both timepoints
geno <- read.csv("data/genomic_data_per_case.csv")
tri  <- geno$Patient[rowSums(geno[-1]) == 6]
samples$tri <- samples$patient %in% tri

dds <- DESeqDataSetFromMatrix(as.matrix(counts[, samples$sample]),
                              column_to_rownames(samples, "sample"), ~ 1)
dds <- dds[rowSums(counts(dds) >= 10) >= 5, ]  # simpler than the paper's filter
sym <- mapIds(org.Hs.eg.db, rownames(dds), "SYMBOL", "ENTREZID")

tidy_res <- function(r) as.data.frame(r) |>
  rownames_to_column("entrez") |>
  mutate(symbol = sym[entrez]) |>
  arrange(padj)

# ---- 1. Pre-treatment: PRCR vs PD (SD kept as its own level, as in the paper) ----
pre <- dds[, dds$time == "Pre" &
             !colnames(dds) %in% c("Pt48_Pre", "Pt67_Pre", "Pt52_Pre", "Pt27_Pre")]  # paper exclusions
design(pre) <- ~ response
pre <- DESeq(pre)
res_pre <- tidy_res(results(pre, name = "response_PRCR_vs_PD"))
write_csv(res_pre, file.path(out, "deseq_pre_PRCR_vs_PD.csv"))

# ---- 2. On vs Pre, paired (patients with both timepoints) ----
paired <- names(which(table(dds$patient) == 2))
onpre  <- dds[, dds$patient %in% paired]
onpre$patient <- factor(onpre$patient)
design(onpre) <- ~ patient + time
onpre <- DESeq(onpre)
res_onpre <- tidy_res(results(onpre, name = "time_On_vs_Pre"))
write_csv(res_onpre, file.path(out, "deseq_on_vs_pre.csv"))

cat("Pre PRCR vs PD, padj < 0.1:", sum(res_pre$padj < 0.1, na.rm = TRUE), "genes\n")
cat("On vs Pre,      padj < 0.1:", sum(res_onpre$padj < 0.1, na.rm = TRUE), "genes\n")
head(res_onpre$symbol, 25) |> unname()

volcano <- function(res, title) {
  ggplot(res, aes(log2FoldChange, -log10(pvalue), colour = padj < 0.1)) +
    geom_point(size = 0.6, alpha = 0.6) +
    geom_text(data = head(res, 15), aes(label = symbol), size = 2.5,
              vjust = -0.6, colour = "black", check_overlap = TRUE) +
    scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick"),
                        na.value = "grey85") +
    labs(title = title, colour = "padj < 0.1") +
    theme_bw()
}
ggsave(file.path(out, "volcano_pre_PRCR_vs_PD.png"),
       volcano(res_pre, "Pre-treatment: PRCR vs PD"), width = 6, height = 5)
ggsave(file.path(out, "volcano_on_vs_pre.png"),
       volcano(res_onpre, "On vs Pre (paired)"), width = 6, height = 5)

# ---- 3. Immune signature scores ----
# Score = mean of per-gene z-scores of VST expression across all samples.
# (Marker-based proxy for CD8 content; swap in proper deconvolution later.)
sigs <- list(
  CD8       = c("CD8A", "CD8B"),
  cytolytic = c("GZMA", "PRF1"),
  IFNG_6    = c("IFNG", "STAT1", "IDO1", "CXCL9", "CXCL10", "HLA-DRA"),  # Ayers 2017
  Tcell_inflamed_18 = c("CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9",
                        "CXCR6", "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3",
                        "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT")
)
vsd <- assay(vst(dds, blind = TRUE))
z   <- t(scale(t(vsd)))
scores <- sapply(sigs, \(g) colMeans(z[na.omit(match(g, sym)), , drop = FALSE]))

score_long <- as.data.frame(scores) |>
  rownames_to_column("sample") |>
  left_join(samples, by = "sample") |>
  pivot_longer(names(sigs), names_to = "signature", values_to = "score")

ggplot(score_long, aes(time, score, fill = response)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1), size = 0.8) +
  facet_wrap(~ signature, scales = "free_y") +
  theme_bw()
ggsave(file.path(out, "signatures_by_response.png"), width = 8, height = 6)

# Per-patient table (pre, on, delta) for the integration step
rna_patient <- score_long |>
  select(patient, response, benefit, cohort, tri, time, signature, score) |>
  pivot_wider(names_from = c(signature, time), values_from = score)
for (s in names(sigs)) {
  rna_patient[[paste0(s, "_delta")]] <- rna_patient[[paste0(s, "_On")]] - rna_patient[[paste0(s, "_Pre")]]
}
write_csv(rna_patient, file.path(out, "rna_scores_per_patient.csv"))

# Paired change in each score, split by response (26 tri-modal patients)
rna_patient |>
  filter(tri) |>
  pivot_longer(ends_with("_delta"), names_to = "signature", values_to = "delta") |>
  group_by(signature, response) |>
  summarise(n = sum(!is.na(delta)), median_delta = median(delta, na.rm = TRUE),
            .groups = "drop") |>
  print(n = Inf)

# Responders (PRCR) vs PD: pre-treatment score and on-treatment change
rna_patient |>
  filter(tri, response != "SD") |>
  pivot_longer(matches("_(Pre|delta)$"), names_to = "feature", values_to = "value") |>
  group_by(feature) |>
  summarise(p = wilcox.test(value ~ droplevels(response))$p.value) |>
  print()
