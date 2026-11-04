# RNA step 4: immune signature / cell-type scores per sample -> per-patient table
# Score = mean z-scored VST expression of marker genes (MCP-counter-style markers,
# Ayers 2017 IFN-g/T-cell-inflamed sets). Melanocytic score = tumour-content proxy.

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
})
out <- "output/rna"; fig <- "figures/rna"
theme_set(theme_bw(base_size = 11))
resp_cols <- c(PD = "#C0392B", SD = "#E1A23B", PRCR = "#2471A3")
vsd     <- readRDS(file.path(out, "vsd.rds"))
samples <- read_csv(file.path(out, "samples.csv"), show_col_types = FALSE) |>
  mutate(time = factor(time, c("Pre", "On")), response = factor(response, c("PD", "SD", "PRCR")))
sym <- rowData(vsd)$symbol

sigs <- list(
  # cell types (MCP-counter-like marker sets)
  T_cells     = c("CD3D", "CD3E", "CD3G", "CD28", "CD6", "TRAT1"),
  CD8_T       = c("CD8A", "CD8B"),
  Cytotoxic   = c("GZMA", "GZMB", "GZMH", "PRF1", "NKG7", "GNLY"),
  NK          = c("KLRD1", "KLRF1", "NCR1", "SH2D1B"),
  B_cells     = c("CD19", "MS4A1", "CD79A", "CR2", "FCRL2"),
  Myeloid     = c("CD14", "CD163", "CD68", "CSF1R", "FCGR3A"),
  Treg        = c("FOXP3", "IL2RA", "CTLA4"),
  # activity programs
  IFNG_6      = c("IFNG", "STAT1", "IDO1", "CXCL9", "CXCL10", "HLA-DRA"),
  Tcell_inflamed_18 = c("CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9",
                        "CXCR6", "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3",
                        "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT"),
  Checkpoints = c("PDCD1", "CD274", "LAG3", "HAVCR2", "TIGIT", "CTLA4"),
  Antigen_pres = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9"),
  # tumour content
  Melanocytic = c("MITF", "PMEL", "TYR", "MLANA", "DCT")
)
z <- t(scale(t(assay(vsd))))
scores <- sapply(sigs, \(g) colMeans(z[na.omit(match(g, sym)), , drop = FALSE]))
write_csv(as_tibble(scores, rownames = "sample"), file.path(out, "signature_scores_per_sample.csv"))

long <- as_tibble(scores, rownames = "sample") |>
  left_join(samples, by = "sample") |>
  pivot_longer(names(sigs), names_to = "signature", values_to = "score") |>
  mutate(signature = factor(signature, names(sigs)))

# ---- Fig 9: paired Pre -> On trajectories per patient, by response ----
key <- c("CD8_T", "Cytotoxic", "IFNG_6", "Tcell_inflamed_18", "Checkpoints", "Melanocytic")
paired <- long |> group_by(patient) |> filter(n_distinct(time) == 2) |> ungroup()
ggplot(filter(paired, signature %in% key), aes(time, score, group = patient, colour = response)) +
  geom_line(alpha = 0.45) + geom_point(size = 0.9) +
  stat_summary(aes(group = response), fun = median, geom = "line", linewidth = 1.4) +
  facet_grid(response ~ signature) +
  scale_colour_manual(values = resp_cols, guide = "none") +
  labs(x = NULL, y = "Signature score (mean z)",
       title = "Per-patient change in immune signatures on nivolumab",
       subtitle = "thin = patient, thick = group median")
ggsave(file.path(fig, "09_signature_trajectories.png"), width = 12, height = 6.5, dpi = 300)

# ---- Per-patient table: Pre, On, delta ----
rna_patient <- long |>
  select(patient, response, benefit, cohort, tri, time, signature, score) |>
  pivot_wider(names_from = c(signature, time), values_from = score, names_sep = "_")
for (s in names(sigs))
  rna_patient[[paste0(s, "_delta")]] <- rna_patient[[paste0(s, "_On")]] - rna_patient[[paste0(s, "_Pre")]]
write_csv(rna_patient, file.path(out, "rna_scores_per_patient.csv"))

# ---- Fig 10: effect sizes, PRCR vs PD (Hodges-Lehmann shift + 95% CI) ----
# Run in all RNA patients and in the 26 tri-modal patients
hl <- function(d) {
  w <- wilcox.test(value ~ droplevels(response), data = d, conf.int = TRUE, exact = FALSE)
  tibble(n_PD = sum(d$response == "PD"), n_PRCR = sum(d$response == "PRCR"),
         est = -w$estimate, lo = -w$conf.int[2], hi = -w$conf.int[1], p = w$p.value)  # PRCR - PD
}
eff <- map_dfr(c(all = FALSE, tri = TRUE), \(only_tri) {
  rna_patient |>
    filter(response != "SD", !only_tri | tri) |>
    pivot_longer(matches("_(Pre|On|delta)$"), names_to = "feature", values_to = "value") |>
    filter(!is.na(value)) |>
    separate(feature, c("signature", "when"), sep = "_(?=[^_]+$)") |>
    group_by(signature, when) |>
    group_modify(~ hl(.x)) |> ungroup()
}, .id = "set") |>
  mutate(signature = factor(signature, rev(names(sigs))), when = factor(when, c("Pre", "On", "delta")),
         padj = p.adjust(p, "BH"), .by = set)
write_csv(eff, file.path(out, "signature_effect_sizes_PRCR_vs_PD.csv"))

ggplot(filter(eff, set == "all"), aes(est, signature, colour = p < 0.05)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.3) +
  facet_wrap(~ when, labeller = as_labeller(c(Pre = "Pre-treatment", On = "On-treatment",
                                              delta = "Change (On - Pre)"))) +
  scale_colour_manual(values = c(`FALSE` = "grey40", `TRUE` = "#C0392B"), name = "p < 0.05") +
  labs(x = "PRCR - PD shift in score (Hodges-Lehmann, 95% CI)", y = NULL,
       title = "Which RNA features separate responders from progressors?",
       subtitle = "All RNA-seq patients (PRCR vs PD; SD excluded)")
ggsave(file.path(fig, "10_effect_sizes_PRCR_vs_PD.png"), width = 9, height = 5, dpi = 300)

# ---- Fig 11: per-patient heatmap of on-treatment change (tri-modal patients) ----
dm <- rna_patient |>
  filter(tri, !is.na(CD8_T_delta)) |>
  arrange(response, desc(Tcell_inflamed_18_delta)) |>
  column_to_rownames("patient")
mat <- t(as.matrix(dm[, paste0(names(sigs), "_delta")]))
rownames(mat) <- names(sigs)
pheatmap(mat, cluster_cols = FALSE, cluster_rows = FALSE,
         annotation_col = dm["response"], gaps_col = cumsum(table(dm$response)),
         annotation_colors = list(response = resp_cols),
         breaks = seq(-1.5, 1.5, length.out = 101),
         color = colorRampPalette(c("#2471A3", "white", "#C0392B"))(100),
         main = "On - Pre change in signature score (tri-modal patients)",
         filename = file.path(fig, "11_heatmap_delta_tri_patients.png"), width = 9, height = 4.5)

# ---- Fig 12: tumour content confound -- melanocytic vs immune score ----
ggplot(long |> select(sample, time, response, signature, score) |>
         pivot_wider(names_from = signature, values_from = score),
       aes(Melanocytic, Tcell_inflamed_18, colour = response, shape = time)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, colour = "black", linewidth = 0.6) +
  scale_colour_manual(values = resp_cols) +
  labs(title = "Immune signal vs tumour (melanocytic) content",
       subtitle = sprintf("Spearman rho = %.2f",
                          cor(scores[, "Melanocytic"], scores[, "Tcell_inflamed_18"], method = "spearman")))
ggsave(file.path(fig, "12_melanocytic_vs_inflamed.png"), width = 6, height = 4.5, dpi = 300)

eff |> filter(set == "all", p < 0.05) |> arrange(p) |> print()
