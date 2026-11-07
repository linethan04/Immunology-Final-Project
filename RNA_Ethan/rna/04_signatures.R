# Step 4: patient-level measures -> figures 01-03 and the table the group joins on.
#  Score = mean z-scored VST expression of marker genes (MCP-counter-style cell-type markers;
#  Ayers 2017 IFN-gamma and T-cell-inflamed sets). Melanocytic score = tumour-content proxy.
source("rna/common.R")
vsd <- readRDS(file.path(cache, "vsd.rds"))
samples <- read_csv(file.path(out, "samples.csv"), show_col_types = FALSE) |> mutate(time = factor(time, c("Pre", "On")))
sym <- rowData(vsd)$symbol

sigs <- list(
  T_cells = c("CD3D", "CD3E", "CD3G", "CD28", "CD6", "TRAT1"),
  CD8_T = c("CD8A", "CD8B"),
  Cytotoxic = c("GZMA", "GZMB", "GZMH", "PRF1", "NKG7", "GNLY"),
  NK = c("KLRD1", "KLRF1", "NCR1", "SH2D1B"),
  B_cells = c("CD19", "MS4A1", "CD79A", "CR2", "FCRL2"),
  Myeloid = c("CD14", "CD163", "CD68", "CSF1R", "FCGR3A"),
  Treg = c("FOXP3", "IL2RA", "CTLA4"),
  IFNG_6 = c("IFNG", "STAT1", "IDO1", "CXCL9", "CXCL10", "HLA-DRA"),
  Tcell_inflamed_18 = c("CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6", "HLA-DQA1",
                        "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT"),
  Checkpoints = c("PDCD1", "CD274", "LAG3", "HAVCR2", "TIGIT", "CTLA4"),
  Antigen_pres = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9"),
  Melanocytic = c("MITF", "PMEL", "TYR", "MLANA", "DCT"))

z <- t(scale(t(assay(vsd))))
scores <- sapply(sigs, \(g) colMeans(z[na.omit(match(g, sym)), , drop = FALSE]))
write_csv(as_tibble(scores, rownames = "sample"), file.path(out, "signature_scores_per_sample.csv"))

long <- as_tibble(scores, rownames = "sample") |> left_join(samples, by = "sample") |>
  pivot_longer(names(sigs), names_to = "signature", values_to = "score") |>
  mutate(signature = factor(signature, names(sigs)))

# Per-patient table: Pre, On, change
rna_patient <- long |> select(patient, response, benefit, cohort, shared, time, signature, score) |>
  pivot_wider(names_from = c(signature, time), values_from = score, names_sep = "_")
for (s in names(sigs)) rna_patient[[paste0(s, "_delta")]] <- rna_patient[[paste0(s, "_On")]] - rna_patient[[paste0(s, "_Pre")]]

# Share of genes changing at least 2-fold between a patient's two biopsies (1 VST unit is about 2-fold)
m <- assay(vsd)
paired <- intersect(rna_patient$patient[!is.na(rna_patient$CD8_T_delta)], samples$patient)
gene_change <- map_dfr(paired, \(p) {
  d <- m[, paste0(p, "_On")] - m[, paste0(p, "_Pre")]
  tibble(patient = p, genes_up = mean(d >= 1), genes_down = mean(d <= -1), genes_changed = genes_up + genes_down)
})
rna_patient <- left_join(rna_patient, gene_change, by = "patient")
write_csv(rna_patient, file.path(out, "rna_patient_features.csv"))
cat("Patients:", nrow(rna_patient), " Paired:", length(paired), " Paired and shared:", sum(rna_patient$shared & !is.na(rna_patient$CD8_T_delta)), "\n")

# Noise floor: same measure between repeat libraries of one biopsy label (the ".2" columns in the count table)
counts <- read.delim(file.path(data_dir, "CountData.BMS038.txt"), check.names = FALSE)[rownames(vsd), ]
reps <- intersect(grep("[.]2$", colnames(counts), value = TRUE), paste0(colnames(counts), ".2"))
rep_vst <- vst(as.matrix(counts[, setdiff(colnames(counts), "HUGO")]), blind = TRUE)
replicate_noise <- tibble(library = reps, genes_changed = map_dbl(reps, \(r) mean(abs(rep_vst[, r] - rep_vst[, sub("[.]2$", "", r)]) >= 1)))
write_csv(replicate_noise, file.path(out, "replicate_noise.csv"))
print(replicate_noise); cat("Median repeat-library share:", median(replicate_noise$genes_changed), "\n")

# Figure 01: baseline immune score (the RNA counterpart of baseline mutation burden)
d <- filter(rna_patient, !is.na(Tcell_inflamed_18_Pre)) |>
  mutate(Response = resp(response), Cohort = recode(cohort, `NIV3-NAIVE` = "No prior ipilimumab", `NIV3-PROG` = "Progressed on ipilimumab"))
p <- ggplot(d, aes(Response, Tcell_inflamed_18_Pre, colour = Response)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(width = .12, seed = 38), size = 2) +
  facet_wrap(~ Cohort) + scale_colour_manual(values = resp_cols) + theme_bw(base_size = 16) + theme(legend.position = "none") +
  labs(title = "Baseline immune activity by clinical response", y = "T-cell-inflamed score (z-score)", x = NULL,
       caption = "Each point is one patient's pre-treatment biopsy. Prior-ipilimumab cohorts shown separately.")
save_plot(p, "01_baseline_immune_by_response")

# Figure 02: gene changes per patient (the RNA counterpart of mutation changes per patient)
d <- filter(rna_patient, !is.na(genes_changed)) |> mutate(Response = resp(response)) |>
  pivot_longer(c(genes_up, genes_down), names_to = "Direction", values_to = "share") |>
  mutate(Direction = factor(Direction, c("genes_down", "genes_up"), c("Down on treatment", "Up on treatment")))
p <- ggplot(d, aes(reorder(patient, genes_changed), share, fill = Direction)) + geom_col() +
  facet_grid(~ Response, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c(DOWN, UP)) + scale_y_continuous(labels = scales::percent) +
  theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 10), legend.position = "bottom") +
  labs(title = "Gene expression changes during nivolumab treatment", y = "Share of genes changing at least 2-fold", x = "Patient", fill = NULL,
       caption = sprintf("%s genes measured; patients sorted by amount of change within each group.", format(nrow(m), big.mark = ",")))
save_plot(p, "02_gene_changes_per_patient", 12, 6)

# Figure 03: tumour and immune change per patient (the RNA counterpart of the neoantigen ranking)
d <- filter(rna_patient, !is.na(CD8_T_delta)) |> mutate(Response = resp(response), ord = Melanocytic_delta) |>
  select(patient, Response, ord, `Tumor score change (melanocytic genes)` = Melanocytic_delta, `CD8 T-cell score change` = CD8_T_delta) |>
  pivot_longer(c(-patient, -Response, -ord), names_to = "measure", values_to = "change") |>
  mutate(measure = factor(measure, c("Tumor score change (melanocytic genes)", "CD8 T-cell score change")))
p <- ggplot(d, aes(change, reorder(patient, ord), colour = Response, shape = Response)) +
  geom_vline(xintercept = 0, linetype = 2) + geom_point(size = 2.8) +
  facet_wrap(~ measure, scales = "free_x") + scale_colour_manual(values = resp_cols) +
  theme_bw(base_size = 14) + theme(axis.text.y = element_text(size = 10), legend.position = "bottom") +
  labs(title = "Tumor and immune change in each patient", x = "On minus Pre (z-score units)", y = "Patient", colour = NULL, shape = NULL,
       caption = "Patients sorted by tumor score change. Left of the dashed line = score fell on treatment.")
save_plot(p, "03_tumor_immune_change_per_patient", 10, 9)

# Supplementary: trajectories, effect sizes, per-patient heatmap, tumour-content check
key <- c("CD8_T", "Cytotoxic", "IFNG_6", "Tcell_inflamed_18", "Checkpoints", "Melanocytic")
tr <- long |> filter(signature %in% key) |> group_by(patient) |> filter(n_distinct(time) == 2) |> ungroup() |> mutate(Response = resp(response))
p <- ggplot(tr, aes(time, score, group = patient, colour = Response)) + geom_line(alpha = 0.45) + geom_point(size = 0.9) +
  stat_summary(aes(group = Response), fun = median, geom = "line", linewidth = 1.4) +
  facet_grid(Response ~ signature) + scale_colour_manual(values = resp_cols, guide = "none") +
  labs(x = NULL, y = "Signature score (mean z)", title = "Per-patient change in signatures", subtitle = "Thin = patient, thick = group median")
save_supp(p, "S08_signature_trajectories", 12, 6.5)

hl <- function(d) {  # Hodges-Lehmann shift, PR/CR minus PD
  w <- wilcox.test(value ~ droplevels(factor(response, c("PD", "PRCR"))), data = d, conf.int = TRUE, exact = FALSE)
  tibble(n_PD = sum(d$response == "PD"), n_PRCR = sum(d$response == "PRCR"), est = -w$estimate, lo = -w$conf.int[2], hi = -w$conf.int[1], p = w$p.value)
}
eff <- map_dfr(c(all = FALSE, shared = TRUE), \(only) {
  rna_patient |> filter(response != "SD", !only | shared) |>
    pivot_longer(matches("_(Pre|On|delta)$"), names_to = "feature", values_to = "value") |> filter(!is.na(value)) |>
    separate(feature, c("signature", "when"), sep = "_(?=[^_]+$)") |>
    group_by(signature, when) |> group_modify(~ hl(.x)) |> ungroup()
}, .id = "set") |> mutate(padj = p.adjust(p, "BH"), .by = set) |>
  mutate(signature = factor(signature, rev(names(sigs))), when = factor(when, c("Pre", "On", "delta")))
write_csv(eff, file.path(out, "signature_effect_sizes_PRCR_vs_PD.csv"))
p <- ggplot(filter(eff, set == "all"), aes(est, signature, colour = p < 0.05)) + geom_vline(xintercept = 0, linetype = 2) +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.3) +
  facet_wrap(~ when, labeller = as_labeller(c(Pre = "Pre-treatment", On = "On-treatment", delta = "Change"))) +
  scale_colour_manual(values = c(`FALSE` = "grey40", `TRUE` = UP), name = "p < 0.05") +
  labs(x = "PR/CR minus PD shift in score (Hodges-Lehmann, 95% CI)", y = NULL, title = "Which RNA scores separate responders from progressors?")
save_supp(p, "S09_effect_sizes", 10, 5.5)

dm <- rna_patient |> filter(shared, !is.na(CD8_T_delta)) |> mutate(Response = resp(response)) |>
  arrange(Response, desc(Tcell_inflamed_18_delta)) |> column_to_rownames("patient")
mat <- t(as.matrix(dm[, paste0(names(sigs), "_delta")])); rownames(mat) <- names(sigs)
pheatmap(mat, cluster_cols = FALSE, cluster_rows = FALSE, annotation_col = dm["Response"],
         gaps_col = cumsum(table(dm$Response)), annotation_colors = list(Response = resp_cols),
         breaks = seq(-1.5, 1.5, length.out = 101), color = colorRampPalette(c(DOWN, "white", UP))(100),
         main = "On minus Pre signature change (patients shared across WES, RNA, TCR)",
         filename = file.path(supp, "S10_heatmap_change_shared_patients.png"), width = 9, height = 4.5)

w <- long |> select(sample, time, response, signature, score) |> pivot_wider(names_from = signature, values_from = score)
p <- ggplot(w, aes(Melanocytic, Tcell_inflamed_18, colour = resp(response), shape = time)) + geom_point(size = 2, alpha = 0.8) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, colour = "black", linewidth = 0.6) + scale_colour_manual(values = resp_cols) +
  labs(colour = "Response", shape = NULL, title = "Immune signal vs tumor content",
       subtitle = sprintf("Spearman rho = %.2f", cor(w$Melanocytic, w$Tcell_inflamed_18, method = "spearman")))
save_supp(p, "S11_tumor_content_vs_immune", 7, 5)
