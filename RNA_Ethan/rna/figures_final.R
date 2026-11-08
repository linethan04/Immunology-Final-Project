# Slide figures 01-05, redrawn from the saved result tables in one consistent, publication-style theme.
# Each slide carries the message, so figures keep short panel titles and no long captions.
source("rna/common.R"); source("rna/theme_pub.R")
rna <- read_csv(file.path(out, "rna_patient_features.csv"), show_col_types = FALSE) |> mutate(Response = resp(response))
boot <- read_csv(file.path(out, "bootstrap_median_ci.csv"), show_col_types = FALSE)
wil <- read_csv(file.path(out, "wilcoxon_PRCR_vs_PD.csv"), show_col_types = FALSE)
pv <- \(m) sprintf("PR/CR vs PD Wilcoxon p = %.2g", wil$p_all[wil$measure == m])

# 01 Baseline immune activity (paper: baseline expression only marginally predictive)
d <- filter(rna, !is.na(Tcell_inflamed_18_Pre)) |>
  mutate(Cohort = recode(cohort, `NIV3-NAIVE` = "Ipilimumab-naive", `NIV3-PROG` = "Progressed on ipilimumab"))
p <- ggplot(d, aes(Response, Tcell_inflamed_18_Pre, colour = Response)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey70") + resp_box() + facet_wrap(~ Cohort) +
  labs(x = NULL, y = "T-cell-inflamed score, pre-treatment\n(18-gene, z-score)", subtitle = pv("baseline_immune_score"))
save_plot(p, "01_baseline_immune_by_response", 9, 5.5)

# 02 Share of genes changing 2-fold per patient
d <- filter(rna, !is.na(genes_changed)) |>
  pivot_longer(c(genes_up, genes_down), names_to = "Direction", values_to = "share") |>
  mutate(Direction = factor(Direction, c("genes_up", "genes_down"), c("Up on treatment", "Down on treatment")))
p <- ggplot(d, aes(reorder(sub("Pt", "", patient), genes_changed), share, fill = Direction)) + geom_col(width = .8) +
  facet_grid(~ Response, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c(UP, DOWN)) + scale_y_continuous(labels = scales::percent, expand = expansion(c(0, .05))) +
  theme(axis.text.x = element_text(size = 10), legend.position = "top") +
  labs(x = "Patient", y = "Genes changing ≥ 2-fold", subtitle = pv("share_genes_changed"))
save_plot(p, "02_gene_changes_per_patient", 11, 5.5)

# 03 Tumor vs CD8 change per patient (dumbbell-free: two aligned dot columns)
d <- filter(rna, !is.na(CD8_T_delta)) |>
  select(patient, Response, ord = Melanocytic_delta, `Tumor (melanocytic) score` = Melanocytic_delta, `CD8 T-cell score` = CD8_T_delta) |>
  pivot_longer(c(-patient, -Response, -ord), names_to = "measure", values_to = "change") |>
  mutate(measure = factor(measure, c("Tumor (melanocytic) score", "CD8 T-cell score")))
p <- ggplot(d, aes(change, reorder(sub("Pt", "", patient), -ord), colour = Response)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") + geom_point(size = 3) +
  facet_wrap(~ measure, scales = "free_x") + scale_colour_manual(values = resp_cols) +
  theme(axis.text.y = element_text(size = 9), legend.position = "top", panel.grid.major.y = element_line(colour = "grey94")) +
  labs(x = "Change on treatment (On − Pre, z-score)", y = "Patient (sorted by tumor score change)")
save_plot(p, "03_tumor_immune_change_per_patient", 9, 8.5)

# 04 Hallmark pathways (paper: immune up, cell cycle down in responders)
programs <- c(HALLMARK_INTERFERON_GAMMA_RESPONSE = "IFN-γ response", HALLMARK_INTERFERON_ALPHA_RESPONSE = "IFN-α response",
              HALLMARK_ALLOGRAFT_REJECTION = "Allograft rejection", HALLMARK_INFLAMMATORY_RESPONSE = "Inflammatory response",
              HALLMARK_E2F_TARGETS = "E2F targets", HALLMARK_G2M_CHECKPOINT = "G2M checkpoint",
              HALLMARK_MYC_TARGETS_V1 = "MYC targets", HALLMARK_MITOTIC_SPINDLE = "Mitotic spindle")
d <- read_csv(file.path(out, "gsea_hallmark.csv"), show_col_types = FALSE) |>
  filter(pathway %in% names(programs), contrast != "All patients") |>
  mutate(label = factor(programs[pathway], rev(programs)),
         program = factor(if_else(pathway %in% names(programs)[1:4], "Immune", "Proliferation"), c("Immune", "Proliferation")),
         column = factor(contrast, c("PR/CR vs PD", "PR/CR", "SD", "PD")),
         panel = factor(if_else(contrast == "PR/CR vs PD", "Pre-treatment\nPR/CR vs PD", "On vs Pre\nwithin each group"),
                        c("Pre-treatment\nPR/CR vs PD", "On vs Pre\nwithin each group")),
         sig = if_else(padj < 0.05, NES, NA_real_))
p <- ggplot(d, aes(column, label, fill = sig)) +
  geom_tile(colour = "white", linewidth = 1.2) + geom_text(aes(label = sprintf("%.1f", NES)), size = 4.8) +
  facet_grid(program ~ panel, scales = "free", space = "free") +
  scale_fill_gradient2(low = DOWN, mid = "white", high = UP, limits = c(-3.8, 3.8), na.value = "grey94", name = "NES") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.title = element_text(), strip.text.y = element_text(angle = 0)) +
  labs(x = NULL, y = NULL, caption = "Normalized enrichment score (fgsea, Hallmark). Grey = FDR ≥ 0.05.")
save_plot(p, "04_pathways_by_response", 9, 5.5)

# 05 On-treatment CD8 is the immune signal that survives the tumor-dilution check
dil <- read_csv(file.path(out, "dilution_check.csv"), show_col_types = FALSE) |>
  mutate(when = if_else(grepl("^CD8 change", model), "CD8 change", "CD8 on treatment"),
         adj = factor(if_else(grepl("adjusted", model), "Adjusted for tumor score", "Unadjusted"), c("Unadjusted", "Adjusted for tumor score")))
pA <- ggplot(filter(rna, !is.na(CD8_T_On)), aes(Response, CD8_T_On, colour = Response)) + resp_box() +
  labs(x = NULL, y = "CD8 T-cell score, on treatment", title = "A  On-treatment CD8", subtitle = pv("cd8_score_on"))
pB <- ggplot(dil, aes(PRCR_effect, when, colour = adj)) + geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(.5), size = .7, linewidth = 1) +
  scale_colour_manual(values = c("grey55", "#F8766D")) + theme(legend.position = "bottom") +
  labs(x = "PR/CR minus PD (95% CI)", y = NULL, title = "B  Dilution check", subtitle = "Linear model, PR/CR vs PD")
p <- pA + pB + plot_layout(widths = c(.8, 1))
save_plot(p, "05_cd8_dilution_check", 11, 5.2)
