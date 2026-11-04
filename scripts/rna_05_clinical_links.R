# RNA step 5: link RNA scores to outcome and to the other modalities
#  - Prediction: AUC of Pre / On / delta scores for PRCR vs PD (aim 3 preview)
#  - Survival: overall survival by on-treatment T-cell-inflamed score
#  - Cross-modality: RNA immune scores vs TCR-seq and mutation burden (aim 2 preview)

suppressPackageStartupMessages({
  library(tidyverse)
  library(pROC)
  library(survival)
  library(patchwork)
})
out <- "output/rna"; fig <- "figures/rna"
theme_set(theme_bw(base_size = 11))
resp_cols <- c(PD = "#C0392B", SD = "#E1A23B", PRCR = "#2471A3")
rna <- read_csv(file.path(out, "rna_scores_per_patient.csv"), show_col_types = FALSE) |>
  mutate(response = factor(response, c("PD", "SD", "PRCR")))

# ---- Fig 13: AUC with DeLong 95% CI ----
feats <- c("CD8_T", "Cytotoxic", "IFNG_6", "Tcell_inflamed_18", "Checkpoints", "Melanocytic")
auc_tab <- expand_grid(signature = feats, when = c("Pre", "On", "delta")) |>
  mutate(res = map2(signature, when, \(s, w) {
    d  <- filter(rna, response != "SD") |> select(response, x = all_of(paste0(s, "_", w))) |> drop_na()
    rc <- roc(d$response, d$x, levels = c("PD", "PRCR"), direction = "auto", quiet = TRUE)
    ci <- ci.auc(rc)
    tibble(n = nrow(d), auc = ci[2], lo = ci[1], hi = ci[3])
  })) |> unnest(res) |>
  mutate(when = factor(when, c("Pre", "On", "delta")), signature = factor(signature, rev(feats)))
write_csv(auc_tab, file.path(out, "auc_PRCR_vs_PD.csv"))

ggplot(auc_tab, aes(auc, signature, colour = when)) +
  geom_vline(xintercept = 0.5, linetype = 2) +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(0.6), size = 0.3) +
  scale_colour_manual(values = c(Pre = "grey55", On = "#8E44AD", delta = "#16A085"),
                      labels = c("Pre-treatment", "On-treatment", "Change")) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "AUC, PRCR vs PD (95% CI, DeLong)", y = NULL, colour = NULL,
       title = "How well does each RNA score classify response?",
       subtitle = "Single-feature AUCs; direction chosen per feature (optimistic)")
ggsave(file.path(fig, "13_auc_by_timepoint.png"), width = 7, height = 4.5, dpi = 300)

# ---- Fig 14: overall survival by on-treatment inflamed score ----
clin <- read.csv("data/bms038_clinical_data.csv") |>
  transmute(patient = PatientID, os_wk = OSWK, death = as.integer(OS_SOR == 0))  # OS_SOR 1 = censored
surv_d <- inner_join(rna, clin, by = "patient") |>
  filter(!is.na(Tcell_inflamed_18_On)) |>
  mutate(group = if_else(Tcell_inflamed_18_On > median(Tcell_inflamed_18_On), "High", "Low"))
fit <- survfit(Surv(os_wk, death) ~ group, surv_d)
cox <- coxph(Surv(os_wk, death) ~ scale(Tcell_inflamed_18_On), surv_d)
hr  <- exp(c(coef(cox), confint(cox)))
lr  <- survdiff(Surv(os_wk, death) ~ group, surv_d)
km  <- broom::tidy(fit) |> mutate(group = sub("group=", "", strata)) |>
  bind_rows(tibble(time = 0, estimate = 1, group = c("High", "Low")))
ggplot(km, aes(time, estimate, colour = group)) +
  geom_step(linewidth = 0.9) +
  geom_point(data = filter(km, n.censor > 0), shape = 3, size = 1.5) +
  scale_colour_manual(values = c(High = "#2471A3", Low = "#C0392B")) +
  labs(x = "Weeks", y = "Overall survival", colour = "On-treatment\nT-cell-inflamed",
       title = "Overall survival by on-treatment T-cell-inflamed score (median split)",
       subtitle = sprintf("n = %d; log-rank p = %.3g; HR per SD = %.2f (95%% CI %.2f-%.2f)",
                          nrow(surv_d), 1 - pchisq(lr$chisq, 1), hr[1], hr[2], hr[3]))
ggsave(file.path(fig, "14_km_os_inflamed_on.png"), width = 7, height = 4.5, dpi = 300)

# ---- Fig 15: RNA vs TCR-seq and mutation burden (cross-modality sanity check) ----
tcr <- read.csv("data/tcr/TCRmetrics.csv", check.names = FALSE) |>
  transmute(patient = sub("^pt", "Pt", Patient),
            tcr_frac_Pre = fraction_productive_of_cells_mass_estimate_Pre,
            tcr_frac_On  = fraction_productive_of_cells_mass_estimate_On,
            clonality_delta = `CLcdr3_On-Pre`)
tmb <- read.csv("data/pre_therapy_nonsynonmous_mutations.csv") |> count(patient = Patient, name = "nonsyn_mut")
x <- rna |> left_join(tcr, by = "patient") |> left_join(tmb, by = "patient")

rho <- \(a, b) { t <- cor.test(a, b, method = "spearman", exact = FALSE)
  sprintf("rho = %.2f, p = %.2g, n = %d", t$estimate, t$p.value, sum(complete.cases(a, b))) }
sc <- \(xv, yv, xl, yl) ggplot(x, aes(.data[[xv]], .data[[yv]], colour = response)) +
  geom_point(size = 2, alpha = 0.85) +
  scale_colour_manual(values = resp_cols, na.translate = FALSE) +
  labs(x = xl, y = yl, subtitle = rho(x[[xv]], x[[yv]]))
p1 <- sc("tcr_frac_Pre", "CD8_T_Pre", "TCR-seq T-cell fraction (Pre)", "RNA CD8 score (Pre)")
p2 <- sc("tcr_frac_On",  "CD8_T_On",  "TCR-seq T-cell fraction (On)",  "RNA CD8 score (On)")
p3 <- sc("clonality_delta", "Cytotoxic_delta", "Change in TCR clonality", "Change in cytotoxic score")
p4 <- sc("nonsyn_mut", "Tcell_inflamed_18_Pre", "Non-synonymous mutations (Pre)", "T-cell-inflamed (Pre)") +
  scale_x_log10()
(p1 + p2) / (p3 + p4) + plot_layout(guides = "collect") +
  plot_annotation(title = "RNA immune scores vs TCR-seq and WES")
ggsave(file.path(fig, "15_rna_vs_tcr_wes.png"), width = 9, height = 7.5, dpi = 300)

write_csv(x, file.path(out, "rna_tcr_wes_merged_preview.csv"))
print(arrange(auc_tab, desc(auc)), n = 18)
print(summary(cox)$coefficients)
