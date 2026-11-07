# Step 5: link RNA to the group's other data (figure 05) and to outcome (supplementary)
source("rna/common.R")
rna <- read_csv(file.path(out, "rna_patient_features.csv"), show_col_types = FALSE) |> mutate(Response = resp(response))

tcr <- read.csv(file.path(data_dir, "tcr/TCRmetrics.csv"), check.names = FALSE) |>
  transmute(patient = sub("^pt", "Pt", Patient),
            tcr_fraction_Pre = fraction_productive_of_cells_mass_estimate_Pre,
            tcr_fraction_On  = fraction_productive_of_cells_mass_estimate_On,
            clonality_delta  = `CLcdr3_On-Pre`)
tmb <- read.csv(file.path(data_dir, "pre_therapy_nonsynonmous_mutations.csv")) |> count(patient = Patient, name = "nonsyn_mutations")
wes <- read.csv("../WES_Kwaku/wes/results/paired_mutation_features.csv") |> transmute(patient = Patient, contraction = Contraction_fraction)
neo <- read.csv("../WES_Kwaku/wes/results/neoantigen_patient_features.csv") |> transmute(patient = Patient, neoantigen_ratio = Observed_expected_ratio)
x <- reduce(list(rna, tcr, tmb, wes, neo), left_join, by = "patient")
write_csv(x, file.path(out, "rna_tcr_wes_merged.csv"))

spearman <- function(a, b, label, d = x) {
  ok <- complete.cases(d[[a]], d[[b]]); t <- cor.test(d[[a]][ok], d[[b]][ok], method = "spearman", exact = FALSE)
  tibble(comparison = label, rho = unname(t$estimate), p = t$p.value, n = sum(ok))
}
cors <- bind_rows(
  spearman("CD8_T_Pre", "tcr_fraction_Pre", "RNA CD8 vs TCR T-cell fraction, Pre"),
  spearman("CD8_T_On", "tcr_fraction_On", "RNA CD8 vs TCR T-cell fraction, On"),
  spearman("Cytotoxic_delta", "clonality_delta", "Change in cytotoxic score vs change in TCR clonality"),
  spearman("Melanocytic_delta", "contraction", "Change in tumor score vs WES contraction"),
  spearman("Melanocytic_delta", "contraction", "Change in tumor score vs WES contraction, excluding PR/CR", filter(x, response != "PRCR")),
  spearman("CD8_T_delta", "contraction", "Change in CD8 score vs WES contraction"),
  spearman("Melanocytic_delta", "neoantigen_ratio", "Change in tumor score vs neoantigen observed/expected"),
  spearman("Tcell_inflamed_18_Pre", "nonsyn_mutations", "Baseline immune score vs mutation burden"))
write_csv(cors, file.path(out, "cross_modality_correlations.csv"))
print(cors)

# Figure 05: does RNA agree with WES and TCR-seq?
sc <- function(xv, yv, xl, yl, title, row) ggplot(x, aes(.data[[xv]], .data[[yv]], colour = Response, shape = Response)) +
  geom_point(size = 3, alpha = 0.9, na.rm = TRUE) + scale_colour_manual(values = resp_cols, na.translate = FALSE) +
  labs(x = xl, y = yl, title = title, subtitle = sprintf("Spearman rho = %.2f, n = %d", cors$rho[row], cors$n[row]), colour = NULL, shape = NULL)
p <- (sc("contraction", "Melanocytic_delta", "WES: share of mutations contracting", "RNA: change in tumor score", "Tumor loss: RNA vs WES", 4) +
      sc("tcr_fraction_On", "CD8_T_On", "TCR-seq: T-cell fraction of cells (On)", "RNA: CD8 T-cell score (On)", "T cells: RNA vs TCR-seq", 2)) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
save_plot(p, "05_rna_vs_wes_tcr", 11, 5.5)

# Supplementary: AUC, PR/CR vs PD
feats <- c("CD8_T", "Cytotoxic", "IFNG_6", "Tcell_inflamed_18", "Checkpoints", "Melanocytic")
auc_tab <- expand_grid(signature = feats, when = c("Pre", "On", "delta")) |>
  mutate(res = map2(signature, when, \(s, w) {
    d  <- filter(rna, response != "SD") |> select(response, v = all_of(paste0(s, "_", w))) |> drop_na()
    ci <- ci.auc(roc(d$response, d$v, levels = c("PD", "PRCR"), direction = "auto", quiet = TRUE))
    tibble(n = nrow(d), auc = ci[2], lo = ci[1], hi = ci[3])
  })) |> unnest(res) |> mutate(when = factor(when, c("Pre", "On", "delta")), signature = factor(signature, rev(feats)))
write_csv(auc_tab, file.path(out, "auc_PRCR_vs_PD.csv"))
p <- ggplot(auc_tab, aes(auc, signature, colour = when)) + geom_vline(xintercept = 0.5, linetype = 2) +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(0.6), size = 0.3) +
  scale_colour_manual(values = c(Pre = "grey55", On = UP, delta = DOWN), labels = c("Pre-treatment", "On-treatment", "Change")) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "AUC, PR/CR vs PD (95% CI, DeLong)", y = NULL, colour = NULL, title = "How well does each RNA score classify response?",
       subtitle = "Single features; direction chosen per feature, so optimistic")
save_supp(p, "S12_auc", 8, 5)

# Supplementary: overall survival by on-treatment T-cell-inflamed score
clin <- read.csv(file.path(data_dir, "bms038_clinical_data.csv")) |>
  transmute(patient = PatientID, os_wk = OSWK, death = as.integer(OS_SOR == 0))  # OS_SOR 1 = censored
s <- inner_join(rna, clin, by = "patient") |> filter(!is.na(Tcell_inflamed_18_On)) |>
  mutate(group = if_else(Tcell_inflamed_18_On > median(Tcell_inflamed_18_On), "High", "Low"))
fit <- survfit(Surv(os_wk, death) ~ group, s)
cox <- coxph(Surv(os_wk, death) ~ scale(Tcell_inflamed_18_On), s)
hr  <- exp(c(coef(cox), confint(cox)))
lr  <- survdiff(Surv(os_wk, death) ~ group, s)
km  <- broom::tidy(fit) |> mutate(group = sub("group=", "", strata)) |> bind_rows(tibble(time = 0, estimate = 1, group = c("High", "Low")))
p <- ggplot(km, aes(time, estimate, colour = group)) + geom_step(linewidth = 0.9) +
  geom_point(data = filter(km, n.censor > 0), shape = 3, size = 1.5) +
  scale_colour_manual(values = c(High = UP, Low = DOWN)) +
  labs(x = "Weeks", y = "Overall survival", colour = "On-treatment\nT-cell-inflamed", title = "Overall survival by on-treatment immune score",
       subtitle = sprintf("n = %d; log-rank p = %.2g; HR per SD = %.2f (95%% CI %.2f-%.2f)", nrow(s), 1 - pchisq(lr$chisq, 1), hr[1], hr[2], hr[3]))
save_supp(p, "S13_survival", 7.5, 5)
