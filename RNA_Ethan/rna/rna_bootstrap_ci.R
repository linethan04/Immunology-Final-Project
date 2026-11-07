# Percentile bootstrap CIs (same method and seed as the WES analysis), Wilcoxon tests,
# and checks on the two main threats: tumour dilution and the choice of patients.
source("rna/common.R")
rna <- read_csv(file.path(out, "rna_patient_features.csv"), show_col_types = FALSE) |>
  mutate(Response = resp(response))

set.seed(2017)
B <- 10000
boot_med  <- function(x) quantile(replicate(B, median(sample(x, replace = TRUE))), c(.025, .975))
boot_diff <- function(a, b) quantile(replicate(B, median(sample(a, replace = TRUE)) - median(sample(b, replace = TRUE))), c(.025, .975))
summarise_measure <- function(col, label) {
  d <- filter(rna, !is.na(.data[[col]])); val <- d[[col]]; grp <- d$Response
  rows <- map_dfr(levels(grp), \(g) { x <- val[grp == g]; ci <- boot_med(x)
    tibble(measure = label, group = g, n = length(x), median = median(x), lo = ci[1], hi = ci[2]) })
  a <- val[grp == "PR/CR"]; b <- val[grp == "PD"]; ci <- boot_diff(a, b)
  bind_rows(rows, tibble(measure = label, group = "PR/CR minus PD", n = NA, median = median(a) - median(b), lo = ci[1], hi = ci[2]))
}
measures <- c(baseline_immune_score = "Tcell_inflamed_18_Pre", baseline_cd8_score = "CD8_T_Pre",
              share_genes_changed = "genes_changed", tumor_score_change = "Melanocytic_delta",
              cd8_score_change = "CD8_T_delta", immune_score_change = "Tcell_inflamed_18_delta",
              ifng_score_change = "IFNG_6_delta", cd8_score_on = "CD8_T_On")
res <- imap_dfr(measures, \(col, label) summarise_measure(col, label))
write_csv(res, file.path(out, "bootstrap_median_ci.csv"))
print(res, n = 40)

wil <- function(d, col) { d <- filter(d, !is.na(.data[[col]]))
  wilcox.test(d[[col]][d$response == "PRCR"], d[[col]][d$response == "PD"], exact = FALSE)$p.value }
tests <- imap_dfr(measures, \(col, label) tibble(measure = label,
  p_all = wil(rna, col), p_shared_patients = wil(filter(rna, shared), col), p_without_Pt3 = wil(filter(rna, patient != "Pt3"), col)))
write_csv(tests, file.path(out, "wilcoxon_PRCR_vs_PD.csv"))
print(tests)

# Dilution check: does the PR/CR effect survive adjusting for the tumour score?
d <- filter(rna, response %in% c("PD", "PRCR")) |> mutate(PRCR = as.integer(response == "PRCR"))
row <- function(f, label) { m <- lm(f, d); ci <- confint(m)["PRCR", ]
  tibble(model = label, PRCR_effect = coef(m)["PRCR"], lo = ci[1], hi = ci[2], p = summary(m)$coefficients["PRCR", 4], n = nobs(m)) }
dilution <- bind_rows(row(CD8_T_delta ~ PRCR, "CD8 change"),
                      row(CD8_T_delta ~ PRCR + Melanocytic_delta, "CD8 change, adjusted for tumor score change"),
                      row(CD8_T_On ~ PRCR, "CD8 on treatment"),
                      row(CD8_T_On ~ PRCR + Melanocytic_On, "CD8 on treatment, adjusted for tumor score on treatment"))
write_csv(dilution, file.path(out, "dilution_check.csv"))
print(dilution)
