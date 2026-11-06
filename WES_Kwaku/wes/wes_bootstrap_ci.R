# Bootstrap 95% CIs for the slide medians (percentile method, 10,000 resamples).
# Run from the repository root after run_wes.R:  Rscript wes/wes_bootstrap_ci.R
out <- 'wes/results/'
set.seed(2017)
B <- 10000
boot_med <- function(x) quantile(replicate(B, median(sample(x, replace = TRUE))), c(.025, .975))
boot_diff <- function(a, b) quantile(replicate(B, median(sample(a, replace = TRUE)) - median(sample(b, replace = TRUE))), c(.025, .975))

summarise <- function(val, grp, label) {
  rows <- lapply(c('PRCR', 'SD', 'PD'), function(g) {
    x <- val[grp == g]; ci <- boot_med(x)
    data.frame(measure = label, group = g, n = length(x), median = median(x), lo = ci[1], hi = ci[2])
  })
  a <- val[grp == 'PRCR']; b <- val[grp == 'PD']; ci <- boot_diff(a, b)
  rbind(do.call(rbind, rows),
        data.frame(measure = label, group = 'PRCR minus PD', n = NA, median = median(a) - median(b), lo = ci[1], hi = ci[2]))
}

b <- read.csv(paste0(out, 'baseline_patient_features.csv'))
b <- b[!is.na(b$NonSynMut), ]
p <- read.csv(paste0(out, 'paired_mutation_features.csv'))
n <- read.csv(paste0(out, 'neoantigen_patient_features.csv'))

res <- rbind(summarise(b$NonSynMut, b$myBOR, 'baseline_nonsyn_mutations'),
             summarise(p$Contraction_fraction, p$Response, 'contraction_fraction'),
             summarise(n$Observed_expected_ratio, n$myBOR, 'neoantigen_obs_exp_ratio'))
rownames(res) <- NULL
write.csv(res, paste0(out, 'bootstrap_median_ci.csv'), row.names = FALSE)
print(res, digits = 3)

# Exploratory rank-sum tests quoted on the slides (responders vs progressors only)
cat('\nWilcoxon PR/CR vs PD, contraction fraction: p =', signif(wilcox.test(p$Contraction_fraction[p$Response == 'PRCR'], p$Contraction_fraction[p$Response == 'PD'])$p.value, 3), '\n')
cat('Wilcoxon PR/CR vs PD, neoantigen obs/exp ratio: p =', signif(wilcox.test(n$Observed_expected_ratio[n$myBOR == 'PRCR'], n$Observed_expected_ratio[n$myBOR == 'PD'])$p.value, 3), '\n')
