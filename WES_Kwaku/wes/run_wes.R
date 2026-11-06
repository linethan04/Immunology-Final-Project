# WES reproduction of Riaz et al., Cell 2017 (nivolumab in melanoma).
# Runs the authors' exome scripts from riazn/bms038_analysis, then makes the
# slide figures and per-patient tables in wes/results/.
# Run from the repository root after wes/setup.R:  Rscript wes/run_wes.R
root <- normalizePath('.', winslash = '/')
repo <- paste0(root, '/external/bms038_analysis/')
out <- paste0(root, '/wes/results/')
dir.create(out, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(repo)) stop('Authors\' repository not found. Run: Rscript wes/setup.R')
library(survival); library(ggplot2)
source(paste0(root, '/external/nrstats/R/graphs_survival.r'))

# Upstream code relies on merge() row order (e.g. df$col[idx] <- ...), so sort
# strings byte-wise, as SQLite does, on every platform.
invisible(Sys.setlocale('LC_COLLATE', 'C'))
# Plain-R stand-in for the two sqldf queries the upstream scripts run
# (select [Key,] COUNT(*) as Name from tbl group by Key). sqldf's compiled
# dependencies are blocked by Windows Smart App Control on some machines.
sqldf <- function(q) {
  m <- regmatches(q, regexec('(?i)select\\s+(?:(\\w+)\\s*,\\s*)?count\\(\\*\\)\\s+as\\s+(\\w+)\\s+from\\s+(\\w+)\\s+group\\s+by\\s+(\\w+)', q, perl = TRUE))[[1]]
  if (length(m) == 0) stop('Unsupported sqldf query: ', q)
  d <- get(m[4], envir = parent.frame())
  tab <- table(d[[m[5]]])
  res <- data.frame(as.vector(tab))
  names(res) <- m[3]
  if (nzchar(m[2])) res <- cbind(setNames(data.frame(names(tab), stringsAsFactors = FALSE), m[2]), res)
  res
}
stext <- function(text, pos = 0, side = 3, ...) mtext(text, side = side, line = pos, cex = .7)
# Compatibility replacement for old nrstats boxplot layout, preserving tests.
tTestandGraph <- function(dat, groupCol, groups, valCol, groupLabels, tit = '', ...) {
  d <- dat[dat[[groupCol]] %in% groups, , drop = FALSE]
  d$Group <- factor(d[[groupCol]], levels = groups, labels = groupLabels)
  a <- d[[valCol]][d$Group == groupLabels[1]]; b <- d[[valCol]][d$Group == groupLabels[2]]
  cap <- sprintf('Welch t-test p=%.3g; Wilcoxon p=%.3g', t.test(a, b)$p.value, wilcox.test(a, b)$p.value)
  ggplot(d, aes(Group, .data[[valCol]])) + geom_boxplot() + geom_point() + theme_bw() + labs(title = tit, caption = cap)
}
# Run an upstream script unchanged on disk: drop workspace clearing and package
# loads we replace above, point PARENTDIR at the local copy, raise PNG resolution.
run_original <- function(file) {
  lines <- readLines(paste0(repo, file))
  lines <- lines[!grepl('^rm\\(list|^library\\(nrstats\\)|^library\\(sqldf\\)', lines)]
  lines[grepl('^PARENTDIR <-', lines)] <- paste0('PARENTDIR <- "', repo, '"')
  lines <- gsub('png(', 'png(width=1800,height=1400,res=180, ', lines, fixed = TRUE)
  eval(parse(text = lines), envir = .GlobalEnv)
}
setwd(repo)
run_original('exome/baseline_analysis.r')
baseline <- ptdat
write.csv(baseline, paste0(out, 'baseline_patient_features.csv'), row.names = FALSE)
run_original('exome/baseline_clonality.r')
write.csv(combdat, paste0(out, 'baseline_clonality_features.csv'), row.names = FALSE)
run_original('clonality_pre_post/clonality.analysis.R')
dir.create(paste0(out, 'original_figures'), showWarnings = FALSE)
file.copy(list.files(paste0(repo, 'output'), pattern = '^(Fig1B|Fig2|SupFig1).*\\.(pdf|png)$', full.names = TRUE), paste0(out, 'original_figures'), overwrite = TRUE)
save_plot <- function(p, name, w = 10, h = 6) {
  ggsave(paste0(out, name, '.png'), p, width = w, height = h, dpi = 300)
  ggsave(paste0(out, name, '.pdf'), p, width = w, height = h)
}
baseline$Response <- factor(baseline$myBOR, levels = c('PRCR', 'SD', 'PD'), labels = c('PR/CR', 'SD', 'PD'))
p <- ggplot(baseline, aes(Response, NonSynMut, color = Response)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(width = .12, seed = 38), size = 2) + facet_wrap(~Cohort) + scale_y_log10() + theme_bw(base_size = 16) + labs(title = 'Baseline mutation counts by clinical response', y = 'Nonsynonymous mutation count (log scale)', x = NULL, caption = 'Available baseline WES cohort; prior-ipilimumab cohorts shown separately.') + theme(legend.position = 'none')
save_plot(p, '01_baseline_mutations_by_response')
changed <- update.ccf(dd, 'pyclone', .95, TRUE)
counts <- as.data.frame.matrix(table(changed$sample, changed$clonal2))
counts$Patient <- rownames(counts); counts$Response <- class[counts$Patient]
counts$Total <- rowSums(counts[, c('negative selection', 'no selection', 'positive selection')])
counts$Contraction_fraction <- counts[['negative selection']] / counts$Total
write.csv(counts, paste0(out, 'paired_mutation_features.csv'), row.names = FALSE)
long <- as.data.frame(table(changed$sample, changed$clonal2)); names(long) <- c('Patient', 'Category', 'Count')
long$Response <- factor(class[as.character(long$Patient)], levels = c('PRCR', 'SD', 'PD'), labels = c('PR/CR', 'SD', 'PD'))
p <- ggplot(long, aes(Patient, Count, fill = Category)) + geom_col(position = 'fill') + facet_grid(~Response, scales = 'free_x', space = 'free_x') + scale_fill_manual(values = c('#277DA8', '#BBBBBB', '#E59845'), labels = c('Contraction', 'Persistence', 'Expansion')) + theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = 'bottom') + labs(title = 'Mutation changes during nivolumab treatment', y = 'Fraction of classified mutations', x = 'Patient', fill = NULL, caption = 'Authors CCF categories and 0.95 threshold; categories alone do not establish selection mechanisms.')
save_plot(p, '02_mutation_changes_per_patient', 12, 6)
counts$Response <- factor(counts$Response, levels = c('PRCR', 'SD', 'PD'), labels = c('PR/CR', 'SD', 'PD'))
p <- ggplot(counts, aes(Response, Contraction_fraction, color = Response)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(width = .12, seed = 38), size = 2) + theme_bw(base_size = 16) + theme(legend.position = 'none') + labs(title = 'Mutation contraction by clinical response', x = NULL, y = 'Fraction of mutations classified as contracting', caption = 'Each point is one patient with paired WES data.')
save_plot(p, '03_contraction_by_response')
neo <- read.csv(paste0(repo, 'neoantigen depletion/Mut_NeoAg_Expected_vs_Observed_BMS038.csv'))
stopifnot(all(neo$E_NeoAgCnt_On > 0))
neo$Observed_expected_ratio <- neo$NeoAgCnt_On / neo$E_NeoAgCnt_On
write.csv(neo, paste0(out, 'neoantigen_patient_features.csv'), row.names = FALSE)
setwd(paste0(repo, 'neoantigen depletion'))
pdf(paste0(out, 'original_figures/neoantigen_observed_expected.pdf'), width = 12, height = 7)
lines <- readLines('ObsvsExpt_NeoAg.R'); eval(parse(text = lines[!grepl('^rm\\(list', lines)])); dev.off()
neo$Response <- factor(neo$myBOR, levels = c('PRCR', 'SD', 'PD'), labels = c('PR/CR', 'SD', 'PD'))
p <- ggplot(neo, aes(reorder(Patient, Observed_expected_ratio), Observed_expected_ratio, color = Response)) + geom_hline(yintercept = 1, linetype = 2) + geom_point(size = 3) + coord_flip() + theme_bw(base_size = 14) + labs(title = 'Observed versus expected on-treatment neoantigens', x = 'Patient', y = 'Observed / expected neoantigen count', caption = 'Below 1: fewer candidates than expected under the authors model; not a baseline-loss fraction.')
save_plot(p, '04_neoantigen_observed_expected', 10, 8)
availability <- read.csv(paste0(repo, 'data/genomic_data_per_case.csv'))
shared <- availability[apply(availability[, -1], 1, function(x) all(!is.na(x) & x == 1)), ]
write.csv(shared, paste0(out, 'all_modality_paired_patient_ids.csv'), row.names = FALSE)
cat('\nCOUNTS: baseline=', nrow(baseline), ' paired WES=', nrow(counts), ' neoantigen=', nrow(neo), ' all-modality pairs=', nrow(shared), '\n')
print(table(baseline$Response)); print(table(counts$Response)); print(table(neo$Response))
writeLines(capture.output(sessionInfo()), paste0(out, 'R_session_info.txt'))
setwd(root)
