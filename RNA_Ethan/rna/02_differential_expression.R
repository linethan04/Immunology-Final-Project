# Step 2: DESeq2 differential expression
#  A) Pre-treatment PR/CR vs PD: is there a baseline signature of response?
#  B) On vs Pre, paired by patient: what does nivolumab change overall?
#  C) On vs Pre within each response group, plus PR/CR change minus PD change
#  D) Benefit (PR/CR + SD) change minus no-benefit change: the paper's grouping
source("rna/common.R")
dds <- readRDS(file.path(cache, "dds.rds"))
vsd <- readRDS(file.path(cache, "vsd.rds"))
sym <- rowData(dds)$symbol
tidy_res <- function(r) as.data.frame(r) |> rownames_to_column("entrez") |>
  mutate(symbol = sym[match(entrez, rownames(dds))]) |> arrange(padj)

# A. Same four sample exclusions as the authors' Table S6A script
pre <- dds[, dds$time == "Pre" & !colnames(dds) %in% c("Pt48_Pre", "Pt67_Pre", "Pt52_Pre", "Pt27_Pre")]
design(pre) <- ~ response
pre <- cached("fit_A", DESeq(pre, quiet = TRUE))
shr <- lfcShrink(pre, coef = "response_PRCR_vs_PD", type = "apeglm", quiet = TRUE)
shr$stat <- results(pre, name = "response_PRCR_vs_PD")$stat
res_A <- tidy_res(shr)

# B. Every patient is their own control
paired <- names(which(table(dds$patient) == 2))
onpre  <- dds[, dds$patient %in% paired]
onpre$patient <- factor(onpre$patient)
design(onpre) <- ~ patient + time
fit_B <- cached("fit_B", DESeq(onpre, quiet = TRUE))
res_B <- tidy_res(results(fit_B, name = "time_On_vs_Pre"))

# C/D. Patients nested in group (DESeq2 vignette design), so On vs Pre is estimated per group
nested_fit <- function(group) {
  x <- onpre; x$grp <- x[[group]]
  x$ind <- as.data.frame(colData(x)) |> mutate(ind = factor(as.integer(droplevels(patient))), .by = grp) |> pull(ind)
  mm <- model.matrix(~ grp + grp:ind + grp:time, colData(x))
  DESeq(x, full = mm[, colSums(mm != 0) > 0], betaPrior = FALSE, quiet = TRUE)
}
fit_C <- cached("fit_C", nested_fit("response"))
res_C <- list(PRCR = results(fit_C, name = "grpPRCR.timeOn"),
              SD   = results(fit_C, name = "grpSD.timeOn"),
              PD   = results(fit_C, name = "grpPD.timeOn"),
              PRCR_vs_PD = results(fit_C, contrast = list("grpPRCR.timeOn", "grpPD.timeOn"))) |> map(tidy_res)
fit_D <- cached("fit_D", nested_fit("benefit"))
res_D <- tidy_res(results(fit_D, contrast = list("grpLB.timeOn", "grpNB.timeOn")))

all_res <- c(list(A_pre_PRCR_vs_PD = res_A, B_on_vs_pre_paired = res_B),
             set_names(res_C, paste0("C_on_vs_pre_", names(res_C))),
             list(D_on_vs_pre_benefit_vs_none = res_D))
iwalk(all_res, \(r, n) write_csv(r, file.path(out, paste0("de_", n, ".csv"))))

de_summary <- imap_dfr(all_res, \(r, n) tibble(contrast = n, genes_tested = sum(!is.na(r$padj)),
                                               fdr_0.1 = sum(r$padj < 0.1, na.rm = TRUE),
                                               fdr_0.2 = sum(r$padj < 0.2, na.rm = TRUE)))
write_csv(de_summary, file.path(out, "de_summary.csv"))
print(de_summary)
print(table(colData(onpre)$response) / 2)

# Supplementary figures
immune_genes <- c("CD8A", "GZMB", "GZMH", "PRF1", "IFNG", "CXCL9", "CXCL10", "CXCL13",
                  "PDCD1", "LAG3", "TIGIT", "CD274", "GBP5", "STAT1", "IDO1", "NKG7")
volcano <- function(res, title) {
  res <- filter(res, !is.na(padj))
  lab <- bind_rows(head(res, 15), filter(res, symbol %in% immune_genes, padj < 0.1)) |> distinct(entrez, .keep_all = TRUE)
  ggplot(res, aes(log2FoldChange, -log10(pvalue))) +
    geom_point(aes(colour = case_when(padj >= 0.1 ~ "ns", log2FoldChange > 0 ~ "Up", TRUE ~ "Down")), size = 0.7, alpha = 0.6) +
    geom_text_repel(data = lab, aes(label = symbol), size = 3, max.overlaps = 30, segment.size = 0.2) +
    scale_colour_manual(values = c(ns = GREY, Up = UP, Down = DOWN), name = "FDR < 0.1") +
    labs(title = title, subtitle = paste(sum(res$padj < 0.1), "genes at FDR < 0.1"), x = "log2 fold change", y = "-log10 p")
}
save_supp(volcano(res_A, "Pre-treatment: PR/CR vs PD"), "S03_volcano_pre_PRCR_vs_PD", 7.5, 5.5)
save_supp(volcano(res_B, "On vs Pre, paired, all patients"), "S04_volcano_on_vs_pre", 7.5, 5.5)

lfc <- inner_join(select(res_C$PRCR, entrez, symbol, PRCR = log2FoldChange),
                  select(res_C$PD, entrez, PD = log2FoldChange), by = "entrez") |>
  left_join(select(res_C$PRCR_vs_PD, entrez, padj_int = padj), by = "entrez") |>
  mutate(sig = !is.na(padj_int) & padj_int < 0.1)
p <- ggplot(lfc, aes(PD, PRCR)) +
  geom_hline(yintercept = 0, colour = "grey60") + geom_vline(xintercept = 0, colour = "grey60") +
  geom_abline(linetype = 2, colour = "grey50") +
  geom_point(aes(colour = sig), size = 0.6, alpha = 0.6) +
  geom_text_repel(data = filter(lfc, symbol %in% immune_genes), aes(label = symbol), size = 3, max.overlaps = 30) +
  scale_colour_manual(values = c(`FALSE` = GREY, `TRUE` = DOWN), labels = c("No", "Yes"), name = "Differs by group\n(FDR < 0.1)") +
  coord_cartesian(xlim = c(-4, 4), ylim = c(-4, 6)) +
  labs(title = "On-treatment change: responders vs progressors", x = "log2 On/Pre in PD", y = "log2 On/Pre in PR/CR")
save_supp(p, "S05_lfc_PRCR_vs_PD", 7.5, 5.5)

top <- head(filter(res_B, !is.na(symbol)), 40)
v   <- assay(vsd)[top$entrez, colnames(onpre)]; rownames(v) <- top$symbol
ann <- as.data.frame(colData(onpre))[, "time", drop = FALSE]; ann$Response <- resp(onpre$response)
ord <- order(ann$time, ann$Response)
pheatmap(t(scale(t(v[, ord]))), cluster_cols = FALSE, show_colnames = FALSE, annotation_col = ann[ord, ],
         breaks = seq(-2.5, 2.5, length.out = 101), color = colorRampPalette(c(DOWN, "white", UP))(100),
         annotation_colors = list(Response = resp_cols, time = c(Pre = "grey80", On = "grey30")),
         main = "Top 40 genes changed on nivolumab (row z-score)", fontsize_row = 7,
         filename = file.path(supp, "S06_heatmap_top_on_vs_pre.png"), width = 8, height = 7)
