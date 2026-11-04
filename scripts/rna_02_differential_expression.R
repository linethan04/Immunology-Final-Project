# RNA step 2: differential expression (DESeq2)
#  A) Pre-treatment PRCR vs PD           -> is there a baseline signature of response?
#  B) On vs Pre, paired by patient       -> what does nivolumab change overall?
#  C) On vs Pre within each response group + interaction (PRCR change - PD change)
#                                         -> is the on-treatment change response-specific?

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(ggrepel)
  library(pheatmap)
})
out <- "output/rna"; fig <- "figures/rna"
theme_set(theme_bw(base_size = 11))
dds <- readRDS(file.path(out, "dds.rds"))
vsd <- readRDS(file.path(out, "vsd.rds"))
sym <- rowData(dds)$symbol

tidy_res <- function(r) as.data.frame(r) |>
  rownames_to_column("entrez") |>
  mutate(symbol = sym[match(entrez, rownames(dds))]) |>
  arrange(padj)

# ---- A. Pre: PRCR vs PD (SD its own level; paper exclusions) ----
pre <- dds[, dds$time == "Pre" &
             !colnames(dds) %in% c("Pt48_Pre", "Pt67_Pre", "Pt52_Pre", "Pt27_Pre")]
design(pre) <- ~ response
pre <- DESeq(pre, quiet = TRUE)
res_pre <- tidy_res(lfcShrink(pre, coef = "response_PRCR_vs_PD", type = "apeglm", quiet = TRUE) |>
                      (\(s) { s$stat <- results(pre, name = "response_PRCR_vs_PD")$stat; s })())
write_csv(res_pre, file.path(out, "de_A_pre_PRCR_vs_PD.csv"))

# ---- B. On vs Pre, paired ----
paired <- names(which(table(dds$patient) == 2))
onpre  <- dds[, dds$patient %in% paired]
onpre$patient <- factor(onpre$patient)
design(onpre) <- ~ patient + time
onpre <- DESeq(onpre, quiet = TRUE)
res_onpre <- tidy_res(results(onpre, name = "time_On_vs_Pre"))
write_csv(res_onpre, file.path(out, "de_B_on_vs_pre_paired.csv"))

# ---- C. Group-specific On vs Pre (patients nested in response; DESeq2 vignette design) ----
cd <- as.data.frame(colData(onpre)) |>
  group_by(response) |>
  mutate(ind = factor(as.integer(droplevels(patient)))) |>
  ungroup()
onpre$ind <- cd$ind
mm <- model.matrix(~ response + response:ind + response:time, colData(onpre))
mm <- mm[, colSums(mm != 0) > 0]          # drop empty patient-index columns
grp <- DESeq(onpre, full = mm, betaPrior = FALSE, quiet = TRUE)
res_C <- list(
  PRCR        = results(grp, name = "responsePRCR.timeOn"),
  SD          = results(grp, name = "responseSD.timeOn"),
  PD          = results(grp, name = "responsePD.timeOn"),
  PRCR_vs_PD  = results(grp, contrast = list("responsePRCR.timeOn", "responsePD.timeOn"))
) |> map(tidy_res)
iwalk(res_C, \(r, n) write_csv(r, file.path(out, paste0("de_C_on_vs_pre_", n, ".csv"))))

n_sig <- function(r) sum(r$padj < 0.1, na.rm = TRUE)
de_summary <- tibble(
  contrast = c("A: Pre PRCR vs PD", "B: On vs Pre (all paired)",
               paste("C: On vs Pre in", c("PRCR", "SD", "PD")), "C: interaction PRCR - PD"),
  n_padj_0.1 = c(n_sig(res_pre), n_sig(res_onpre), map_int(res_C, n_sig)[c(1:3, 4)]))
write_csv(de_summary, file.path(out, "de_summary.csv"))
print(de_summary)
cat("Paired patients per group:\n"); print(table(cd$response) / 2)

# ---- Figures ----
immune_genes <- c("CD8A", "GZMB", "GZMH", "PRF1", "IFNG", "CXCL9", "CXCL10", "CXCL13",
                  "PDCD1", "LAG3", "TIGIT", "CD274", "GBP5", "STAT1", "IDO1", "NKG7")

volcano <- function(res, title, n_lab = 15) {
  res <- filter(res, !is.na(padj))
  lab <- bind_rows(head(res, n_lab), filter(res, symbol %in% immune_genes, padj < 0.1)) |>
    distinct(entrez, .keep_all = TRUE)
  ggplot(res, aes(log2FoldChange, -log10(pvalue))) +
    geom_point(aes(colour = case_when(padj >= 0.1 ~ "ns",
                                      log2FoldChange > 0 ~ "up", TRUE ~ "down")),
               size = 0.7, alpha = 0.6) +
    geom_text_repel(data = lab, aes(label = symbol), size = 2.8, max.overlaps = 30,
                    segment.size = 0.2) +
    scale_colour_manual(values = c(ns = "grey75", up = "#C0392B", down = "#2471A3"),
                        name = "padj < 0.1") +
    labs(title = title, subtitle = paste(n_sig(res), "genes at padj < 0.1"),
         x = "log2 fold change", y = "-log10 p")
}
ggsave(file.path(fig, "03_volcano_pre_PRCR_vs_PD.png"),
       volcano(res_pre, "A. Pre-treatment: PRCR vs PD (shrunken LFC)"),
       width = 6.5, height = 5, dpi = 300)
ggsave(file.path(fig, "04_volcano_on_vs_pre.png"),
       volcano(res_onpre, "B. On- vs pre-treatment (paired, all patients)"),
       width = 6.5, height = 5, dpi = 300)

# Response-specific change: per-gene LFC in PRCR vs in PD
lfc <- inner_join(select(res_C$PRCR, entrez, symbol, PRCR = log2FoldChange),
                  select(res_C$PD,   entrez, PD   = log2FoldChange), by = "entrez") |>
  left_join(select(res_C$PRCR_vs_PD, entrez, padj_int = padj), by = "entrez") |>
  mutate(sig = !is.na(padj_int) & padj_int < 0.1)
ggplot(lfc, aes(PD, PRCR)) +
  geom_hline(yintercept = 0, colour = "grey60") + geom_vline(xintercept = 0, colour = "grey60") +
  geom_abline(linetype = 2, colour = "grey50") +
  geom_point(aes(colour = sig), size = 0.6, alpha = 0.6) +
  geom_text_repel(data = filter(lfc, symbol %in% immune_genes), aes(label = symbol),
                  size = 2.8, max.overlaps = 30) +
  scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#8E44AD"),
                      name = "interaction\npadj < 0.1") +
  coord_cartesian(xlim = c(-4, 4), ylim = c(-4, 6)) +
  labs(title = "C. On-treatment change: responders vs progressors",
       x = "log2FC On/Pre in PD patients", y = "log2FC On/Pre in PRCR patients")
ggsave(file.path(fig, "05_lfc_PRCR_vs_PD_change.png"), width = 6.5, height = 5, dpi = 300)

# Heatmap: top 40 on-vs-pre genes, paired samples, ordered by time then response
top <- head(filter(res_onpre, !is.na(symbol)), 40)
v   <- assay(vsd)[top$entrez, colnames(onpre)]
rownames(v) <- top$symbol
ann <- as.data.frame(colData(onpre))[, c("time", "response")]
ord <- order(ann$time, ann$response)
pheatmap(t(scale(t(v[, ord]))), cluster_cols = FALSE, show_colnames = FALSE,
         annotation_col = ann[ord, ], breaks = seq(-2.5, 2.5, length.out = 101),
         color = colorRampPalette(c("#2471A3", "white", "#C0392B"))(100),
         annotation_colors = list(response = c(PD = "#C0392B", SD = "#E1A23B", PRCR = "#2471A3"),
                                  time = c(Pre = "grey80", On = "grey30")),
         main = "Top 40 genes changed on nivolumab (row z-score)", fontsize_row = 7,
         filename = file.path(fig, "06_heatmap_top_on_vs_pre.png"), width = 8, height = 7)
