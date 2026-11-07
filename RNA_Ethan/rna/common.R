# Shared paths, colors and helpers. Run everything from the RNA_Ethan folder.
suppressPackageStartupMessages({
  library(DESeq2); library(tidyverse); library(ggrepel); library(pheatmap)
  library(patchwork); library(fgsea); library(msigdbr); library(pROC); library(survival)
})
select <- dplyr::select  # Bioconductor packages mask it
data_dir <- "../data"
out  <- "rna/results"
supp <- file.path(out, "supplementary_figures")
cache <- file.path(out, "cache")
for (d in c(supp, cache)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Same response colors and labels as the WES figures, so the talks read as one
resp_cols <- c(`PR/CR` = "#F8766D", SD = "#00BA38", PD = "#619CFF")
resp <- function(x) factor(x, levels = c("PRCR", "SD", "PD"), labels = c("PR/CR", "SD", "PD"))
UP <- "#E59845"; DOWN <- "#277DA8"; GREY <- "#BBBBBB"

theme_set(theme_bw(base_size = 14))
save_plot <- function(p, name, w = 10, h = 6) {
  ggsave(file.path(out, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(out, paste0(name, ".pdf")), p, width = w, height = h)
}
save_supp <- function(p, name, w = 8, h = 5)
  ggsave(file.path(supp, paste0(name, ".png")), p, width = w, height = h, dpi = 300)

# Run a slow fit once, then reuse it
cached <- function(name, expr) {
  f <- file.path(cache, paste0(name, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  x <- expr; saveRDS(x, f); x
}
