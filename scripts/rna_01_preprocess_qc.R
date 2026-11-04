# RNA step 1: load counts + metadata, filter, VST, QC figures
# Run from project root. Outputs reused by rna_02..rna_05.

suppressPackageStartupMessages({
  library(DESeq2)
  library(org.Hs.eg.db)
  library(tidyverse)
})
out <- "output/rna"; fig <- "figures/rna"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 11))
resp_cols <- c(PD = "#C0392B", SD = "#E1A23B", PRCR = "#2471A3")

counts <- read.delim("data/CountData.BMS038.txt", check.names = FALSE)  # rows = Entrez IDs

geno <- read.csv("data/genomic_data_per_case.csv")
tri  <- geno$Patient[rowSums(geno[-1]) == 6]  # WES + TCR + RNA at both timepoints

samples <- read.csv("data/SampleTableCorrected.9.19.16.csv") |>
  transmute(sample   = Sample,
            patient  = sub("_.*", "", Sample),  # PatientID column mislabels Pt7_Pre as Pt72
            time     = factor(PreOn, c("Pre", "On")),
            response = factor(Response, c("PD", "SD", "PRCR")),
            benefit  = Response.1,
            cohort   = Cohort) |>
  filter(!is.na(response), sample %in% colnames(counts)) |>
  mutate(tri = patient %in% tri)

dds <- DESeqDataSetFromMatrix(as.matrix(counts[, samples$sample]),
                              column_to_rownames(samples, "sample"), ~ 1)
dds <- dds[rowSums(counts(dds) >= 10) >= 5, ]
rowData(dds)$symbol <- mapIds(org.Hs.eg.db, rownames(dds), "SYMBOL", "ENTREZID")
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)

saveRDS(dds, file.path(out, "dds.rds"))
saveRDS(vsd, file.path(out, "vsd.rds"))
write_csv(samples, file.path(out, "samples.csv"))

cat(ncol(dds), "samples,", nrow(dds), "genes after filtering\n")
print(with(samples, table(response, time)))

# ---- QC 1: library size per sample ----
qc <- as.data.frame(colData(dds)) |>
  rownames_to_column("sample") |>
  mutate(lib_M = colSums(counts(dds)) / 1e6)

ggplot(qc, aes(reorder(sample, lib_M), lib_M, fill = response)) +
  geom_col() +
  facet_grid(~ time, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = resp_cols) +
  labs(x = NULL, y = "Assigned reads (millions)", title = "Library size per sample") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(file.path(fig, "01_qc_library_size.png"), width = 8, height = 3.5, dpi = 300)

# ---- QC 2: PCA on top-500 variable genes ----
pca <- plotPCA(vsd, intgroup = c("time", "response"), ntop = 500, returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
ggplot(pca, aes(PC1, PC2, colour = response, shape = time)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_colour_manual(values = resp_cols) +
  labs(x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)"),
       title = "PCA of VST expression (top 500 variable genes)")
ggsave(file.path(fig, "02_qc_pca.png"), width = 6, height = 4.5, dpi = 300)

# What drives PC1? Correlate with marker genes (immune vs melanocytic)
pc   <- prcomp(t(assay(vsd)[head(order(-matrixStats::rowVars(assay(vsd))), 500), ]))
sym  <- rowData(vsd)$symbol
expr <- function(g) assay(vsd)[match(g, sym), ]
cat("cor(PC1, CD8A) =", round(cor(pc$x[, 1], expr("CD8A")), 2),
    "| cor(PC1, MITF) =", round(cor(pc$x[, 1], expr("MITF")), 2), "\n")
