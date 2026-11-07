# Step 1: read counts and sample table, filter genes, normalize (VST), QC figures
source("rna/common.R")
suppressPackageStartupMessages(library(org.Hs.eg.db))

counts <- read.delim(file.path(data_dir, "CountData.BMS038.txt"), check.names = FALSE)  # rows = Entrez IDs
geno   <- read.csv(file.path(data_dir, "genomic_data_per_case.csv"))
shared <- geno$Patient[rowSums(geno[-1]) == 6]  # WES + TCR + RNA at both biopsies

samples <- read.csv(file.path(data_dir, "SampleTableCorrected.9.19.16.csv")) |>
  transmute(sample   = Sample,
            patient  = sub("_.*", "", Sample),  # PatientID column mislabels Pt7_Pre as Pt72
            time     = factor(PreOn, c("Pre", "On")),
            response = factor(Response, c("PD", "SD", "PRCR")),
            benefit  = factor(Response.1, c("NB", "LB")),
            cohort   = Cohort) |>
  filter(!is.na(response), sample %in% colnames(counts)) |>
  mutate(shared = patient %in% shared)

dds <- DESeqDataSetFromMatrix(as.matrix(counts[, samples$sample]),
                              column_to_rownames(samples, "sample"), ~ 1)
dds <- dds[rowSums(counts(dds) >= 10) >= 5, ]
rowData(dds)$symbol <- mapIds(org.Hs.eg.db, rownames(dds), "SYMBOL", "ENTREZID")
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)
saveRDS(dds, file.path(cache, "dds.rds"))
saveRDS(vsd, file.path(cache, "vsd.rds"))
write_csv(samples, file.path(out, "samples.csv"))

cat(ncol(dds), "samples,", nrow(dds), "genes after filtering\n")
print(with(samples, table(response, time)))

# QC: reads per sample
qc <- as.data.frame(colData(dds)) |> rownames_to_column("sample") |>
  mutate(reads_M = colSums(counts(dds)) / 1e6, Response = resp(response))
p <- ggplot(qc, aes(reorder(sample, reads_M), reads_M, fill = Response)) + geom_col() +
  facet_grid(~ time, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = resp_cols) +
  labs(x = NULL, y = "Assigned reads (millions)", title = "Library size per sample") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_supp(p, "S01_library_size", 9, 4)

# QC: PCA on the 500 most variable genes
pca <- plotPCA(vsd, intgroup = c("time", "response"), ntop = 500, returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, colour = resp(response), shape = time)) +
  geom_point(size = 2.5, alpha = 0.85) + scale_colour_manual(values = resp_cols) +
  labs(x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)"), colour = "Response", shape = NULL,
       title = "PCA of normalized expression (top 500 genes)")
save_supp(p, "S02_pca", 7, 5)
