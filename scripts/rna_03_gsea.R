# RNA step 3: pathway-level view -- fgsea on MSigDB Hallmark sets
# Genes ranked by DESeq2 Wald statistic for each contrast from rna_02.

suppressPackageStartupMessages({
  library(tidyverse)
  library(fgsea)
  library(msigdbr)
})
out <- "output/rna"; fig <- "figures/rna"
theme_set(theme_bw(base_size = 11))
set.seed(1)

h <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark <- split(h$gene_symbol, h$gs_name)

contrasts <- c(`Pre: PRCR vs PD`     = "de_A_pre_PRCR_vs_PD",
               `On vs Pre: all`      = "de_B_on_vs_pre_paired",
               `On vs Pre: PRCR`     = "de_C_on_vs_pre_PRCR",
               `On vs Pre: SD`       = "de_C_on_vs_pre_SD",
               `On vs Pre: PD`       = "de_C_on_vs_pre_PD")

gsea <- imap_dfr(contrasts, \(f, nm) {
  r <- read_csv(file.path(out, paste0(f, ".csv")), show_col_types = FALSE) |>
    filter(!is.na(stat), !is.na(symbol)) |>
    slice_max(abs(stat), n = 1, by = symbol, with_ties = FALSE)
  fgsea(hallmark, setNames(r$stat, r$symbol), minSize = 15, maxSize = 500) |>
    as_tibble() |> mutate(contrast = nm)
})
write_csv(select(gsea, -leadingEdge), file.path(out, "gsea_hallmark.csv"))

# Dot plot: 20 pathways with the smallest padj in any contrast
keep <- gsea |> summarise(p = min(padj), .by = pathway) |> slice_min(p, n = 20) |> pull(pathway)  # top 20
d <- gsea |>
  filter(pathway %in% keep) |>
  mutate(pathway = str_remove(pathway, "HALLMARK_") |> str_replace_all("_", " ") |> str_to_sentence(),
         contrast = factor(contrast, names(contrasts)))
lvl <- d |> filter(contrast == "On vs Pre: PRCR") |> arrange(NES) |> pull(pathway)
d$pathway <- factor(d$pathway, unique(c(lvl, d$pathway)))

ggplot(d, aes(contrast, pathway, colour = NES, size = -log10(padj))) +
  geom_point() +
  geom_point(data = filter(d, padj < 0.05), shape = 1, colour = "black", stroke = 0.6) +
  scale_colour_gradient2(low = "#2471A3", mid = "white", high = "#C0392B") +
  scale_size_continuous(range = c(0.5, 5)) +
  labs(x = NULL, y = NULL, title = "Hallmark pathway enrichment (fgsea, top 20 pathways)",
       subtitle = "Outlined = padj < 0.05; NES > 0 = higher in PRCR / On-treatment") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(fig, "07_gsea_hallmark_dotplot.png"), width = 8.5, height = 6.5, dpi = 300)

# Enrichment curve for the headline pathway in responders
r <- read_csv(file.path(out, "de_C_on_vs_pre_PRCR.csv"), show_col_types = FALSE) |>
  filter(!is.na(stat), !is.na(symbol)) |>
  slice_max(abs(stat), n = 1, by = symbol, with_ties = FALSE)
p <- plotEnrichment(hallmark$HALLMARK_INTERFERON_GAMMA_RESPONSE, setNames(r$stat, r$symbol)) +
  labs(title = "IFN-gamma response: On vs Pre in responders (PRCR)",
       x = "Gene rank (up on treatment -> down)", y = "Enrichment score")
ggsave(file.path(fig, "08_gsea_ifng_PRCR.png"), p, width = 6, height = 4, dpi = 300)

gsea |> filter(padj < 0.05) |> arrange(contrast, desc(NES)) |>
  select(contrast, pathway, NES, padj) |> print(n = 60)
