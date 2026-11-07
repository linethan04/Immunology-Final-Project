# Step 3: pathway view. fgsea on MSigDB Hallmark sets, genes ranked by the DESeq2 Wald statistic.
source("rna/common.R")
set.seed(1)
h <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark <- split(h$gene_symbol, h$gs_name)
contrasts <- c(`PR/CR vs PD` = "A_pre_PRCR_vs_PD", `All patients` = "B_on_vs_pre_paired",
               `PR/CR` = "C_on_vs_pre_PRCR", SD = "C_on_vs_pre_SD", PD = "C_on_vs_pre_PD")
gsea <- imap_dfr(contrasts, \(f, nm) {
  r <- read_csv(file.path(out, paste0("de_", f, ".csv")), show_col_types = FALSE) |>
    filter(!is.na(stat), !is.na(symbol)) |> slice_max(abs(stat), n = 1, by = symbol, with_ties = FALSE)
  as_tibble(fgsea(hallmark, setNames(r$stat, r$symbol), minSize = 15, maxSize = 500)) |> mutate(contrast = nm)
})
write_csv(select(gsea, -leadingEdge), file.path(out, "gsea_hallmark.csv"))

# Figure 04: the two programs the paper reports, immune up and proliferation down
programs <- c(HALLMARK_INTERFERON_GAMMA_RESPONSE = "Interferon-gamma response", HALLMARK_INTERFERON_ALPHA_RESPONSE = "Interferon-alpha response",
              HALLMARK_ALLOGRAFT_REJECTION = "Allograft rejection (T-cell killing)", HALLMARK_INFLAMMATORY_RESPONSE = "Inflammatory response",
              HALLMARK_E2F_TARGETS = "E2F targets", HALLMARK_G2M_CHECKPOINT = "G2M checkpoint",
              HALLMARK_MYC_TARGETS_V1 = "MYC targets", HALLMARK_MITOTIC_SPINDLE = "Mitotic spindle")
d <- gsea |> filter(pathway %in% names(programs), contrast != "All patients") |>
  mutate(label = factor(programs[pathway], rev(programs)),
         program = factor(if_else(pathway %in% names(programs)[1:4], "Immune", "Proliferation"), c("Immune", "Proliferation")),
         column = factor(contrast, c("PR/CR vs PD", "PR/CR", "SD", "PD")),
         panel = factor(if_else(contrast == "PR/CR vs PD", "Before treatment", "On vs Pre, by response"), c("Before treatment", "On vs Pre, by response")),
         sig_nes = if_else(padj < 0.05, NES, NA_real_))
p <- ggplot(d, aes(column, label, fill = sig_nes)) +
  geom_tile(colour = "white", linewidth = 1) + geom_text(aes(label = sprintf("%.1f", NES)), size = 4.6) +
  facet_grid(program ~ panel, scales = "free", space = "free") +
  scale_fill_gradient2(low = DOWN, mid = "white", high = UP, midpoint = 0, limits = c(-3.8, 3.8), na.value = "grey93", name = "NES") +
  labs(x = NULL, y = NULL, title = "Hallmark pathways before and during nivolumab",
       caption = "Number = normalized enrichment score. Orange = higher in PR/CR (left) or on treatment (right); grey = FDR >= 0.05.") +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey92"))
save_plot(p, "04_pathways_by_response", 10, 6)
gsea |> filter(pathway %in% names(programs)) |> select(contrast, pathway, NES, padj) |> arrange(pathway, contrast) |> print(n = 50)

# Supplementary: top 20 pathways across all contrasts
keep <- gsea |> summarise(p = min(padj), .by = pathway) |> slice_min(p, n = 20, with_ties = FALSE) |> pull(pathway)
d2 <- gsea |> filter(pathway %in% keep) |>
  mutate(pathway = str_to_sentence(str_replace_all(str_remove(pathway, "HALLMARK_"), "_", " ")),
         contrast = factor(contrast, names(contrasts)))
lvl <- d2 |> filter(contrast == "PR/CR") |> arrange(NES) |> pull(pathway)
d2$pathway <- factor(d2$pathway, unique(c(lvl, d2$pathway)))
p <- ggplot(d2, aes(contrast, pathway, colour = NES, size = -log10(padj))) + geom_point() +
  geom_point(data = filter(d2, padj < 0.05), shape = 1, colour = "black", stroke = 0.6) +
  scale_colour_gradient2(low = DOWN, mid = "white", high = UP) + scale_size_continuous(range = c(0.5, 5)) +
  labs(x = NULL, y = NULL, title = "Hallmark enrichment, top 20 pathways", subtitle = "Outlined = FDR < 0.05") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.text.y = element_text(size = 10))
save_supp(p, "S07_gsea_top20", 9, 7)
