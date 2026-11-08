# Step 6: RNA vs Krithika's TCR clonal expansion (tcr-tracking/outputs). Figure 06 + table.
# Paper analogue: Riaz Fig. 5 counts clones that expand on treatment and relates them to response.
# Here we ask whether those expansions show up in the RNA of the same biopsy.
source("rna/common.R"); source("rna/theme_pub.R")
tcr_out <- "../tcr-tracking/outputs"
exp <- read_csv(file.path(tcr_out, "expansion_summary_aa_vj.csv"), show_col_types = FALSE) |>
  transmute(patient = paste0("Pt", patient), n_expanded, n_expanded_novel,
            expanded_rar = rar_n_expanded_mean, frac_expanded = frac_on_templates_in_expanded)
div <- read_csv(file.path(tcr_out, "diversity_summary_aa.csv"), show_col_types = FALSE) |>
  transmute(patient = paste0("Pt", patient), clon_rar = rar_fixed_clonality_delta)
x <- read_csv(file.path(out, "rna_tcr_wes_merged.csv"), show_col_types = FALSE) |> select(-Response) |>
  inner_join(exp, by = "patient") |> left_join(div, by = "patient") |>
  filter(!is.na(CD8_T_delta)) |> mutate(Response = resp(response))
cat("Patients with paired RNA and TCR expansion:", nrow(x), "\n")

rna_feats <- c(CD8_T_delta = "CD8 T-cell score change", Cytotoxic_delta = "Cytotoxic score change",
               IFNG_6_delta = "IFN-γ score change", Tcell_inflamed_18_delta = "T-cell-inflamed change",
               Melanocytic_delta = "Tumor score change", CD8_T_On = "CD8 T-cell score, On")
tcr_feats <- c(expanded_rar = "Expanded clones (depth-matched)", frac_expanded = "Share of On T cells in expanded clones",
               clon_rar = "Clonality change", tcr_fraction_On = "T-cell fraction (On)")
cors <- expand_grid(rna = names(rna_feats), tcr = names(tcr_feats)) |>
  mutate(res = map2(rna, tcr, \(a, b) {
    t <- cor.test(x[[a]], x[[b]], method = "spearman", exact = FALSE)
    tibble(rho = unname(t$estimate), p = t$p.value, n = sum(complete.cases(x[[a]], x[[b]])))
  })) |> unnest(res) |> mutate(padj = p.adjust(p, "BH"), rna = rna_feats[rna], tcr = tcr_feats[tcr])
write_csv(cors, file.path(out, "rna_vs_tcr_expansion.csv")); print(cors, n = 30)

# Expansion by response (paper-style), and agreement with RNA
wx <- wilcox.test(expanded_rar ~ response, data = filter(x, response != "SD"), exact = FALSE)$p.value
grp <- x |> summarise(n = n(), med_expanded = median(expanded_rar, na.rm = TRUE), .by = Response) |> mutate(p_PRCR_vs_PD = wx)
write_csv(grp, file.path(out, "tcr_expansion_by_response.csv")); print(grp)

lab <- \(a, b) with(filter(cors, rna == rna_feats[a], tcr == tcr_feats[b]), sprintf("Spearman ρ = %.2f, p = %.2g, n = %d", rho, p, n))
pA <- ggplot(x, aes(Response, expanded_rar, colour = Response)) + resp_box() +
  labs(x = NULL, y = "Expanded TCR clones\n(depth-matched)", title = "A  TCR expansion by response",
       subtitle = sprintf("PR/CR vs PD Wilcoxon p = %.2g", wx))
sc <- function(xv, yv, xl, yl, title) ggplot(x, aes(.data[[xv]], .data[[yv]], colour = Response)) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, colour = "grey30", fill = "grey88", linewidth = .6, formula = y ~ x) +
  geom_point(size = 3, alpha = .9) + scale_colour_manual(values = resp_cols) + scale_x_log10() +
  labs(x = xl, y = yl, title = title, subtitle = lab(yv, xv))
t8 <- cor.test(x$CD8_T_On, x$tcr_fraction_On, method = "spearman", exact = FALSE)
pB <- ggplot(x, aes(tcr_fraction_On, CD8_T_On, colour = Response)) +
  geom_smooth(aes(group = 1), method = "lm", colour = "grey30", fill = "grey88", linewidth = .6, formula = y ~ x) +
  geom_point(size = 3, alpha = .9) + scale_colour_manual(values = resp_cols) +
  labs(x = "TCR-seq: T-cell fraction of cells (On)", y = "RNA: CD8 T-cell score (On)", title = "B  How many T cells: agree",
       subtitle = sprintf("Spearman ρ = %.2f, p = %.1g, n = %d", t8$estimate, t8$p.value, sum(complete.cases(x$CD8_T_On, x$tcr_fraction_On))))
pC <- sc("expanded_rar", "CD8_T_delta", "Expanded TCR clones (log scale)", "RNA: CD8 score change", "C  Which clones expand: weak")
p <- (pA + pB + pC) + plot_layout(guides = "collect", widths = c(.8, 1, 1)) & theme(legend.position = "bottom")
save_plot(p, "06_rna_vs_tcr_expansion", 13, 5)
