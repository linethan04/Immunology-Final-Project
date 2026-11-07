# One-time install of any missing packages. Data come from ../data (shared repo folder).
cran <- c("tidyverse", "ggrepel", "pheatmap", "patchwork", "pROC", "survival", "broom", "rmarkdown", "BiocManager")
bioc <- c("DESeq2", "apeglm", "org.Hs.eg.db", "fgsea", "genefilter")
msig <- "msigdbr"
need <- setdiff(c(cran, msig), rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
need <- setdiff(bioc, rownames(installed.packages()))
if (length(need)) BiocManager::install(need, update = FALSE)
