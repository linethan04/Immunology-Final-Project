# Run the full RNA workflow from the RNA_Ethan folder. First run takes about 20 minutes (DESeq2);
# later runs reuse the fits in rna/results/cache.
steps <- c("01_preprocess_qc", "02_differential_expression", "03_gsea", "04_signatures",
           "05_clinical_links", "rna_bootstrap_ci", "paper_comparison")
for (s in steps) { message("== ", s); source(file.path("rna", paste0(s, ".R")), local = new.env()) }
writeLines(capture.output(sessionInfo()), "rna/results/R_session_info.txt")
