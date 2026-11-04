# Run the full RNA-seq workflow in order (from project root).
# 01 preprocess/QC -> 02 DESeq2 (~10 min) -> 03 GSEA -> 04 signatures -> 05 clinical/cross-modality
# Tables -> output/rna/, figures -> figures/rna/
for (s in sprintf("scripts/rna_0%d_%s.R", 1:5,
                  c("preprocess_qc", "differential_expression", "gsea",
                    "immune_signatures", "clinical_links"))) {
  message("== ", s)
  source(s, local = new.env())
}
