# Presenter script: RNA-seq arm (RNA_seq_nivolumab.pptx)

About 1–1.5 min per slide, roughly 20–25 min total. **Bold** = the one sentence to land on each slide.
Figure files are in `figures/rna/`; numbers come from `output/rna/`.

---

## Part 0: Introduction

**Slide 1 · Title**
"We re-analyzed the CheckMate 038 data from Riaz et al. (2017). Patients with advanced melanoma were biopsied before nivolumab and again about four weeks in. I led the RNA-seq arm, Krithika the TCR-seq and Kwaku the exome and integration. The coloured dots are our response groups, used in every figure: red = progressive disease (PD), gold = stable disease (SD), blue = partial/complete response (PRCR)."

**Slide 2 · Why biopsy during treatment?**
- Nivolumab blocks PD-1, an inhibitory receptor on exhausted T cells. If tumor-reactive T cells are present, blocking PD-1 lets them kill.
- Baseline biomarkers (PD-L1, mutation burden) only weakly predict who benefits. **The on-treatment biopsy lets us see the immune response while it is happening.**
- Dataset: 103 RNA-seq samples after QC, 42 patients with both biopsies, and 25 patients with all three data types plus a response call.

**Slide 3 · Aims**
Aim 1 characterizes each modality, Aim 2 asks whether events co-occur in the same patients, and Aim 3 asks about prediction. **RNA-seq provides the immune and tumor-state measurement that the other two modalities get compared against.**

---

## Part 1: Workflow and code

**Slide 4 · Section divider.**

**Slide 5 · Pipeline**
Five scripts, run in order by `scripts/rna_00_run_all.R`:
1. `rna_01_preprocess_qc.R`: read raw counts (Entrez IDs × samples); fix the Pt7/Pt72 label by taking the patient from the sample name; keep genes with ≥10 reads in ≥5 samples; apply the variance-stabilizing transform (VST) so samples can be compared on one scale.
2. `rna_02_differential_expression.R`: DESeq2 negative-binomial models for three questions (next slide).
3. `rna_03_gsea.R`: rank all genes by the DESeq2 Wald statistic and ask whether whole Hallmark pathways sit at the top or bottom of the list (fgsea).
4. `rna_04_immune_signatures.R`: turn about 20,000 genes into 12 interpretable scores per sample: cell types (T, CD8, NK, B, myeloid, Treg), programs (IFN-γ, T-cell-inflamed, checkpoints, antigen presentation) and a melanocytic score as a tumor-content proxy.
5. `rna_05_clinical_links.R`: link those scores to response (AUC), survival (Kaplan–Meier/Cox) and the TCR/WES data.
**The hand-off file `rna_scores_per_patient.csv` is what Kwaku joins into the integrated table.**

**Slide 6 · Statistical design choices**
- *Paired design* `~ patient + time`: every patient gets their own intercept, so we only test within-patient change. This matters because tumors differ enormously between people.
- *Nested design* `~ response + response:patient + response:time`: the same pairing, but the On-vs-Pre effect is estimated separately per response group. Then we test PRCR change minus PD change, which is a formal interaction test. **This tests whether the change is specific to responders, rather than just whether anything changes.**
- *Signature scores*: z-score each gene across samples, then average the genes in a set. This is simple and transparent, and it is robust to any one gene being noisy.
- *Effect sizes*: with 9–24 patients per group, p-values alone mislead. We report the Hodges–Lehmann shift (the median of all PRCR−PD differences) with a 95% CI, and AUC with DeLong CIs.

---

## Part 2: Results

**Slide 7 · Section divider.**

**Slide 8 · QC (figs 01, 02)**
- *Library size*: bars are reads per sample, coloured by response. They are evenly sized and the colours are mixed, with no group systematically shallower.
- *PCA*: each point is a sample, using the 500 most variable genes. Responders and progressors overlap, so **no technical batch lines up with our comparison.** The spread along PC1 comes from a few outlier samples with unusual tumor/stroma composition.

**Slide 9 · Pre-treatment PRCR vs PD (fig 03)**
- Volcano plot: x = log2 fold change (right = higher in responders), y = significance. Coloured = FDR < 0.1.
- 223 genes pass, but **the strongest hits are skin/keratinocyte genes (SPRR2A, SERPINB3, LCE family), likely from biopsy sampling rather than tumor biology.** Immune genes such as CXCL13 lean toward responders.
- Conclusion: the baseline biopsy gives at most a weak predictive signal, as the original paper also found.

**Slide 10 · On vs Pre, all patients (fig 04)**
- The same plot type, now comparing each patient's On biopsy to their own Pre biopsy.
- 44 genes at FDR < 0.1, led by **PDCD1 (PD-1 itself), CXCL9/10 (IFN-γ-induced chemokines that recruit T cells), GZMH/PRF1 (killing machinery), IFNG, LAG3**. This is the pharmacodynamic signature of PD-1 blockade.

**Slide 11 · Responders vs progressors change (fig 05)**
- Each dot is a gene. x = its On/Pre change in PD patients, y = the change in PRCR patients. Points on the diagonal change equally in both groups. Purple = significant interaction.
- PRCR: 3,582 genes change. PD: **zero**. SD: 28.
- The immune genes sit above the diagonal: they go up in responders. **The big drops in responders (bottom of the plot) are melanoma lineage genes (SOX10, MLANA, PRAME, PAX3), so tumor cells are disappearing from the biopsy.**

**Slide 12 · Heatmap (fig 06)**
- Rows = the top 40 On-vs-Pre genes, columns = samples (Pre on the left, On on the right, coloured by response). Red = above that gene's average.
- The On half is redder, and the PRCR On samples are the reddest. The genes cluster into one co-regulated T-cell/IFN-γ program.

**Slide 13 · Pathways (fig 07)**
- Dot colour = normalized enrichment score (red = up on treatment / in PRCR), dot size = significance, outline = FDR < 0.05.
- **Responders: interferon-γ/α, allograft rejection and inflammatory pathways go up; E2F, MYC and G2M proliferation pathways crash.** These are two sides of one event: immune attack plus loss of dividing tumor cells.
- SD shows the strongest interferon induction of any group but little proliferation loss. PD shows modest immune changes only.

**Slide 14 · Patient trajectories (fig 09)**
- Each thin line is one patient's score at Pre and On; the thick line is the group median. Rows = response group.
- PRCR: CD8 and cytotoxic scores rise, and the melanocytic score collapses (median falls from about 0 to −1.6).
- SD: the immune scores rise as much as or more than in PRCR, but the tumor barely changes.
- PD: flat.
- **Interpretation: immune activation is necessary but not sufficient; what separates PRCR from SD is whether the tumor is actually cleared.**

**Slide 15 · Effect sizes (fig 10)**
- Each row is a score; the dot is the PRCR−PD difference with a 95% CI. Red = p < 0.05.
- Pre-treatment: every interval crosses zero.
- On-treatment: 10 of 12 scores are higher in responders (CD8 shift +1.4 SD, FDR 0.006).
- Change: the clearest separator is melanocytic loss.
- **This is the central finding: the information is in the on-treatment biopsy, not the baseline.**

**Slide 16 · Per-patient heatmap (fig 11)**
- Columns = the 25 tri-modal patients grouped by response, rows = the On−Pre change in each score.
- Some PD patients (e.g. Pt28, Pt31) gain immune signal without tumor loss. **These discordant patients are where TCR clonal expansion and neoantigen loss (Aim 2) can explain what RNA alone can't.**

**Slide 17 · AUC (fig 13)**
- AUC = the probability that a random responder scores higher than a random progressor (0.5 = coin flip).
- Pre-treatment scores: AUC 0.60–0.68, with CIs reaching 0.5.
- On-treatment CD8: 0.87 (0.75–1.00).
- Melanocytic change: 0.96, but this is partly circular, because tumor shrinkage is how response is defined.
- Caveat: these are single features with direction chosen after the fact, so they are optimistic. The real test is the cross-validated multi-modal model in Aim 3.

**Slide 18 · Survival (fig 14)**
- Kaplan–Meier overall survival, split at the median on-treatment T-cell-inflamed score.
- The high group lives longer (HR 0.74 per SD, 95% CI 0.49–1.10, log-rank p = 0.09). **The trend points the right way but isn't significant at n = 48.**
- Caveat (landmark bias): only patients alive at week 4 could be biopsied, so survival really starts at the biopsy.

**Slide 19 · Cross-modality (fig 15)**
- Top row: **the RNA CD8 score correlates with the T-cell fraction measured independently by TCR-seq (ρ = 0.81 Pre, 0.88 On)**, which validates the RNA scores.
- Bottom-left: change in TCR clonality (expansion of a few clones) vs change in cytotoxic score. There is a weak positive trend (ρ = 0.27), which the TCR analysis can refine.
- Bottom-right: mutation burden doesn't predict baseline inflammation (ρ = 0.14), so WES and RNA carry different information. That is good for a combined model.

**Slide 20 · Tumor-content caveat (fig 12)**
- Each point is a sample: x = melanocytic (tumor) score, y = T-cell-inflamed score. They are weakly anti-correlated (ρ = −0.24).
- In bulk RNA, when tumor cells die, every non-tumor signal rises proportionally. **Part of the immune "increase" in responders may be dilution from tumor loss.** The next step is to adjust for purity (ESTIMATE or proper deconvolution) before the final integration.

**Slide 21 · Conclusions and next steps**
1. Baseline RNA barely predicts response.
2. In responders, nivolumab triggers a coordinated IFN-γ / CD8 / cytotoxic program while tumor and proliferation genes collapse.
3. PD tumors stay flat. SD tumors activate immunity but don't clear the tumor.
4. On-treatment CD8 (AUC 0.87) is the best single RNA feature, and it agrees with TCR-seq.
Next: purity adjustment and deconvolution, merging with TCR/WES per patient, a nested-CV classifier, and finalizing the response definitions (Pt23 = NE, SD).

---

### Likely questions
- *Why does your Pre DE count (223) differ from the paper?* We used a simpler gene filter and current DESeq2 without IHW. The direction and the top genes are similar.
- *Why not CIBERSORT?* `immunedeconv` isn't installed yet. Marker scores are a transparent first pass, and they validate against TCR-seq.
- *Isn't SD a response?* Clinically it's "benefit" for some definitions. We kept it separate because the RNA shows it is biologically intermediate.
