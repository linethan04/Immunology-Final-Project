# Background and explanation notes: Riaz et al. 2017 reanalysis (RNA section)

# Modality order when all three pairs are compared

**Keep WES → TCR → RNA.** Comparing every pair makes this order work even better, because each speaker can compare their data with what has already been shown:

| Order | Speaker | Ends with this comparison |
|---|---|---|
| 1. WES | Kwaku | none yet (nothing earlier to compare with) |
| 2. TCR | Krithika | **WES vs TCR**: did patients whose tumors lost mutations also gain T-cell clones? This is the link the paper itself draws. |
| 3. RNA | Ethan | **WES vs RNA** and **TCR vs RNA**: RNA checks both of the earlier results in the tissue. |

This also follows the paper's own order (DNA, then T cells, then gene expression) and the biology (targets, then attackers, then outcome). Put a single 3-pair comparison slide at the very end, after all three of us have presented. Kwaku and Krithika fill in WES vs TCR, and Ethan presents TCR vs RNA. WES vs RNA can be Ethan's or Kwaku's; the numbers are already available (below).

---

# 1. The background, from scratch

**Cells, DNA, RNA, protein.** Every cell carries DNA, the instruction manual. To use one instruction (a **gene**), the cell copies it into **RNA**, and the RNA is translated into a **protein** that does the work. The DNA is the same in almost every cell. Which genes are being copied into RNA, and how much, is different for each cell type. That is why RNA tells you what a cell is *doing*, while DNA tells you what it *could* do.

**Melanoma.** A skin cancer that starts in **melanocytes**, the pigment cells that make melanin. Melanomas usually still switch on melanocyte genes (MLANA, PMEL, TYR, DCT, MITF), so those genes act as a tag for "tumor cells are here." Melanomas carry a very large number of **mutations** (DNA spelling errors), largely from UV damage. That matters for the immune system.

**T cells and the T-cell receptor (TCR).** T cells are immune cells that kill infected or abnormal cells. Each T cell has one receptor, the **TCR**, that recognizes one specific target. A **clone** (or **clonotype**) is a family of T cells descended from one parent, all sharing the same TCR. When a T cell finds its target it multiplies, and that is called **clonal expansion**. **CD8 T cells** are the killer type.

**Neoantigens.** A mutated tumor protein can be chopped up and displayed on the cell surface. If that fragment looks foreign, it is a **neoantigen**, a target a T cell can recognize. More mutations usually means more potential neoantigens. **Tumor mutational burden (TMB)** is the count of mutations.

**PD-1, the brake.** T cells carry **PD-1**, an "off switch." Tumors display **PD-L1**, which presses that switch, so T cells that have found the tumor get shut down. **Nivolumab** is an antibody that blocks PD-1 and takes the brake off. This is **checkpoint blockade immunotherapy**: it doesn't attack the tumor directly; it frees the patient's own T cells to do it. **Ipilimumab** blocks a different brake, **CTLA-4**.

**CheckMate 038.** A clinical trial of nivolumab in advanced melanoma run by Bristol-Myers Squibb (the data files are labelled "BMS038"). Its special feature is that it was designed to take tumor **biopsies** (small tissue samples) **before treatment and again about 4 weeks into treatment**. Most trials only have the "before" sample. It enrolled two groups:
- **Ipilimumab-naive**: never had ipilimumab.
- **Ipilimumab-progressed**: already treated with ipilimumab, and the cancer grew anyway.

The released clinical table has 73 patients.

**Responders vs non-responders.** Response is judged on scans using standard criteria (RECIST):
- **CR**, complete response: the tumor disappears.
- **PR**, partial response: shrinks by at least 30%.
- **SD**, stable disease: neither shrinks enough nor grows enough.
- **PD**, progressive disease: grows by at least 20%, or new tumors appear.

We group PR and CR as **"PR/CR" (responders)** and compare them mainly with **PD (non-responders)**. SD is in between. The paper sometimes used **"benefit" (CR/PR/SD) vs "no benefit" (PD)** instead, which is why some numbers differ by definition.

**The research question.** When nivolumab works, what changes inside the tumor, and can we see it in the tumor's DNA, its T cells and its gene expression? A second question: does anything measured *before* treatment predict who will respond?

---

# 2. Each modality: what it measures, how, and what the paper found

## WES: whole-exome sequencing (Kwaku)

- **What it is.** Reads the DNA sequence of the **exome**, the ~2% of the genome that codes for proteins. Tumor DNA is compared with normal DNA (from blood) to find **somatic mutations**, the ones present only in the tumor.
- **Key methods:**
  - **TMB**: mutation count.
  - **Neoantigen prediction**: will a mutated fragment fit the patient's display molecules (HLA)?
  - **Clonality**, estimated with a tool called PyClone: what fraction of tumor cells carries each mutation? A mutation in every tumor cell is "clonal"; one in a subset is "subclonal."
- **Contraction vs expansion.** Comparing the two biopsies, a mutation **contracts** if the share of cells carrying it drops on treatment, which means the cells carrying it were killed. It **expands** if that share grows.
- **What WES is for.** Shows which targets exist, and whether the tumor cells carrying them disappear on treatment.
- **Paper's finding.** Responders' tumors showed much more mutation contraction, and the lost mutations were enriched for neoantigens. That is evidence the immune system eliminated the cells displaying those targets. Baseline TMB was only modestly associated with response.

## TCR-seq: T-cell receptor sequencing (Krithika)

- **What it is.** Sequences the most variable part of each T-cell receptor in the biopsy, the **CDR3** region of the TCR beta chain (using the Adaptive ImmunoSEQ assay). Each unique CDR3 sequence is one clonotype, and its read count shows how big that clone is.
- **Key measures:**
  - **Richness**: how many different clones there are.
  - **Clonality**: 1 − evenness (evenness is based on Shannon entropy). High clonality means a few clones dominate, as when a few T-cell families are multiplying against a target.
  - **T-cell fraction**: the share of cells in the biopsy that are T cells.
  - **Expanded clones**: clones significantly bigger in the on-treatment biopsy than before. Krithika also used **depth-matched** (rarefied) counts, so a biopsy sequenced more deeply doesn't look like it has more expansions.
  - **Tracking**: which clones were already present before treatment and which are new.
- **What TCR-seq is for.** Shows which T cells responded: is the attack driven by a few multiplying clones?
- **Paper's finding.** The T-cell repertoire changed on treatment, and the diversity changes differed by cohort. The paper also reports that clonal expansion went along with mutation contraction, which is the WES-TCR link.

## RNA-seq: RNA sequencing (Ethan)

- **What it is.** Captures the RNA from the whole biopsy and counts how many RNA fragments ("reads") come from each gene. This is **bulk RNA-seq**: the biopsy is a smoothie of tumor cells, T cells, other immune cells and support tissue, all blended, so each gene's count is a mixture.
- **Key methods:**
  - **Normalization**: correcting for how deeply each sample was sequenced (variance-stabilizing transformation, VST).
  - **Differential expression (DESeq2)**: for each gene, is it higher or lower between biopsies or groups beyond what noise explains?
  - **Signature scores**: average several marker genes into one number, such as a CD8 score from CD8A/CD8B, a tumor score from the melanocyte genes, or the 18-gene T-cell-inflamed score.
  - **Gene set enrichment (fgsea)**: do whole pathways (sets of genes that work together) move together?
- **What RNA-seq is for.** Shows what the tissue is doing: how much immune activity there is, and how much tumor is left. It is the only one of the three that sees both the immune side and the tumor side at once.
- **Paper's finding.** In responders, immune genes went up on treatment, while cell-cycle (division) and melanin/melanocyte genes went down. Baseline expression was only marginally predictive.

---

# 3. Statistical terms used on the slides

- **Median**: the middle value. It is used instead of the average because a few extreme patients would distort an average.
- **z-score**: a gene's value expressed as "how many standard deviations above or below that gene's average across all samples." Averaging z-scores across genes makes a signature score, where 0 is average, positive is higher than typical and negative is lower.
- **2-fold change**: the level doubles or halves.
- **p-value**: if there were truly no difference, how surprising would data like this be? Small (below 0.05) means "probably a real difference." It is not the size of the effect.
- **Wilcoxon (rank-sum) test**: compares two groups by ranking the values, so it is robust to outliers and small groups.
- **95% confidence interval (CI)**: the range of plausible values for the true effect. If it **crosses zero**, "no difference" is still plausible.
- **Bootstrap**: estimating that range by resampling the patients thousands of times.
- **Spearman ρ (rho)**: a correlation from −1 to +1, based on ranks. +1 means they rise together perfectly, 0 means no relationship, −1 means one rises as the other falls.
- **FDR (false discovery rate)**: when you test 20,000 genes, some come up "significant" by luck. FDR < 0.1 means at most ~10% of the genes you call changed are expected to be false.
- **NES (normalized enrichment score)**: how strongly a whole pathway shifts up (positive) or down (negative).
- **Linear model adjustment**: estimating the responder effect while holding another variable, here the tumor score, constant.

---

# 4. The RNA figures in depth

The RNA section runs as a story: before treatment there is barely a signal → responders change more consistently on treatment → the change that marks a response is tumor loss → immune activation happens in everyone, but proliferation collapses only in responders → what's left in responders is CD8-rich even after correcting for tumor loss.

**Slide 4: workflow (no figure).** Counts, then normalization, then signature scores, then comparing each patient's two biopsies. Comparing a patient with themselves matters because melanomas differ hugely between people. 49 patients have a pre-treatment biopsy, and 42 have both (9 PR/CR, 15 SD, 18 PD).

**Slide 5: baseline immune activity (recreates the paper).**
- **What it shows.** One dot per patient: the T-cell-inflamed score from the pre-treatment biopsy (18 genes from Ayers et al. 2017, measuring how much IFN-γ-driven T-cell activity a tumor has). Boxes show the median and middle half. The two panels split the two cohorts.
- **Result.** PR/CR median 0.19, SD −0.47, PD −0.34. PR/CR minus PD is +0.53, but the 95% CI (−0.73 to 1.33) crosses zero, p = 0.20. The ipilimumab-progressed responders all start "hot" (immune-active); in ipilimumab-naive patients the groups overlap.
- **Meaning.** One snapshot before treatment can't sort patients well. The paper agrees ("marginally predictive"). What we changed: one fixed published score instead of genome-wide testing, and the cohort split shown directly.

**Slide 6: how many genes change (builds on the paper).**
- **What it shows.** One bar per patient: the % of genes that went up (orange) or down (blue) at least 2-fold between biopsies, grouped by response.
- **Result.** Median 8.6% (PR/CR), 4.7% (SD), 4.5% (PD), p = 0.048. Two sequencing runs of the same biopsy differ by ~3%, so that is the noise floor. DESeq2 genes shifting *consistently* across a group: 3,582 in PR/CR, 28 in SD, 0 in PD.
- **Meaning.** Responders change in a shared direction; progressors change, if at all, each in their own way.
- **Link to the paper.** The paper reported counts of changed genes, but we could not reproduce their exact numbers: before treatment we get 397 genes vs their 189 with their design. The biological direction does agree.

**Slide 7: tumor loss vs CD8 gain per patient (new).**
- **What it shows.** Each row is a patient, sorted by the change in the tumor (melanocyte-gene) score. The left column is the tumor score change, the right column is the CD8 score change. Left of the dashed line means the score fell.
- **Result.** Tumor score: PR/CR −1.07 vs PD −0.01, p = 0.0002. CD8: PR/CR +0.47, SD +0.31, PD 0.00, p = 0.07. Exceptions: Pt94 and Pt101 respond with little tumor-score change.
- **Meaning.** The clearest RNA mark of response is the tumor disappearing from the biopsy.
- **Caveat.** Response is *defined* by tumor shrinkage, so this partly restates the outcome. It describes response rather than predicting it.

**Slide 8: pathways (recreates the paper).**
- **What it shows.** A heatmap of NES. The left column compares responders vs progressors *before* treatment. The right columns show *on vs pre* within each group. The top block is immune pathways (IFN-γ, IFN-α, allograft rejection, which is essentially T-cell killing, and inflammatory response). The bottom block is proliferation pathways (E2F targets, G2M checkpoint, MYC targets, mitotic spindle), which are the cell-division machinery. Grey means not significant.
- **Result.** Immune pathways rise in every group, including PD (IFN-γ +2.7 PR/CR, +3.5 SD, +2.2 PD). Proliferation collapses only in responders (E2F −3.5) and rises in PD (+1.3).
- **Meaning.** Nivolumab triggers immune signalling broadly, but only in responders does the tumor stop dividing. Immune activation alone isn't enough.
- **Link to the paper.** Same conclusion as the paper (immune up, cell cycle down in responders); adding the SD and PD columns is what makes "not enough" visible.

**Slide 9: the dilution check (new).**
- **The problem.** Bulk RNA is a mixture. If tumor cells die, T cells become a bigger share of what's left even if no new T cells arrive, like removing fruit from a smoothie so the remaining ingredients' share rises.
- **Panel A.** On-treatment CD8 score by response (p = 0.00035).
- **Panel B.** The PR/CR − PD effect with and without adjusting for the tumor score. The CD8 *rise* goes from +0.67 to +0.21, and its CI crosses zero, so it can't be separated from dilution. The on-treatment CD8 *level* goes from +1.30 to +1.21 and holds.
- **Meaning.** By week 4, responders' biopsies are genuinely CD8-rich, beyond what tumor loss alone explains.

**Slide 10: TCR vs RNA comparison (new analysis, 24 patients with both).**
- **Panel A.** Expanded clones by response (depth-matched): PR/CR 12, SD 4, PD 9, p = 0.75. This mirrors the paper's expansion-by-response comparison and is not significant in this subset.
- **Panel B.** RNA CD8 score vs TCR T-cell fraction: ρ = 0.86. Two completely different assays agree on *how many* T cells are there, which also validates the RNA CD8 score.
- **Panel C.** Expanded clones vs RNA CD8 change: ρ = 0.27, not significant. No expansion or clonality measure tracks any RNA score after correction.
- **Meaning.** RNA measures the size of the immune response; only TCR-seq sees its clonal makeup, meaning which specific T-cell families are multiplying. They add to each other rather than repeat each other.

---

# 5. RNA results in context: the full picture

| Comparison | What it asks | Result |
|---|---|---|
| **WES vs TCR** (Krithika/Kwaku; the paper's link) | Do tumors that lose mutations also gain expanding T-cell clones? | Theirs to fill in |
| **WES vs RNA** (numbers already in the RNA results) | Does losing mutations (WES) match losing tumor gene expression (RNA)? | RNA tumor score change vs WES contraction ρ = −0.69 (n = 25), still −0.50 without responders. Two independent assays see the same tumor loss. Caveat: both partly reflect how much tumor is in the biopsy. Baseline immune score vs mutation burden ρ = 0.14: nearly unrelated, so they carry different information. |
| **TCR vs RNA** (Ethan) | Do they agree about the T-cell response? | They agree on *how many* T cells (ρ = 0.86); RNA doesn't capture clonal *expansion* (ρ ≤ 0.31). |

**Answer to the research question, combining all three.** When nivolumab works, tumor cells carrying mutations are eliminated (WES contraction). The tumor's own gene programs, melanocyte identity and cell division, shut down (RNA). The biopsy fills with CD8 T cells (RNA + TCR agree), with specific clones expanding, which only TCR-seq can show. Immune signalling switches on even in non-responders, so what separates responders is the tumor actually dying. Before treatment, no single measure from any of the three sorts patients well.
