# Presenter script

Matches `RNA_10slides_Ethan.pptx` (also in its speaker notes). Numbers come from `rna/results/`.

## My RNA section: about 4.5 minutes (slides 4, 5, 6, 7, 8, 9)

### Slide 4: RNA-seq turns each biopsy into gene and program scores, compared within each patient (~30 s)

My part is the RNA. We start from the authors' read counts, 19,578 genes after filtering. We normalize for sequencing depth, then average marker genes into scores, for example CD8 T cells, interferon-gamma, and a tumor score built from melanoma identity genes. Then we compare each patient's two biopsies, so every patient is their own control. That matters because melanomas differ so much from person to person that a between-patient comparison would mostly measure those differences. 49 patients have a pre-treatment biopsy, and 42 have both.

### Slide 5: Before treatment, immune activity only weakly predicts who will respond (~50 s)

First question: does immune activity before treatment predict response? Each dot is one patient's T-cell-inflamed score, a published 18-gene measure of how immune-active a tumor is. Responders sit a little higher, but the difference from progressors, 0.53, has a confidence interval that crosses zero. The one clear pattern is on the right: responders who had already failed ipilimumab all start with 'hot' tumors. That matches the paper, which also called baseline expression only marginally predictive. In the patients who had never had ipilimumab, responders and progressors overlap almost completely. So a single snapshot before treatment is not enough; to understand response, we need to look at what changes during treatment.

### Slide 6: On treatment, responders' genes shift more, and far more consistently, than progressors' (~50 s)

Next, the two biopsies. Each bar is a patient: the share of genes that went up, in orange, or down, in blue, at least two-fold. Responders change more, a median of 8.6% of genes against 4.5% in progressors. But the bigger difference is consistency. When we ask which genes shift together across a group, we find 3,582 in responders, 28 in stable disease and 0 in progressors. For scale, two sequencing runs of the same biopsy already differ by about 3 percent of genes, so the raw amount of change is modest. Responders move in a shared direction; progressors, when they change, change in their own way.

### Slide 7: Tumor loss, more than immune gain, marks a response (~50 s)

So which genes are changing? Each row is a patient, sorted by the change in a tumor score made from five melanoma identity genes, like MLANA and tyrosinase. Responders, in red, cluster at the top with the biggest drops: a median of −1.07 against −0.01 in progressors, p about 0.0002. On the right, the CD8 T-cell score rises in responders, but it also rises in stable disease, and the responder-versus-progressor test is only p about 0.07. It isn't universal: patients 94 and 101 respond with almost no change in the tumor score. But overall, the clearest RNA sign of response is the tumor disappearing from the biopsy.

### Slide 8: Immune pathways rise in every group; only responders shut down proliferation (~45 s)

Pathways show the same thing across thousands of genes. Orange means up on treatment, blue down, grey not significant. The top block is immune pathways: interferon-gamma rises in every group, even progressors at +2.2. The bottom block is proliferation: it collapses only in responders, E2F targets at −3.5, while progressors' tumors keep dividing. This recreates the paper's result, immune up and cell cycle down in responders, and adding the other groups shows that immune activation alone is not enough.

### Slide 9: The on-treatment CD8 level is the immune signal that survives the tumor-dilution check (~45 s)

One catch with bulk RNA: it is a mixture. When tumor cells die, immune cells become a bigger share of what's left, even if no new T cells arrive. So we adjusted for the tumor score. The responder CD8 rise drops from +0.67 to +0.21, so it can't be separated from dilution. But the on-treatment CD8 level holds, +1.30 before and +1.21 after. So by week four, responders' tumors are gone and what remains is full of CD8 T cells. That's the RNA story, and it leads into how RNA compares with the TCR data.

## Shared slides: introduction, modalities, comparison (not in the 5-minute count)

### Slide 2: Riaz et al. biopsied melanomas before and during nivolumab to see what changes when PD-1 blockade works (~40 s)

Nivolumab is an antibody that blocks PD-1, the brake on T cells. In the CheckMate 038 trial, Riaz and colleagues took a biopsy from patients with advanced melanoma before the first dose and again about four weeks in, and sequenced the tumor DNA, the T-cell receptors and the RNA. Half the cohort had never had ipilimumab; the other half had already progressed on it. Their headline: in responders the tumor loses mutations, T-cell clones expand, and immune genes switch on. We re-analyzed their public data to test each of those claims.

### Slide 3: We follow the biology in order: the targets, the T cells that respond, then what the tissue is doing (~30 s)

We present three data types in this order. Kwaku starts with whole-exome sequencing: which mutations, and so which targets, the tumor has, and whether they disappear. Krithika follows with T-cell receptor sequencing: which T-cell clones respond and expand. I finish with RNA, which reads out the whole biopsy at once, both the immune activity and the tumor that is left, so it ties the other two together and leads into our comparison.

### Slide 10: RNA and TCR agree on how many T cells are present, but only TCR-seq sees which clones expand (~35 s)

To compare modalities: RNA and TCR sequencing agree almost perfectly on how many T cells are in a biopsy, rho 0.86. But Krithika's count of expanding clones barely tracks any RNA score, for example rho 0.27 with CD8 change, and no expansion or clonality measure survives multiple-testing correction against any RNA score. So RNA measures the size of the immune response, and TCR-seq measures its clonal makeup: complementary, not redundant.
