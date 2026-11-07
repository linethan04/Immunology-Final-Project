# Presenter script: RNA-seq results

Matches `RNA_11slides_Ethan.pptx` (also in its speaker notes). Total about 7 minutes. Numbers come from `rna/results/`.

## Slide 1: Does gene expression change when nivolumab works? (~20 s)

Hi, I'm Ethan, and I led the bulk RNA sequencing part of our project. Kwaku asked whether tumors lose their mutations when nivolumab works. I ask the matching question for gene expression: what changes in the biopsy, and is that change specific to patients who respond?

## Slide 2: Nivolumab releases the brake; RNA shows what follows (~35 s)

Kwaku covered how mutations create targets. RNA picks up the story from there. T cells recognize the tumor, but PD-1 acts as a brake. Nivolumab blocks PD-1. If that works, active T cells should switch on interferon-gamma and killing genes, and cancer cells should die, so melanoma genes should fade from the biopsy. DNA tells us which targets exist. RNA tells us what the cells in the biopsy are doing. So RNA can measure both the immune response and its effect on the tumor.

## Slide 3: The study sampled each tumor twice (~30 s)

This is the same study design Kwaku showed: one biopsy before the first dose, a second about four weeks in, and scans that later sort patients into responders, stable disease and progressive disease. For RNA, 49 patients with a response call have a pre-treatment biopsy and 42 have both biopsies. Those 42 are my main comparison. Having both biopsies lets every patient act as their own control, which matters because tumors differ so much between people.

## Slide 4: RNA-seq turns each biopsy into gene and program scores (~40 s)

Here is the workflow. We start from the authors' table of read counts: how many RNA reads came from each gene in each biopsy. After removing rarely detected genes, 19,578 remain. Second, we normalize, because some samples were simply sequenced more deeply. Third, we average marker genes into scores for programs such as CD8 T cells, interferon-gamma and melanoma identity. Fourth, we compare the two biopsies, using DESeq2 for single genes and pathway analysis for groups of genes. Take patient 30, a responder: 39 percent of genes changed at least two-fold, the CD8 score rose by 2.5, and the melanoma score fell by 1.5.

## Slide 5: Baseline immune activity is a weak predictor (~40 s)

First: does immune activity before treatment predict response? This is the RNA version of Kwaku's mutation-burden question. Each dot is a patient's T-cell-inflamed score, a published 18-gene measure of how immune-active a tumor is. Medians are 0.19 for responders, −0.47 for stable disease and −0.34 for progressors. The difference between responders and progressors is 0.53, with a 95 percent bootstrap interval from −0.73 to 1.33, so it crosses zero. One cohort detail matches the paper: the four responders who had already progressed on ipilimumab all sit at or above average, the 'hot tumor' pattern the authors describe, but responders without prior ipilimumab do not. The authors also called baseline expression only marginally predictive.

## Slide 6: Each gene is tracked between the two biopsies (~45 s)

Next, the two biopsies. Each bar is one patient, and its height is the share of genes that changed at least two-fold: orange went up, blue went down. Responders change most, with a median of 8.6%, against 4.7% in stable disease and 4.5% in progressors. That gap is modest, and repeat libraries of the same biopsy already differ by about 3 percent. Some progressors, like patients 78 and 103, change a lot. The bigger difference is consistency. When DESeq2 asks which genes shift together across a group, it finds 3,582 in responders, 28 in stable disease and 0 in progressors. A formal test finds 1,721 genes whose change differs between responders and progressors. Responders move in a shared direction; progressors move, if at all, in their own directions.

## Slide 7: Responders lose tumor signal; immune gains are less specific (~45 s)

So which genes? Here each row is a patient, sorted by the change in a score built from five melanoma identity genes. Responders pile up at the bottom left: median change −1.07 against −0.01 in progressors. The difference is −1.05, with an interval of −1.77 to −0.23 and a rank-sum p of about 0.0002. The top single genes in responders tell the same story: SOX10, MLANA, PRAME and PAX3 all fall. On the right is the CD8 T-cell score. It rises in responders, median +0.47, but also in stable disease, +0.31, and the responder-versus-progressor test gives p of only about 0.07. Immune activation is clearly happening, but tumor loss is what marks a response. It isn't universal: patients 94 and 101 respond with little change in the tumor score.

## Slide 8: Immune pathways rise broadly; proliferation falls in responders (~40 s)

Pathways make the same point across thousands of genes. Orange means higher, blue lower, and grey is not significant. In the top block, immune pathways rise on treatment in every group, and stable disease is strongest, with interferon-gamma at +3.5. In the bottom block, proliferation falls hard in responders, with E2F targets at −3.5, while most of those programs rise in progressors, at +1.3. That matches the paper, which reported immune genes up and cell-cycle and melanin genes down in responders. The left column shows immune pathways already lean higher in responders before treatment, again as the paper describes. Enrichment p-values run small because genes in a pathway move together, so I read the direction and pattern here more than the exact significance.

## Slide 9: RNA agrees with WES and TCR sequencing (~40 s)

Do the RNA results agree with data measured a completely different way? On the left, each point is a patient with both RNA and exome data. The more a patient's mutations contract in Kwaku's analysis, the more the RNA tumor score falls: rho of −0.69 across 25 patients, and still −0.50 with responders removed. On the right, the RNA CD8 score tracks the T-cell fraction measured directly by TCR sequencing in Krithika's data: rho of 0.88. That validates the score. Both left-panel measures partly reflect how much tumor is in the biopsy, so they aren't fully independent. Baseline immune score and mutation burden barely correlate, rho of 0.14, so RNA and WES carry different information for a combined model.

## Slide 10: What could make genes look changed? (~60 s)

Before the takeaways, what could fool this? First and most important is dilution. Bulk RNA is a mixture, so when tumor cells die, every other cell type becomes a bigger share. When I adjust the responder CD8 rise for tumor loss, it drops from +0.67 to +0.21 and the interval crosses zero, so that rise can't be separated from dilution. The on-treatment CD8 level holds up, at +1.30 before and +1.21 after adjustment. Second, response is defined by tumor shrinkage, so falling melanoma genes partly restate the outcome rather than predict it. Third, I couldn't reproduce the paper's gene counts. With their designs and a current DESeq2, I get 397 versus their 189 before treatment and 157 versus 475 on treatment. Their responder comparison, 2,670 genes, has no public script, and grouping stable disease with responders gives me only 7; responders versus progressors gives 3,479. The biological direction agrees. The paper also removed patient 3 as an outlier; I kept that patient, and the main results hold either way. Finally, groups are small, 9 responders and 18 progressors, at one time point.

## Slide 11: Takeaways and next steps (~35 s)

Three takeaways. First, what tracks response on RNA is tumor loss: melanoma and proliferation programs fall in responders, and that agrees with Kwaku's mutation contraction. Second, immune activation alone isn't enough. Immune pathways rise in every group, including patients whose tumors don't shrink, and the baseline score doesn't separate responders. The on-treatment CD8 level does, with a difference of +1.66, and it survives the tumor-content adjustment. Next, I'd use a proper deconvolution method to estimate cell types directly, repeat this on the 26 patients our team shares, and combine RNA, WES and TCR in a cross-validated model. Thank you.
