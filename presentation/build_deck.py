"""Build presentation/RNA_seq_nivolumab.pptx from figures/rna/*.png.
Pure standard library (writes the OOXML parts directly). Run from project root:
    python3 presentation/build_deck.py
"""
import struct, zipfile
from pathlib import Path
from xml.sax.saxutils import escape

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "figures" / "rna"
OUT = ROOT / "presentation" / "RNA_seq_nivolumab.pptx"

EMU = 914400
W, H = 13.333, 7.5
# palette: deep ink navy dominant, responder blue / progressor red / SD gold as data accents
INK, PAPER, MUTED, CARD = "16213E", "FFFFFF", "5B6475", "EEF2F7"
BLUE, RED, GOLD, TEAL, ICE = "2471A3", "C0392B", "E1A23B", "16A085", "BFD7EA"
HEAD, BODY = "Cambria", "Calibri"
e = lambda inch: int(round(inch * EMU))


def png_size(p):
    with open(p, "rb") as f:
        f.read(16)
        return struct.unpack(">II", f.read(8))


class Slide:
    def __init__(self, bg=PAPER):
        self.bg, self.shapes, self.images, self.nid = bg, [], [], 2

    def _id(self):
        self.nid += 1
        return self.nid

    def rect(self, x, y, w, h, fill, geom="rect", line=None, name="Shape"):
        ln = f'<a:ln w="12700"><a:solidFill><a:srgbClr val="{line}"/></a:solidFill></a:ln>' if line else "<a:ln><a:noFill/></a:ln>"
        adj = '<a:avLst><a:gd name="adj" fmla="val 8000"/></a:avLst>' if geom == "roundRect" else "<a:avLst/>"
        self.shapes.append(
            f'<p:sp><p:nvSpPr><p:cNvPr id="{self._id()}" name="{name}"/><p:cNvSpPr/><p:nvPr/></p:nvSpPr>'
            f'<p:spPr><a:xfrm><a:off x="{e(x)}" y="{e(y)}"/><a:ext cx="{e(w)}" cy="{e(h)}"/></a:xfrm>'
            f'<a:prstGeom prst="{geom}">{adj}</a:prstGeom><a:solidFill><a:srgbClr val="{fill}"/></a:solidFill>{ln}</p:spPr>'
            f'<p:txBody><a:bodyPr/><a:lstStyle/><a:p><a:endParaRPr lang="en-US"/></a:p></p:txBody></p:sp>')

    def text(self, x, y, w, h, paras, anchor="t", name="Text"):
        """paras: list of dicts {t | runs, size, color, bold, italic, font, bullet, align, after}"""
        ps = []
        for p in paras:
            runs = p.get("runs") or [{"t": p["t"]}]
            ppr = f' algn="{p.get("align", "l")}"'
            inner = f'<a:spcAft><a:spcPts val="{int(p.get("after", 6) * 100)}"/></a:spcAft>'
            if p.get("bullet"):
                ppr += ' marL="285750" indent="-285750"'
                inner += f'<a:buClr><a:srgbClr val="{p.get("bucolor", BLUE)}"/></a:buClr><a:buFont typeface="Arial"/><a:buChar char="&#8226;"/>'
            else:
                inner += "<a:buNone/>"
            rs = ""
            for r in runs:
                sz = int(r.get("size", p.get("size", 16)) * 100)
                b = ' b="1"' if r.get("bold", p.get("bold")) else ""
                i = ' i="1"' if r.get("italic", p.get("italic")) else ""
                col = r.get("color", p.get("color", INK))
                font = r.get("font", p.get("font", BODY))
                rs += (f'<a:r><a:rPr lang="en-US" sz="{sz}"{b}{i} dirty="0"><a:solidFill><a:srgbClr val="{col}"/></a:solidFill>'
                       f'<a:latin typeface="{font}"/><a:cs typeface="{font}"/></a:rPr><a:t>{escape(r["t"])}</a:t></a:r>')
            ps.append(f"<a:p><a:pPr{ppr}>{inner}</a:pPr>{rs}</a:p>")
        self.shapes.append(
            f'<p:sp><p:nvSpPr><p:cNvPr id="{self._id()}" name="{name}"/><p:cNvSpPr txBox="1"/><p:nvPr/></p:nvSpPr>'
            f'<p:spPr><a:xfrm><a:off x="{e(x)}" y="{e(y)}"/><a:ext cx="{e(w)}" cy="{e(h)}"/></a:xfrm>'
            f'<a:prstGeom prst="rect"><a:avLst/></a:prstGeom><a:noFill/></p:spPr>'
            f'<p:txBody><a:bodyPr wrap="square" lIns="0" tIns="0" rIns="0" bIns="0" anchor="{anchor}"><a:noAutofit/></a:bodyPr>'
            f'<a:lstStyle/>{"".join(ps)}</p:txBody></p:sp>')

    def image(self, path, x, y, w, h, align="c", alt=""):
        """Fit image inside box (x, y, w, h) preserving aspect; returns placed box."""
        pw, ph = png_size(path)
        s = min(w / pw, h / ph)
        iw, ih = pw * s, ph * s
        ix = x + (w - iw) / 2 if align == "c" else x
        iy = y + (h - ih) / 2
        rid = f"rId{len(self.images) + 2}"
        self.images.append((rid, Path(path)))
        self.shapes.append(
            f'<p:pic><p:nvPicPr><p:cNvPr id="{self._id()}" name="{escape(Path(path).stem)}" descr="{escape(alt)}"/>'
            f'<p:cNvPicPr><a:picLocks noChangeAspect="1"/></p:cNvPicPr><p:nvPr/></p:nvPicPr>'
            f'<p:blipFill><a:blip r:embed="{rid}"/><a:stretch><a:fillRect/></a:stretch></p:blipFill>'
            f'<p:spPr><a:xfrm><a:off x="{e(ix)}" y="{e(iy)}"/><a:ext cx="{e(iw)}" cy="{e(ih)}"/></a:xfrm>'
            f'<a:prstGeom prst="rect"><a:avLst/></a:prstGeom></p:spPr></p:pic>')
        return ix, iy, iw, ih

    def xml(self):
        return (f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                f'<p:sld {NS}><p:cSld><p:bg><p:bgPr><a:solidFill><a:srgbClr val="{self.bg}"/></a:solidFill>'
                f'<a:effectLst/></p:bgPr></p:bg><p:spTree><p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr>'
                f'<p:grpSpPr/>{"".join(self.shapes)}</p:spTree></p:cSld><p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:sld>')


NS = ('xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" '
      'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" '
      'xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main"')

# ---------------------------------------------------------------- slide helpers
slides = []
n_content = [0]


def content(title, kicker=None):
    s = Slide()
    s.text(0.6, 0.38, 12.1, 0.8, [{"t": title, "size": 30, "bold": True, "font": HEAD, "color": INK}],
           anchor="b", name="Title")
    if kicker:
        s.text(0.6, 1.22, 12.1, 0.4, [{"t": kicker, "size": 15, "italic": True, "color": MUTED}], name="Kicker")
    n_content[0] += 1
    s.text(12.1, 7.0, 0.8, 0.3, [{"t": str(len(slides) + 1), "size": 10, "color": MUTED, "align": "r"}], name="Page")
    slides.append(s)
    return s


def divider(num, title, sub):
    s = Slide(INK)
    s.text(0.9, 2.3, 3, 1.2, [{"t": num, "size": 72, "bold": True, "font": HEAD, "color": GOLD}], name="Number")
    s.text(0.9, 3.6, 11, 1.0, [{"t": title, "size": 40, "bold": True, "font": HEAD, "color": PAPER}], name="Title")
    s.text(0.9, 4.6, 11, 0.8, [{"t": sub, "size": 18, "color": ICE}], name="Subtitle")
    slides.append(s)


def takeaways(s, x, y, w, h, head, items, color=BLUE):
    s.rect(x, y, w, h, CARD, "roundRect", name="Takeaway card")
    paras = [{"t": head, "size": 17, "bold": True, "font": HEAD, "color": color, "after": 10}]
    for it in items:
        if isinstance(it, list):
            paras.append({"runs": it, "size": 15, "bullet": True, "bucolor": color, "after": 8})
        else:
            paras.append({"t": it, "size": 15, "bullet": True, "bucolor": color, "after": 8})
    s.text(x + 0.3, y + 0.3, w - 0.6, h - 0.5, paras, name="Takeaways")


def fig_slide(title, kicker, fig, head, items, fig_w=7.9, color=BLUE, alt=""):
    s = content(title, kicker)
    top = 1.85
    s.image(FIG / fig, 0.6, top, fig_w, 7.0 - top - 0.15, alt=alt or title)
    x = 0.6 + fig_w + 0.35
    takeaways(s, x, top, W - x - 0.6, 7.0 - top - 0.15, head, items, color)
    return s


B = lambda t: {"t": t, "bold": True}

# ---------------------------------------------------------------- 1. title
s = Slide(INK)
s.text(0.9, 1.0, 11.5, 0.5, [{"t": "IMMUNOLOGY FINAL PROJECT  ·  RNA-SEQ ARM", "size": 14, "bold": True, "color": GOLD}])
s.text(0.9, 1.7, 11.5, 2.2, [{"t": "What changes in the tumor when nivolumab works?", "size": 44, "bold": True,
                               "font": HEAD, "color": PAPER}])
s.text(0.9, 4.1, 11.5, 1.0, [{"t": "Bulk RNA-seq re-analysis of Riaz et al., Cell 2017 (CheckMate 038): "
                               "pre- vs on-treatment biopsies in advanced melanoma", "size": 20, "color": ICE}])
s.text(0.9, 5.8, 11.5, 0.9, [{"t": "Ethan (RNA-seq)  ·  Krithika (TCR-seq)  ·  Kwaku (WES + integration)",
                               "size": 16, "color": PAPER}])
for i, c in enumerate([RED, GOLD, BLUE]):
    s.rect(0.9 + i * 0.45, 5.35, 0.3, 0.3, c, "ellipse", name="Response dot")
slides.append(s)

# ---------------------------------------------------------------- 2. background
s = content("Why biopsy patients during treatment?",
            "Baseline biomarkers for anti-PD-1 are weak; the on-treatment tumor may show the response as it happens")
s.text(0.6, 1.95, 6.4, 4.9, [
    {"runs": [B("Nivolumab"), {"t": " blocks PD-1 on T cells, releasing a brake on anti-tumor T-cell killing."}], "bullet": True, "after": 12},
    {"runs": [B("CheckMate 038"), {"t": " took biopsies before treatment (Pre) and at ~week 4 (On), with RNA-seq, TCR-seq and exome sequencing."}], "bullet": True, "after": 12},
    {"runs": [B("Response"), {"t": " (RECIST): PRCR = partial/complete response, SD = stable disease, PD = progression."}], "bullet": True, "after": 12},
    {"runs": [B("Our question: "), {"t": "do responders show immune activation and tumor loss in the same patients, and can RNA help predict response?"}], "bullet": True},
], name="Background")
stats = [("103", "RNA-seq samples\nafter QC filter"), ("42", "patients with paired\nPre + On biopsies"),
         ("25", "tri-modal patients\nwith a response call"), ("19,578", "genes kept\n(≥10 counts in ≥5 samples)")]
for i, (num, lab) in enumerate(stats):
    cx, cy = 7.5 + (i % 2) * 2.7, 1.95 + (i // 2) * 2.45
    s.rect(cx, cy, 2.45, 2.2, CARD, "roundRect", name="Stat card")
    s.text(cx + 0.2, cy + 0.25, 2.05, 0.9, [{"t": num, "size": 40, "bold": True, "font": HEAD, "color": BLUE}])
    s.text(cx + 0.2, cy + 1.2, 2.05, 0.9, [{"t": lab, "size": 13, "color": MUTED}])

# ---------------------------------------------------------------- 3. aims
s = content("Three project aims, and what RNA-seq contributes", "RNA-seq supplies the immune and tumor-state readout for every patient")
aims = [("1", "Characterize each modality", "Which genes and pathways differ between responders and progressors, before and during treatment?", RED),
        ("2", "Integrate per patient", "Do immune activation (RNA), clonal expansion (TCR) and neoantigen loss (WES) happen in the same patients?", GOLD),
        ("3", "Predict response", "Which features best separate responders, with honest small-n uncertainty?", BLUE)]
rna = ["DE genes, pathways, signature scores", "per-patient Pre / On / delta table", "AUCs of every score, ready for the model"]
for i, (n, t, d, c) in enumerate(aims):
    x = 0.6 + i * 4.1
    s.rect(x, 2.0, 3.8, 4.7, CARD, "roundRect", name="Aim card")
    s.rect(x + 0.35, 2.35, 0.8, 0.8, c, "ellipse", name="Aim badge")
    s.text(x + 0.35, 2.35, 0.8, 0.8, [{"t": n, "size": 28, "bold": True, "font": HEAD, "color": PAPER, "align": "ctr"}], anchor="ctr")
    s.text(x + 0.35, 3.4, 3.1, 0.9, [{"t": t, "size": 20, "bold": True, "font": HEAD}])
    s.text(x + 0.35, 4.3, 3.1, 1.4, [{"t": d, "size": 15, "color": MUTED}])
    s.text(x + 0.35, 5.75, 3.1, 0.8, [{"runs": [{"t": "RNA: ", "bold": True, "color": c}, {"t": rna[i]}], "size": 14}])

# ---------------------------------------------------------------- 4. divider
divider("01", "Workflow and code", "Five R scripts, run in order with scripts/rna_00_run_all.R")

# ---------------------------------------------------------------- 5. pipeline
s = content("The RNA-seq pipeline", "Each script reads the previous step's output; tables go to output/rna/, figures to figures/rna/")
steps = [("rna_01", "Preprocess + QC", "Counts → filter → VST\nLibrary size, PCA"),
         ("rna_02", "Differential\nexpression", "DESeq2: Pre PRCR vs PD;\npaired On vs Pre;\nper-group change"),
         ("rna_03", "Pathways", "fgsea on 50 MSigDB\nHallmark sets"),
         ("rna_04", "Immune\nsignatures", "12 marker scores per\nsample → per-patient\nPre / On / delta"),
         ("rna_05", "Clinical +\ncross-modal", "AUC, survival,\nRNA vs TCR / WES")]
bw, gap = 2.18, 0.3
for i, (code, t, d) in enumerate(steps):
    x = 0.6 + i * (bw + gap)
    s.rect(x, 2.1, bw, 1.6, INK, "roundRect", name="Step")
    s.text(x + 0.15, 2.2, bw - 0.3, 0.4, [{"t": code, "size": 12, "bold": True, "color": GOLD, "align": "ctr"}])
    s.text(x + 0.15, 2.6, bw - 0.3, 1.0, [{"t": t, "size": 17, "bold": True, "font": HEAD, "color": PAPER, "align": "ctr"}], anchor="ctr")
    s.text(x + 0.1, 3.95, bw - 0.2, 1.4, [{"t": d, "size": 13, "color": MUTED, "align": "ctr"}])
    if i < 4:
        s.rect(x + bw + 0.04, 2.75, gap - 0.08, 0.3, BLUE, "rightArrow", name="Arrow")
s.rect(0.6, 5.55, 12.1, 1.25, CARD, "roundRect", name="Output card")
s.text(0.9, 5.7, 11.5, 1.0, [
    {"runs": [B("Hand-off to the team: "), {"t": "output/rna/rna_scores_per_patient.csv"}], "size": 15, "after": 4},
    {"t": "One row per patient: Pre, On and On − Pre values for every signature, plus response, cohort and a tri-modal flag. "
          "This joins directly to Kwaku's integration table.", "size": 14, "color": MUTED}])

# ---------------------------------------------------------------- 6. stats design
s = content("Key statistical design choices", "Each choice protects against a specific failure mode of small, paired data")
cards = [("Paired design", "~ patient + time", "Each patient is their own control, which removes between-patient baseline differences."),
         ("Response-specific change", "~ response + response:patient + response:time",
          "Estimates the On vs Pre change separately in PRCR, SD and PD, then tests PRCR − PD (an interaction)."),
         ("Signature scores", "mean z-score of VST genes",
          "Robust, interpretable per-sample scores for cell types and programs (CD8, cytotoxic, IFN-γ, melanocytic)."),
         ("Effect sizes, not just p", "Hodges–Lehmann shift + 95% CI; AUC (DeLong CI)",
          "With n ≈ 9–24 per group, CIs show how much is uncertain; nothing is tuned on the outcome.")]
for i, (t, code, d) in enumerate(cards):
    x, y = 0.6 + (i % 2) * 6.15, 1.95 + (i // 2) * 2.5
    s.rect(x, y, 5.95, 2.25, CARD, "roundRect", name="Design card")
    s.text(x + 0.3, y + 0.22, 5.4, 0.45, [{"t": t, "size": 18, "bold": True, "font": HEAD, "color": BLUE}])
    s.text(x + 0.3, y + 0.72, 5.4, 0.4, [{"t": code, "size": 13, "font": "Courier New", "color": INK}])
    s.text(x + 0.3, y + 1.18, 5.4, 1.0, [{"t": d, "size": 14, "color": MUTED}])

# ---------------------------------------------------------------- 7. divider
divider("02", "Results", "QC → differential expression → pathways → signatures → outcome")

# ---------------------------------------------------------------- 8. QC
s = content("Data quality is consistent across groups", "No batch structure by response or timepoint, so the downstream comparisons are fair")
s.image(FIG / "01_qc_library_size.png", 0.6, 1.85, 6.6, 2.75, alt="Library size per sample")
s.image(FIG / "02_qc_pca.png", 7.4, 1.85, 5.3, 4.95, alt="PCA of VST expression")
takeaways(s, 0.6, 4.8, 6.6, 2.05, "What to notice", [
    "Library sizes are comparable in PD, SD and PRCR at both timepoints.",
    "The PCA shows no separation by response, so the main variation is tumor/immune composition, not artifacts."], color=TEAL)

# ---------------------------------------------------------------- 9-11 DE
fig_slide("Before treatment, responders look only modestly different",
          "Aim 1: pre-treatment PRCR vs PD (paper exclusions applied)", "03_volcano_pre_PRCR_vs_PD.png",
          "What to notice", ["223 genes at FDR < 0.1, mostly with small shrunken effects.",
                             "The top hits are epidermal/keratinocyte genes (SPRR2A, SERPINB3, LCE), which likely reflect skin in some biopsies rather than tumor biology.",
                             "Immune genes (e.g. CXCL13) lean higher in responders, but there is no clean baseline biomarker, consistent with the paper.",
                             "Motivates looking at the on-treatment biopsy."])
fig_slide("Nivolumab switches on T-cell and interferon genes",
          "Aim 1: paired On vs Pre across all 42 patients with both biopsies", "04_volcano_on_vs_pre.png",
          "What to notice", ["44 genes change at FDR < 0.1.",
                             "The top hits are the PD-1 axis and T-cell effector genes: PDCD1, CXCL9, GZMH, PRF1, IFNG, LAG3, GBP5.",
                             "This is the expected immune response to PD-1 blockade, averaged over everyone."], color=RED)
fig_slide("The on-treatment change is driven by responders",
          "Aim 1: per-gene On/Pre change in PRCR (y) vs PD (x); purple = interaction FDR < 0.1",
          "05_lfc_PRCR_vs_PD_change.png", "What to notice", [
              [B("PRCR: 3,582"), {"t": " genes change. "}, B("PD: 0."), {"t": " SD: 28."}],
              "Immune genes rise mainly in responders (points above the diagonal).",
              "The strongest drops in responders are melanoma lineage genes (SOX10, MLANA, PRAME, PAX3): the tumor is being cleared.",
              "1,721 genes differ in PRCR vs PD change."], fig_w=7.0, color=GOLD)
fig_slide("Top on-treatment genes, sample by sample",
          "Top 40 On vs Pre genes, row z-scores; samples ordered by timepoint then response", "06_heatmap_top_on_vs_pre.png",
          "What to notice", ["The On block is redder, especially for PRCR samples.",
                             "Some PD patients also gain immune genes, so activation alone is not enough.",
                             "The genes form a coherent T-cell/IFN-γ program, not scattered hits."], fig_w=7.0)

# ---------------------------------------------------------------- 13 GSEA
fig_slide("Pathways: interferon up, proliferation down in responders",
          "fgsea on Hallmark sets; colour = direction (NES), size = significance", "07_gsea_hallmark_dotplot.png",
          "What to notice", ["IFN-γ / IFN-α, allograft rejection and inflammatory pathways rise On-treatment (IFN-γ NES 2.7 in PRCR).",
                             "E2F, MYC and G2M (tumor proliferation) fall sharply only in PRCR.",
                             "PD patients show weaker immune gains and no proliferation loss.",
                             "SD patients show strong interferon induction: an intermediate group."], fig_w=7.6, color=RED)

# ---------------------------------------------------------------- 14 trajectories
s = content("Every patient's trajectory, Pre → On",
            "Signature scores per patient (thin lines) with group medians (thick); rows = response group")
s.image(FIG / "09_signature_trajectories.png", 0.6, 1.75, 9.0, 5.15, alt="Signature trajectories")
takeaways(s, 9.85, 1.75, 2.88, 5.15, "What to notice", [
    [B("PRCR: "), {"t": "immune scores rise; tumor score collapses."}],
    [B("SD: "), {"t": "immune scores rise; tumor barely moves."}],
    [B("PD: "), {"t": "flat on every axis."}],
    "Immune activation is necessary but not sufficient."])

# ---------------------------------------------------------------- 15 effect sizes
s = content("On-treatment RNA separates responders; baseline does not",
            "PRCR − PD shift in each score with 95% CI (Hodges–Lehmann); red = p < 0.05")
s.image(FIG / "10_effect_sizes_PRCR_vs_PD.png", 0.6, 1.8, 8.6, 5.05, alt="Effect sizes")
takeaways(s, 9.45, 1.8, 3.28, 5.05, "What to notice", [
    "Pre: every CI crosses 0.",
    "On: 10 of 12 scores higher in PRCR (CD8 shift +1.4, FDR 0.006).",
    "Change: the clearest signal is loss of melanocytic genes (FDR 0.006)."])

# ---------------------------------------------------------------- 16 heatmap per patient
s = content("Patient-level view for integration (Aim 2)",
            "On − Pre change in each score for the 25 tri-modal patients; these values feed the merged table")
s.image(FIG / "11_heatmap_delta_tri_patients.png", 0.6, 1.8, 12.1, 3.9, alt="Delta heatmap")
s.rect(0.6, 5.85, 12.1, 1.0, CARD, "roundRect")
s.text(0.9, 5.98, 11.5, 0.8, [{"runs": [B("Read across a column: "), {"t": "responders such as Pt30 and Pt3 gain every immune program and lose the tumor signature. "
    "Some PD patients (Pt28, Pt31) gain immune signal without tumor loss. These are the cases where TCR and neoantigen data should explain the difference."}], "size": 14}])

# ---------------------------------------------------------------- 17 AUC
fig_slide("Which RNA feature best classifies response? (Aim 3)",
          "Single-feature AUC, PRCR vs PD, with DeLong 95% CI", "13_auc_by_timepoint.png",
          "What to notice", ["Pre-treatment scores: AUC 0.60 to 0.68, with CIs touching 0.5.",
                             "On-treatment CD8: AUC 0.87 (0.75 to 1.00).",
                             "Melanocytic change: AUC 0.96, but partly circular, since tumor shrinkage is the response.",
                             "Next: combine with TCR and WES features in a cross-validated model."], fig_w=7.4, color=TEAL)

# ---------------------------------------------------------------- 18 survival
fig_slide("A higher on-treatment inflamed score trends toward longer survival",
          "Overall survival, median split on the 18-gene T-cell-inflamed score (On biopsy)", "14_km_os_inflamed_on.png",
          "What to notice", ["The curves separate early and stay apart.",
                             "HR per SD = 0.74 (0.49 to 1.10); log-rank p = 0.09, which is not significant with n = 48.",
                             "Caveat: the score is measured at week 4, so patients must survive to the biopsy (landmark bias)."],
          fig_w=7.6, color=BLUE)

# ---------------------------------------------------------------- 19 cross modality
fig_slide("RNA agrees with TCR-seq; mutation burden is a separate axis",
          "Preview of Aim 2: RNA scores vs Krithika's TCR metrics and Kwaku's mutation counts", "15_rna_vs_tcr_wes.png",
          "What to notice", [[B("Validation: "), {"t": "the RNA CD8 score tracks the TCR-seq T-cell fraction (ρ = 0.81 Pre, 0.88 On)."}],
                             "Change in TCR clonality vs change in cytotoxic score: a weak trend (ρ = 0.27).",
                             "Mutation burden does not predict baseline inflammation (ρ = 0.14), so it adds independent information."],
          fig_w=7.0, color=GOLD)

# ---------------------------------------------------------------- 20 confound
fig_slide("Caveat: tumor content confounds immune scores",
          "Bulk RNA mixes tumor and immune cells; more tumor means less room for immune signal",
          "12_melanocytic_vs_inflamed.png",
          "Why it matters", ["Immune and melanocytic scores are weakly anti-correlated across samples (ρ = −0.24).",
                             "Part of the On-treatment immune 'gain' in responders reflects tumor loss.",
                             "Fix: adjust for tumor purity (e.g. ESTIMATE or deconvolution) before integration."],
          fig_w=7.0, color=RED)

# ---------------------------------------------------------------- 21 conclusions
s = Slide(INK)
s.text(0.9, 0.6, 11.5, 0.9, [{"t": "Conclusions and next steps", "size": 36, "bold": True, "font": HEAD, "color": PAPER}])
cols = [("What RNA-seq shows", GOLD, [
            "Baseline expression barely separates responders.",
            "In responders, nivolumab drives an IFN-γ / CD8 / cytotoxic program while tumor and proliferation genes collapse.",
            "PD tumors stay flat; SD tumors activate immunity without clearing the tumor.",
            "On-treatment CD8 score reaches AUC 0.87 and tracks the TCR-seq T-cell fraction."]),
        ("Next steps", ICE, [
            "Purity-adjust scores; proper deconvolution (CIBERSORTx / quanTIseq).",
            "Merge RNA deltas with TCR expansion and neoantigen loss (Aim 2).",
            "Cross-validated multi-modal classifier with nested CV (Aim 3).",
            "Settle the response definitions: Pt23 (NE) and SD as responder or not."])]
for i, (h, c, items) in enumerate(cols):
    x = 0.9 + i * 6.0
    paras = [{"t": h, "size": 20, "bold": True, "font": HEAD, "color": c, "after": 12}] + \
            [{"t": t, "size": 17, "color": PAPER, "bullet": True, "bucolor": c, "after": 14} for t in items]
    s.text(x, 1.8, 5.5, 5.0, paras)
slides.append(s)

# ---------------------------------------------------------------- package
THEME = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<a:theme xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" name="RNA deck"><a:themeElements>
<a:clrScheme name="RNA deck"><a:dk1><a:srgbClr val="{INK}"/></a:dk1><a:lt1><a:srgbClr val="FFFFFF"/></a:lt1>
<a:dk2><a:srgbClr val="{MUTED}"/></a:dk2><a:lt2><a:srgbClr val="{CARD}"/></a:lt2>
<a:accent1><a:srgbClr val="{BLUE}"/></a:accent1><a:accent2><a:srgbClr val="{RED}"/></a:accent2>
<a:accent3><a:srgbClr val="{GOLD}"/></a:accent3><a:accent4><a:srgbClr val="{TEAL}"/></a:accent4>
<a:accent5><a:srgbClr val="8E44AD"/></a:accent5><a:accent6><a:srgbClr val="{ICE}"/></a:accent6>
<a:hlink><a:srgbClr val="{BLUE}"/></a:hlink><a:folHlink><a:srgbClr val="8E44AD"/></a:folHlink></a:clrScheme>
<a:fontScheme name="RNA deck"><a:majorFont><a:latin typeface="{HEAD}"/><a:ea typeface=""/><a:cs typeface=""/></a:majorFont>
<a:minorFont><a:latin typeface="{BODY}"/><a:ea typeface=""/><a:cs typeface=""/></a:minorFont></a:fontScheme>
<a:fmtScheme name="RNA deck"><a:fillStyleLst><a:solidFill><a:schemeClr val="phClr"/></a:solidFill><a:solidFill><a:schemeClr val="phClr"/></a:solidFill><a:solidFill><a:schemeClr val="phClr"/></a:solidFill></a:fillStyleLst>
<a:lnStyleLst><a:ln w="6350"><a:solidFill><a:schemeClr val="phClr"/></a:solidFill></a:ln><a:ln w="12700"><a:solidFill><a:schemeClr val="phClr"/></a:solidFill></a:ln><a:ln w="19050"><a:solidFill><a:schemeClr val="phClr"/></a:solidFill></a:ln></a:lnStyleLst>
<a:effectStyleLst><a:effectStyle><a:effectLst/></a:effectStyle><a:effectStyle><a:effectLst/></a:effectStyle><a:effectStyle><a:effectLst/></a:effectStyle></a:effectStyleLst>
<a:bgFillStyleLst><a:solidFill><a:schemeClr val="phClr"/></a:solidFill><a:solidFill><a:schemeClr val="phClr"/></a:solidFill><a:solidFill><a:schemeClr val="phClr"/></a:solidFill></a:bgFillStyleLst>
</a:fmtScheme></a:themeElements><a:objectDefaults/><a:extraClrSchemeLst/></a:theme>'''

EMPTY_TREE = '<p:spTree><p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr><p:grpSpPr/></p:spTree>'
MASTER = (f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><p:sldMaster {NS}><p:cSld><p:bg><p:bgRef idx="1001"><a:schemeClr val="bg1"/></p:bgRef></p:bg>{EMPTY_TREE}</p:cSld>'
          '<p:clrMap bg1="lt1" tx1="dk1" bg2="lt2" tx2="dk2" accent1="accent1" accent2="accent2" accent3="accent3" accent4="accent4" accent5="accent5" accent6="accent6" hlink="hlink" folHlink="folHlink"/>'
          '<p:sldLayoutIdLst><p:sldLayoutId id="2147483649" r:id="rId1"/></p:sldLayoutIdLst></p:sldMaster>')
LAYOUT = (f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><p:sldLayout {NS} type="blank" preserve="1"><p:cSld name="Blank">{EMPTY_TREE}</p:cSld>'
          '<p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:sldLayout>')
REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PREL = "http://schemas.openxmlformats.org/package/2006/relationships"


def rels(items):
    return (f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="{PREL}">' +
            "".join(f'<Relationship Id="{i}" Type="{REL}/{t}" Target="{tg}"/>' for i, t, tg in items) + "</Relationships>")


n = len(slides)
pres = (f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><p:presentation {NS} saveSubsetFonts="1">'
        '<p:sldMasterIdLst><p:sldMasterId id="2147483648" r:id="rId1"/></p:sldMasterIdLst><p:sldIdLst>' +
        "".join(f'<p:sldId id="{256 + i}" r:id="rId{i + 3}"/>' for i in range(n)) +
        f'</p:sldIdLst><p:sldSz cx="{e(W)}" cy="{e(H)}"/><p:notesSz cx="6858000" cy="9144000"/></p:presentation>')
ctypes = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
          '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
          '<Default Extension="xml" ContentType="application/xml"/><Default Extension="png" ContentType="image/png"/>'
          '<Override PartName="/ppt/presentation.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.presentation.main+xml"/>'
          '<Override PartName="/ppt/slideMasters/slideMaster1.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slideMaster+xml"/>'
          '<Override PartName="/ppt/slideLayouts/slideLayout1.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slideLayout+xml"/>'
          '<Override PartName="/ppt/theme/theme1.xml" ContentType="application/vnd.openxmlformats-officedocument.theme+xml"/>'
          '<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>'
          '<Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>' +
          "".join(f'<Override PartName="/ppt/slides/slide{i + 1}.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slide+xml"/>' for i in range(n)) +
          "</Types>")
core = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" '
        'xmlns:dc="http://purl.org/dc/elements/1.1/"><dc:title>RNA-seq analysis: nivolumab in melanoma</dc:title><dc:creator>Ethan</dc:creator></cp:coreProperties>')
app = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties">'
       f'<Application>build_deck.py</Application><Slides>{n}</Slides></Properties>')

OUT.parent.mkdir(exist_ok=True)
media = {}
with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED) as z:
    z.writestr("[Content_Types].xml", ctypes)
    z.writestr("_rels/.rels", rels([("rId1", "officeDocument", "ppt/presentation.xml")]).replace(
        "</Relationships>",
        '<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>'
        f'<Relationship Id="rId3" Type="{REL}/extended-properties" Target="docProps/app.xml"/></Relationships>'))
    z.writestr("docProps/core.xml", core)
    z.writestr("docProps/app.xml", app)
    z.writestr("ppt/presentation.xml", pres)
    z.writestr("ppt/_rels/presentation.xml.rels", rels(
        [("rId1", "slideMaster", "slideMasters/slideMaster1.xml"), ("rId2", "theme", "theme/theme1.xml")] +
        [(f"rId{i + 3}", "slide", f"slides/slide{i + 1}.xml") for i in range(n)]))
    z.writestr("ppt/theme/theme1.xml", THEME)
    z.writestr("ppt/slideMasters/slideMaster1.xml", MASTER)
    z.writestr("ppt/slideMasters/_rels/slideMaster1.xml.rels",
               rels([("rId1", "slideLayout", "../slideLayouts/slideLayout1.xml"), ("rId2", "theme", "../theme/theme1.xml")]))
    z.writestr("ppt/slideLayouts/slideLayout1.xml", LAYOUT)
    z.writestr("ppt/slideLayouts/_rels/slideLayout1.xml.rels", rels([("rId1", "slideMaster", "../slideMasters/slideMaster1.xml")]))
    for i, s in enumerate(slides, 1):
        z.writestr(f"ppt/slides/slide{i}.xml", s.xml())
        items = [("rId1", "slideLayout", "../slideLayouts/slideLayout1.xml")]
        for rid, p in s.images:
            if p not in media:
                media[p] = f"image{len(media) + 1}.png"
                z.write(p, f"ppt/media/{media[p]}")
            items.append((rid, "image", f"../media/{media[p]}"))
        z.writestr(f"ppt/slides/_rels/slide{i}.xml.rels", rels(items))
print(f"wrote {OUT.relative_to(ROOT)} ({n} slides, {len(media)} images)")
