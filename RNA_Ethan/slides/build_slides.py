"""Build the RNA slides and presenter script from rna/results, on the theme of Kwaku's WES deck.
Standard library only. Run from the RNA_Ethan folder: python3 slides/build_slides.py
Writes slides/RNA_11slides_Ethan.pptx and slides/presentation_script.md.
"""
import csv, re, struct, zipfile
from pathlib import Path
from xml.sax.saxutils import escape

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "rna" / "results"
TEMPLATE = ROOT.parent / "WES_Kwaku" / "slides" / "WES_10slides_Kwaku.pptx"
OUT = ROOT / "slides" / "RNA_10slides_Ethan.pptx"
SCRIPT = ROOT / "slides" / "presentation_script.md"
AUTHOR = "Ethan"


# ---------- numbers: everything quoted comes from the result tables ----------
def table(name):
    return list(csv.DictReader(open(RES / name, encoding="utf8")))

def num(x, d=2, sign=False):
    x = float(x)
    if round(x, d) == 0: x = 0.0  # no "-0.00"
    s = f"{x:+.{d}f}" if sign else f"{x:.{d}f}"
    return s.replace("+0." + "0" * d, "0." + "0" * d).replace("-", "−")

boot = {(r["measure"], r["group"]): r for r in table("bootstrap_median_ci.csv")}
def med(m, g, d=2, sign=False): return num(boot[m, g]["median"], d, sign)
def ci(m, g, d=2, sign=False):
    r = boot[m, g]
    return f"{num(r['lo'], d, sign)} to {num(r['hi'], d, sign)}"
def n(m, g): return boot[m, g]["n"]
def pct(m, g, d=1): return f"{float(boot[m, g]['median']) * 100:.{d}f}%"
def pct_ci(m, g):
    r = boot[m, g]
    return f"{float(r['lo']) * 100:.1f}–{float(r['hi']) * 100:.1f}%"

wil = {r["measure"]: r for r in table("wilcoxon_PRCR_vs_PD.csv")}
def pv(m, col="p_all"):
    p = float(wil[m][col])
    return f"{p:.1g}" if p < 0.001 else (f"{p:.3f}" if p < 0.01 else f"{p:.2f}")

de = {r["contrast"]: r for r in table("de_summary.csv")}
n_prcr = f"{int(de['C_on_vs_pre_PRCR']['fdr_0.1']):,}"
n_sd = de["C_on_vs_pre_SD"]["fdr_0.1"]
n_pd = de["C_on_vs_pre_PD"]["fdr_0.1"]
n_int = f"{int(de['C_on_vs_pre_PRCR_vs_PD']['fdr_0.1']):,}"
n_on_pre = de["B_on_vs_pre_paired"]["fdr_0.1"]

cor = {r["comparison"]: r for r in table("cross_modality_correlations.csv")}
def rho(k): return num(cor[k]["rho"], 2)
RW, RT = "Change in tumor score vs WES contraction", "RNA CD8 vs TCR T-cell fraction, On"
RW2 = RW + ", excluding PR/CR"

dil = {r["model"]: r for r in table("dilution_check.csv")}
d1, d2 = dil["CD8 change"], dil["CD8 change, adjusted for tumor score change"]
d3, d4 = dil["CD8 on treatment"], dil["CD8 on treatment, adjusted for tumor score on treatment"]

pap = {r["comparison"]: r for r in table("paper_comparison.csv")}
pA, pC = pap["Pre: PR/CR vs PD"], pap["On vs Pre, all patients"]
pD, pE = pap["On vs Pre change: benefit vs no benefit"], pap["On vs Pre change: PR/CR vs PD"]
n_pap_e = f"{int(float(pE['our_pipeline_fdr0.2'])):,}"

noise = sorted(float(r["genes_changed"]) for r in table("replicate_noise.csv"))
noise_med = (noise[len(noise) // 2 - 1] + noise[len(noise) // 2]) / 2 * 100
noise_max = noise[-1] * 100

pts = {r["patient"]: r for r in table("rna_patient_features.csv")}
pt30 = pts["Pt30"]
genes = sum(1 for _ in open(RES / "de_C_on_vs_pre_PD.csv")) - 1
n_pre = sum(r["Tcell_inflamed_18_Pre"] not in ("", "NA") for r in pts.values())
n_pair = sum(r["CD8_T_delta"] not in ("", "NA") for r in pts.values())
gsea = {(r["pathway"], r["contrast"]): r for r in table("gsea_hallmark.csv")}
def nes(p, c): return num(gsea["HALLMARK_" + p, c]["NES"], 1, True)


# ---------- drawing: same shapes and scheme colours as the WES deck ----------
E = 914400
def emu(v): return int(round(v * E))

class Slide:
    def __init__(self, layout, title, notes, title_ph=True):
        self.layout, self.notes, self.title_text = layout, notes, title
        self.parts, self.images = [], []
        self.nid = 1
        if title_ph:
            color = "bg1" if layout == 4 else "tx1"
            self.parts.append(
                f'<p:sp><p:nvSpPr><p:cNvPr id="{self.id()}" name="Title"/><p:cNvSpPr/><p:nvPr><p:ph idx="100" type="title" hasCustomPrompt="1"/></p:nvPr></p:nvSpPr>'
                f'<p:spPr><a:xfrm><a:off x="{emu(.6)}" y="{emu(.4)}"/><a:ext cx="{emu(12.1)}" cy="{emu(.95)}"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom><a:noFill/></p:spPr>'
                f'<p:txBody><a:bodyPr wrap="square" lIns="0" tIns="0" rIns="0" bIns="0" rtlCol="0" anchor="ctr"/><a:lstStyle/><a:p><a:pPr algn="l" indent="0" marL="0"><a:buNone/></a:pPr>'
                f'{run(title, 30 if len(title) < 50 else 24, color, True, True)}</a:p></p:txBody></p:sp>')

    def id(self):
        self.nid += 1
        return self.nid

    def shape(self, x, y, w, h, fill, geom="roundRect", name="Shape", adj=None):
        av = f'<a:avLst><a:gd name="adj" fmla="val {adj}"/></a:avLst>' if adj else "<a:avLst/>"
        clr = f'<a:srgbClr val="{fill}"/>' if len(fill) == 6 and fill.isupper() else f'<a:schemeClr val="{fill}"/>'
        self.parts.append(
            f'<p:sp><p:nvSpPr><p:cNvPr id="{self.id()}" name="{name}"/><p:cNvSpPr/><p:nvPr/></p:nvSpPr><p:spPr><a:xfrm><a:off x="{emu(x)}" y="{emu(y)}"/><a:ext cx="{emu(w)}" cy="{emu(h)}"/></a:xfrm>'
            f'<a:prstGeom prst="{geom}">{av}</a:prstGeom><a:solidFill>{clr}</a:solidFill><a:ln><a:noFill/></a:ln></p:spPr></p:sp>')

    def dot(self, x, y, d, fill, name="Dot"):
        self.shape(x, y, d, d, fill, "ellipse", name)

    def line(self, x, y, w):
        self.parts.append(
            f'<p:sp><p:nvSpPr><p:cNvPr id="{self.id()}" name="Timeline line"/><p:cNvSpPr/><p:nvPr/></p:nvSpPr><p:spPr><a:xfrm><a:off x="{emu(x)}" y="{emu(y)}"/><a:ext cx="{emu(w)}" cy="0"/></a:xfrm>'
            f'<a:prstGeom prst="line"><a:avLst/></a:prstGeom><a:noFill/><a:ln w="25400"><a:solidFill><a:schemeClr val="accent6"/></a:solidFill></a:ln></p:spPr></p:sp>')

    def text(self, x, y, w, h, paras, anchor="t", name="Text"):
        """paras: list of paragraphs; a paragraph is a list of runs (t, size, colour, bold, heading font, italic)."""
        ps = []
        for p in paras:
            opts = {}
            if isinstance(p, dict):
                opts, p = p, p["runs"]
            ppr = (f'<a:pPr marL="228600" indent="-228600"><a:spcAft><a:spcPts val="600"/></a:spcAft><a:buFont typeface="Arial"/><a:buChar char="&#8226;"/></a:pPr>' if opts.get("bullet")
                   else f'<a:pPr algn="{opts.get("align", "l")}" indent="0" marL="0"><a:buNone/></a:pPr>')
            ps.append(f'<a:p>{ppr}{"".join(run(*r) for r in p)}</a:p>')
        self.parts.append(
            f'<p:sp><p:nvSpPr><p:cNvPr id="{self.id()}" name="{name}"/><p:cNvSpPr txBox="1"/><p:nvPr/></p:nvSpPr><p:spPr><a:xfrm><a:off x="{emu(x)}" y="{emu(y)}"/><a:ext cx="{emu(w)}" cy="{emu(h)}"/></a:xfrm>'
            f'<a:prstGeom prst="rect"><a:avLst/></a:prstGeom><a:noFill/></p:spPr><p:txBody><a:bodyPr wrap="square" lIns="0" tIns="0" rIns="0" bIns="0" rtlCol="0" anchor="{anchor}"/><a:lstStyle/>{"".join(ps)}</p:txBody></p:sp>')

    def image(self, name, x, y, w, h, alt):
        path = RES / name
        pw, ph = png_size(path)
        s = min(w / pw, h / ph)
        iw, ih = pw * s, ph * s
        self.images.append(path)
        rid = f"rId{len(self.images) + 2}"
        self.parts.append(
            f'<p:pic><p:nvPicPr><p:cNvPr id="{self.id()}" name="{escape(path.stem)}" descr="{escape(alt, {chr(34): "&quot;"})}"/><p:cNvPicPr><a:picLocks noChangeAspect="1"/></p:cNvPicPr><p:nvPr/></p:nvPicPr>'
            f'<p:blipFill><a:blip r:embed="{rid}"/><a:stretch><a:fillRect/></a:stretch></p:blipFill><p:spPr><a:xfrm><a:off x="{emu(x)}" y="{emu(y)}"/><a:ext cx="{emu(iw)}" cy="{emu(ih)}"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom></p:spPr></p:pic>')

    def number(self):
        self.parts.append(
            f'<p:sp><p:nvSpPr><p:cNvPr id="{self.id()}" name="Slide Number"/><p:cNvSpPr><a:spLocks noGrp="1"/></p:cNvSpPr><p:nvPr><p:ph type="sldNum" sz="quarter" idx="4294967295"/></p:nvPr></p:nvSpPr>'
            f'<p:spPr><a:xfrm><a:off x="11155680" y="6355080"/><a:ext cx="457200" cy="320040"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom></p:spPr><p:txBody><a:bodyPr/><a:lstStyle><a:lvl1pPr><a:defRPr sz="1000"><a:solidFill><a:schemeClr val="accent5"/></a:solidFill></a:defRPr></a:lvl1pPr></a:lstStyle>'
            '<a:p><a:pPr algn="r"/><a:fld id="{F7021451-1387-4CA6-816F-3879F97B5CBC}" type="slidenum"><a:rPr b="0" lang="en-US"/><a:t>1</a:t></a:fld><a:endParaRPr lang="en-US"/></a:p></p:txBody></p:sp>')

    def xml(self):
        return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><p:sld xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main">'
                '<p:cSld><p:spTree><p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr><p:grpSpPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="0"/><a:chOff x="0" y="0"/><a:chExt cx="0" cy="0"/></a:xfrm></p:grpSpPr>'
                + "".join(self.parts) + '</p:spTree></p:cSld><p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:sld>')


def run(t, size=15, color="tx2", bold=False, heading=False, italic=False):
    font = '<a:latin typeface="Cambria"/><a:cs typeface="Cambria"/>' if heading else ""
    attrs = (' b="1"' if bold else "") + (' i="1"' if italic else "")
    return (f'<a:r><a:rPr lang="en-US" sz="{int(size * 100)}"{attrs} dirty="0"><a:solidFill><a:schemeClr val="{color}"/></a:solidFill>{font}</a:rPr><a:t>{escape(t)}</a:t></a:r>')

def png_size(p):
    with open(p, "rb") as f:
        f.read(16)
        return struct.unpack(">II", f.read(8))

GROUPS = [("accent1", "PR/CR"), ("accent2", "Stable"), ("accent3", "Progressive")]
KEYS = ["PR/CR", "SD", "PD"]

def stat_rows(s, x, y, rows, val_w=1.45, step=.95):
    """Dot, big value and a label with a CI: the right-hand column on the WES result slides."""
    for i, (color, value, label, sub) in enumerate(rows):
        yy = y + i * step
        s.dot(x, yy + .25, .3, color, f"Dot {i + 1}")
        s.text(x + .45, yy, val_w, .8, [[(value, 30, "tx1", True, True)]], "ctr", f"Value {i + 1}")
        s.text(x + .45 + val_w + .05, yy, 3.6 - .5 - val_w, .8, [[(label, 14)], [(sub, 12, "accent5")]], "ctr", f"Label {i + 1}")


# ---------- slides ----------
slides = []
tx = {r["tcr"] + "|" + r["rna"]: r for r in table("rna_vs_tcr_expansion.csv")}
tg = {r["Response"]: r for r in table("tcr_expansion_by_response.csv")}
EXP = "Expanded clones (depth-matched)"

def fig_slide(title, img, alt, bullets, paper, notes, mine=True):
    s = Slide(3, title, notes); s.mine = mine
    s.image(img, .6, 1.55, 8.5, 4.95, alt)
    s.text(9.4, 1.6, 3.3, 2.9, [{"bullet": True, "runs": [(b, 15)]} for b in bullets], name="Bullets")
    tag, body = paper
    s.shape(9.4, 4.55, 3.3, 1.95, "bg2", adj=6000, name="Paper card")
    s.text(9.6, 4.7, 2.95, 1.75, [[(tag, 13, "accent1", True)], [(body, 12)]], name="Paper note")
    s.number(); slides.append(s); return s

# 1 Title
s = Slide(2, "", "", title_ph=False); s.mine = False
s.parts.append(
    f'<p:sp><p:nvSpPr><p:cNvPr id="{s.id()}" name="Title"/><p:cNvSpPr/><p:nvPr><p:ph idx="100" type="title" hasCustomPrompt="1"/></p:nvPr></p:nvSpPr><p:spPr><a:xfrm><a:off x="{emu(.8)}" y="{emu(2.1)}"/><a:ext cx="{emu(11.7)}" cy="{emu(1.9)}"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom><a:noFill/></p:spPr>'
    f'<p:txBody><a:bodyPr wrap="square" rtlCol="0"/><a:lstStyle/><a:p><a:pPr algn="l" indent="0" marL="0"><a:buNone/></a:pPr>{run("How tumors and T cells change when nivolumab works", 40, "bg1", True, True)}</a:p></p:txBody></p:sp>')
s.title_text = "How tumors and T cells change when nivolumab works"
s.parts.append(
    f'<p:sp><p:nvSpPr><p:cNvPr id="{s.id()}" name="Subtitle"/><p:cNvSpPr/><p:nvPr><p:ph idx="101" type="body" hasCustomPrompt="1"/></p:nvPr></p:nvSpPr><p:spPr><a:xfrm><a:off x="{emu(.8)}" y="{emu(4.25)}"/><a:ext cx="{emu(11.7)}" cy="{emu(1.6)}"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom><a:noFill/></p:spPr>'
    f'<p:txBody><a:bodyPr wrap="square" rtlCol="0"/><a:lstStyle/>'
    f'<a:p><a:pPr algn="l" indent="0" marL="0"><a:buNone/></a:pPr>{run("A reanalysis of Riaz et al., Cell 2017: WES, TCR-seq and RNA-seq", 18, "bg1", True)}</a:p>'
    f'<a:p><a:pPr algn="l" indent="0" marL="0"><a:buNone/></a:pPr>{run("Kwaku (WES)  ·  Krithika (TCR)  ·  Ethan (RNA)", 18, "accent6")}</a:p></p:txBody></p:sp>')
for x, c in zip((.8, 1.25, 1.7), ("accent1", "accent2", "accent3")):
    s.dot(x, 1.2, .32, c, "Motif dot")
slides.append(s)

# 2 Introduction
s = Slide(3, "Riaz et al. biopsied melanomas before and during nivolumab to see what changes when PD-1 blockade works",
          "Nivolumab is an antibody that blocks PD-1, the brake on T cells. In the CheckMate 038 trial, Riaz and colleagues took a biopsy from patients with advanced melanoma before the first dose and again about four weeks in, "
          "and sequenced the tumor DNA, the T-cell receptors and the RNA. Half the cohort had never had ipilimumab; the other half had already progressed on it. "
          "Their headline: in responders the tumor loses mutations, T-cell clones expand, and immune genes switch on. We re-analyzed their public data to test each of those claims. [~40 s]")
s.mine = False
s.line(2.1, 2.05, 9.0)
for i, (head, body, c) in enumerate([("Biopsy 1", "Before the first dose", "tx2"), ("Nivolumab", "Anti-PD-1: releases the T-cell brake", "accent4"),
                                     ("Biopsy 2", "About 4 weeks in", "tx2"), ("Response", "PR/CR shrank · SD stable · PD grew", "accent1")]):
    x = .75 + i * 3
    s.dot(x + .95, 1.65, .8, c, f"Timeline dot {i + 1}")
    s.text(x + .95, 1.65, .8, .8, [{"runs": [(str(i + 1), 22, "bg1", True)], "align": "ctr"}], "ctr", f"Timeline number {i + 1}")
    s.text(x, 2.6, 2.7, .45, [{"runs": [(head, 18, "tx1", True)], "align": "ctr"}], "ctr", f"Timeline head {i + 1}")
    s.text(x, 3.05, 2.7, .7, [{"runs": [(body, 14)], "align": "ctr"}], name=f"Timeline body {i + 1}")
s.text(.6, 3.95, 12, .45, [[("CheckMate 038 · advanced melanoma · ipilimumab-naive and ipilimumab-progressed cohorts", 15, "accent5", False, False, True)]], "ctr", "Cohort")
for i, (head, body, c) in enumerate([("Tumor DNA", "Responders' tumors lose mutations and neoantigens", "accent1"),
                                     ("T-cell receptors", "T-cell clones expand in responders", "accent2"),
                                     ("RNA", "Immune genes rise; cell-cycle and melanin genes fall", "accent3")]):
    x = .6 + i * 4.1
    s.shape(x, 4.6, 3.85, 1.75, "bg2", adj=6000, name=f"Finding card {i + 1}")
    s.dot(x + .3, 4.85, .35, c, f"Finding dot {i + 1}")
    s.text(x + .8, 4.8, 2.9, .45, [[(head, 17, "tx1", True)]], "ctr", f"Finding head {i + 1}")
    s.text(x + .3, 5.35, 3.3, .9, [[(body, 15)]], name=f"Finding body {i + 1}")
s.number(); slides.append(s)

# 3 Modalities
s = Slide(3, "We follow the biology in order: the targets, the T cells that respond, then what the tissue is doing",
          "We present three data types in this order. Kwaku starts with whole-exome sequencing: which mutations, and so which targets, the tumor has, and whether they disappear. "
          "Krithika follows with T-cell receptor sequencing: which T-cell clones respond and expand. I finish with RNA, which reads out the whole biopsy at once, both the immune activity and the tumor that is left, "
          "so it ties the other two together and leads into our comparison. [~30 s]")
s.mine = False
for i, (head, who, body, c) in enumerate([("WES", "Kwaku", "Which mutations and neoantigens exist, and are they lost on treatment?", "accent1"),
                                          ("TCR-seq", "Krithika", "Which T-cell clones respond and expand?", "accent2"),
                                          ("RNA-seq", "Ethan", "What is the biopsy doing: immune activity and remaining tumor?", "accent3")]):
    x = .6 + i * 4.1
    s.shape(x, 1.7, 3.85, 2.75, "bg2", adj=5000, name=f"Modality card {i + 1}")
    s.dot(x + .35, 2.0, .7, c, f"Modality dot {i + 1}")
    s.text(x + .35, 2.0, .7, .7, [{"runs": [(str(i + 1), 22, "bg1", True)], "align": "ctr"}], "ctr", f"Modality number {i + 1}")
    s.text(x + 1.25, 2.0, 2.5, .7, [[(head, 24, "tx1", True, True)], [(who, 14, "accent5")]], "ctr", f"Modality head {i + 1}")
    s.text(x + .35, 3.0, 3.2, 2.1, [[(body, 17)]], name=f"Modality body {i + 1}")
    if i < 2:
        s.text(x + 3.85, 3.2, .25, .5, [{"runs": [("→", 22, "accent6")], "align": "ctr"}], "ctr", f"Arrow {i + 1}")
s.text(.6, 4.8, 12, .8, [[("Why this order? ", 16, "tx1", True), ("DNA sets the targets, TCR shows who attacks them, and RNA shows the outcome in the tissue, linking both.", 16)]], name="Order")
s.number(); slides.append(s)

# 4 RNA workflow
s = Slide(3, "RNA-seq turns each biopsy into gene and program scores, compared within each patient",
          f"My part is the RNA. We start from the authors' read counts, {genes:,} genes after filtering. We normalize for sequencing depth, then average marker genes into scores, for example CD8 T cells, "
          f"interferon-gamma, and a tumor score built from melanoma identity genes. Then we compare each patient's two biopsies, so every patient is their own control. That matters because melanomas differ so much from person to person that a between-patient comparison would mostly measure those differences. "
          f"{n_pre} patients have a pre-treatment biopsy, and {n_pair} have both. [~30 s]")
s.mine = True
steps = [("Count genes", f"Authors' read counts; {genes:,} genes pass the filter"), ("Normalize", "Correct for sequencing depth (DESeq2 VST)"),
         ("Score programs", "Average marker genes: CD8, IFN-γ, tumor (melanocytic)…"), ("Compare biopsies", "On vs Pre per patient: DESeq2 genes, fgsea pathways")]
for i, (head, body) in enumerate(steps):
    x = .6 + i * 3.08
    s.shape(x, 1.65, 2.85, 2.6, "bg2", adj=5000, name=f"Step card {i + 1}")
    s.dot(x + .3, 1.95, .6, "accent1" if i == 3 else "tx2", f"Step dot {i + 1}")
    s.text(x + .3, 1.95, .6, .6, [{"runs": [(str(i + 1), 18, "bg1", True)], "align": "ctr"}], "ctr", f"Step number {i + 1}")
    s.text(x + .3, 2.7, 2.4, .5, [[(head, 18, "tx1", True)]], "ctr", f"Step head {i + 1}")
    s.text(x + .3, 3.2, 2.35, 1.0, [[(body, 14)]], name=f"Step body {i + 1}")
for i, (v, lab, c) in enumerate([(str(n_pre), "patients with a pre-treatment biopsy", "tx1"), (str(n_pair), "with both biopsies (main comparison)", "accent1"),
                                 (f"{n('tumor_score_change', 'PR/CR')} · {n('tumor_score_change', 'SD')} · {n('tumor_score_change', 'PD')}", "paired PR/CR · SD · PD", "tx1")]):
    x = .6 + i * 4.1
    s.text(x, 4.75, 1.9, 1, [[(v, 40, c, True, True)]], "ctr", f"Cohort number {i + 1}")
    s.text(x + 1.95, 4.75, 2, 1, [[(lab, 14)]], "ctr", f"Cohort label {i + 1}")
s.number(); slides.append(s)

# 5 Baseline
m = "baseline_immune_score"
fig_slide("Before treatment, immune activity only weakly predicts who will respond", "01_baseline_immune_by_response.png",
          "Boxplots of pre-treatment T-cell-inflamed score by response, split by prior ipilimumab",
          [f"Medians: PR/CR {med(m, 'PR/CR')}, SD {med(m, 'SD')}, PD {med(m, 'PD')}",
           f"PR/CR − PD {med(m, 'PR/CR minus PD', 2, True)}; 95% CI {ci(m, 'PR/CR minus PD', 2, True)} crosses zero",
           "Ipilimumab-progressed responders start 'hot'"],
          ("Recreates the paper", "Paper: baseline expression only marginally predictive, split by ipilimumab cohort. Our change: one fixed 18-gene score (Ayers 2017) instead of genome-wide testing."),
          f"First question: does immune activity before treatment predict response? Each dot is one patient's T-cell-inflamed score, a published 18-gene measure of how immune-active a tumor is. "
          f"Responders sit a little higher, but the difference from progressors, {med(m, 'PR/CR minus PD', 2)}, has a confidence interval that crosses zero. "
          f"The one clear pattern is on the right: responders who had already failed ipilimumab all start with 'hot' tumors. That matches the paper, which also called baseline expression only marginally predictive. "
          f"In the patients who had never had ipilimumab, responders and progressors overlap almost completely. So a single snapshot before treatment is not enough; to understand response, we need to look at what changes during treatment. [~50 s]")

# 6 Genes changing
m = "share_genes_changed"
fig_slide("On treatment, responders' genes shift more, and far more consistently, than progressors'", "02_gene_changes_per_patient.png",
          "Stacked bars of the share of genes up or down at least two-fold per patient, by response",
          [f"Median genes changing 2-fold: {pct(m, 'PR/CR')} PR/CR vs {pct(m, 'PD')} PD",
           f"Genes shifting together (DESeq2, FDR < 0.1): {n_prcr} PR/CR, {n_sd} SD, {n_pd} PD"],
          ("Builds on the paper", f"Paper counted On vs Pre DE genes; we could not reproduce its exact counts (before treatment: {pA['authors_design_this_R']} ours vs {pA['paper_q0.2']} paper). New: per-patient view and a repeat-library noise floor (~{noise_med:.0f}%)."),
          f"Next, the two biopsies. Each bar is a patient: the share of genes that went up, in orange, or down, in blue, at least two-fold. Responders change more, a median of {pct(m, 'PR/CR')} of genes against {pct(m, 'PD')} in progressors. "
          f"But the bigger difference is consistency. When we ask which genes shift together across a group, we find {n_prcr} in responders, {n_sd} in stable disease and {n_pd} in progressors. "
          f"For scale, two sequencing runs of the same biopsy already differ by about 3 percent of genes, so the raw amount of change is modest. Responders move in a shared direction; progressors, when they change, change in their own way. [~50 s]")

# 7 Tumor vs CD8
m, c8 = "tumor_score_change", "cd8_score_change"
fig_slide("Tumor loss, more than immune gain, marks a response", "03_tumor_immune_change_per_patient.png",
          "Per-patient dot plots of tumor score change and CD8 score change, sorted by tumor score change",
          [f"Tumor score: PR/CR {med(m, 'PR/CR')} vs PD {med(m, 'PD')} (p ≈ {pv(m)})",
           f"CD8 score: PR/CR {med(c8, 'PR/CR', 2, True)}, SD {med(c8, 'SD', 2, True)}, PD {med(c8, 'PD', 2)} (p ≈ {pv(c8)})",
           "Exceptions: Pt94 and Pt101"],
          ("New beyond the paper", "Paper showed melanin pathway genes falling in responders. We turn this into a per-patient tumor score and line it up against CD8 change."),
          f"So which genes are changing? Each row is a patient, sorted by the change in a tumor score made from five melanoma identity genes, like MLANA and tyrosinase. "
          f"Responders, in red, cluster at the top with the biggest drops: a median of {med(m, 'PR/CR')} against {med(m, 'PD')} in progressors, p about {pv(m)}. "
          f"On the right, the CD8 T-cell score rises in responders, but it also rises in stable disease, and the responder-versus-progressor test is only p about {pv(c8)}. "
          f"It isn't universal: patients 94 and 101 respond with almost no change in the tumor score. But overall, the clearest RNA sign of response is the tumor disappearing from the biopsy. [~50 s]")

# 8 Pathways
fig_slide("Immune pathways rise in every group; only responders shut down proliferation", "04_pathways_by_response.png",
          "Heatmap of normalized enrichment scores for immune and proliferation Hallmark pathways",
          [f"IFN-γ: {nes('INTERFERON_GAMMA_RESPONSE', 'PR/CR')} PR/CR, {nes('INTERFERON_GAMMA_RESPONSE', 'SD')} SD, {nes('INTERFERON_GAMMA_RESPONSE', 'PD')} PD",
           f"E2F targets: {nes('E2F_TARGETS', 'PR/CR')} PR/CR vs {nes('E2F_TARGETS', 'PD')} PD"],
          ("Recreates the paper", "Paper: immune genes up and cell-cycle genes down in responders (GSEA). Our change: Hallmark sets with fgsea, plus SD and PD columns so the contrast between groups is visible."),
          f"Pathways show the same thing across thousands of genes. Orange means up on treatment, blue down, grey not significant. "
          f"The top block is immune pathways: interferon-gamma rises in every group, even progressors at {nes('INTERFERON_GAMMA_RESPONSE', 'PD')}. "
          f"The bottom block is proliferation: it collapses only in responders, E2F targets at {nes('E2F_TARGETS', 'PR/CR')}, while progressors' tumors keep dividing. "
          f"This recreates the paper's result, immune up and cell cycle down in responders, and adding the other groups shows that immune activation alone is not enough. [~45 s]")

# 9 CD8 dilution
fig_slide("The on-treatment CD8 level is the immune signal that survives the tumor-dilution check", "05_cd8_dilution_check.png",
          "Boxplot of on-treatment CD8 score by response and forest plot of PR/CR minus PD effects with and without tumor-score adjustment",
          [f"CD8 rise: {num(d1['PRCR_effect'], 2, True)} → {num(d2['PRCR_effect'], 2, True)} after adjustment",
           f"On-treatment CD8: {num(d3['PRCR_effect'], 2, True)} → {num(d4['PRCR_effect'], 2, True)}, still clear of zero"],
          ("New beyond the paper", "Bulk RNA is a mixture: when tumor cells die, other cells' share grows. We add an adjustment for tumor content."),
          f"One catch with bulk RNA: it is a mixture. When tumor cells die, immune cells become a bigger share of what's left, even if no new T cells arrive. "
          f"So we adjusted for the tumor score. The responder CD8 rise drops from {num(d1['PRCR_effect'], 2, True)} to {num(d2['PRCR_effect'], 2, True)}, so it can't be separated from dilution. "
          f"But the on-treatment CD8 level holds, {num(d3['PRCR_effect'], 2, True)} before and {num(d4['PRCR_effect'], 2, True)} after. "
          f"So by week four, responders' tumors are gone and what remains is full of CD8 T cells. That's the RNA story, and it leads into how RNA compares with the TCR data. [~45 s]")

# 10 Comparison of modalities
r_exp = tx[EXP + "|CD8 T-cell score change"]
rt24 = num(tx["T-cell fraction (On)|CD8 T-cell score, On"]["rho"])
s = Slide(3, "RNA and TCR agree on how many T cells are present, but only TCR-seq sees which clones expand",
          f"To compare modalities: RNA and TCR sequencing agree almost perfectly on how many T cells are in a biopsy, rho {rt24}. "
          f"But Krithika's count of expanding clones barely tracks any RNA score, for example rho {num(r_exp['rho'])} with CD8 change, and no expansion or clonality measure survives multiple-testing correction against any RNA score. "
          f"So RNA measures the size of the immune response, and TCR-seq measures its clonal makeup: complementary, not redundant. [~35 s]")
s.mine = False
s.image("06_rna_vs_tcr_expansion.png", .6, 1.5, 8.7, 3.4, "Three panels: TCR expansion by response, RNA CD8 vs TCR T-cell fraction, and RNA CD8 change vs expanded clones")
s.text(9.6, 1.5, 3.1, .4, [[("Modality comparisons", 16, "tx1", True)]], "ctr", "Matrix header")
for i, (pair, body, filled) in enumerate([("WES vs TCR", "Kwaku & Krithika", False), ("WES vs RNA", "Kwaku", False),
                                          ("TCR vs RNA", f"T-cell amount ρ {rt24}; expansion vs CD8 change ρ {num(r_exp['rho'])}", True)]):
    y = 2.0 + i * 1.0
    s.shape(9.6, y, 3.1, .9, "accent3" if filled else "bg2", adj=8000, name=f"Pair card {i + 1}")
    s.text(9.8, y + .08, 2.8, .75, [[(pair, 14, "bg1" if filled else "tx1", True)], [(body, 11, "bg1" if filled else "accent5")]], "ctr", f"Pair text {i + 1}")
s.text(.6, 5.15, 12.1, 1.3, [
    {"bullet": True, "runs": [(f"Expanded clones by response (depth-matched median): PR/CR {float(tg['PR/CR']['med_expanded']):.0f}, SD {float(tg['SD']['med_expanded']):.0f}, PD {float(tg['PD']['med_expanded']):.0f}; PR/CR vs PD p = {float(tg['PR/CR']['p_PRCR_vs_PD']):.2f} (n = {tg['PR/CR']['n']} vs {tg['PD']['n']})", 14)]},
    {"bullet": True, "runs": [("Panel A recreates the paper's expansion-by-response comparison; B and C are new links between modalities", 14)]}], name="Comparison bullets")
s.number(); slides.append(s)


# ---------- package: reuse the WES deck's theme, master and layouts ----------
NOTES_NS = 'xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main"'
REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"

def notes_xml(text, i):
    return (f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><p:notes {NOTES_NS}><p:cSld><p:spTree><p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr><p:grpSpPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="0"/><a:chOff x="0" y="0"/><a:chExt cx="0" cy="0"/></a:xfrm></p:grpSpPr>'
            '<p:sp><p:nvSpPr><p:cNvPr id="2" name="Slide Image Placeholder 1"/><p:cNvSpPr><a:spLocks noGrp="1" noRot="1" noChangeAspect="1"/></p:cNvSpPr><p:nvPr><p:ph type="sldImg"/></p:nvPr></p:nvSpPr><p:spPr/></p:sp>'
            f'<p:sp><p:nvSpPr><p:cNvPr id="3" name="Notes Placeholder 2"/><p:cNvSpPr><a:spLocks noGrp="1"/></p:cNvSpPr><p:nvPr><p:ph type="body" idx="1"/></p:nvPr></p:nvSpPr><p:spPr/><p:txBody><a:bodyPr/><a:lstStyle/><a:p><a:r><a:rPr lang="en-US" dirty="0"/><a:t>{escape(text)}</a:t></a:r></a:p></p:txBody></p:sp>'
            '</p:spTree></p:cSld><p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:notes>')

src = zipfile.ZipFile(TEMPLATE)
names = src.namelist()
skip = re.compile(r"ppt/(slides|notesSlides|media)/")
parts = {nm: src.read(nm) for nm in names if not skip.match(nm) and not nm.endswith("/")}

for nm in list(parts):
    if nm.startswith("ppt/slideLayouts/slideLayout") and nm.endswith(".xml"):
        parts[nm] = parts[nm].decode("utf8").replace("Analysis: Kwaku Darkwah (WES)", f"Analysis: {AUTHOR} (RNA)").encode("utf8")
theme = "ppt/theme/theme1.xml"
parts[theme] = parts[theme].decode("utf8").replace("Nivolumab WES", "Nivolumab RNA").encode("utf8")

# content types
ct = parts["[Content_Types].xml"].decode("utf8")
ct = re.sub(r'<Override PartName="/ppt/(slides|notesSlides)/[^>]*/>', "", ct)
extra = "".join(f'<Override PartName="/ppt/slides/slide{i}.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slide+xml"/>'
                f'<Override PartName="/ppt/notesSlides/notesSlide{i}.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.notesSlide+xml"/>' for i in range(1, len(slides) + 1))
parts["[Content_Types].xml"] = ct.replace("</Types>", extra + "</Types>").encode("utf8")

# presentation parts: slides get rId101+, everything else keeps its id
rels = parts["ppt/_rels/presentation.xml.rels"].decode("utf8")
rels = re.sub(r'<Relationship [^>]*relationships/slide"[^>]*/>', "", rels)
new = "".join(f'<Relationship Id="rId{100 + i}" Type="{REL}/slide" Target="slides/slide{i}.xml"/>' for i in range(1, len(slides) + 1))
parts["ppt/_rels/presentation.xml.rels"] = rels.replace("</Relationships>", new + "</Relationships>").encode("utf8")
pres = parts["ppt/presentation.xml"].decode("utf8")
lst = "".join(f'<p:sldId id="{255 + i}" r:id="rId{100 + i}"/>' for i in range(1, len(slides) + 1))
parts["ppt/presentation.xml"] = re.sub(r"<p:sldIdLst>.*?</p:sldIdLst>", f"<p:sldIdLst>{lst}</p:sldIdLst>", pres, flags=re.S).encode("utf8")

# document properties
core = parts["docProps/core.xml"].decode("utf8")
core = re.sub(r"<dc:title>.*?</dc:title>", "<dc:title>RNA-seq results: gene expression under nivolumab</dc:title>", core)
core = re.sub(r"<dc:subject>.*?</dc:subject>", "<dc:subject>Reanalysis of Riaz et al., Cell 2017</dc:subject>", core)
core = core.replace("Kwaku Darkwah", AUTHOR)
parts["docProps/core.xml"] = core.encode("utf8")
parts["docProps/app.xml"] = (f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties"><Application>Microsoft Office PowerPoint</Application><Slides>{len(slides)}</Slides><Notes>{len(slides)}</Notes></Properties>').encode("utf8")

media = {}
for i, sl in enumerate(slides, 1):
    r = [f'<Relationship Id="rId1" Type="{REL}/slideLayout" Target="../slideLayouts/slideLayout{sl.layout}.xml"/>',
         f'<Relationship Id="rId2" Type="{REL}/notesSlide" Target="../notesSlides/notesSlide{i}.xml"/>']
    for j, img in enumerate(sl.images, 3):
        media[img.name] = img.read_bytes()
        r.append(f'<Relationship Id="rId{j}" Type="{REL}/image" Target="../media/{img.name}"/>')
    parts[f"ppt/slides/slide{i}.xml"] = sl.xml().encode("utf8")
    parts[f"ppt/slides/_rels/slide{i}.xml.rels"] = (f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">{"".join(r)}</Relationships>').encode("utf8")
    parts[f"ppt/notesSlides/notesSlide{i}.xml"] = notes_xml(sl.notes, i).encode("utf8")
    parts[f"ppt/notesSlides/_rels/notesSlide{i}.xml.rels"] = (f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="{REL}/notesMaster" Target="../notesMasters/notesMaster1.xml"/><Relationship Id="rId2" Type="{REL}/slide" Target="../slides/slide{i}.xml"/></Relationships>').encode("utf8")
for k, v in media.items():
    parts[f"ppt/media/{k}"] = v

with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED) as z:
    z.writestr("[Content_Types].xml", parts.pop("[Content_Types].xml"))
    for k, v in parts.items():
        z.writestr(k, v)

# ---------- presenter script, from the same notes ----------
def secs(sl): return sum(int(x) for x in re.findall(r"\[~(\d+) s\]", sl.notes))
mine = sum(secs(sl) for sl in slides if sl.mine)
lines = ["# Presenter script\n", f"Matches `{OUT.name}` (also in its speaker notes). Numbers come from `rna/results/`.\n",
         f"## My RNA section: about {mine / 60:.1f} minutes (slides {', '.join(str(i) for i, sl in enumerate(slides, 1) if sl.mine)})\n"]
def add(i, sl):
    text, t = re.match(r"(.*) \[(~\d+ s)\]$", sl.notes, re.S).groups()
    lines.append(f"### Slide {i}: {sl.title_text} ({t})\n\n{text}\n")
for i, sl in enumerate(slides, 1):
    if sl.mine: add(i, sl)
lines.append("## Shared slides: introduction, modalities, comparison (not in the 5-minute count)\n")
for i, sl in enumerate(slides, 1):
    if not sl.mine and sl.notes: add(i, sl)
SCRIPT.write_text("\n".join(lines), encoding="utf8")
print(f"Wrote {OUT.name} ({len(slides)} slides; my part {mine} s) and {SCRIPT.name}")
