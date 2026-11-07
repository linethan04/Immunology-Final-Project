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
OUT = ROOT / "slides" / "RNA_11slides_Ethan.pptx"
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
                f'{run(title, 32, color, True, True)}</a:p></p:txBody></p:sp>')

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

# 1 Title
s = Slide(2, "", "", title_ph=False)
s.notes = ("Hi, I'm Ethan, and I led the bulk RNA sequencing part of our project. Kwaku asked whether tumors lose their mutations when nivolumab works. "
           "I ask the matching question for gene expression: what changes in the biopsy, and is that change specific to patients who respond? [~20 s]")
s.parts.append(
    f'<p:sp><p:nvSpPr><p:cNvPr id="{s.id()}" name="Title"/><p:cNvSpPr/><p:nvPr><p:ph idx="100" type="title" hasCustomPrompt="1"/></p:nvPr></p:nvSpPr><p:spPr><a:xfrm><a:off x="{emu(.8)}" y="{emu(2.1)}"/><a:ext cx="{emu(11.7)}" cy="{emu(1.9)}"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom><a:noFill/></p:spPr>'
    f'<p:txBody><a:bodyPr wrap="square" rtlCol="0"/><a:lstStyle/><a:p><a:pPr algn="l" indent="0" marL="0"><a:buNone/></a:pPr>{run("Does gene expression change when nivolumab works?", 44, "bg1", True, True)}</a:p></p:txBody></p:sp>')
s.title_text = "Does gene expression change when nivolumab works?"
s.parts.append(
    f'<p:sp><p:nvSpPr><p:cNvPr id="{s.id()}" name="Subtitle"/><p:cNvSpPr/><p:nvPr><p:ph idx="101" type="body" hasCustomPrompt="1"/></p:nvPr></p:nvSpPr><p:spPr><a:xfrm><a:off x="{emu(.8)}" y="{emu(4.25)}"/><a:ext cx="{emu(11.7)}" cy="{emu(1.6)}"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom><a:noFill/></p:spPr>'
    f'<p:txBody><a:bodyPr wrap="square" rtlCol="0"/><a:lstStyle/>'
    f'<a:p><a:pPr algn="l" indent="0" marL="0"><a:buNone/></a:pPr>{run("Bulk RNA sequencing (RNA-seq) results", 18, "bg1", True)}</a:p>'
    f'<a:p><a:pPr algn="l" indent="0" marL="0"><a:buNone/></a:pPr>{run(f"{AUTHOR}, RNA lead  ·  with Krithika (TCR) and Kwaku (WES)", 18, "accent6")}</a:p>'
    f'<a:p><a:pPr algn="l" indent="0" marL="0"><a:buNone/></a:pPr>{run("Reanalysis of Riaz et al., Cell 2017", 18, "accent6")}</a:p></p:txBody></p:sp>')
for x, c in zip((.8, 1.25, 1.7), ("accent1", "accent2", "accent3")):
    s.dot(x, 1.2, .32, c, "Motif dot")
slides.append(s)

# 2 Biology
s = Slide(3, "Nivolumab releases the brake; RNA shows what follows",
          "Kwaku covered how mutations create targets. RNA picks up the story from there. T cells recognize the tumor, but PD-1 acts as a brake. Nivolumab blocks PD-1. "
          "If that works, active T cells should switch on interferon-gamma and killing genes, and cancer cells should die, so melanoma genes should fade from the biopsy. "
          "DNA tells us which targets exist. RNA tells us what the cells in the biopsy are doing. So RNA can measure both the immune response and its effect on the tumor. [~35 s]")
chain = [("T cells", "Recognize the tumor, but PD-1 brakes them", "tx2"),
         ("Nivolumab", "Blocks PD-1 and releases the brake", "accent4"),
         ("Immune genes", "Active T cells switch on IFN-γ and killing genes", "accent3"),
         ("Tumor genes", "As cancer cells die, melanoma genes fade", "accent1")]
for i, (head, body, c) in enumerate(chain):
    x = .6 + i * 3.12
    s.shape(x, 1.6, 2.75, 2.25, "bg2", adj=6486, name=f"Chain card {i + 1}")
    s.dot(x + .3, 1.9, .4, c, f"Chain dot {i + 1}")
    s.text(x + .85, 1.85, 1.8, .5, [[(head, 18, "tx1", True)]], "ctr", f"Chain head {i + 1}")
    s.text(x + .3, 2.55, 2.15, 1.2, [[(body, 15)]], name=f"Chain body {i + 1}")
    if i < 3:
        s.text(x + 2.75, 2.45, .37, .5, [{"runs": [("→", 24, "accent6")], "align": "ctr"}], "ctr", f"Chain arrow {i + 1}")
s.shape(6.84, 4.25, 5.87, .5, "accent1", adj=20000, name="RNA bracket")
s.text(6.84, 4.25, 5.87, .5, [{"runs": [("RNA measures these two, before and during treatment", 14, "bg1", True)], "align": "ctr"}], "ctr", "RNA bracket label")
s.text(.6, 5.1, 11.6, 1.2, [[("Why RNA? ", 18, "tx1", True), ("DNA shows which targets exist. RNA shows what the cells are doing. If T cells attack, immune genes should rise and melanoma genes should fall in patients who respond.", 18, "tx1")]], name="Why RNA")
s.number(); slides.append(s)

# 3 Design
s = Slide(3, "The study sampled each tumor twice",
          f"This is the same study design Kwaku showed: one biopsy before the first dose, a second about four weeks in, and scans that later sort patients into responders, stable disease and progressive disease. "
          f"For RNA, {n_pre} patients with a response call have a pre-treatment biopsy and {n_pair} have both biopsies. Those {n_pair} are my main comparison. "
          f"Having both biopsies lets every patient act as their own control, which matters because tumors differ so much between people. [~30 s]")
s.line(2.1, 2.15, 9.0)
tl = [("Biopsy 1", "Taken before the first dose", "tx2"), ("Nivolumab", "Anti-PD-1 antibody that releases the brake on T cells", "accent4"),
      ("Biopsy 2", "Taken about 4 weeks into treatment", "tx2"), ("Response", "Scans sort patients: PR/CR shrank, SD stable, PD grew", "accent1")]
for i, (head, body, c) in enumerate(tl):
    x = .75 + i * 3
    s.dot(x + .95, 1.75, .8, c, f"Timeline dot {i + 1}")
    s.text(x + .95, 1.75, .8, .8, [{"runs": [(str(i + 1), 22, "bg1", True)], "align": "ctr"}], "ctr", f"Timeline number {i + 1}")
    s.text(x, 2.75, 2.7, .45, [{"runs": [(head, 18, "tx1", True)], "align": "ctr"}], "ctr", f"Timeline head {i + 1}")
    s.text(x, 3.2, 2.7, .9, [{"runs": [(body, 14)], "align": "ctr"}], name=f"Timeline body {i + 1}")
coh = [(str(n_pre), "patients with RNA before treatment", "tx1"), (str(n_pair), "with RNA at both biopsies: my main comparison", "accent1"),
       ("26", "with WES, RNA and TCR at both biopsies (team's shared set)", "tx1")]
for i, (v, lab, c) in enumerate(coh):
    x = .6 + i * 4.1
    s.text(x, 4.9, 1.3, 1, [[(v, 48, c, True, True)]], "ctr", f"Cohort number {i + 1}")
    s.text(x + 1.35, 4.9, 2.55, 1, [[(lab, 14)]], "ctr", f"Cohort label {i + 1}")
s.number(); slides.append(s)

# 4 Method
s = Slide(3, "RNA-seq turns each biopsy into gene and program scores",
          f"Here is the workflow. We start from the authors' table of read counts: how many RNA reads came from each gene in each biopsy. After removing rarely detected genes, {genes:,} remain. "
          f"Second, we normalize, because some samples were simply sequenced more deeply. Third, we average marker genes into scores for programs such as CD8 T cells, interferon-gamma and melanoma identity. "
          f"Fourth, we compare the two biopsies, using DESeq2 for single genes and pathway analysis for groups of genes. Take patient 30, a responder: {float(pt30['genes_changed']) * 100:.0f} percent of genes changed at least two-fold, "
          f"the CD8 score rose by {float(pt30['CD8_T_delta']):.1f}, and the melanoma score fell by {abs(float(pt30['Melanocytic_delta'])):.1f}. [~40 s]")
steps = [("Count genes", f"Count RNA reads per gene. {genes:,} genes pass the filter.", "tx2"),
         ("Normalize", "Adjust for sequencing depth so that samples are comparable (DESeq2, VST).", "tx2"),
         ("Score programs", "Average marker genes into 12 scores: T cells, IFN-γ, melanoma and more.", "tx2"),
         ("Change on treatment", "Compare biopsies per patient: DESeq2 for genes, fgsea for pathways.", "accent1")]
for i, (head, body, c) in enumerate(steps):
    x = .6 + i * 3.08
    s.shape(x, 1.55, 2.85, 2.85, "bg2", adj=5000, name=f"Step card {i + 1}")
    s.dot(x + .3, 1.85, .6, c, f"Step dot {i + 1}")
    s.text(x + .3, 1.85, .6, .6, [{"runs": [(str(i + 1), 18, "bg1", True)], "align": "ctr"}], "ctr", f"Step number {i + 1}")
    s.text(x + .3, 2.6, 2.4, .5, [[(head, 18, "tx1", True)]], "ctr", f"Step head {i + 1}")
    s.text(x + .3, 3.1, 2.35, 1.2, [[(body, 14)]], name=f"Step body {i + 1}")
s.text(.6, 4.6, 8, .4, [[("Example: patient 30, a responder", 16, "tx1", True)]], "ctr", "Example header")
ex = [(f"{genes:,}", "genes measured in both biopsies", "tx1"), (f"{float(pt30['genes_changed']) * 100:.0f}%", "of genes changed 2-fold or more", "tx1"),
      (num(pt30["CD8_T_delta"], 1, True), f"CD8 score ({num(pt30['CD8_T_Pre'], 1)} to {num(pt30['CD8_T_On'], 1, True)})", "tx1"),
      (num(pt30["Melanocytic_delta"], 1, True), f"melanoma score ({num(pt30['Melanocytic_Pre'], 1, True)} to {num(pt30['Melanocytic_On'], 1)})", "accent1")]
for i, (v, lab, c) in enumerate(ex):
    x = .6 + i * 3.08
    s.text(x, 5.0, 2.4, .7, [[(v, 36, c, True, True)]], "ctr", f"Example value {i + 1}")
    s.text(x, 5.7, 2.4, .6, [[(lab, 14)]], name=f"Example label {i + 1}")
    if i < 3:
        s.text(x + 2.5, 5.05, .45, .55, [{"runs": [("→", 24, "accent6")], "align": "ctr"}], "ctr", f"Example arrow {i + 1}")
s.number(); slides.append(s)

# 5 Baseline
m = "baseline_immune_score"
s = Slide(3, "Baseline immune activity is a weak predictor",
          f"First: does immune activity before treatment predict response? This is the RNA version of Kwaku's mutation-burden question. Each dot is a patient's T-cell-inflamed score, a published 18-gene measure of how immune-active a tumor is. "
          f"Medians are {med(m, 'PR/CR')} for responders, {med(m, 'SD')} for stable disease and {med(m, 'PD')} for progressors. The difference between responders and progressors is {med(m, 'PR/CR minus PD', 2)}, with a 95 percent bootstrap interval from {ci(m, 'PR/CR minus PD')}, so it crosses zero. "
          f"One cohort detail matches the paper: the four responders who had already progressed on ipilimumab all sit at or above average, the 'hot tumor' pattern the authors describe, but responders without prior ipilimumab do not. The authors also called baseline expression only marginally predictive. [~40 s]")
s.image("01_baseline_immune_by_response.png", .6, 1.6, 8.0, 4.8, "Boxplots of pre-treatment T-cell-inflamed score by response, split by prior ipilimumab")
s.text(9.1, 1.65, 3.6, .6, [[("Median T-cell-inflamed score", 16, "tx1", True)], [("Both cohorts combined", 12, "accent5")]], name="Median header")
stat_rows(s, 9.1, 2.35, [(c, med(m, k), f"{lab}  (n = {n(m, k)})", f"95% CI {ci(m, k)}") for (c, lab), k in zip(GROUPS, KEYS)])
s.text(9.1, 5.3, 3.6, 1.2, [[(f"PR/CR − PD: {med(m, 'PR/CR minus PD', 2, True)} (95% CI {ci(m, 'PR/CR minus PD', 2, True)}). The interval crosses zero.", 14, "tx2", False, False, True)]], name="Baseline takeaway")
s.number(); slides.append(s)

# 6 Genes
m = "share_genes_changed"
s = Slide(3, "Each gene is tracked between the two biopsies",
          f"Next, the two biopsies. Each bar is one patient, and its height is the share of genes that changed at least two-fold: orange went up, blue went down. Responders change most, with a median of {pct(m, 'PR/CR')}, against {pct(m, 'SD')} in stable disease and {pct(m, 'PD')} in progressors. "
          f"That gap is modest, and repeat libraries of the same biopsy already differ by about {noise_med:.0f} percent. Some progressors, like patients 78 and 103, change a lot. The bigger difference is consistency. "
          f"When DESeq2 asks which genes shift together across a group, it finds {n_prcr} in responders, {n_sd} in stable disease and {n_pd} in progressors. A formal test finds {n_int} genes whose change differs between responders and progressors. Responders move in a shared direction; progressors move, if at all, in their own directions. [~45 s]")
s.image("02_gene_changes_per_patient.png", .6, 1.6, 8.6, 4.3, "Stacked bars of the share of genes up or down at least two-fold in each patient, by response")
s.text(9.55, 1.65, 3.2, .6, [[("Genes shifting together", 16, "tx1", True)], [("DESeq2 On vs Pre within group, FDR < 0.1", 12, "accent5")]], name="DE header")
stat_rows(s, 9.55, 2.3, [("accent1", n_prcr, "PR/CR", f"{n('share_genes_changed', 'PR/CR')} patients"), ("accent2", n_sd, "Stable", f"{n('share_genes_changed', 'SD')} patients"), ("accent3", n_pd, "Progressive", f"{n('share_genes_changed', 'PD')} patients")], val_w=1.3, step=.82)
s.text(9.55, 4.85, 3.2, 1.3, [[(f"Median share changing 2-fold: {pct(m, 'PR/CR')}, {pct(m, 'SD')} and {pct(m, 'PD')}. Repeat libraries differ by about {noise_med:.0f}%, so consistency matters more than size.", 12, "tx2", False, False, True)]], name="DE note")
s.number(); slides.append(s)

# 7 Tumor and immune
m, c8 = "tumor_score_change", "cd8_score_change"
s = Slide(3, "Responders lose tumor signal; immune gains are less specific",
          f"So which genes? Here each row is a patient, sorted by the change in a score built from five melanoma identity genes. Responders pile up at the bottom left: median change {med(m, 'PR/CR')} against {med(m, 'PD')} in progressors. The difference is {med(m, 'PR/CR minus PD')}, with an interval of {ci(m, 'PR/CR minus PD')} and a rank-sum p of about {pv(m)}. "
          f"The top single genes in responders tell the same story: SOX10, MLANA, PRAME and PAX3 all fall. On the right is the CD8 T-cell score. It rises in responders, median {med(c8, 'PR/CR', 2, True)}, but also in stable disease, {med(c8, 'SD', 2, True)}, and the responder-versus-progressor test gives p of only about {pv(c8)}. "
          f"Immune activation is clearly happening, but tumor loss is what marks a response. It isn't universal: patients 94 and 101 respond with little change in the tumor score. [~45 s]")
s.text(.6, 1.6, 5.4, .4, [[("Change in tumor score on treatment (median)", 16, "tx1", True)]], "ctr", "Tumor header")
for i, (v, lab, sub, c) in enumerate([(med(m, "PR/CR"), "PR/CR", f"95% CI {ci(m, 'PR/CR')}", "accent1"), (med(m, "PD"), "Progressive", f"95% CI {ci(m, 'PD')}", "accent3")]):
    x = .6 + i * 2.8
    s.shape(x, 2.15, 2.5, 1.6, "bg2", adj=6486, name=f"Tumor card {i + 1}")
    s.text(x + .25, 2.2, 2.1, .85, [[(v, 44, c, True, True)]], "ctr", f"Tumor value {i + 1}")
    s.text(x + .25, 3.05, 2.1, .6, [[(lab, 16)], [(sub, 12, "accent5")]], name=f"Tumor label {i + 1}")
s.text(.6, 4.0, 5.3, 2.4, [
    {"bullet": True, "runs": [(f"PR/CR − PD: {med(m, 'PR/CR minus PD')} (95% CI {ci(m, 'PR/CR minus PD')}); Wilcoxon p ≈ {pv(m)} (n = {n(m, 'PR/CR')} vs {n(m, 'PD')})", 15)]},
    {"bullet": True, "runs": [(f"CD8 change: PR/CR {med(c8, 'PR/CR', 2, True)}, SD {med(c8, 'SD', 2, True)}, PD {med(c8, 'PD', 2)}; Wilcoxon p ≈ {pv(c8)}", 15)]},
    {"bullet": True, "runs": [("Not universal: Pt94 and Pt101 respond with little tumor score change", 15)]}], name="Tumor bullets")
s.image("03_tumor_immune_change_per_patient.png", 6.2, 1.5, 6.5, 5.2, "Per-patient dot plots of tumor score change and CD8 score change, sorted by tumor score change")
s.number(); slides.append(s)

# 8 Pathways
s = Slide(3, "Immune pathways rise broadly; proliferation falls in responders",
          f"Pathways make the same point across thousands of genes. Orange means higher, blue lower, and grey is not significant. In the top block, immune pathways rise on treatment in every group, and stable disease is strongest, with interferon-gamma at {nes('INTERFERON_GAMMA_RESPONSE', 'SD')}. "
          f"In the bottom block, proliferation falls hard in responders, with E2F targets at {nes('E2F_TARGETS', 'PR/CR')}, while most of those programs rise in progressors, at {nes('E2F_TARGETS', 'PD')}. That matches the paper, which reported immune genes up and cell-cycle and melanin genes down in responders. "
          f"The left column shows immune pathways already lean higher in responders before treatment, again as the paper describes. Enrichment p-values run small because genes in a pathway move together, so I read the direction and pattern here more than the exact significance. [~40 s]")
s.image("04_pathways_by_response.png", .6, 1.6, 8.0, 4.8, "Heatmap of normalized enrichment scores for immune and proliferation Hallmark pathways, before treatment and on treatment by response")
s.text(9.1, 1.65, 3.6, .4, [[("Reading the pattern", 16, "tx1", True)]], "ctr", "Pathway header")
for i, (c, head, body) in enumerate([
        ("accent4", "Immune: up in all groups", f"IFN-γ NES {nes('INTERFERON_GAMMA_RESPONSE', 'PR/CR')} PR/CR, {nes('INTERFERON_GAMMA_RESPONSE', 'SD')} SD, {nes('INTERFERON_GAMMA_RESPONSE', 'PD')} PD"),
        ("accent1", "Proliferation: splits groups", f"E2F NES {nes('E2F_TARGETS', 'PR/CR')} PR/CR, {nes('E2F_TARGETS', 'SD')} SD (n.s.), {nes('E2F_TARGETS', 'PD')} PD"),
        ("tx2", "Before treatment", "Immune pathways lean higher in PR/CR, but scores overlap (slide 5)")]):
    yy = 2.25 + i * 1.25
    s.dot(9.1, yy + .05, .3, c, f"Pathway dot {i + 1}")
    s.text(9.55, yy, 3.2, 1.15, [[(head, 15, "tx1", True)], [(body, 14)]], name=f"Pathway text {i + 1}")
s.number(); slides.append(s)

# 9 Cross-modality
s = Slide(3, "RNA agrees with WES and TCR sequencing",
          f"Do the RNA results agree with data measured a completely different way? On the left, each point is a patient with both RNA and exome data. The more a patient's mutations contract in Kwaku's analysis, the more the RNA tumor score falls: rho of {rho(RW)} across {cor[RW]['n']} patients, and still {rho(RW2)} with responders removed. "
          f"On the right, the RNA CD8 score tracks the T-cell fraction measured directly by TCR sequencing in Krithika's data: rho of {rho(RT)}. That validates the score. Both left-panel measures partly reflect how much tumor is in the biopsy, so they aren't fully independent. "
          f"Baseline immune score and mutation burden barely correlate, rho of {rho('Baseline immune score vs mutation burden')}, so RNA and WES carry different information for a combined model. [~40 s]")
s.image("05_rna_vs_wes_tcr.png", .6, 1.6, 8.6, 4.3, "Two scatter plots: change in RNA tumor score vs WES contraction, and RNA CD8 score vs TCR-seq T-cell fraction")
s.text(9.55, 1.65, 3.2, .6, [[("Spearman correlation", 16, "tx1", True)], [("Same patients, independent assays", 12, "accent5")]], name="Correlation header")
for i, (c, v, lab, sub) in enumerate([("accent1", rho(RW), "Tumor loss vs WES contraction", f"n = {cor[RW]['n']}"), ("accent3", rho(RT), "CD8 vs TCR T-cell fraction", f"n = {cor[RT]['n']}")]):
    yy = 2.35 + i * 1.3
    s.text(9.55, yy, 3.2, .65, [[(v, 36, c if c == "accent1" else "tx1", True, True)]], "ctr", f"Correlation value {i + 1}")
    s.text(9.55, yy + .65, 3.2, .55, [[(lab, 14)], [(sub, 12, "accent5")]], name=f"Correlation label {i + 1}")
s.text(9.55, 5.0, 3.2, 1.1, [[(f"Baseline immune score vs mutation burden: rho {rho('Baseline immune score vs mutation burden')}. They carry different information.", 12, "tx2", False, False, True)]], name="Correlation note")
s.number(); slides.append(s)

# 10 Caveats
s = Slide(3, "What could make genes look changed?",
          f"Before the takeaways, what could fool this? First and most important is dilution. Bulk RNA is a mixture, so when tumor cells die, every other cell type becomes a bigger share. When I adjust the responder CD8 rise for tumor loss, it drops from {num(d1['PRCR_effect'], 2, True)} to {num(d2['PRCR_effect'], 2, True)} and the interval crosses zero, so that rise can't be separated from dilution. The on-treatment CD8 level holds up, at {num(d3['PRCR_effect'], 2, True)} before and {num(d4['PRCR_effect'], 2, True)} after adjustment. "
          f"Second, response is defined by tumor shrinkage, so falling melanoma genes partly restate the outcome rather than predict it. Third, I couldn't reproduce the paper's gene counts. With their designs and a current DESeq2, I get {pA['authors_design_this_R']} versus their {pA['paper_q0.2']} before treatment and {pC['authors_design_this_R']} versus {pC['paper_q0.2']} on treatment. Their responder comparison, {int(float(pD['paper_q0.2'])):,} genes, has no public script, and grouping stable disease with responders gives me only {pD['our_pipeline_fdr0.2']}; responders versus progressors gives {n_pap_e}. The biological direction agrees. The paper also removed patient 3 as an outlier; I kept that patient, and the main results hold either way. "
          f"Finally, groups are small, {n('tumor_score_change', 'PR/CR')} responders and {n('tumor_score_change', 'PD')} progressors, at one time point. [~60 s]")
cav = [("Tumor loss dilutes other signals", f"Fewer cancer cells inflate other cell types' share. Adjusting for tumor loss cuts the CD8 rise from {num(d1['PRCR_effect'], 2, True)} to {num(d2['PRCR_effect'], 2, True)}. On-treatment CD8 holds ({num(d3['PRCR_effect'], 2, True)} to {num(d4['PRCR_effect'], 2, True)})."),
       ("Tumor score is close to the outcome", "Response is defined by tumor shrinkage, so falling melanoma genes partly restate it. The score describes response; it is not yet a predictor."),
       ("Paper's gene counts not reproduced", f"Authors' designs in current DESeq2: {pA['authors_design_this_R']} vs {pA['paper_q0.2']} genes (Pre), {pC['authors_design_this_R']} vs {pC['paper_q0.2']} (On vs Pre). Direction agrees. No public script for their {int(float(pD['paper_q0.2'])):,}."),
       ("Small groups, one snapshot", f"{n('tumor_score_change', 'PR/CR')} responders and {n('tumor_score_change', 'PD')} progressors with paired RNA, sampled once at about week 4. Main results hold without Pt3, whom the paper excluded.")]
for i, (head, body) in enumerate(cav):
    x, y = .6 + (i % 2) * 6.15, 1.55 + (i // 2) * 2.5
    s.shape(x, y, 5.95, 2.3, "bg2", adj=5000, name=f"Caveat card {i + 1}")
    s.dot(x + .3, y + .3, .4, "accent4", f"Caveat dot {i + 1}")
    s.text(x + .85, y + .25, 4.9, .5, [[(head, 19, "tx1", True)]], "ctr", f"Caveat head {i + 1}")
    s.text(x + .3, y + .9, 5.35, 1.3, [[(body, 15)]], name=f"Caveat body {i + 1}")
s.number(); slides.append(s)

# 11 Takeaways
c8on = "cd8_score_on"
s = Slide(4, "Takeaways and next steps",
          f"Three takeaways. First, what tracks response on RNA is tumor loss: melanoma and proliferation programs fall in responders, and that agrees with Kwaku's mutation contraction. Second, immune activation alone isn't enough. Immune pathways rise in every group, including patients whose tumors don't shrink, and the baseline score doesn't separate responders. "
          f"The on-treatment CD8 level does, with a difference of {med(c8on, 'PR/CR minus PD', 2, True)}, and it survives the tumor-content adjustment. Next, I'd use a proper deconvolution method to estimate cell types directly, repeat this on the 26 patients our team shares, and combine RNA, WES and TCR in a cross-validated model. Thank you. [~35 s]")
take = [("Tumor loss stands out", "Melanoma and proliferation programs fall in responders, in step with Kwaku's mutation contraction.", "accent1"),
        ("Immunity is not enough", f"Immune pathways rise in all groups. Baseline scores overlap, but on-treatment CD8 separates PR/CR from PD ({med(c8on, 'PR/CR minus PD', 2, True)}).", "accent4"),
        ("Next steps", "Deconvolve cell types, repeat on the 26 patients our team shares, and combine RNA, WES and TCR in a cross-validated model.", "accent3")]
for i, (head, body, c) in enumerate(take):
    x = .6 + i * 4.125
    s.shape(x, 1.75, 3.85, 4.4, "tx2", adj=4000, name=f"Takeaway card {i + 1}")
    s.dot(x + .35, 2.1, .55, c, f"Takeaway dot {i + 1}")
    s.text(x + .35, 2.1, .55, .55, [{"runs": [(str(i + 1), 18, "bg1", True)], "align": "ctr"}], "ctr", f"Takeaway number {i + 1}")
    s.text(x + .35, 2.9, 3.25, .55, [[(head, 21, "bg1", True, True)]], "ctr", f"Takeaway head {i + 1}")
    s.text(x + .35, 3.55, 3.15, 2.4, [[(body, 16, "bg2")]], name=f"Takeaway body {i + 1}")
s.text(.6, 6.45, 6, .5, [[("Thank you", 18, "accent6", False, False, True)]], "ctr", "Thanks")
slides.append(s)


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
secs = sum(int(x) for sl in slides for x in re.findall(r"\[~(\d+) s\]", sl.notes))
lines = [f"# Presenter script: RNA-seq results\n",
         f"Matches `RNA_11slides_Ethan.pptx` (also in its speaker notes). Total about {round(secs / 60)} minutes. Numbers come from `rna/results/`.\n"]
for i, sl in enumerate(slides, 1):
    text, t = re.match(r"(.*) \[(~\d+ s)\]$", sl.notes, re.S).groups()
    lines.append(f"## Slide {i}: {sl.title_text} ({t})\n\n{text}\n")
SCRIPT.write_text("\n".join(lines), encoding="utf8")
print(f"Wrote {OUT.name} ({len(slides)} slides, {secs} s) and {SCRIPT.name}")
