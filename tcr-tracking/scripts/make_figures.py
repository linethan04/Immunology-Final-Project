"""Make slide-ready figures from the CSVs in outputs/ (no raw data needed).

Usage: python scripts/make_figures.py        ->  figures/fig1..fig4 (.png)
Colors: blue = benefit (CR/PR/SD), orange = no benefit (PD); markers differ too, so colour is not the only cue.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT, FIG = Path("outputs"), Path("figures")
FIG.mkdir(exist_ok=True)
C = {"benefit": "#2a78d6", "no_benefit": "#eb6834"}
M = {"benefit": "o", "no_benefit": "^"}
LAB = {"benefit": "Benefit (CR/PR/SD)", "no_benefit": "No benefit (PD)"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
plt.rcParams.update({"font.size": 12, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "axes.axisbelow": True,
                     "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb"})
COH = {"NIV3-NAIVE": "Ipi-naive", "NIV3-PROG": "Ipi-progressed"}


def legend(ax, **kw):
    h = [plt.Line2D([], [], marker=M[k], color=C[k], ls="", ms=8, label=LAB[k]) for k in C]
    ax.legend(handles=h, frameon=False, **kw)


def save(fig, name):
    fig.tight_layout()
    fig.savefig(FIG / f"{name}.png", dpi=200)
    plt.close(fig)
    print("wrote", FIG / f"{name}.png")


div = pd.read_csv(OUT / "diversity_summary_aa.csv")
tests = pd.read_csv(OUT / "diversity_group_tests_aa.csv")
trk = pd.read_csv(OUT / "tracking_aa_vj.csv")
exp = pd.read_csv(OUT / "expansion_summary_aa_vj.csv")

# 1. raw richness change tracks sequencing depth; rarefied does not
fig, axs = plt.subplots(1, 3, figsize=(14, 4.6), sharey=True)
for ax, v, t in zip(axs, ["raw", "rar_pair", "rar_fixed"], ["Raw", "Equal depth within pair", "Fixed 3,500 templates"]):
    d = div.dropna(subset=[f"{v}_richness_fc"])
    for k in C:
        s = d[d.benefit == k]
        ax.scatter(s.log2_depth_ratio_on_over_pre, np.log2(s[f"{v}_richness_fc"]), c=C[k], marker=M[k], s=45,
                   edgecolor="#fcfcfb", linewidth=1)
    ax.axhline(0, color=MUTED, lw=1)
    ax.set_title(t, loc="left", color=INK)
    ax.set_xlabel("log2 depth ratio (on / pre templates)")
axs[0].set_ylabel("log2 richness fold change (on / pre)")
legend(axs[0], loc="upper left", fontsize=10)
save(fig, "fig1_depth_confound")

# 2. richness change by response group, three depth treatments
fig, axs = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
rng = np.random.default_rng(0)
for ax, coh in zip(axs, COH):
    pos, ticks = 0, []
    for v, vl in [("raw", "Raw"), ("rar_pair", "Equal depth\nin pair"), ("rar_fixed", "Fixed\n3,500")]:
        p = tests[(tests.version == v) & (tests.metric == "richness_fc") & (tests.cohort == coh)].p_mannwhitney.iloc[0]
        for i, k in enumerate(C):
            s = div[(div.Cohort == coh) & (div.benefit == k)][f"{v}_richness_fc"].dropna()
            x = pos + i * 0.7 + rng.uniform(-0.12, 0.12, len(s))
            ax.scatter(x, np.log2(s), c=C[k], marker=M[k], s=50, edgecolor="#fcfcfb", linewidth=1)
            ax.hlines(np.log2(s.median()), pos + i * 0.7 - 0.25, pos + i * 0.7 + 0.25, color=INK, lw=2)
        ax.text(pos + 0.35, 1.02, f"p = {p:.3f}", transform=ax.get_xaxis_transform(), ha="center", color=MUTED, fontsize=10)
        ticks.append((pos + 0.35, vl))
        pos += 2
    ax.axhline(0, color=MUTED, lw=1)
    ax.set_xticks([t[0] for t in ticks], [t[1] for t in ticks])
    ax.set_title(COH[coh], loc="left", pad=22)
axs[0].set_ylabel("log2 richness fold change (on / pre)")
legend(axs[0], loc="lower left", fontsize=10)
save(fig, "fig2_richness_by_response")

# 3. replacement: fraction of on-treatment templates from baseline clones
t = trk.sort_values("raw_frac_on_templates_from_baseline").reset_index(drop=True)
fig, ax = plt.subplots(figsize=(11, 4.8))
x = np.arange(len(t))
ax.vlines(x, t.raw_frac_on_templates_from_baseline, t.rar_frac_on_templates_from_baseline, color=GRID, lw=3)
ax.scatter(x, t.raw_frac_on_templates_from_baseline, s=55, c="#52514e", marker="s", label="Raw", zorder=3)
for k in C:
    s = t[t.benefit == k]
    ax.scatter(s.index, s.rar_frac_on_templates_from_baseline, s=60, c=C[k], marker=M[k], zorder=4,
               edgecolor="#fcfcfb", label=f"Equal depth: {LAB[k]}")
ax.set_xticks(x, [f"Pt{p}" for p in t.patient], rotation=90, fontsize=9)
ax.set_ylabel("Fraction of on-treatment templates\nfrom baseline clones (aa+V+J)")
ax.set_ylim(0, 1)
ax.legend(frameon=False, fontsize=10, loc="upper left", ncol=3)
save(fig, "fig3_replacement")

# 4. expansion: pre-existing vs novel expanded clones, and the depth-matched count
e = exp.sort_values("n_expanded").reset_index(drop=True)
fig, axs = plt.subplots(1, 2, figsize=(14, 4.8), gridspec_kw={"width_ratios": [1.6, 1]})
x = np.arange(len(e))
axs[0].bar(x, e.n_expanded_preexisting, color="#2a78d6", label="Seen at baseline", width=0.75, edgecolor="#fcfcfb", linewidth=1.2)
axs[0].bar(x, e.n_expanded_novel, bottom=e.n_expanded_preexisting, color="#eda100", label="Novel (absent at baseline)", width=0.75, edgecolor="#fcfcfb", linewidth=1.2)
axs[0].set_xticks(x, [f"Pt{p}" for p in e.patient], rotation=90, fontsize=9)
axs[0].set_ylabel("Significantly expanded clonotypes (FDR < 0.05)")
axs[0].set_title("Full depth", loc="left")
axs[0].legend(frameon=False, fontsize=10, loc="upper left")
axs[1].scatter(e.n_expanded.clip(lower=0.5), e.rar_n_expanded_mean.clip(lower=0.5), c="#52514e", s=45)
axs[1].set_xscale("log"); axs[1].set_yscale("log")
axs[1].plot([0.5, 1000], [0.5, 1000], color=MUTED, lw=1, ls="--")
axs[1].text(0.6, 600, "dashed line: equal", color=MUTED, fontsize=10)
axs[1].set_xlabel("Expanded clonotypes, full depth (log)")
axs[1].set_ylabel("Expanded clonotypes, both samples\nat 3,500 templates (log)")
axs[1].set_title("Depth matched", loc="left")
save(fig, "fig4_expansion")
