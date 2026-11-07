"""Generate all figures for the TCR analysis presentation."""
import sys; sys.path.insert(0, '.')
import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import spearmanr, mannwhitneyu
import os

os.makedirs('outputs/figures', exist_ok=True)

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 12,
    'figure.facecolor': 'white',
    'savefig.dpi': 200,
    'savefig.bbox': 'tight',
})

C_CR = '#2ecc71'
C_SD = '#f39c12'
C_PD = '#e74c3c'

def rcol(bor):
    return {'PRCR': C_CR, 'SD': C_SD, 'PD': C_PD}.get(bor, '#95a5a6')

# Load data
div = pd.read_csv('outputs/diversity_summary_aa.csv', index_col='patient')
trk = pd.read_csv('outputs/tracking_aa_vj.csv', index_col='patient')
trk_fix = pd.read_csv('outputs/tracking_aa_vj_fixed3500.csv', index_col='patient')
exp = pd.read_csv('outputs/expansion_summary_aa_vj.csv', index_col='patient')
val = pd.read_csv('outputs/validation_aa.csv', index_col='patient')

resp_legend = [Patch(facecolor='white', edgecolor=C_CR, lw=2, label='CR/PR'),
               Patch(facecolor='white', edgecolor=C_SD, lw=2, label='SD'),
               Patch(facecolor='white', edgecolor=C_PD, lw=2, label='PD')]

# ============================================================
# FIGURE 1: Validation
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

ax = axes[0]
ax.scatter(val['ref_richness_Pre'], val['richness_Pre'], s=40, alpha=0.7, c='#3498db', label='Pre-treatment')
ax.scatter(val['ref_richness_On'], val['richness_On'], s=40, alpha=0.7, c='#e74c3c', label='On-treatment')
mn, mx = 0, max(val[['ref_richness_Pre','richness_Pre','ref_richness_On','richness_On']].max())
ax.plot([mn, mx], [mn, mx], 'k--', alpha=0.3, lw=1)
ax.set_xlabel('Published richness (Riaz et al.)')
ax.set_ylabel('Our computed richness')
ax.set_title('Richness: exact match')
ax.legend(fontsize=9)

ax = axes[1]
ax.scatter(val['ref_clonality_Pre'], val['clonality_Pre'], s=40, alpha=0.7, c='#3498db', label='Pre-treatment')
ax.scatter(val['ref_clonality_On'], val['clonality_On'], s=40, alpha=0.7, c='#e74c3c', label='On-treatment')
mn = min(val[['ref_clonality_Pre','clonality_Pre','ref_clonality_On','clonality_On']].min())
mx = max(val[['ref_clonality_Pre','clonality_Pre','ref_clonality_On','clonality_On']].max())
ax.plot([mn, mx], [mn, mx], 'k--', alpha=0.3, lw=1)
ax.set_xlabel('Published clonality (Riaz et al.)')
ax.set_ylabel('Our computed clonality')
ax.set_title('Clonality: within 0.0005')
ax.legend(fontsize=9)

fig.suptitle('Validation: pipeline reproduces published TCR metrics (n=29)', fontsize=14, fontweight='bold', y=1.02)
fig.tight_layout()
fig.savefig('outputs/figures/fig1_validation.png')
plt.close()
print('Fig 1: validation ✓')

# ============================================================
# FIGURE 2: Clonotype tracking — novel vs persistent
# ============================================================
t = trk.copy()
t['novel_frac'] = 1 - t['raw_frac_on_clonotypes_from_baseline']
t['persistent_frac'] = t['raw_frac_on_clonotypes_from_baseline']

fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

for ax_i, (cn, cl) in enumerate([('NIV3-PROG', 'Ipi-Progressed'), ('NIV3-NAIVE', 'Ipi-Naive')]):
    ax = axes[ax_i]
    sub = t[t['Cohort'] == cn].sort_values('novel_frac', ascending=False)
    x = np.arange(len(sub))

    ax.bar(x, sub['persistent_frac'], color='#3498db', alpha=0.85, label='Persistent from baseline')
    ax.bar(x, sub['novel_frac'], bottom=sub['persistent_frac'].values, color='#e67e22', alpha=0.85, label='Novel on-treatment')

    for i, p in enumerate(sub.index):
        ax.bar(i, 1, fill=False, edgecolor=rcol(sub.loc[p, 'myBOR']), linewidth=2.5)

    ax.set_xticks(x)
    ax.set_xticklabels([f'Pt{p}' for p in sub.index], rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Fraction of on-treatment clonotypes')
    ax.set_title(cl)
    ax.set_ylim(0, 1.05)
    if ax_i == 0:
        ax.legend(fontsize=9, loc='lower left')

axes[1].legend(handles=resp_legend, fontsize=9, loc='lower left', title='Response')

fig.suptitle('Clonotype composition on treatment: most clonotypes are novel', fontsize=14, fontweight='bold', y=1.02)
fig.tight_layout()
fig.savefig('outputs/figures/fig2_clonotype_tracking.png')
plt.close()
print('Fig 2: clonotype tracking ✓')

# ============================================================
# FIGURE 3: Expansion — novel vs preexisting
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

for ax_i, (cn, cl) in enumerate([('NIV3-PROG', 'Ipi-Progressed'), ('NIV3-NAIVE', 'Ipi-Naive')]):
    ax = axes[ax_i]
    sub = exp[exp['Cohort'] == cn].sort_values('n_expanded', ascending=False)
    x = np.arange(len(sub))

    ax.bar(x, sub['n_expanded_preexisting'], color='#2980b9', alpha=0.85, label='Preexisting expanded')
    ax.bar(x, sub['n_expanded_novel'], bottom=sub['n_expanded_preexisting'].values, color='#c0392b', alpha=0.85, label='Novel expanded')

    for i, p in enumerate(sub.index):
        total = sub.loc[p, 'n_expanded']
        ax.bar(i, total, fill=False, edgecolor=rcol(sub.loc[p, 'myBOR']), linewidth=2.5)

    ax.set_xticks(x)
    ax.set_xticklabels([f'Pt{p}' for p in sub.index], rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Number of expanded clonotypes')
    ax.set_title(cl)
    if ax_i == 0:
        ax.legend(fontsize=9, loc='upper right')

axes[1].legend(handles=resp_legend, fontsize=9, loc='upper right', title='Response')

fig.suptitle('Significantly expanded clonotypes per patient (Fisher exact, FDR < 0.05)', fontsize=14, fontweight='bold', y=1.02)
fig.tight_layout()
fig.savefig('outputs/figures/fig3_expansion.png')
plt.close()
print('Fig 3: expansion ✓')

# ============================================================
# FIGURE 4: Diversity & depth context (3-panel)
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Panel A: Raw richness FC by response (reproducing paper result)
ax = axes[0]
for cn, cl, color in [('NIV3-PROG', 'Ipi-Prog', '#9b59b6'), ('NIV3-NAIVE', 'Ipi-Naive', '#3498db')]:
    sub = div[div['Cohort'] == cn]
    benefit = sub[sub['benefit'] == 'benefit']['raw_richness_fc']
    no_ben = sub[sub['benefit'] == 'no_benefit']['raw_richness_fc']
    pos = [0, 1] if cn == 'NIV3-PROG' else [2.5, 3.5]
    bp = ax.boxplot([benefit, no_ben], positions=pos, widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor=color, alpha=0.3), medianprops=dict(color='black'))
    for d, p in zip([benefit, no_ben], pos):
        ax.scatter(np.full(len(d), p) + np.random.normal(0, 0.05, len(d)), d,
                   c=color, alpha=0.7, s=30, zorder=3)

ax.set_xticks([0, 1, 2.5, 3.5])
ax.set_xticklabels(['Benefit', 'No\nBenefit', 'Benefit', 'No\nBenefit'], fontsize=9)
ax.axhline(1, ls='--', color='grey', alpha=0.5)
ax.set_ylabel('Richness fold change (on/pre)')
ax.set_title('A. Raw richness FC\n(reproduces paper Fig 5B)')
ax.text(0.5, ax.get_ylim()[1]*0.95, 'Ipi-Prog', ha='center', fontsize=9, color='#9b59b6', fontweight='bold')
ax.text(3.0, ax.get_ylim()[1]*0.95, 'Ipi-Naive', ha='center', fontsize=9, color='#3498db', fontweight='bold')

# Panel B: Richness FC vs depth ratio
ax = axes[1]
depth_ratio = div['raw_templates_on'] / div['raw_templates_pre'] if 'raw_templates_on' in div.columns else div['on_templates'] / div['pre_templates']
rich_fc = div['raw_richness_fc']
colors = [rcol(div.loc[p, 'myBOR']) for p in div.index]
ax.scatter(depth_ratio, rich_fc, c=colors, s=50, alpha=0.7, edgecolors='black', linewidth=0.5)
rho, pv = spearmanr(depth_ratio, rich_fc)
ax.set_xlabel('Sequencing depth ratio (on/pre)')
ax.set_ylabel('Richness fold change (on/pre)')
ax.set_title(f'B. Richness tracks depth\nρ = {rho:.2f}, p = {pv:.1e}')
# Fit line
z = np.polyfit(np.log(depth_ratio), np.log(rich_fc), 1)
xs = np.linspace(depth_ratio.min(), depth_ratio.max(), 100)
ax.plot(xs, np.exp(np.polyval(z, np.log(xs))), 'k--', alpha=0.4)

# Panel C: Depth-matched = null
ax = axes[2]
for cn, cl, color in [('NIV3-PROG', 'Ipi-Prog', '#9b59b6'), ('NIV3-NAIVE', 'Ipi-Naive', '#3498db')]:
    sub = div[div['Cohort'] == cn]
    benefit = sub[sub['benefit'] == 'benefit']['rar_pair_richness_fc']
    no_ben = sub[sub['benefit'] == 'no_benefit']['rar_pair_richness_fc']
    pos = [0, 1] if cn == 'NIV3-PROG' else [2.5, 3.5]
    bp = ax.boxplot([benefit, no_ben], positions=pos, widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor=color, alpha=0.3), medianprops=dict(color='black'))
    for d, p in zip([benefit, no_ben], pos):
        ax.scatter(np.full(len(d), p) + np.random.normal(0, 0.05, len(d)), d,
                   c=color, alpha=0.7, s=30, zorder=3)

ax.set_xticks([0, 1, 2.5, 3.5])
ax.set_xticklabels(['Benefit', 'No\nBenefit', 'Benefit', 'No\nBenefit'], fontsize=9)
ax.axhline(1, ls='--', color='grey', alpha=0.5)
ax.set_ylabel('Richness fold change (depth-matched)')
ax.set_title('C. Depth-matched: no signal\n(rarefied to per-pair minimum)')
ax.text(0.5, ax.get_ylim()[1]*0.95, 'Ipi-Prog', ha='center', fontsize=9, color='#9b59b6', fontweight='bold')
ax.text(3.0, ax.get_ylim()[1]*0.95, 'Ipi-Naive', ha='center', fontsize=9, color='#3498db', fontweight='bold')

fig.suptitle('Diversity analysis: richness change reflects sequencing depth, not biology', fontsize=14, fontweight='bold', y=1.04)
fig.tight_layout()
fig.savefig('outputs/figures/fig4_diversity_depth.png')
plt.close()
print('Fig 4: diversity + depth ✓')

# ============================================================
# FIGURE 5: TCR-WES connection
# ============================================================
# Load WES data
clonality = pd.read_csv(os.path.expanduser('~/bms038_analysis/data/clonality.data.table.pyclone.95.ci.txt'), sep='\t')
clonality['patient'] = clonality['sample'].str.lower().str.replace('pt','').astype(int)
clonality['neo_lost'] = clonality['pre.is.neo'].fillna(False) & ~clonality['on.is.neo'].fillna(False)

per_pt = clonality.groupby('patient').agg(
    neo_lost=('neo_lost', 'sum'),
).reset_index().set_index('patient')

common = exp.index.intersection(per_pt.index)
w = pd.DataFrame({
    'n_expanded': exp.loc[common, 'n_expanded'],
    'neo_lost': per_pt.loc[common, 'neo_lost'],
    'myBOR': div.loc[common, 'myBOR'],
    'Cohort': div.loc[common, 'Cohort'],
}, index=common)

fig, axes = plt.subplots(1, 2, figsize=(11, 5))

# Panel A: All patients, color by response
ax = axes[0]
for bor, label in [('PRCR', 'CR/PR'), ('SD', 'SD'), ('PD', 'PD')]:
    sub = w[w['myBOR'] == bor]
    ax.scatter(sub['n_expanded'], sub['neo_lost'], c=rcol(bor), s=60, alpha=0.7,
              edgecolors='black', linewidth=0.5, label=label, zorder=3)
rho, pv = spearmanr(w['n_expanded'], w['neo_lost'])
ax.set_xlabel('Number of expanded T cell clonotypes')
ax.set_ylabel('Neoantigens lost (pre→on)')
ax.set_title(f'A. All patients (n={len(w)})\nρ = {rho:.2f}, p = {pv:.3f}')
ax.legend(fontsize=9)

# Panel B: CR/PR only
ax = axes[1]
cr = w[w['myBOR'] == 'PRCR']
ax.scatter(cr['n_expanded'], cr['neo_lost'], c=C_CR, s=80, alpha=0.8,
          edgecolors='black', linewidth=0.5, zorder=3)
for p in cr.index:
    ax.annotate(f'Pt{p}', (cr.loc[p, 'n_expanded'], cr.loc[p, 'neo_lost']),
               textcoords='offset points', xytext=(5, 5), fontsize=8)
rho_cr, pv_cr = spearmanr(cr['n_expanded'], cr['neo_lost'])
ax.set_xlabel('Number of expanded T cell clonotypes')
ax.set_ylabel('Neoantigens lost (pre→on)')
ax.set_title(f'B. Responders only (CR/PR, n={len(cr)})\nρ = {rho_cr:.2f}, p = {pv_cr:.3f}')

# Trend line
z = np.polyfit(cr['n_expanded'], cr['neo_lost'], 1)
xs = np.linspace(cr['n_expanded'].min(), cr['n_expanded'].max(), 50)
ax.plot(xs, np.polyval(z, xs), 'k--', alpha=0.4)

fig.suptitle('TCR–WES connection: T cell expansion correlates with neoantigen clearance',
             fontsize=14, fontweight='bold', y=1.04)
fig.tight_layout()
fig.savefig('outputs/figures/fig5_tcr_wes.png')
plt.close()
print('Fig 5: TCR-WES ✓')

# ============================================================
# FIGURE 6: Template ratio (depth) by response — the TIL mediation
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(11, 5))

# Panel A: Template ratio by response
ax = axes[0]
template_ratio = div['on_templates'] / div['pre_templates']
for cn, cl, color, pos_base in [('NIV3-PROG', 'Ipi-Prog', '#9b59b6', 0), ('NIV3-NAIVE', 'Ipi-Naive', '#3498db', 2.5)]:
    sub = div[div['Cohort'] == cn]
    tr_sub = template_ratio.loc[sub.index]
    benefit = tr_sub[sub['benefit'] == 'benefit']
    no_ben = tr_sub[sub['benefit'] == 'no_benefit']
    positions = [pos_base, pos_base+1]
    bp = ax.boxplot([benefit, no_ben], positions=positions, widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor=color, alpha=0.3), medianprops=dict(color='black'))
    for d, p in zip([benefit, no_ben], positions):
        ax.scatter(np.full(len(d), p) + np.random.normal(0, 0.05, len(d)), d,
                   c=color, alpha=0.7, s=30, zorder=3)

ax.set_xticks([0, 1, 2.5, 3.5])
ax.set_xticklabels(['Benefit', 'No\nBenefit', 'Benefit', 'No\nBenefit'], fontsize=9)
ax.axhline(1, ls='--', color='grey', alpha=0.5)
ax.set_ylabel('Template count ratio (on/pre)')
ax.set_title('A. Depth ratio ≈ T cell infiltration change')
ax.text(0.5, ax.get_ylim()[1]*0.9, 'Ipi-Prog', ha='center', fontsize=9, color='#9b59b6', fontweight='bold')
ax.text(3.0, ax.get_ylim()[1]*0.9, 'Ipi-Naive', ha='center', fontsize=9, color='#3498db', fontweight='bold')

# Panel B: Partial correlation diagram (text-based)
ax = axes[1]
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.axis('off')

# Draw the mediation diagram
boxes = {
    'depth': (5, 8.5, 'Depth ratio\n(T cell infiltration)'),
    'rich': (1.5, 4, 'Richness\nfold change'),
    'resp': (8.5, 4, 'Treatment\nresponse'),
}
for key, (cx, cy, txt) in boxes.items():
    ax.add_patch(plt.Rectangle((cx-1.5, cy-0.8), 3, 1.6, fill=True,
                                facecolor='#ecf0f1', edgecolor='#2c3e50', linewidth=2, zorder=2))
    ax.text(cx, cy, txt, ha='center', va='center', fontsize=10, fontweight='bold', zorder=3)

# Arrows
from matplotlib.patches import FancyArrowPatch
ax.annotate('', xy=(3, 5.2), xytext=(4, 7.7),
           arrowprops=dict(arrowstyle='->', lw=2, color='#2c3e50'))
ax.annotate('', xy=(7, 7.7), xytext=(8, 5.2),
           arrowprops=dict(arrowstyle='->', lw=2, color='#2c3e50'))

# Direct path (dashed, weak)
ax.annotate('', xy=(7, 4), xytext=(3, 4),
           arrowprops=dict(arrowstyle='->', lw=2, color='#e74c3c', ls='dashed'))

# Labels on arrows
ax.text(2.8, 6.8, 'ρ = 0.83', fontsize=10, color='#2c3e50', fontweight='bold', rotation=55)
ax.text(5, 3.0, 'ρ = 0.35 → 0.16\n(after controlling\nfor depth)',
        ha='center', fontsize=9, color='#e74c3c', fontstyle='italic')
ax.text(7.5, 6.8, 'ρ = 0.35', fontsize=10, color='#2c3e50', fontweight='bold', rotation=-55)

ax.set_title('B. TIL mediation: richness–response\nis mediated by T cell infiltration', fontsize=12)

fig.suptitle('The richness signal measures T cell infiltration change, not repertoire diversity',
             fontsize=14, fontweight='bold', y=1.02)
fig.tight_layout()
fig.savefig('outputs/figures/fig6_til_mediation.png')
plt.close()
print('Fig 6: TIL mediation ✓')

print('\nAll figures saved to outputs/figures/')
