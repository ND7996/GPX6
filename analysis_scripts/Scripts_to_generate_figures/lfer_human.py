# =================== HUMAN CODE WITH UNIFIED DISTANCE SCALE ===================
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from pathlib import Path
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# =================== JCIM PUBLICATION STYLE ===================
ACS_PATH = r"./analysis_scripts/Scripts_to_generate_figures/Figures"
if ACS_PATH not in sys.path:
    sys.path.append(ACS_PATH)

from acsfonts import set_jcim_style, jcim_figure, save_jcim_figure

set_jcim_style()

# =================== DATA LOADING ===================
input_csv = Path(r"./analysis_scripts\human_HBONDS.csv")
df = pd.read_csv(input_csv)

# =================== OUTPUT FOLDER ===================
outdir = Path(r"./analysis_scripts/Scripts_to_generate_figures/Figures/Results\LFER")
outdir.mkdir(parents=True, exist_ok=True)

for c in ['dg_star', 'dg_star_error', 'dg0', 'dg0_error']:
    df[c] = pd.to_numeric(df[c], errors='coerce')

df['level_num'] = df['level'].str.replace('Level', '').astype(int)

# =================== SEM HELPER ===================
def sem(series):
    n = series.dropna().shape[0]
    return np.nan if n < 2 else series.std(ddof=1) / np.sqrt(n)

# =================== Î”Î”G CALCULATION ===================
ref = df[df['mutation'] == 'humansec'].set_index('level')
if len(ref) == 0:
    raise ValueError("No 'humansec' reference found!")
fallback = ref.iloc[0]

all_data = []
for level in df['level'].unique():
    sub = df[df['level'] == level].copy()
    r = ref.loc[level] if level in ref.index else fallback
    sub['Î”Î”Gâ€¡'] = sub['dg_star'] - r['dg_star']
    sub['Î”Î”Gâ°'] = sub['dg0'] - r['dg0']
    sub['Î”Î”Gâ€¡_err'] = np.sqrt(sub['dg_star_error']**2 + r['dg_star_error']**2)
    sub['Î”Î”Gâ°_err'] = np.sqrt(sub['dg0_error']**2 + r['dg0_error']**2)
    all_data.append(sub)

dfp = pd.concat(all_data)
dfp = dfp[~dfp['mutation'].isin(['humansec', 'humancys'])].copy()

# =================== DISTANCE COLUMN ===================
dist_col = None
possible_names = [
    'distance_to_Sec49', 'dist_to_Sec49', 'distance_Sec49', 'dist_Sec49',
    'Distance_to_Sec49', 'dist_to_Sec/Cys49', 'dist_to_Cys49',
    'distance_to_Cys49', 'dist_to_Cys49'
]
for col in possible_names:
    if col in dfp.columns:
        dist_col = col
        break
print(f"Found distance column: '{dist_col}'" if dist_col else "No distance column found")

if dist_col:
    dfp[dist_col] = pd.to_numeric(dfp[dist_col], errors='coerce')

# =================== PER-MUTATION MEANS ===================
means_per_mut = dfp.groupby('mutation').agg(
    Î”Î”G0_mean=('Î”Î”Gâ°', 'mean'),
    Î”Î”G0_sem=('Î”Î”Gâ°', sem),
    Î”Î”Gddagger_mean=('Î”Î”Gâ€¡', 'mean'),
    Î”Î”Gddagger_sem=('Î”Î”Gâ€¡', sem),
    Î”Î”Gddagger_err_mean=('Î”Î”Gâ€¡_err', 'mean'),
).reset_index()

if dist_col:
    means_per_mut[dist_col] = dfp.groupby('mutation')[dist_col].mean().values

# =================== LINEAR REGRESSION (using ALL raw data points) ===================
mask_reg = dfp[['Î”Î”Gâ°', 'Î”Î”Gâ€¡']].notna().all(axis=1)
x_all = dfp.loc[mask_reg, 'Î”Î”Gâ°'].values
y_all = dfp.loc[mask_reg, 'Î”Î”Gâ€¡'].values

slope, intercept, r_value, p_value, std_err = stats.linregress(x_all, y_all)

x_fit = np.linspace(x_all.min() - 0.5, x_all.max() + 0.5, 400)
y_fit = slope * x_fit + intercept

n_pts = len(x_all)

# Count mutations
n_mutations = len(means_per_mut)
print(f"\nLFER regression:")
print(f"  Î± (slope)   = {slope:.3f} Â± {std_err:.3f}")
print(f"  intercept   = {intercept:.3f}")
print(f"  RÂ²          = {r_value**2:.3f}")
print(f"  p-value     = {p_value:.3e}")
print(f"  n (points)  = {n_pts}")
print(f"  Number of mutations: {n_mutations}")

# =================== UNIFIED DISTANCE SCALE (SAME FOR HUMAN AND MOUSE) ===================
# Set these to the same values for both plots!
UNIFIED_DISTANCE_MIN = 3.0  # Adjust based on your data
UNIFIED_DISTANCE_MAX = 25.0  # Adjust based on your data

# ===========================================================================
# FIGURE 1 â€” MAIN LFER
# ===========================================================================
fig = plt.figure(figsize=(14, 6))
gs = fig.add_gridspec(1, 2, width_ratios=[2.5, 1.2], wspace=0.05)

ax_main = fig.add_subplot(gs[0])
ax_legend = fig.add_subplot(gs[1])

# Colourmap setup - USE RdBu_r WITH UNIFIED SCALE
if dist_col:
    mask_valid = dfp[['Î”Î”Gâ°', 'Î”Î”Gâ€¡', dist_col]].notna().all(axis=1)
    dfp_valid = dfp[mask_valid].copy()
    # USE UNIFIED MIN/MAX FOR DISTANCE
    norm = plt.Normalize(vmin=UNIFIED_DISTANCE_MIN, vmax=UNIFIED_DISTANCE_MAX)
    cmap_v = plt.colormaps['RdBu_r']
else:
    dfp_valid = dfp[dfp[['Î”Î”Gâ°', 'Î”Î”Gâ€¡']].notna().all(axis=1)].copy()

# Raw per-level scatter - grey dots (ALL data points)
ax_main.scatter(dfp_valid['Î”Î”Gâ°'], dfp_valid['Î”Î”Gâ€¡'],
                color='gray', edgecolors='none', s=15, alpha=0.35, zorder=2)

# Regression line
ax_main.plot(x_fit, y_fit, color='#d62728', lw=2, zorder=3, alpha=0.85)

# Per-mutation mean markers - colored by distance with UNIFIED scale
marker_size = 55
if dist_col and dist_col in means_per_mut.columns:
    sc = ax_main.scatter(means_per_mut['Î”Î”G0_mean'], means_per_mut['Î”Î”Gddagger_mean'],
                         c=means_per_mut[dist_col], cmap='RdBu_r', norm=norm,
                         s=marker_size, edgecolors='black', linewidth=0.8, alpha=0.95, zorder=6)
else:
    sc = ax_main.scatter(means_per_mut['Î”Î”G0_mean'], means_per_mut['Î”Î”Gddagger_mean'],
                         color='#4c72b0', s=marker_size, edgecolors='black', linewidth=0.8, alpha=0.95, zorder=6)

# Error bars
ax_main.errorbar(means_per_mut['Î”Î”G0_mean'], means_per_mut['Î”Î”Gddagger_mean'],
                 xerr=means_per_mut['Î”Î”G0_sem'], yerr=means_per_mut['Î”Î”Gddagger_err_mean'],
                 fmt='none', elinewidth=1, capsize=3, capthick=1,
                 color='black', alpha=0.65, zorder=5)

# RÂ² annotation
ax_main.text(0.04, 0.96, f'$R^2$ = {r_value**2:.2f}',
             transform=ax_main.transAxes, va='top', ha='left',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                       edgecolor='black', linewidth=0.8, alpha=0.95))

ax_main.set_xlabel(r'$\Delta\Delta G^{\circ}$ (kcal mol$^{-1}$)')
ax_main.set_ylabel(r'$\Delta\Delta G^{\ddagger}$ (kcal mol$^{-1}$)')

# Axis limits matching mouse
ax_main.set_xlim(-5, 15)
ax_main.set_ylim(-5, 12.5)

# Colorbar with UNIFIED scale
if dist_col:
    cbar = plt.colorbar(sc, ax=ax_main, pad=0.02, aspect=25, shrink=0.85)
    cbar.set_label('Distance to residue 49 (Ã…)', rotation=270, labelpad=15)

# =================== RIGHT PANEL - COMPLETELY OFF ===================
ax_legend.axis('off')

# Main legend
legend_elements = [
    plt.scatter([], [], color='gray', edgecolors='none', s=15, alpha=0.40, label='Per-level data'),
    plt.scatter([], [], color='white', edgecolors='black', s=marker_size, linewidth=0.8,
                label='Mutation mean $\\pm$ SEM'),
    plt.Line2D([0], [0], color='#d62728', lw=2, label='Linear fit'),
]
ax_main.legend(handles=legend_elements, loc='lower right', framealpha=0.9,
               edgecolor='black', frameon=True, handlelength=1.5, borderpad=0.5, fontsize=8)

# Save
fig.savefig(outdir / "human_Figure1_LFER.png", dpi=600, bbox_inches='tight')
fig.savefig(outdir / "human_Figure1_LFER.pdf", dpi=600, bbox_inches='tight')
plt.show()
print(f"Saved: {outdir / 'human_Figure1_LFER.png'}")

# [REST OF HUMAN FIGURES 2-5 remain the same...]

