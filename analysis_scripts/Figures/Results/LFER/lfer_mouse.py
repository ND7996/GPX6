# =================== MOUSE CODE WITH UNIFIED DISTANCE SCALE ===================
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# =================== JCIM PUBLICATION STYLE ===================
ACS_PATH = r"D:\PhD_Thesis\analysis\FINAL_PUBLICATION_FIGURES"
if ACS_PATH not in sys.path:
    sys.path.append(ACS_PATH)

from acsfonts import set_jcim_style, jcim_figure, save_jcim_figure
set_jcim_style()

# =================== DATA LOADING ===================
input_csv = Path(r"D:\PhD_Thesis\analysis\dipole\mouse_HBONDS.csv")
df = pd.read_csv(input_csv)

# =================== OUTPUT FOLDER ===================
outdir = Path(r"D:\PhD_Thesis\analysis\FINAL_PUBLICATION_FIGURES\Results\LFER")
outdir.mkdir(parents=True, exist_ok=True)

for c in ['dg_star', 'dg_star_error', 'dg0', 'dg0_error']:
    df[c] = pd.to_numeric(df[c], errors='coerce')

df['level_num'] = df['level'].str.replace('Level', '').astype(int)

# =================== SEM HELPER ===================
def sem(series):
    n = series.dropna().shape[0]
    return np.nan if n < 2 else series.std(ddof=1) / np.sqrt(n)

# =================== ΔΔG CALCULATION ===================
ref = df[df['mutation'] == 'C49U'].set_index('level')
fallback = ref.iloc[0]

all_data = []
for level in df['level'].unique():
    sub = df[df['level'] == level].copy()
    r = ref.loc[level] if level in ref.index else fallback

    sub['ΔΔG‡'] = sub['dg_star'] - r['dg_star']
    sub['ΔΔG⁰'] = sub['dg0'] - r['dg0']
    sub['ΔΔG‡_err'] = np.sqrt(sub['dg_star_error']**2 + r['dg_star_error']**2)
    sub['ΔΔG⁰_err'] = np.sqrt(sub['dg0_error']**2 + r['dg0_error']**2)

    all_data.append(sub)

dfp = pd.concat(all_data)
dfp = dfp[dfp['mutation'] != 'C49U'].copy()

# =================== DISTANCE COLUMN ===================
dist_col = None
possible_names = [
    'distance_to_Cys49', 'dist_to_Cys49', 'distance_Cys49',
    'dist_C49', 'Distance_to_Cys49', 'dist_to_Sec/Cys49'
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
    ΔΔG0_mean=('ΔΔG⁰', 'mean'),
    ΔΔG0_sem=('ΔΔG⁰', sem),
    ΔΔGddagger_mean=('ΔΔG‡', 'mean'),
    ΔΔGddagger_sem=('ΔΔG‡', sem),
    ΔΔGddagger_err_mean=('ΔΔG‡_err', 'mean'),
).reset_index()

if dist_col:
    means_per_mut[dist_col] = dfp.groupby('mutation')[dist_col].mean().values

# =================== LINEAR REGRESSION ===================
mask_reg = dfp[['ΔΔG⁰', 'ΔΔG‡']].notna().all(axis=1)
x_all = dfp.loc[mask_reg, 'ΔΔG⁰'].values
y_all = dfp.loc[mask_reg, 'ΔΔG‡'].values

slope, intercept, r_value, p_value, std_err = stats.linregress(x_all, y_all)

# =================== FIX: use fixed axis limits, not data range ===================
x_fit = np.linspace(-5, 15, 400)
y_fit = slope * x_fit + intercept

n_mutations = len(means_per_mut)

# =================== UNIFIED DISTANCE SCALE (SAME FOR HUMAN AND MOUSE) ===================
UNIFIED_DISTANCE_MIN = 3.0
UNIFIED_DISTANCE_MAX = 25.0

# =================== FIGURE ===================
fig = plt.figure(figsize=(14, 6))
gs = fig.add_gridspec(1, 2, width_ratios=[2.5, 1.2], wspace=0.05)

ax_main = fig.add_subplot(gs[0])
ax_legend = fig.add_subplot(gs[1])

# =================== COLOR MAP WITH UNIFIED SCALE ===================
if dist_col:
    mask_valid = dfp[['ΔΔG⁰', 'ΔΔG‡', dist_col]].notna().all(axis=1)
    dfp_valid = dfp[mask_valid].copy()
    norm = plt.Normalize(vmin=UNIFIED_DISTANCE_MIN, vmax=UNIFIED_DISTANCE_MAX)
    cmap_v = plt.colormaps['RdBu_r']
else:
    dfp_valid = dfp[dfp[['ΔΔG⁰', 'ΔΔG‡']].notna().all(axis=1)].copy()

# =================== MAIN SCATTER ===================
ax_main.scatter(
    dfp_valid['ΔΔG⁰'],
    dfp_valid['ΔΔG‡'],
    color='gray',
    s=15,
    alpha=0.35,
    zorder=2
)

ax_main.plot(x_fit, y_fit, color='#d62728', lw=2, zorder=3)

marker_size = 55

if dist_col and dist_col in means_per_mut.columns:
    sc = ax_main.scatter(
        means_per_mut['ΔΔG0_mean'],
        means_per_mut['ΔΔGddagger_mean'],
        c=means_per_mut[dist_col],
        cmap=cmap_v,
        norm=norm,
        s=marker_size,
        edgecolors='black',
        linewidth=0.8,
        zorder=6
    )
else:
    sc = ax_main.scatter(
        means_per_mut['ΔΔG0_mean'],
        means_per_mut['ΔΔGddagger_mean'],
        color='#4c72b0',
        s=marker_size,
        edgecolors='black',
        linewidth=0.8,
        zorder=6
    )

ax_main.errorbar(
    means_per_mut['ΔΔG0_mean'],
    means_per_mut['ΔΔGddagger_mean'],
    xerr=means_per_mut['ΔΔG0_sem'],
    yerr=means_per_mut['ΔΔGddagger_err_mean'],
    fmt='none',
    ecolor='black',
    alpha=0.6,
    capsize=3,
    zorder=5
)

ax_main.text(
    0.04, 0.96,
    f'$R^2$ = {r_value**2:.2f}',
    transform=ax_main.transAxes,
    va='top',
    bbox=dict(facecolor='white', edgecolor='black', alpha=0.9)
)

ax_main.set_xlabel(r'$\Delta\Delta G^{\circ}$ (kcal mol$^{-1}$)')
ax_main.set_ylabel(r'$\Delta\Delta G^{\ddagger}$ (kcal mol$^{-1}$)')

# =================== FIXED SCALE ===================
ax_main.set_xlim(-5, 15)
ax_main.set_ylim(-5, 12.5)

# =================== COLORBAR WITH UNIFIED SCALE ===================
if dist_col:
    cbar = plt.colorbar(sc, ax=ax_main, pad=0.02, aspect=25, shrink=0.85)
    cbar.set_label('Distance to residue 49 (Å)', rotation=270, labelpad=15)

# =================== RIGHT PANEL - EMPTY ===================
ax_legend.axis('off')

# =================== MAIN LEGEND ===================
legend_elements = [
    plt.scatter([], [], color='gray', edgecolors='none', s=15, alpha=0.40, label='Per-level data'),
    plt.scatter([], [], color='white', edgecolors='black', s=marker_size, linewidth=0.8,
                label='Mutation mean $\\pm$ SEM'),
    plt.Line2D([0], [0], color='#d62728', lw=2, label='Linear fit'),
]
ax_main.legend(handles=legend_elements, loc='lower right', framealpha=0.9,
               edgecolor='black', frameon=True, handlelength=1.5, borderpad=0.5, fontsize=8)

# =================== SAVE ===================
fig.savefig(outdir / "mouse_Figure1_LFER.png", dpi=600, bbox_inches='tight')
fig.savefig(outdir / "mouse_Figure1_LFER.pdf", dpi=600, bbox_inches='tight')

plt.show()

print(f"Saved: {outdir / 'mouse_Figure1_LFER.png'}")
print(f"\nLFER regression:")
print(f"  α (slope)   = {slope:.3f} ± {std_err:.3f}")
print(f"  R²          = {r_value**2:.3f}")
print(f"  p-value     = {p_value:.3e}")
print(f"  n (points)  = {len(x_all)}")
print(f"  Number of mutations: {n_mutations}")