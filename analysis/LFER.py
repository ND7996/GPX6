import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import seaborn as sns

# =================== PUBLICATION STYLE ===================
plt.rcParams.update({
    'font.size': 14,
    'axes.labelsize': 18,
    'axes.titlesize': 20,
    'xtick.labelsize': 15,
    'ytick.labelsize': 15,
    'legend.fontsize': 16,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})
sns.set_style("whitegrid")

# =================== DATA LOADING ===================
df = pd.read_csv(r"D:\PhD_Thesis\GPX6\analysis\mouse_mutants_FINAL_with_distances.csv")

for c in ['dg_star', 'dg_star_error', 'dg0', 'dg0_error']:
    df[c] = pd.to_numeric(df[c], errors='coerce')

df['level_num'] = df['level'].str.replace('Level', '').astype(int)

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
    all_data.append(sub)

dfp = pd.concat(all_data)
dfp = dfp[dfp['mutation'] != 'C49U'].copy()

# =================== AUTO-FIND DISTANCE COLUMN ===================
dist_col = None
possible_names = ['distance_to_Cys49', 'dist_to_Cys49', 'distance_Cys49', 'dist_C49', 'Distance_to_Cys49']
for col in possible_names:
    if col in dfp.columns:
        dist_col = col
        break
print(f"Found distance column: '{dist_col}'" if dist_col else "No distance column found")

# =================== OUTPUT FOLDER ===================
outdir = Path(r"D:\PhD_Thesis\GPX6\analysis\Figures_LFER")
outdir.mkdir(parents=True, exist_ok=True)

# =================== 1. MAIN LFER – ONLY R² IN THE BOX ===================
fig = plt.figure(figsize=(10.5, 7.8))
gs = fig.add_gridspec(1, 1, left=0.12, right=0.82, bottom=0.14, top=0.88)
ax = fig.add_subplot(gs[0])

sc = ax.scatter(dfp['ΔΔG⁰'], dfp['ΔΔG‡'],
                c=dfp['level_num'], cmap='viridis_r',
                s=16, edgecolors='black', linewidth=0.4, alpha=0.94, zorder=3)

ax.errorbar(dfp['ΔΔG⁰'], dfp['ΔΔG‡'], yerr=dfp['ΔΔG‡_err'],
            fmt='none', elinewidth=0.9, capsize=3, capthick=0.9,
            color='gray', alpha=0.6, zorder=1)

mask = dfp[['ΔΔG⁰', 'ΔΔG‡']].notna().all(axis=1)
slope, intercept, r_value, _, std_err = stats.linregress(dfp.loc[mask, 'ΔΔG⁰'], dfp.loc[mask, 'ΔΔG‡'])

x_fit = np.linspace(-3.5, 16, 200)
y_fit = slope * x_fit + intercept
ax.plot(x_fit, y_fit, color='#d62728', lw=3.8, zorder=5)

# ONLY R² — NOTHING ELSE
text = f'R² = {r_value**2:.3f}'
ax.text(0.04, 0.96, text, transform=ax.transAxes, fontsize=20, fontweight='bold',
        va='top', ha='left',
        bbox=dict(boxstyle='round,pad=0.8', facecolor='white', edgecolor='black', linewidth=1.5, alpha=0.98))

ax.set_xlabel(r'$\Delta\Delta$G° (kcal mol⁻¹)')
ax.set_ylabel(r'$\Delta\Delta$G‡ (kcal mol⁻¹)')
ax.set_title('Linear Free-Energy Relationship', fontsize=22, pad=20)
ax.set_xlim(-3.5, 16)
ax.set_ylim(-2, 22)

cax = fig.add_axes([0.85, 0.14, 0.025, 0.74])
cbar = fig.colorbar(sc, cax=cax)
cbar.set_label('Evolutionary Level', rotation=270, labelpad=28, fontsize=17)

plt.savefig(outdir / "1_LFER_final.png", dpi=600, facecolor='white')
plt.savefig(outdir / "1_LFER_final.pdf")
plt.show()

print(f"LFER → R² = {r_value**2:.3f} | n = {mask.sum()} (Φ removed from plot as requested)")

# =================== 2. MEAN EFFECT PER LEVEL ===================
plt.figure(figsize=(9,6))
means = dfp.groupby('level_num')['ΔΔG‡'].mean()
errors = dfp.groupby('level_num')['ΔΔG‡_err'].mean()
means.plot(kind='bar', yerr=errors, capsize=6, color='#4c72b0', edgecolor='black', alpha=0.9)
plt.ylabel(r'Mean $\Delta\Delta$G‡ (kcal mol⁻¹)')
plt.xlabel('Evolutionary Level')
plt.title('Average Transition-State Destabilization per Level')
plt.tight_layout()
plt.savefig(outdir / "2_Mean_by_level.png", dpi=600)
plt.show()

# =================== 3. DISTANCE VS EFFECT ===================
if dist_col:
    plt.figure(figsize=(9,6))
    plt.scatter(dfp[dist_col], dfp['ΔΔG‡'],
                c=dfp['level_num'], cmap='viridis_r',
                s=16, edgecolors='black', linewidth=0.4, alpha=0.94)
    plt.colorbar(label='Evolutionary Level')
    plt.xlabel('Distance to Cys49 (Å)')
    plt.ylabel(r'$\Delta\Delta$G‡ (kcal mol⁻¹)')
    plt.title('Structural Distance vs Functional Impact')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(outdir / "3_Distance_vs_effect.png", dpi=600)
    plt.show()

# =================== 4. HEATMAP ===================
pivot = dfp.pivot_table(values='ΔΔG‡', index='level_num', columns='mutation', aggfunc='mean')
plt.figure(figsize=(16,8))
sns.heatmap(pivot, cmap='RdBu_r', center=0, linewidths=0.5,
            cbar_kws={'label': r'$\Delta\Delta$G‡ (kcal mol⁻¹)'})
plt.title('Mutational Effects Across Evolutionary Levels', fontsize=20, pad=20)
plt.tight_layout()
plt.savefig(outdir / "4_Heatmap.png", dpi=600)
plt.show()

# =================== DONE ===================
print("\nAll figures generated successfully!")
print("Only R² is shown on the main LFER plot — no Φ, no slope, no error bars in text.")
print(f"Output folder: {outdir}")