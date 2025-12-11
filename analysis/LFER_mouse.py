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
    'legend.fontsize': 11,
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
possible_names = ['distance_to_Cys49', 'dist_to_Cys49', 'distance_Cys49', 'dist_C49', 
                  'Distance_to_Cys49', 'dist_to_Sec/Cys49']
for col in possible_names:
    if col in dfp.columns:
        dist_col = col
        break
print(f"Found distance column: '{dist_col}'" if dist_col else "No distance column found")

# Convert distance to numeric
if dist_col:
    dfp[dist_col] = pd.to_numeric(dfp[dist_col], errors='coerce')

# =================== OUTPUT FOLDER ===================
outdir = Path(r"D:\PhD_Thesis\GPX6\analysis\Figures_LFER")
outdir.mkdir(parents=True, exist_ok=True)

# =================== 1. MAIN LFER – COLOR BY DISTANCE WITH CLEAN LEGEND ===================
fig = plt.figure(figsize=(16, 8))
gs = fig.add_gridspec(1, 2, left=0.08, right=0.98, bottom=0.12, top=0.90, 
                      width_ratios=[2.5, 1], wspace=0.35)
ax_main = fig.add_subplot(gs[0])
ax_legend = fig.add_subplot(gs[1])

# Filter data with valid distances
mask_valid = dfp[['ΔΔG⁰', 'ΔΔG‡', dist_col]].notna().all(axis=1)
dfp_valid = dfp[mask_valid].copy()

# Get unique mutations with their average distance
mutation_distances = dfp_valid.groupby('mutation')[dist_col].mean().sort_values()

# Create scatter plot colored by distance
sc = ax_main.scatter(dfp_valid['ΔΔG⁰'], dfp_valid['ΔΔG‡'],
                     c=dfp_valid[dist_col], cmap='viridis',
                     s=40, edgecolors='black', linewidth=0.6, alpha=0.85, zorder=3)

# Add error bars
ax_main.errorbar(dfp_valid['ΔΔG⁰'], dfp_valid['ΔΔG‡'], yerr=dfp_valid['ΔΔG‡_err'],
                 fmt='none', elinewidth=0.9, capsize=3, capthick=0.9,
                 color='gray', alpha=0.4, zorder=1)

# Calculate regression
mask_reg = dfp[['ΔΔG⁰', 'ΔΔG‡']].notna().all(axis=1)
slope, intercept, r_value, _, std_err = stats.linregress(
    dfp.loc[mask_reg, 'ΔΔG⁰'], dfp.loc[mask_reg, 'ΔΔG‡'])

x_fit = np.linspace(-3.5, 16, 200)
y_fit = slope * x_fit + intercept
ax_main.plot(x_fit, y_fit, color='#d62728', lw=3.5, zorder=5, alpha=0.7)

# R² box
text = f'R² = {r_value**2:.3f}'
ax_main.text(0.04, 0.96, text, transform=ax_main.transAxes, fontsize=20, fontweight='bold',
             va='top', ha='left',
             bbox=dict(boxstyle='round,pad=0.8', facecolor='white', 
                      edgecolor='black', linewidth=1.5, alpha=0.98))

ax_main.set_xlabel(r'$\Delta\Delta$G° (kcal/mol)')
ax_main.set_ylabel(r'$\Delta\Delta$G‡ (kcal/mol)')
ax_main.set_title('Linear Free-Energy Relationship', fontsize=22, pad=20)
ax_main.set_xlim(-3.5, 16)
ax_main.set_ylim(-2, 22)

# Create colorbar
cbar = plt.colorbar(sc, ax=ax_main, pad=0.02, aspect=30)
cbar.set_label('Distance to Residue 49 (Å)', rotation=270, labelpad=28, fontsize=16)

# =================== LEGEND PANEL ===================
ax_legend.axis('off')
ax_legend.set_xlim(0, 1)
ax_legend.set_ylim(0, 1)

# Title for legend
ax_legend.text(0.5, 0.98, 'Mutations by Distance', 
               fontsize=16, fontweight='bold', ha='center', va='top')

# Create color-coded list of mutations
cmap = plt.cm.viridis
norm = plt.Normalize(vmin=dfp_valid[dist_col].min(), vmax=dfp_valid[dist_col].max())

y_pos = 0.92
line_height = 0.032

for mutation, distance in mutation_distances.items():
    color = cmap(norm(distance))
    
    # Color box
    ax_legend.add_patch(plt.Rectangle((0.05, y_pos - 0.015), 0.08, 0.025, 
                                      facecolor=color, edgecolor='black', linewidth=0.8))
    
    # Mutation and distance text
    ax_legend.text(0.16, y_pos, f'{mutation}', fontsize=11, va='center', ha='left', 
                   fontweight='bold')
    ax_legend.text(0.85, y_pos, f'{distance:.1f} Å', fontsize=10, va='center', ha='right',
                   color='gray')
    
    y_pos -= line_height
    
    if y_pos < 0.05:  # Start a new column if needed
        break

plt.savefig(outdir / "1_LFER_distance_colored_mouse.png", dpi=600, facecolor='white')
plt.savefig(outdir / "1_LFER_distance_colored_mouse.pdf")
plt.show()

print(f"LFER → R² = {r_value**2:.3f} | n = {mask_reg.sum()}")
print(f"Distance range: {dfp_valid[dist_col].min():.1f} - {dfp_valid[dist_col].max():.1f} Å")
print(f"\nMutations ordered by distance:")
for mutation, distance in mutation_distances.items():
    print(f"  {mutation}: {distance:.1f} Å")

# =================== 2. MEAN EFFECT PER LEVEL ===================
plt.figure(figsize=(9,6))
means = dfp.groupby('level_num')['ΔΔG‡'].mean()
errors = dfp.groupby('level_num')['ΔΔG‡_err'].mean()
means.plot(kind='bar', yerr=errors, capsize=6, color='#4c72b0', edgecolor='black', alpha=0.9)
plt.ylabel(r'Mean $\Delta\Delta$G‡ (kcal/mol)')
plt.xlabel('Evolutionary Level')
plt.title('Average Transition-State Destabilization per Level')
plt.tight_layout()
plt.savefig(outdir / "2_Mean_by_level.png", dpi=600)
plt.show()

# =================== 3. DISTANCE VS EFFECT ===================
if dist_col:
    fig, ax = plt.subplots(figsize=(10,7))
    
    # Scatter plot
    sc = ax.scatter(dfp_valid[dist_col], dfp_valid['ΔΔG‡'],
                    c=dfp_valid[dist_col], cmap='viridis',
                    s=40, edgecolors='black', linewidth=0.6, alpha=0.85)
    
    cbar = plt.colorbar(sc, label='Distance to Residue 49 (Å)')
    ax.set_xlabel('Distance to Residue 49 (Å)', fontsize=16)
    ax.set_ylabel(r'$\Delta\Delta$G‡ (kcal/mol)', fontsize=16)
    ax.set_title('Structural Distance vs Functional Impact', fontsize=18)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(outdir / "3_Distance_vs_effect.png", dpi=600)
    plt.show()

# =================== 4. HEATMAP ===================
pivot = dfp.pivot_table(values='ΔΔG‡', index='level_num', columns='mutation', aggfunc='mean')
plt.figure(figsize=(16,8))
sns.heatmap(pivot, cmap='RdBu_r', center=0, linewidths=0.5,
            cbar_kws={'label': r'$\Delta\Delta$G‡ (kcal/mol)'})
plt.title('Mutational Effects Across Evolutionary Levels', fontsize=20, pad=20)
plt.tight_layout()
plt.savefig(outdir / "4_Heatmap.png", dpi=600)
plt.show()

# =================== DONE ===================
print("\nAll figures generated successfully!")
print(f"Main LFER plot with clean color-coded legend")
print(f"Output folder: {outdir}")