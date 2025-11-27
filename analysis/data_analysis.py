import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit
from pathlib import Path
import seaborn as sns
from matplotlib.patches import Ellipse
from sklearn.decomposition import PCA

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
    sub['ΔΔG⁰_err'] = np.sqrt(sub['dg0_error']**2 + r['dg0_error']**2)
    all_data.append(sub)

dfp = pd.concat(all_data)
dfp = dfp[dfp['mutation'] != 'C49U'].copy()

# Auto-find distance column
dist_col = None
for col in ['distance_to_Cys49', 'dist_to_Cys49', 'distance_Cys49', 'dist_C49', 'Distance_to_Cys49']:
    if col in dfp.columns:
        dist_col = col
        break

outdir = Path(r"D:\PhD_Thesis\GPX6\analysis\Figures_LFER_Extended")
outdir.mkdir(parents=True, exist_ok=True)

# =================== ANALYSIS 1: BRØNSTED COEFFICIENT (α) PER LEVEL ===================
print("\n" + "="*70)
print("ANALYSIS 1: BRØNSTED COEFFICIENT (α) BY EVOLUTIONARY LEVEL")
print("="*70)

n_levels = len(dfp['level_num'].unique())
n_cols = 3
n_rows = int(np.ceil(n_levels / n_cols))
fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 6*n_rows))
axes = axes.flatten()

bronsted_results = []
for idx, level in enumerate(sorted(dfp['level_num'].unique())):
    sub = dfp[dfp['level_num'] == level].copy()
    mask = sub[['ΔΔG⁰', 'ΔΔG‡']].notna().all(axis=1)
    
    if mask.sum() > 2:
        slope, intercept, r, _, se = stats.linregress(sub.loc[mask, 'ΔΔG⁰'], 
                                                        sub.loc[mask, 'ΔΔG‡'])
        
        axes[idx].scatter(sub['ΔΔG⁰'], sub['ΔΔG‡'], s=80, alpha=0.7, 
                         edgecolors='black', linewidth=1)
        x_fit = np.linspace(sub['ΔΔG⁰'].min()-1, sub['ΔΔG⁰'].max()+1, 100)
        axes[idx].plot(x_fit, slope*x_fit + intercept, 'r-', lw=2.5)
        
        axes[idx].set_title(f'Level {level}: α = {slope:.3f}±{se:.3f}\nR² = {r**2:.3f}', 
                           fontsize=16)
        axes[idx].set_xlabel(r'$\Delta\Delta$G° (kcal mol⁻¹)')
        axes[idx].set_ylabel(r'$\Delta\Delta$G‡ (kcal mol⁻¹)')
        axes[idx].grid(alpha=0.3)
        
        bronsted_results.append({
            'level': level,
            'alpha': slope,
            'alpha_error': se,
            'R²': r**2,
            'n_mutations': mask.sum()
        })
        
        print(f"Level {level}: α = {slope:.4f} ± {se:.4f} | R² = {r**2:.4f} | n = {mask.sum()}")

# Hide unused subplots
for idx in range(n_levels, len(axes)):
    axes[idx].axis('off')

plt.tight_layout()
plt.savefig(outdir / "5_Bronsted_by_level.png", dpi=600)
plt.show()

br_df = pd.DataFrame(bronsted_results)

# Plot Brønsted coefficient evolution
plt.figure(figsize=(10, 6))
plt.errorbar(br_df['level'], br_df['alpha'], yerr=br_df['alpha_error'], 
             fmt='o-', markersize=12, linewidth=2.5, capsize=8, capthick=2)
plt.axhline(0.5, color='red', linestyle='--', linewidth=2, label='α = 0.5 (Hammond)')
plt.xlabel('Evolutionary Level')
plt.ylabel('Brønsted Coefficient (α)')
plt.title('Transition State Position Across Evolution')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(outdir / "6_Bronsted_evolution.png", dpi=600)
plt.show()

# =================== ANALYSIS 2: MUTATION-SPECIFIC COUPLING CONSTANTS ===================
print("\n" + "="*70)
print("ANALYSIS 2: MUTATION-SPECIFIC COUPLING (Φ-VALUES)")
print("="*70)

mutations_with_data = dfp.groupby('mutation').filter(lambda x: len(x) >= 3)['mutation'].unique()

phi_results = []
n_muts = len(mutations_with_data)
n_cols_mut = 4
n_rows_mut = int(np.ceil(n_muts / n_cols_mut))
fig, axes = plt.subplots(n_rows_mut, n_cols_mut, figsize=(20, 5*n_rows_mut))
axes = axes.flatten()

for idx, mut in enumerate(mutations_with_data):
    sub = dfp[dfp['mutation'] == mut].copy()
    mask = sub[['ΔΔG⁰', 'ΔΔG‡']].notna().all(axis=1)
    
    if mask.sum() >= 3:
        slope, intercept, r, p, se = stats.linregress(sub.loc[mask, 'ΔΔG⁰'], 
                                                        sub.loc[mask, 'ΔΔG‡'])
        
        axes[idx].errorbar(sub['ΔΔG⁰'], sub['ΔΔG‡'], 
                          xerr=sub['ΔΔG⁰_err'], yerr=sub['ΔΔG‡_err'],
                          fmt='o', markersize=8, capsize=5, alpha=0.7)
        
        x_fit = np.linspace(sub['ΔΔG⁰'].min()-0.5, sub['ΔΔG⁰'].max()+0.5, 100)
        axes[idx].plot(x_fit, slope*x_fit + intercept, 'r-', lw=2)
        
        axes[idx].set_title(f'{mut}: Φ = {slope:.3f}±{se:.3f}', fontsize=14)
        axes[idx].set_xlabel(r'$\Delta\Delta$G°')
        axes[idx].set_ylabel(r'$\Delta\Delta$G‡')
        axes[idx].grid(alpha=0.3)
        
        phi_results.append({
            'mutation': mut,
            'phi': slope,
            'phi_error': se,
            'R²': r**2,
            'p_value': p
        })
        
        print(f"{mut:10s}: Φ = {slope:.4f} ± {se:.4f} | R² = {r**2:.4f} | p = {p:.4e}")

# Hide unused subplots
for idx in range(n_muts, len(axes)):
    axes[idx].axis('off')

plt.tight_layout()
plt.savefig(outdir / "7_Phi_values_per_mutation.png", dpi=600)
plt.show()

phi_df = pd.DataFrame(phi_results)
if dist_col and not phi_df.empty:
    phi_df = phi_df.merge(dfp[['mutation', dist_col]].drop_duplicates(), on='mutation', how='left')
    
    plt.figure(figsize=(10, 6))
    plt.scatter(phi_df[dist_col], phi_df['phi'], s=120, alpha=0.7, 
               edgecolors='black', linewidth=1.5)
    for _, row in phi_df.iterrows():
        plt.annotate(row['mutation'], (row[dist_col], row['phi']), 
                    fontsize=9, ha='right')
    plt.xlabel('Distance to Cys49 (Å)')
    plt.ylabel('Φ-value (Coupling Constant)')
    plt.title('Structural Coupling vs Distance')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(outdir / "8_Phi_vs_distance.png", dpi=600)
    plt.show()

# =================== ANALYSIS 3: RESIDUAL ANALYSIS & OUTLIER DETECTION ===================
print("\n" + "="*70)
print("ANALYSIS 3: RESIDUAL ANALYSIS & OUTLIERS")
print("="*70)

mask_all = dfp[['ΔΔG⁰', 'ΔΔG‡']].notna().all(axis=1)
slope_all, intercept_all, _, _, _ = stats.linregress(dfp.loc[mask_all, 'ΔΔG⁰'], 
                                                       dfp.loc[mask_all, 'ΔΔG‡'])

dfp.loc[mask_all, 'predicted'] = slope_all * dfp.loc[mask_all, 'ΔΔG⁰'] + intercept_all
dfp.loc[mask_all, 'residual'] = dfp.loc[mask_all, 'ΔΔG‡'] - dfp.loc[mask_all, 'predicted']
dfp.loc[mask_all, 'std_residual'] = dfp.loc[mask_all, 'residual'] / dfp.loc[mask_all, 'residual'].std()

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Residuals vs fitted
axes[0].scatter(dfp.loc[mask_all, 'predicted'], dfp.loc[mask_all, 'residual'],
               c=dfp.loc[mask_all, 'level_num'], cmap='viridis_r', s=60, alpha=0.7)
axes[0].axhline(0, color='red', linestyle='--', linewidth=2)
axes[0].set_xlabel(r'Predicted $\Delta\Delta$G‡')
axes[0].set_ylabel('Residual')
axes[0].set_title('Residuals vs Fitted Values')
axes[0].grid(alpha=0.3)

# Q-Q plot
stats.probplot(dfp.loc[mask_all, 'residual'].dropna(), dist="norm", plot=axes[1])
axes[1].set_title('Q-Q Plot (Normality Check)')
axes[1].grid(alpha=0.3)

# Histogram of residuals
axes[2].hist(dfp.loc[mask_all, 'residual'].dropna(), bins=20, edgecolor='black', alpha=0.7)
axes[2].set_xlabel('Residual')
axes[2].set_ylabel('Frequency')
axes[2].set_title('Distribution of Residuals')
axes[2].grid(alpha=0.3)

plt.tight_layout()
plt.savefig(outdir / "9_Residual_analysis.png", dpi=600)
plt.show()

# Identify outliers
outliers = dfp[np.abs(dfp['std_residual']) > 2].copy()
print(f"\nOutliers (|standardized residual| > 2): {len(outliers)}")
if not outliers.empty:
    print(outliers[['mutation', 'level', 'ΔΔG⁰', 'ΔΔG‡', 'residual', 'std_residual']].to_string())

# =================== ANALYSIS 4: ENERGY LANDSCAPE CLUSTERING ===================
print("\n" + "="*70)
print("ANALYSIS 4: MUTATION CLUSTERING IN ENERGY SPACE")
print("="*70)

from sklearn.cluster import KMeans

features = dfp[['ΔΔG⁰', 'ΔΔG‡']].dropna()
features_indexed = dfp[['ΔΔG⁰', 'ΔΔG‡']].dropna()
kmeans = KMeans(n_clusters=3, random_state=42)
clusters = kmeans.fit_predict(features)

# Create figure with annotations
fig, ax = plt.subplots(figsize=(12, 8))

# Define cluster colors and labels
cluster_colors = ['#FF6B6B', '#4ECDC4', '#95E1D3']
cluster_labels = ['High Barrier\n(Slow Mutations)', 
                 'Moderate Impact\n(Balanced)', 
                 'Low Barrier\n(Fast Mutations)']

# Sort clusters by mean ΔΔG‡ to assign meaningful labels
cluster_means = []
for i in range(3):
    cluster_mask = clusters == i
    mean_ddg_star = features.iloc[cluster_mask]['ΔΔG‡'].mean()
    cluster_means.append((i, mean_ddg_star))
cluster_means.sort(key=lambda x: x[1], reverse=True)

# Map old cluster IDs to new sorted IDs
cluster_mapping = {old_id: new_id for new_id, (old_id, _) in enumerate(cluster_means)}
sorted_clusters = np.array([cluster_mapping[c] for c in clusters])

# Plot each cluster
for i in range(3):
    mask = sorted_clusters == i
    ax.scatter(features.iloc[mask]['ΔΔG⁰'], 
              features.iloc[mask]['ΔΔG‡'],
              c=cluster_colors[i], 
              s=100, 
              alpha=0.7, 
              edgecolors='black',
              linewidth=1,
              label=cluster_labels[i])

# Plot cluster centers (use 'star' or '*' marker instead of emoji)
sorted_centers = np.array([kmeans.cluster_centers_[old_id] for old_id, _ in cluster_means])
ax.scatter(sorted_centers[:, 0], sorted_centers[:, 1],
          s=600, c='gold', marker='*', 
          edgecolors='black', linewidths=3,
          label='Cluster Centers', zorder=10)

# Add LFER line for reference
mask_all = dfp[['ΔΔG⁰', 'ΔΔG‡']].notna().all(axis=1)
slope_all, intercept_all, _, _, _ = stats.linregress(dfp.loc[mask_all, 'ΔΔG⁰'], 
                                                       dfp.loc[mask_all, 'ΔΔG‡'])
x_fit = np.linspace(features['ΔΔG⁰'].min()-1, features['ΔΔG⁰'].max()+1, 100)
ax.plot(x_fit, slope_all*x_fit + intercept_all, 
       'k--', lw=2.5, alpha=0.5, label='Global LFER')

# Add interpretation zones
ax.axhline(0, color='gray', linestyle=':', linewidth=1.5, alpha=0.5)
ax.axvline(0, color='gray', linestyle=':', linewidth=1.5, alpha=0.5)

ax.set_xlabel(r'$\Delta\Delta$G° (kcal mol⁻¹)', fontsize=18)
ax.set_ylabel(r'$\Delta\Delta$G‡ (kcal mol⁻¹)', fontsize=18)
ax.set_title('Mutation Classification by Energetic Profile', fontsize=20, pad=20)
ax.legend(loc='upper left', fontsize=13, framealpha=0.95)
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(outdir / "10_Energy_clustering.png", dpi=600)
plt.show()

# Print cluster statistics
print("\nCluster Interpretations:")
for i in range(3):
    cluster_data = features.iloc[sorted_clusters == i]
    n_points = len(cluster_data)
    mean_dg0 = cluster_data['ΔΔG⁰'].mean()
    mean_dg_star = cluster_data['ΔΔG‡'].mean()
    
    print(f"\n{cluster_labels[i].replace(chr(10), ' ')}:")
    print(f"  n = {n_points} mutations")
    print(f"  Mean ΔΔG° = {mean_dg0:.2f} kcal/mol")
    print(f"  Mean ΔΔG‡ = {mean_dg_star:.2f} kcal/mol")
    print(f"  Interpretation: ", end="")
    
    if mean_dg_star > 10:
        print("Strong kinetic barriers - mutations severely impair catalysis")
    elif mean_dg_star > 5:
        print("Moderate kinetic impact - partial loss of catalytic efficiency")
    else:
        print("Minimal kinetic perturbation - enzyme maintains activity")

# =================== ANALYSIS 5: MUTATION MATRIX & TRENDS ===================
print("\n" + "="*70)
print("ANALYSIS 5: EVOLUTIONARY TRENDS IN MUTATION EFFECTS")
print("="*70)

pivot_dg0 = dfp.pivot_table(values='ΔΔG⁰', index='mutation', columns='level_num', aggfunc='mean')
pivot_dg_star = dfp.pivot_table(values='ΔΔG‡', index='mutation', columns='level_num', aggfunc='mean')

fig, axes = plt.subplots(2, 1, figsize=(14, 16))

# Plot 1: Transition State Effects
sns.heatmap(pivot_dg_star, cmap='YlOrRd', linewidths=0.8, ax=axes[0],
           cbar_kws={'label': r'$\Delta\Delta$G‡ (kcal mol⁻¹)'}, 
           vmin=0, annot=False)
axes[0].set_title(r'Transition State Destabilization ($\Delta\Delta$G‡) Across Evolution', 
                 fontsize=18, pad=15)
axes[0].set_ylabel('Mutation', fontsize=16)
axes[0].set_xlabel('Evolutionary Level', fontsize=16)

# Plot 2: Ground State Effects
sns.heatmap(pivot_dg0, cmap='coolwarm', center=0, linewidths=0.8, ax=axes[1],
           cbar_kws={'label': r'$\Delta\Delta$G° (kcal mol⁻¹)'}, annot=False)
axes[1].set_title(r'Ground State Stabilization ($\Delta\Delta$G°) Across Evolution', 
                 fontsize=18, pad=15)
axes[1].set_ylabel('Mutation', fontsize=16)
axes[1].set_xlabel('Evolutionary Level', fontsize=16)

plt.tight_layout()
plt.savefig(outdir / "11_Dual_heatmaps.png", dpi=600)
plt.show()

# Calculate trends within each level
print("\nEvolutionary Trends Analysis:")
print("-" * 70)

for level in sorted(dfp['level_num'].unique()):
    level_data = dfp[dfp['level_num'] == level]
    
    print(f"\n🧬 Level {level}:")
    print(f"   Number of mutations: {level_data['mutation'].nunique()}")
    
    # Transition state statistics
    mean_ddg_star = level_data['ΔΔG‡'].mean()
    std_ddg_star = level_data['ΔΔG‡'].std()
    max_ddg_star = level_data['ΔΔG‡'].max()
    print(f"   ΔΔG‡: Mean = {mean_ddg_star:.2f} ± {std_ddg_star:.2f} kcal/mol")
    print(f"        Range: 0 to {max_ddg_star:.2f} kcal/mol")
    
    # Ground state statistics
    mean_ddg0 = level_data['ΔΔG⁰'].mean()
    std_ddg0 = level_data['ΔΔG⁰'].std()
    print(f"   ΔΔG°:  Mean = {mean_ddg0:.2f} ± {std_ddg0:.2f} kcal/mol")
    
    # Catalytic efficiency interpretation
    if mean_ddg_star < 5:
        efficiency = "HIGH - Mutations well-tolerated"
    elif mean_ddg_star < 10:
        efficiency = "MODERATE - Significant sensitivity"
    else:
        efficiency = "LOW - High mutational vulnerability"
    print(f"   Catalytic Robustness: {efficiency}")

# Analyze cross-level trends
print("\n" + "="*70)
print("CROSS-EVOLUTIONARY ANALYSIS")
print("="*70)

level_summary = dfp.groupby('level_num').agg({
    'ΔΔG‡': ['mean', 'std', 'min', 'max'],
    'ΔΔG⁰': ['mean', 'std'],
    'mutation': 'count'
}).round(2)

print("\nLevel-by-Level Summary:")
print(level_summary.to_string())

# Test for increasing/decreasing trend
levels = sorted(dfp['level_num'].unique())
mean_barriers = [dfp[dfp['level_num']==l]['ΔΔG‡'].mean() for l in levels]

from scipy.stats import spearmanr
corr, p_value = spearmanr(levels, mean_barriers)
print(f"\nTrend Test (Spearman correlation):")
print(f"  ρ = {corr:.3f}, p = {p_value:.4f}")

if p_value < 0.05:
    if corr > 0:
        trend = "INCREASING - Later levels show higher mutational sensitivity"
    else:
        trend = "DECREASING - Later levels show greater mutational robustness"
else:
    trend = "NO SIGNIFICANT TREND - Mutational effects stable across evolution"
    
print(f"  Interpretation: {trend}")

# =================== SUMMARY REPORT ===================
print("\n" + "="*70)
print("SUMMARY STATISTICS")
print("="*70)

print(f"\nGlobal LFER:")
mask = dfp[['ΔΔG⁰', 'ΔΔG‡']].notna().all(axis=1)
slope, intercept, r, p, se = stats.linregress(dfp.loc[mask, 'ΔΔG⁰'], dfp.loc[mask, 'ΔΔG‡'])
print(f"  α = {slope:.4f} ± {se:.4f}")
print(f"  R² = {r**2:.4f}")
print(f"  p-value = {p:.2e}")
print(f"  n = {mask.sum()} mutations")

print(f"\nBrønsted coefficient range: {br_df['alpha'].min():.3f} to {br_df['alpha'].max():.3f}")
print(f"Mean α across levels: {br_df['alpha'].mean():.3f} ± {br_df['alpha'].std():.3f}")

if not phi_df.empty:
    print(f"\nΦ-value range: {phi_df['phi'].min():.3f} to {phi_df['phi'].max():.3f}")
    print(f"Mean Φ across mutations: {phi_df['phi'].mean():.3f} ± {phi_df['phi'].std():.3f}")

print(f"\n{'='*70}")
print(f"All figures saved to: {outdir}")
print(f"{'='*70}")