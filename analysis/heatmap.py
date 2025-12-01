import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm

# Read the CSV file
file_path = r'D:\PhD_Thesis\GPX6\analysis\mouse_mutants_FINAL_with_distances.csv'
df = pd.read_csv(file_path)

# Clean column names (strip whitespace)
df.columns = df.columns.str.strip()

print("Columns in the dataset:")
print(df.columns.tolist())
print("\nFirst few rows:")
print(df.head())

# Get WT values from Level0 mousecys
wt_row = df[(df['level'] == 'Level0') & (df['mutation'] == 'mousecys')].iloc[0]
wt_dg_star = wt_row['dg_star']
wt_dg0 = wt_row['dg0']

print(f"\nWT dg_star: {wt_dg_star}")
print(f"WT dg0: {wt_dg0}")

# Remove all mousecys rows for the heatmap (since they're all reference)
df_mutants = df[df['mutation'] != 'mousecys'].copy()

# Prepare data for heatmap - pivot tables
pivot_dg_star = df_mutants.pivot(index='mutation', columns='level', values='dg_star')
pivot_dg0 = df_mutants.pivot(index='mutation', columns='level', values='dg0')

# Sort columns by level number
level_order = sorted(pivot_dg_star.columns, key=lambda x: int(x.replace('Level', '')))
pivot_dg_star = pivot_dg_star[level_order]
pivot_dg0 = pivot_dg0[level_order]

# Create figure with two subplots
fig, axes = plt.subplots(1, 2, figsize=(20, 12))

# Define colormap centered on WT values
cmap = sns.diverging_palette(250, 10, as_cmap=True)  # Blue to Red

# Plot dg_star heatmap
vmin_star = pivot_dg_star.min().min()
vmax_star = pivot_dg_star.max().max()
norm_dg_star = TwoSlopeNorm(vmin=vmin_star, vcenter=wt_dg_star, vmax=vmax_star)

sns.heatmap(pivot_dg_star, annot=False, cmap=cmap, norm=norm_dg_star,
            cbar_kws={'label': 'ΔG* value'}, ax=axes[0], linewidths=0.5)
axes[0].set_title(f'ΔG* Heatmap (WT = {wt_dg_star:.2f})', fontsize=16, fontweight='bold')
axes[0].set_xlabel('Level', fontsize=12)
axes[0].set_ylabel('Mutation', fontsize=12)
axes[0].tick_params(axis='x', rotation=45)

# Plot dg0 heatmap
vmin_dg0 = pivot_dg0.min().min()
vmax_dg0 = pivot_dg0.max().max()
norm_dg0 = TwoSlopeNorm(vmin=vmin_dg0, vcenter=wt_dg0, vmax=vmax_dg0)

sns.heatmap(pivot_dg0, annot=False, cmap=cmap, norm=norm_dg0,
            cbar_kws={'label': 'ΔG0 value'}, ax=axes[1], linewidths=0.5)
axes[1].set_title(f'ΔG0 Heatmap (WT = {wt_dg0:.2f})', fontsize=16, fontweight='bold')
axes[1].set_xlabel('Level', fontsize=12)
axes[1].set_ylabel('Mutation', fontsize=12)
axes[1].tick_params(axis='x', rotation=45)

plt.tight_layout()
plt.savefig('mutant_heatmap.png', dpi=300, bbox_inches='tight')
print("\nHeatmap saved as 'mutant_heatmap.png'")
plt.show()

# Print summary statistics
print("\n=== Summary Statistics ===")
print("\nMutations with dg_star > WT (destabilizing):")
high_dg_star = df_mutants[df_mutants['dg_star'] > wt_dg_star][['level', 'mutation', 'dg_star']].copy()
high_dg_star['diff'] = high_dg_star['dg_star'] - wt_dg_star
print(high_dg_star.sort_values('diff', ascending=False).head(10))

print("\nMutations with dg_star < WT (stabilizing):")
low_dg_star = df_mutants[df_mutants['dg_star'] < wt_dg_star][['level', 'mutation', 'dg_star']].copy()
low_dg_star['diff'] = wt_dg_star - low_dg_star['dg_star']
print(low_dg_star.sort_values('diff', ascending=False).head(10))