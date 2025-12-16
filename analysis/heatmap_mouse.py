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

# Get WT values from Level0 mousecys (the reference)
wt_row = df[(df['level'] == 'Level0') & (df['mutation'] == 'mousecys')].iloc[0]
wt_dg_star = wt_row['dg_star']
wt_dg0 = wt_row['dg0']

print(f"\nWT (mousecys) dg_star: {wt_dg_star}")
print(f"WT (mousecys) dg0: {wt_dg0}")

# Remove only mousecys rows for the heatmap (keep C49U)
df_mutants = df[df['mutation'] != 'mousecys'].copy()

print(f"\nWorking with {len(df_mutants)} mutant rows (excluding mousecys, keeping C49U)")

# Check for duplicates
duplicates = df_mutants[df_mutants.duplicated(subset=['mutation', 'level'], keep=False)]
if not duplicates.empty:
    print("\n=== WARNING: Duplicate mutation-level combinations found ===")
    print(f"Number of duplicate rows: {len(duplicates)}")
    print(duplicates[['mutation', 'level', 'dg_star', 'dg0']].sort_values(['mutation', 'level']))
    print("\nAveraging duplicate values...")

# Handle duplicates by averaging values for each mutation-level combination
df_mutants_agg = df_mutants.groupby(['mutation', 'level'], as_index=False).agg({
    'dg_star': 'mean',
    'dg0': 'mean',
    'dist_to_Sec/Cys49': 'first',  # Distance should be the same for all levels
    'dist_to_Gln83': 'first'
})

print(f"\nAfter aggregation: {len(df_mutants_agg)} unique mutation-level combinations")

# Prepare data for heatmap - pivot tables
pivot_dg_star = df_mutants_agg.pivot(index='mutation', columns='level', values='dg_star')
pivot_dg0 = df_mutants_agg.pivot(index='mutation', columns='level', values='dg0')

# Sort columns by level number
level_order = sorted(pivot_dg_star.columns, key=lambda x: int(x.replace('Level', '')))
pivot_dg_star = pivot_dg_star[level_order]
pivot_dg0 = pivot_dg0[level_order]

# Get distance information for sorting mutations
mutation_distances = df_mutants_agg.groupby('mutation').agg({
    'dist_to_Sec/Cys49': 'first'
}).reset_index()

# Sort mutations by distance (handling NaN values by putting them first)
mutation_distances['dist_to_Sec/Cys49'] = mutation_distances['dist_to_Sec/Cys49'].fillna(-1)
mutation_order = mutation_distances.sort_values('dist_to_Sec/Cys49')['mutation'].tolist()

# Reorder pivot tables by distance
pivot_dg_star = pivot_dg_star.reindex(mutation_order)
pivot_dg0 = pivot_dg0.reindex(mutation_order)

print(f"\nHeatmap dimensions: {pivot_dg_star.shape[0]} mutations x {pivot_dg_star.shape[1]} levels")

for mut in mutation_order:
    dist = mutation_distances[mutation_distances['mutation'] == mut]['dist_to_Sec/Cys49'].values[0]
    if dist == -1:
        print(f"  {mut}: N/A (at active site)")
    else:
        print(f"  {mut}: {dist:.2f} Å")

# Create figure with two subplots
fig, axes = plt.subplots(1, 2, figsize=(24, 10))

# Define colormap centered on WT values
cmap = sns.diverging_palette(250, 10, as_cmap=True)  # Blue to Red

# Plot dg_star heatmap
vmin_star = pivot_dg_star.min().min()
vmax_star = pivot_dg_star.max().max()
norm_dg_star = TwoSlopeNorm(vmin=vmin_star, vcenter=wt_dg_star, vmax=vmax_star)

sns.heatmap(pivot_dg_star, annot=False, cmap=cmap, norm=norm_dg_star,
            cbar_kws={'label': 'ΔG* value'}, ax=axes[0], linewidths=0.5)
axes[0].set_title(f'ΔG* Heatmap (WT mousecys = {wt_dg_star:.2f})', 
                  fontsize=16, fontweight='bold')
axes[0].set_xlabel('Level', fontsize=12)
axes[0].set_ylabel('Mutation (ordered by distance)', fontsize=12)
axes[0].tick_params(axis='x', rotation=45)

# Plot dg0 heatmap
vmin_dg0 = pivot_dg0.min().min()
vmax_dg0 = pivot_dg0.max().max()
norm_dg0 = TwoSlopeNorm(vmin=vmin_dg0, vcenter=wt_dg0, vmax=vmax_dg0)

sns.heatmap(pivot_dg0, annot=False, cmap=cmap, norm=norm_dg0,
            cbar_kws={'label': 'ΔG0 value'}, ax=axes[1], linewidths=0.5)
axes[1].set_title(f'ΔG0 Heatmap (WT mousecys = {wt_dg0:.2f})', 
                  fontsize=16, fontweight='bold')
axes[1].set_xlabel('Level', fontsize=12)
axes[1].set_ylabel('Mutation (ordered by distance)', fontsize=12)
axes[1].tick_params(axis='x', rotation=45)

plt.tight_layout()
plt.savefig('mutant_heatmap_mouse_by_distance.png', dpi=300, bbox_inches='tight')
print("\nHeatmap saved as 'mutant_heatmap_mouse_by_distance.png'")
plt.show()

# Print summary statistics (using aggregated data)
print("\n=== Summary Statistics ===")
print("\nTop 10 mutations with dg_star > WT (most destabilizing):")
high_dg_star = df_mutants_agg[df_mutants_agg['dg_star'] > wt_dg_star][['level', 'mutation', 'dg_star', 'dist_to_Sec/Cys49']].copy()
if len(high_dg_star) > 0:
    high_dg_star['diff'] = high_dg_star['dg_star'] - wt_dg_star
    print(high_dg_star.sort_values('diff', ascending=False).head(10))
else:
    print("None found")

print("\nTop 10 mutations with dg_star < WT (most stabilizing):")
low_dg_star = df_mutants_agg[df_mutants_agg['dg_star'] < wt_dg_star][['level', 'mutation', 'dg_star', 'dist_to_Sec/Cys49']].copy()
if len(low_dg_star) > 0:
    low_dg_star['diff'] = wt_dg_star - low_dg_star['dg_star']
    print(low_dg_star.sort_values('diff', ascending=False).head(10))
else:
    print("None found")

print("\nTop 10 mutations with dg0 < WT (most negative/favorable):")
low_dg0 = df_mutants_agg[df_mutants_agg['dg0'] < wt_dg0][['level', 'mutation', 'dg0', 'dist_to_Sec/Cys49']].copy()
if len(low_dg0) > 0:
    low_dg0['diff'] = wt_dg0 - low_dg0['dg0']
    print(low_dg0.sort_values('diff', ascending=False).head(10))
else:
    print("None found")

# Additional comparison: C49U (selenocysteine) vs mousecys
print("\n=== C49U (Selenocysteine) vs Mousecys (Cysteine) Comparison ===")
c49u_data = df[df['mutation'] == 'C49U'][['level', 'dg_star', 'dg0']]
if len(c49u_data) > 0:
    print(f"C49U dg_star: {c49u_data.iloc[0]['dg_star']:.2f} (diff from WT: {c49u_data.iloc[0]['dg_star'] - wt_dg_star:.2f})")
    print(f"C49U dg0: {c49u_data.iloc[0]['dg0']:.2f} (diff from WT: {c49u_data.iloc[0]['dg0'] - wt_dg0:.2f})")