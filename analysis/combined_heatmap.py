import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm

# Set high-quality rendering parameters
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['pdf.fonttype'] = 42  # TrueType fonts for PLOS
plt.rcParams['ps.fonttype'] = 42

# Read both CSV files
mouse_file_path = r'D:\PhD_Thesis\GPX6\analysis\mouse_mutants_FINAL_with_distances.csv'
human_file_path = r'D:\PhD_Thesis\GPX6\analysis\human_mutants_FINAL_with_distances.csv'

df_mouse = pd.read_csv(mouse_file_path)
df_human = pd.read_csv(human_file_path)

# Clean column names
df_mouse.columns = df_mouse.columns.str.strip()
df_human.columns = df_human.columns.str.strip()

print("="*80)
print("MOUSE DATA")
print("="*80)
print(f"Columns: {df_mouse.columns.tolist()}")
print(f"\nFirst few rows:\n{df_mouse.head()}")

print("\n" + "="*80)
print("HUMAN DATA")
print("="*80)
print(f"Columns: {df_human.columns.tolist()}")
print(f"\nFirst few rows:\n{df_human.head()}")

# Get WT values
wt_mouse_row = df_mouse[(df_mouse['level'] == 'Level0') & (df_mouse['mutation'] == 'mousecys')].iloc[0]
wt_mouse_dg_star = wt_mouse_row['dg_star']
wt_mouse_dg0 = wt_mouse_row['dg0']

wt_human_row = df_human[(df_human['level'] == 'Level0') & (df_human['mutation'] == 'humansec')].iloc[0]
wt_human_dg_star = wt_human_row['dg_star']
wt_human_dg0 = wt_human_row['dg0']

print(f"\nWT Mouse (mousecys) - dg_star: {wt_mouse_dg_star:.2f}, dg0: {wt_mouse_dg0:.2f}")
print(f"WT Human (humansec) - dg_star: {wt_human_dg_star:.2f}, dg0: {wt_human_dg0:.2f}")

# Remove WT rows for heatmaps
df_mouse_mutants = df_mouse[df_mouse['mutation'] != 'mousecys'].copy()
df_human_mutants = df_human[df_human['mutation'] != 'humansec'].copy()

print(f"\nMouse mutants: {len(df_mouse_mutants)} rows")
print(f"Human mutants: {len(df_human_mutants)} rows")

# Handle duplicates by averaging
def process_mutants(df_mutants, dataset_name):
    duplicates = df_mutants[df_mutants.duplicated(subset=['mutation', 'level'], keep=False)]
    if not duplicates.empty:
        print(f"\n{dataset_name}: {len(duplicates)} duplicate rows found, averaging...")
    
    df_agg = df_mutants.groupby(['mutation', 'level'], as_index=False).agg({
        'dg_star': 'mean',
        'dg0': 'mean',
        'dist_to_Sec/Cys49': 'first',
        'dist_to_Gln83': 'first'
    })
    
    print(f"{dataset_name}: {len(df_agg)} unique mutation-level combinations")
    return df_agg

df_mouse_agg = process_mutants(df_mouse_mutants, "MOUSE")
df_human_agg = process_mutants(df_human_mutants, "HUMAN")

# Create pivot tables
def create_pivots(df_agg):
    pivot_dg_star = df_agg.pivot(index='mutation', columns='level', values='dg_star')
    pivot_dg0 = df_agg.pivot(index='mutation', columns='level', values='dg0')
    
    level_order = sorted(pivot_dg_star.columns, key=lambda x: int(x.replace('Level', '')))
    pivot_dg_star = pivot_dg_star[level_order]
    pivot_dg0 = pivot_dg0[level_order]
    
    mutation_distances = df_agg.groupby('mutation').agg({
        'dist_to_Sec/Cys49': 'first'
    }).reset_index()
    
    mutation_distances['dist_to_Sec/Cys49'] = mutation_distances['dist_to_Sec/Cys49'].fillna(-1)
    mutation_order = mutation_distances.sort_values('dist_to_Sec/Cys49')['mutation'].tolist()
    
    pivot_dg_star = pivot_dg_star.reindex(mutation_order)
    pivot_dg0 = pivot_dg0.reindex(mutation_order)
    
    return pivot_dg_star, pivot_dg0, mutation_order, mutation_distances

mouse_pivot_dg_star, mouse_pivot_dg0, mouse_mutation_order, mouse_distances = create_pivots(df_mouse_agg)
human_pivot_dg_star, human_pivot_dg0, human_mutation_order, human_distances = create_pivots(df_human_agg)

print(f"\nMouse: {mouse_pivot_dg_star.shape[0]} mutations x {mouse_pivot_dg_star.shape[1]} levels")
print(f"Human: {human_pivot_dg_star.shape[0]} mutations x {human_pivot_dg_star.shape[1]} levels")

# Colormap
cmap = sns.diverging_palette(250, 10, as_cmap=True)

# ============================================================================
# LARGE READABLE HEATMAPS - ΔG* Comparison
# Generate at 300 DPI (instead of 600) for better memory management
# while maintaining excellent quality and readability
# ============================================================================
n_mutations_mouse = mouse_pivot_dg_star.shape[0]
n_mutations_human = human_pivot_dg_star.shape[0]
n_levels_mouse = mouse_pivot_dg_star.shape[1]
n_levels_human = human_pivot_dg_star.shape[1]

# Calculate LARGE dimensions for readability
# Height: 0.5 inches per mutation (generous spacing)
height_needed = max(n_mutations_mouse, n_mutations_human) * 0.5
fig_height = max(12, min(height_needed, 25))  # Between 12-25 inches

# Width: generous space per level
width_per_panel = max(n_levels_mouse, n_levels_human) * 0.7
fig_width = max(20, width_per_panel * 2.5)  # Minimum 20", accommodate both panels

print(f"\nCreating ΔG* figure: {fig_width:.1f} x {fig_height:.1f} inches")
print(f"At 300 DPI: {int(fig_width*300)} x {int(fig_height*300)} pixels")

fig1, axes1 = plt.subplots(1, 2, figsize=(fig_width, fig_height), dpi=100)

# --- Mouse dg_star ---
vmin_star_mouse = mouse_pivot_dg_star.min().min()
vmax_star_mouse = mouse_pivot_dg_star.max().max()
norm_dg_star_mouse = TwoSlopeNorm(vmin=vmin_star_mouse, vcenter=wt_mouse_dg_star, vmax=vmax_star_mouse)

sns.heatmap(mouse_pivot_dg_star, annot=False, cmap=cmap, norm=norm_dg_star_mouse,
            cbar_kws={'label': 'ΔG* (kcal/mol)', 'shrink': 0.65}, 
            ax=axes1[0], linewidths=0.8, linecolor='lightgray',
            square=False)

axes1[0].set_title(f'Mouse ΔG*\n(WT = {wt_mouse_dg_star:.2f} kcal/mol)', 
                   fontsize=24, fontweight='bold', pad=20)
axes1[0].set_xlabel('Mutation Level', fontsize=18, fontweight='bold', labelpad=12)
axes1[0].set_ylabel('Mutation (ordered by distance)', fontsize=18, fontweight='bold', labelpad=12)
axes1[0].tick_params(axis='x', rotation=45, labelsize=14, pad=6)
axes1[0].tick_params(axis='y', rotation=0, labelsize=14, pad=4, length=4)

# Colorbar formatting
cbar_mouse = axes1[0].collections[0].colorbar
cbar_mouse.set_label('ΔG* (kcal/mol)', fontsize=16, fontweight='bold', labelpad=12)
cbar_mouse.ax.tick_params(labelsize=13)

# --- Human dg_star ---
vmin_star_human = human_pivot_dg_star.min().min()
vmax_star_human = human_pivot_dg_star.max().max()
norm_dg_star_human = TwoSlopeNorm(vmin=vmin_star_human, vcenter=wt_human_dg_star, vmax=vmax_star_human)

sns.heatmap(human_pivot_dg_star, annot=False, cmap=cmap, norm=norm_dg_star_human,
            cbar_kws={'label': 'ΔG* (kcal/mol)', 'shrink': 0.65}, 
            ax=axes1[1], linewidths=0.8, linecolor='lightgray',
            square=False)

axes1[1].set_title(f'Human ΔG*\n(WT = {wt_human_dg_star:.2f} kcal/mol)', 
                   fontsize=24, fontweight='bold', pad=20)
axes1[1].set_xlabel('Mutation Level', fontsize=18, fontweight='bold', labelpad=12)
axes1[1].set_ylabel('Mutation (ordered by distance)', fontsize=18, fontweight='bold', labelpad=12)
axes1[1].tick_params(axis='x', rotation=45, labelsize=14, pad=6)
axes1[1].tick_params(axis='y', rotation=0, labelsize=14, pad=4, length=4)

# Colorbar formatting
cbar_human = axes1[1].collections[0].colorbar
cbar_human.set_label('ΔG* (kcal/mol)', fontsize=16, fontweight='bold', labelpad=12)
cbar_human.ax.tick_params(labelsize=13)

# Adjust layout
plt.subplots_adjust(left=0.10, right=0.92, top=0.95, bottom=0.08, wspace=0.35)

# Save at 300 DPI - excellent quality, manageable file size
output_path_star = r'D:\PhD_Thesis\GPX6\analysis\comparison_dg_star_mouse_vs_human.png'
plt.savefig(output_path_star, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
print(f"\nΔG* saved: {output_path_star}")
print("Resolution: 300 DPI (PLOS ONE compliant, excellent for combined figures)")

plt.close()

# ============================================================================
# LARGE READABLE HEATMAPS - ΔG0 Comparison
# ============================================================================
print(f"\nCreating ΔG0 figure: {fig_width:.1f} x {fig_height:.1f} inches")

fig2, axes2 = plt.subplots(1, 2, figsize=(fig_width, fig_height), dpi=100)

# --- Mouse dg0 ---
vmin_dg0_mouse = mouse_pivot_dg0.min().min()
vmax_dg0_mouse = mouse_pivot_dg0.max().max()
norm_dg0_mouse = TwoSlopeNorm(vmin=vmin_dg0_mouse, vcenter=wt_mouse_dg0, vmax=vmax_dg0_mouse)

sns.heatmap(mouse_pivot_dg0, annot=False, cmap=cmap, norm=norm_dg0_mouse,
            cbar_kws={'label': 'ΔG0 (kcal/mol)', 'shrink': 0.65}, 
            ax=axes2[0], linewidths=0.8, linecolor='lightgray',
            square=False)

axes2[0].set_title(f'Mouse ΔG0\n(WT = {wt_mouse_dg0:.2f} kcal/mol)', 
                   fontsize=24, fontweight='bold', pad=20)
axes2[0].set_xlabel('Mutation Level', fontsize=18, fontweight='bold', labelpad=12)
axes2[0].set_ylabel('Mutation (ordered by distance)', fontsize=18, fontweight='bold', labelpad=12)
axes2[0].tick_params(axis='x', rotation=45, labelsize=14, pad=6)
axes2[0].tick_params(axis='y', rotation=0, labelsize=14, pad=4, length=4)

# Colorbar formatting
cbar_mouse_dg0 = axes2[0].collections[0].colorbar
cbar_mouse_dg0.set_label('ΔG0 (kcal/mol)', fontsize=16, fontweight='bold', labelpad=12)
cbar_mouse_dg0.ax.tick_params(labelsize=13)

# --- Human dg0 ---
vmin_dg0_human = human_pivot_dg0.min().min()
vmax_dg0_human = human_pivot_dg0.max().max()
norm_dg0_human = TwoSlopeNorm(vmin=vmin_dg0_human, vcenter=wt_human_dg0, vmax=vmax_dg0_human)

sns.heatmap(human_pivot_dg0, annot=False, cmap=cmap, norm=norm_dg0_human,
            cbar_kws={'label': 'ΔG0 (kcal/mol)', 'shrink': 0.65}, 
            ax=axes2[1], linewidths=0.8, linecolor='lightgray',
            square=False)

axes2[1].set_title(f'Human ΔG0\n(WT = {wt_human_dg0:.2f} kcal/mol)', 
                   fontsize=24, fontweight='bold', pad=20)
axes2[1].set_xlabel('Mutation Level', fontsize=18, fontweight='bold', labelpad=12)
axes2[1].set_ylabel('Mutation (ordered by distance)', fontsize=18, fontweight='bold', labelpad=12)
axes2[1].tick_params(axis='x', rotation=45, labelsize=14, pad=6)
axes2[1].tick_params(axis='y', rotation=0, labelsize=14, pad=4, length=4)

# Colorbar formatting
cbar_human_dg0 = axes2[1].collections[0].colorbar
cbar_human_dg0.set_label('ΔG0 (kcal/mol)', fontsize=16, fontweight='bold', labelpad=12)
cbar_human_dg0.ax.tick_params(labelsize=13)

# Adjust layout
plt.subplots_adjust(left=0.10, right=0.92, top=0.95, bottom=0.08, wspace=0.35)

# Save at 300 DPI
output_path_dg0 = r'D:\PhD_Thesis\GPX6\analysis\comparison_dg0_mouse_vs_human.png'
plt.savefig(output_path_dg0, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
print(f"\nΔG0 saved: {output_path_dg0}")

plt.close()

print("\n" + "="*80)
print("LARGE READABLE HEATMAPS CREATED")
print("="*80)
print(f"Figure dimensions: {fig_width:.1f} x {fig_height:.1f} inches")
print("Resolution: 300 DPI")
print("- Excellent quality for publication AND combining")
print("- Larger fonts and spacing than before")
print("- Memory-efficient for combining into multi-panel figures")
print("- Still meets PLOS ONE requirements (300-600 DPI)")
print("="*80)

# Summary statistics
print("\n" + "="*80)
print("SUMMARY STATISTICS")
print("="*80)

print("\n" + "-"*80)
print("MOUSE - Top 10 mutations with dg_star > WT (most destabilizing):")
print("-"*80)
high_dg_star_mouse = df_mouse_agg[df_mouse_agg['dg_star'] > wt_mouse_dg_star][['level', 'mutation', 'dg_star', 'dist_to_Sec/Cys49']].copy()
if len(high_dg_star_mouse) > 0:
    high_dg_star_mouse['diff'] = high_dg_star_mouse['dg_star'] - wt_mouse_dg_star
    print(high_dg_star_mouse.sort_values('diff', ascending=False).head(10))
else:
    print("None found")

print("\n" + "-"*80)
print("HUMAN - Top 10 mutations with dg_star > WT (most destabilizing):")
print("-"*80)
high_dg_star_human = df_human_agg[df_human_agg['dg_star'] > wt_human_dg_star][['level', 'mutation', 'dg_star', 'dist_to_Sec/Cys49']].copy()
if len(high_dg_star_human) > 0:
    high_dg_star_human['diff'] = high_dg_star_human['dg_star'] - wt_human_dg_star
    print(high_dg_star_human.sort_values('diff', ascending=False).head(10))
else:
    print("None found")

print("\n" + "="*80)
print("Analysis complete!")
print("="*80)
print("\nNOTE: These heatmaps are at 300 DPI instead of 600 DPI.")
print("This makes them:")
print("  - Still high quality (exceeds PLOS minimum)")
print("  - LARGER in physical size (better readability)")
print("  - BIGGER fonts (more readable)")
print("  - Manageable for combining into multi-panel figures")
print("  - No memory errors when combining!")