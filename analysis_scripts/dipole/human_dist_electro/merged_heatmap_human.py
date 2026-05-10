import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle

# ==============================================================================
# READ ALL THREE CSV FILES
# ==============================================================================
print("\n" + "="*70)
print("LOADING HUMAN DATA FROM THREE SOURCES")
print("="*70)

# 1. Human activation energy data
human_df = pd.read_csv('human_data.csv')
print(f"✓ Human ΔG* data: {len(human_df)} rows")

# 2. Distance measurements
dist_df = pd.read_csv('human_distances.csv')
print(f"✓ Distance data: {len(dist_df)} rows")

# 3. H-bond data
hbond_df = pd.read_csv('human_hbonds.csv')
print(f"✓ H-bond data: {len(hbond_df)} rows")

# ==============================================================================
# EXTRACT RESIDUE NUMBERS FROM MUTATION NAMES
# ==============================================================================
def extract_residue_number(mutation):
    """Extract residue number from mutation name like 'Y48F' -> '48'"""
    import re
    match = re.search(r'(\d+)', mutation)
    return match.group(1) if match else None

# Add residue numbers to human data
human_df['residue'] = human_df['mutation'].apply(extract_residue_number)

# Calculate average ΔG* across all levels for each mutation
avg_dg = human_df[~human_df['mutation'].isin(['humansec', 'humancys'])].groupby('residue').agg({
    'dg_star': 'mean',
    'dg0': 'mean',
    'mutation': 'first'
}).reset_index()

print(f"\n✓ Calculated average ΔG* for {len(avg_dg)} mutations")

# ==============================================================================
# PROCESS DISTANCE DATA
# ==============================================================================
# Pivot distance data to get one row per mutation
dist_pivot = dist_df.pivot_table(
    index=['Mutation_Residue', 'Mutation_ResName'],
    columns='Active_Site_Residue',
    values='CA_to_CA_Distance_Angstrom'
).reset_index()

dist_pivot.columns.name = None
dist_pivot.columns = ['Residue', 'ResName', 'Dist_49', 'Dist_83', 'Dist_198']
dist_pivot['Residue'] = dist_pivot['Residue'].astype(str)

print(f"✓ Processed distances for {len(dist_pivot)} mutations")

# ==============================================================================
# PROCESS H-BOND DATA
# ==============================================================================
# Remove self-interactions (SEC49 to SEC49)
hbond_df = hbond_df[~((hbond_df['Donor_Residue'] == hbond_df['Acceptor_Residue']) & 
                       (hbond_df['Donor_Residue'] == '49'))]

# Count H-bonds involving each residue
hbond_counts = {}
active_site_hbonds = {}

for _, row in hbond_df.iterrows():
    donor = str(row['Donor_Residue'])
    acceptor = str(row['Acceptor_Residue'])
    
    hbond_counts[donor] = hbond_counts.get(donor, 0) + 1
    hbond_counts[acceptor] = hbond_counts.get(acceptor, 0) + 1
    
    # Check if involves active site (49, 83, 198)
    if donor in ['49', '83', '198'] or acceptor in ['49', '83', '198']:
        if donor in ['49', '83', '198']:
            active_site_hbonds[acceptor] = active_site_hbonds.get(acceptor, 0) + 1
        else:
            active_site_hbonds[donor] = active_site_hbonds.get(donor, 0) + 1

print(f"✓ Counted H-bonds for {len(hbond_counts)} residues")

# ==============================================================================
# MERGE ALL DATA
# ==============================================================================
# Merge activation energy with distances
merged_df = avg_dg.merge(dist_pivot, left_on='residue', right_on='Residue', how='inner')

# Add H-bond counts
merged_df['Total_HBonds'] = merged_df['residue'].map(hbond_counts).fillna(0).astype(int)
merged_df['Active_Site_HBonds'] = merged_df['residue'].map(active_site_hbonds).fillna(0).astype(int)

# Calculate minimum distance
merged_df['Min_Dist'] = merged_df[['Dist_49', 'Dist_83', 'Dist_198']].min(axis=1)

# Sort by minimum distance
merged_df = merged_df.sort_values('Min_Dist')

print(f"\n✓ MERGED DATASET: {len(merged_df)} mutations with complete data")
print("\nColumns:", list(merged_df.columns))

# ==============================================================================
# CREATE MERGED HEATMAP
# ==============================================================================
print("\n" + "="*70)
print("CREATING COMPREHENSIVE MERGED HEATMAP")
print("="*70)

fig, axes = plt.subplots(1, 7, figsize=(22, 12), 
                         gridspec_kw={'width_ratios': [1.2, 1, 1, 1, 1, 0.6, 0.6]})

# Prepare labels
residue_labels = [f"{row['ResName']}{row['residue']}" for _, row in merged_df.iterrows()]

# ==============================================================================
# Column 1: ΔG‡ (Human) - SAME COLOR SCALE AS YOUR HEATMAP
# ==============================================================================
data_dg = merged_df[['dg_star']].values
im1 = axes[0].imshow(data_dg, aspect='auto', cmap='RdBu_r', vmin=12, vmax=28)
axes[0].set_yticks(range(len(residue_labels)))
axes[0].set_yticklabels(residue_labels, fontsize=11, fontweight='bold')
axes[0].set_xticks([0])
axes[0].set_xticklabels(['ΔG‡\nHuman\n(kcal/mol)'], fontsize=11, fontweight='bold')
axes[0].set_title('Activation\nBarrier\n(Averaged)', fontsize=13, fontweight='bold', pad=10)

for i in range(len(data_dg)):
    value = data_dg[i,0]
    if value < 16:
        text_color = 'darkblue'
    elif value > 20:
        text_color = 'darkred'
    else:
        text_color = 'black'
    axes[0].text(0, i, f'{value:.1f}', 
                ha='center', va='center', fontsize=10, fontweight='bold', color=text_color)

cbar1 = plt.colorbar(im1, ax=axes[0], fraction=0.046, pad=0.04)
cbar1.set_label('ΔG‡ (kcal/mol)', fontsize=10, fontweight='bold')

# ==============================================================================
# Column 2: ΔG₀ (Ground State)
# ==============================================================================
data_dg0 = merged_df[['dg0']].values
im2 = axes[1].imshow(data_dg0, aspect='auto', cmap='RdBu', vmin=-65, vmax=-52)
axes[1].set_yticks([])
axes[1].set_xticks([0])
axes[1].set_xticklabels(['ΔG₀\nHuman\n(kcal/mol)'], fontsize=11, fontweight='bold')
axes[1].set_title('Ground State\nEnergy\n(Averaged)', fontsize=13, fontweight='bold', pad=10)

for i in range(len(data_dg0)):
    value = data_dg0[i,0]
    axes[1].text(0, i, f'{value:.1f}', 
                ha='center', va='center', fontsize=9, fontweight='bold', color='white')

cbar2 = plt.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.04)
cbar2.set_label('ΔG₀ (kcal/mol)', fontsize=10, fontweight='bold')

# ==============================================================================
# Column 3: Distance to SEC49
# ==============================================================================
data_49 = merged_df[['Dist_49']].values
im3 = axes[2].imshow(data_49, aspect='auto', cmap='viridis_r', vmin=0, vmax=25)
axes[2].set_yticks([])
axes[2].set_xticks([0])
axes[2].set_xticklabels(['Dist to\nSEC49\n(Å)'], fontsize=11, fontweight='bold')
axes[2].set_title('Distance to\nSEC49', fontsize=13, fontweight='bold', pad=10)

for i in range(len(data_49)):
    if not np.isnan(data_49[i,0]):
        color = 'white' if data_49[i,0] > 12 else 'black'
        axes[2].text(0, i, f'{data_49[i,0]:.1f}', 
                    ha='center', va='center', fontsize=9, fontweight='bold', color=color)

cbar3 = plt.colorbar(im3, ax=axes[2], fraction=0.046, pad=0.04)
cbar3.set_label('Distance (Å)', fontsize=10, fontweight='bold')

# ==============================================================================
# Column 4: Distance to GLN83
# ==============================================================================
data_83 = merged_df[['Dist_83']].values
im4 = axes[3].imshow(data_83, aspect='auto', cmap='viridis_r', vmin=0, vmax=25)
axes[3].set_yticks([])
axes[3].set_xticks([0])
axes[3].set_xticklabels(['Dist to\nGLN83\n(Å)'], fontsize=11, fontweight='bold')
axes[3].set_title('Distance to\nGLN83', fontsize=13, fontweight='bold', pad=10)

for i in range(len(data_83)):
    if not np.isnan(data_83[i,0]):
        color = 'white' if data_83[i,0] > 12 else 'black'
        axes[3].text(0, i, f'{data_83[i,0]:.1f}', 
                    ha='center', va='center', fontsize=9, fontweight='bold', color=color)

cbar4 = plt.colorbar(im4, ax=axes[3], fraction=0.046, pad=0.04)
cbar4.set_label('Distance (Å)', fontsize=10, fontweight='bold')

# ==============================================================================
# Column 5: Distance to RES198
# ==============================================================================
data_198 = merged_df[['Dist_198']].values
im5 = axes[4].imshow(data_198, aspect='auto', cmap='viridis_r', vmin=0, vmax=25)
axes[4].set_yticks([])
axes[4].set_xticks([0])
axes[4].set_xticklabels(['Dist to\nRES198\n(Å)'], fontsize=11, fontweight='bold')
axes[4].set_title('Distance to\nRES198', fontsize=13, fontweight='bold', pad=10)

for i in range(len(data_198)):
    if not np.isnan(data_198[i,0]):
        color = 'white' if data_198[i,0] > 12 else 'black'
        axes[4].text(0, i, f'{data_198[i,0]:.1f}', 
                    ha='center', va='center', fontsize=9, fontweight='bold', color=color)

cbar5 = plt.colorbar(im5, ax=axes[4], fraction=0.046, pad=0.04)
cbar5.set_label('Distance (Å)', fontsize=10, fontweight='bold')

# ==============================================================================
# Column 6: Total H-bonds
# ==============================================================================
data_hb = merged_df[['Total_HBonds']].values
im6 = axes[5].imshow(data_hb, aspect='auto', cmap='Greens', vmin=0, vmax=6)
axes[5].set_yticks([])
axes[5].set_xticks([0])
axes[5].set_xticklabels(['Total\nH-bonds'], fontsize=10, fontweight='bold')
axes[5].set_title('H-bonds\n(Total)', fontsize=13, fontweight='bold', pad=10)

for i in range(len(data_hb)):
    value = int(data_hb[i,0])
    color = 'white' if value > 3 else 'black'
    axes[5].text(0, i, f'{value}', 
                ha='center', va='center', fontsize=10, fontweight='bold', color=color)

# ==============================================================================
# Column 7: Active Site H-bonds
# ==============================================================================
data_hb_active = merged_df[['Active_Site_HBonds']].values
im7 = axes[6].imshow(data_hb_active, aspect='auto', cmap='Oranges', vmin=0, vmax=6)
axes[6].set_yticks([])
axes[6].set_xticks([0])
axes[6].set_xticklabels(['To Active\nSite'], fontsize=10, fontweight='bold')
axes[6].set_title('H-bonds to\nActive Site', fontsize=13, fontweight='bold', pad=10)

for i in range(len(data_hb_active)):
    value = int(data_hb_active[i,0])
    color = 'white' if value > 3 else 'black'
    axes[6].text(0, i, f'{value}', 
                ha='center', va='center', fontsize=10, fontweight='bold', color=color)

# ==============================================================================
# TITLE
# ==============================================================================
plt.suptitle('Comprehensive Human Enzyme Analysis: Activation Barriers, Distances, and H-bonding\n' +
             '(Sorted by distance | Red = High ΔG‡ [BAD] | Blue = Low ΔG‡ [GOOD] | ΔG* averaged across all mutation levels)', 
             fontsize=15, fontweight='bold', y=0.98)

plt.tight_layout()
plt.savefig('final_merged_heatmap_human.png', dpi=300, bbox_inches='tight')
print("\n✓ Saved: final_merged_heatmap_human.png")

# ==============================================================================
# CORRELATION ANALYSIS
# ==============================================================================
print("\n" + "="*70)
print("CORRELATION ANALYSIS (HUMAN)")
print("="*70)

corr_dg_49 = merged_df['dg_star'].corr(merged_df['Dist_49'])
corr_dg_83 = merged_df['dg_star'].corr(merged_df['Dist_83'])
corr_dg_198 = merged_df['dg_star'].corr(merged_df['Dist_198'])
corr_dg_hb = merged_df['dg_star'].corr(merged_df['Total_HBonds'])
corr_dg_active_hb = merged_df['dg_star'].corr(merged_df['Active_Site_HBonds'])

print(f"\nΔG‡ vs Distance to SEC49:      r = {corr_dg_49:.3f}")
print(f"\nΔG‡ vs Distance to GLN83:      r = {corr_dg_83:.3f}")
print(f"ΔG‡ vs Distance to RES198:     r = {corr_dg_198:.3f}")
print(f"ΔG‡ vs Total H-bonds:          r = {corr_dg_hb:.3f}")
print(f"ΔG‡ vs Active Site H-bonds:    r = {corr_dg_active_hb:.3f}")

# ==============================================================================
# KEY FINDINGS
# ==============================================================================
print("\n" + "="*70)
print("KEY FINDINGS (HUMAN)")
print("="*70)

print("\n🔴 HIGHEST ΔG‡ (Worst for catalysis - RED in heatmap):")
worst = merged_df.nlargest(5, 'dg_star')
for _, row in worst.iterrows():
    print(f"  {row['ResName']}{row['residue']}: ΔG‡={row['dg_star']:.1f} kcal/mol | " +
          f"Min Dist={row['Min_Dist']:.1f} Å | H-bonds to active={int(row['Active_Site_HBonds'])}")

print("\n🔵 LOWEST ΔG‡ (Best for catalysis - BLUE in heatmap):")
best = merged_df.nsmallest(5, 'dg_star')
for _, row in best.iterrows():
    print(f"  {row['ResName']}{row['residue']}: ΔG‡={row['dg_star']:.1f} kcal/mol | " +
          f"Min Dist={row['Min_Dist']:.1f} Å | H-bonds to active={int(row['Active_Site_HBonds'])}")

print("\n💚 MOST H-bonded to active site:")
top_hb = merged_df.nlargest(5, 'Active_Site_HBonds')
for _, row in top_hb.iterrows():
    if row['Active_Site_HBonds'] > 0:
        print(f"  {row['ResName']}{row['residue']}: {int(row['Active_Site_HBonds'])} H-bonds | " +
              f"ΔG‡={row['dg_star']:.1f} kcal/mol | Min Dist={row['Min_Dist']:.1f} Å")

# ==============================================================================
# SAVE MERGED DATA
# ==============================================================================
merged_df.to_csv('merged_analysis_human.csv', index=False)
print(f"\n✓ Saved merged data to: merged_analysis_human.csv")

print("\n" + "="*70)
print("ANALYSIS COMPLETE!")
print("="*70)
print("\nInterpretation:")
print("• RED (high ΔG‡) = mutations that INCREASE barrier → WORSE for catalysis")
print("• BLUE (low ΔG‡) = mutations that DECREASE barrier → BETTER for catalysis")
print("• Close residues have stronger effects on activation energy")
print("• H-bonds to active site indicate electrostatic coupling")
print("• Human WT baseline = 16.40 kcal/mol")
print("="*70 + "\n")

plt.show()
