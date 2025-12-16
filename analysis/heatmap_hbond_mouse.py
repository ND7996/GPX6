import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle

# ==============================================================================
# FILE PATHS
# ==============================================================================
DISTANCE_CSV = r"D:\PhD_Thesis\GPX6\analysis\mutation_distances.csv"
HBOND_CSV = r"D:\PhD_Thesis\GPX6\analysis\mutation_hbonds.csv"

# ==============================================================================
# READ CSV FILES
# ==============================================================================
print("\n" + "="*70)
print("LOADING DATA FROM CSV FILES")
print("="*70)

# Read distance data
try:
    dist_df = pd.read_csv(DISTANCE_CSV)
    print(f"✓ Loaded distances: {len(dist_df)} rows")
    print(f"  Columns: {list(dist_df.columns)}")
except Exception as e:
    print(f"✗ Error loading {DISTANCE_CSV}: {e}")
    exit(1)

# Read H-bond data
try:
    hbond_df = pd.read_csv(HBOND_CSV)
    print(f"✓ Loaded H-bonds: {len(hbond_df)} rows")
    print(f"  Columns: {list(hbond_df.columns)}")
except Exception as e:
    print(f"✗ Error loading {HBOND_CSV}: {e}")
    hbond_df = pd.DataFrame()  # Empty dataframe if no H-bonds

# ==============================================================================
# PROCESS DISTANCE DATA
# ==============================================================================
print("\n" + "="*70)
print("PROCESSING DISTANCE DATA")
print("="*70)

# Get unique mutations
unique_mutations = dist_df['Mutation_Residue'].unique()
print(f"Found {len(unique_mutations)} unique mutation sites")

# Create summary dataframe - one row per mutation
mutation_summary = []

for mut_resi in unique_mutations:
    mut_data = dist_df[dist_df['Mutation_Residue'] == mut_resi]
    
    # Get mutation info
    mut_resn = mut_data['Mutation_ResName'].iloc[0]
    
    # Get distances to each active site residue
    distances = {}
    for _, row in mut_data.iterrows():
        active_site = row['Active_Site_Residue']
        dist = row['CA_to_CA_Distance_Angstrom']
        distances[f'Dist_to_{active_site}'] = dist
    
    # Calculate average distance to active site
    avg_dist = np.mean([v for v in distances.values()])
    min_dist = np.min([v for v in distances.values()])
    
    summary_row = {
        'Mutation_Residue': mut_resi,
        'Mutation_ResName': mut_resn,
        'Min_Distance': min_dist,
        'Avg_Distance': avg_dist,
    }
    summary_row.update(distances)
    mutation_summary.append(summary_row)

mut_df = pd.DataFrame(mutation_summary)
mut_df = mut_df.sort_values('Min_Distance')

print(f"Created summary for {len(mut_df)} mutations")

# ==============================================================================
# PROCESS H-BOND DATA
# ==============================================================================
print("\n" + "="*70)
print("PROCESSING H-BOND DATA")
print("="*70)

if not hbond_df.empty:
    # Count H-bonds involving each residue
    hbond_counts = {}
    active_site_hbonds = {}  # H-bonds with active site (49, 83, 196)
    
    for _, row in hbond_df.iterrows():
        donor = str(row['Donor_Residue'])
        acceptor = str(row['Acceptor_Residue'])
        
        hbond_counts[donor] = hbond_counts.get(donor, 0) + 1
        hbond_counts[acceptor] = hbond_counts.get(acceptor, 0) + 1
        
        # Check if involves active site
        if donor in ['49', '83', '196'] or acceptor in ['49', '83', '196']:
            if donor in ['49', '83', '196']:
                active_site_hbonds[acceptor] = active_site_hbonds.get(acceptor, 0) + 1
            else:
                active_site_hbonds[donor] = active_site_hbonds.get(donor, 0) + 1
    
    # Add H-bond info to mutation dataframe
    mut_df['Total_HBonds'] = mut_df['Mutation_Residue'].astype(str).map(hbond_counts).fillna(0)
    mut_df['Active_Site_HBonds'] = mut_df['Mutation_Residue'].astype(str).map(active_site_hbonds).fillna(0)
    
    print(f"Found {len(hbond_counts)} residues with H-bonds")
    print(f"Active site H-bonds: {len(active_site_hbonds)} residues")
else:
    mut_df['Total_HBonds'] = 0
    mut_df['Active_Site_HBonds'] = 0
    print("No H-bond data available")

# ==============================================================================
# CREATE DISTANCE HEATMAP
# ==============================================================================
print("\n" + "="*70)
print("CREATING DISTANCE HEATMAP")
print("="*70)

fig, axes = plt.subplots(1, 3, figsize=(16, 10))

# Prepare data
residue_labels = [f"{row['Mutation_ResName']}{row['Mutation_Residue']}" 
                  for _, row in mut_df.iterrows()]

# Column 1: Distance to residue 49
if 'Dist_to_49' in mut_df.columns:
    data1 = mut_df[['Dist_to_49']].values
    im1 = axes[0].imshow(data1, aspect='auto', cmap='viridis_r', vmin=0, vmax=25)
    axes[0].set_yticks(range(len(residue_labels)))
    axes[0].set_yticklabels(residue_labels, fontsize=10)
    axes[0].set_xticks([0])
    axes[0].set_xticklabels(['CYS49'], fontsize=11, fontweight='bold')
    axes[0].set_title('Distance to\nCYS49 (Å)', fontsize=12, fontweight='bold', pad=10)
    
    for i in range(len(data1)):
        color = 'white' if data1[i,0] > 12 else 'black'
        axes[0].text(0, i, f'{data1[i,0]:.1f}', 
                    ha='center', va='center', fontsize=9, fontweight='bold', color=color)
    
    plt.colorbar(im1, ax=axes[0], fraction=0.046, pad=0.04, label='Distance (Å)')

# Column 2: Distance to residue 83
if 'Dist_to_83' in mut_df.columns:
    data2 = mut_df[['Dist_to_83']].values
    im2 = axes[1].imshow(data2, aspect='auto', cmap='viridis_r', vmin=0, vmax=25)
    axes[1].set_yticks([])
    axes[1].set_xticks([0])
    axes[1].set_xticklabels(['GLN83'], fontsize=11, fontweight='bold')
    axes[1].set_title('Distance to\nGLN83 (Å)', fontsize=12, fontweight='bold', pad=10)
    
    for i in range(len(data2)):
        color = 'white' if data2[i,0] > 12 else 'black'
        axes[1].text(0, i, f'{data2[i,0]:.1f}', 
                    ha='center', va='center', fontsize=9, fontweight='bold', color=color)
    
    plt.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.04, label='Distance (Å)')

# Column 3: Distance to residue 196
if 'Dist_to_196' in mut_df.columns:
    data3 = mut_df[['Dist_to_196']].values
    im3 = axes[2].imshow(data3, aspect='auto', cmap='viridis_r', vmin=0, vmax=25)
    axes[2].set_yticks([])
    axes[2].set_xticks([0])
    axes[2].set_xticklabels(['RES196'], fontsize=11, fontweight='bold')
    axes[2].set_title('Distance to\nRES196 (Å)', fontsize=12, fontweight='bold', pad=10)
    
    for i in range(len(data3)):
        color = 'white' if data3[i,0] > 12 else 'black'
        axes[2].text(0, i, f'{data3[i,0]:.1f}', 
                    ha='center', va='center', fontsize=9, fontweight='bold', color=color)
    
    plt.colorbar(im3, ax=axes[2], fraction=0.046, pad=0.04, label='Distance (Å)')

plt.suptitle('Mutation Site Distances to Active Site Residues', 
             fontsize=14, fontweight='bold', y=0.98)

plt.tight_layout()
plt.savefig('mutation_distances_heatmap.png', dpi=300, bbox_inches='tight')
print("✓ Saved: mutation_distances_heatmap.png")

# ==============================================================================
# CREATE COMPREHENSIVE HEATMAP WITH H-BONDS
# ==============================================================================
print("\n" + "="*70)
print("CREATING COMPREHENSIVE HEATMAP")
print("="*70)

fig2, axes2 = plt.subplots(1, 4, figsize=(18, 10), 
                           gridspec_kw={'width_ratios': [1, 1, 0.5, 0.5]})

# Column 1: Minimum distance to active site
data_min = mut_df[['Min_Distance']].values
im1 = axes2[0].imshow(data_min, aspect='auto', cmap='viridis_r', vmin=0, vmax=25)
axes2[0].set_yticks(range(len(residue_labels)))
axes2[0].set_yticklabels(residue_labels, fontsize=10)
axes2[0].set_xticks([0])
axes2[0].set_xticklabels(['Min Dist\n(Å)'], fontsize=11, fontweight='bold')
axes2[0].set_title('Closest Distance\nto Active Site', fontsize=12, fontweight='bold', pad=10)

for i in range(len(data_min)):
    color = 'white' if data_min[i,0] > 12 else 'black'
    axes2[0].text(0, i, f'{data_min[i,0]:.1f}', 
                ha='center', va='center', fontsize=9, fontweight='bold', color=color)

plt.colorbar(im1, ax=axes2[0], fraction=0.046, pad=0.04, label='Distance (Å)')

# Column 2: Average distance to active site
data_avg = mut_df[['Avg_Distance']].values
im2 = axes2[1].imshow(data_avg, aspect='auto', cmap='viridis_r', vmin=0, vmax=25)
axes2[1].set_yticks([])
axes2[1].set_xticks([0])
axes2[1].set_xticklabels(['Avg Dist\n(Å)'], fontsize=11, fontweight='bold')
axes2[1].set_title('Average Distance\nto Active Site', fontsize=12, fontweight='bold', pad=10)

for i in range(len(data_avg)):
    color = 'white' if data_avg[i,0] > 12 else 'black'
    axes2[1].text(0, i, f'{data_avg[i,0]:.1f}', 
                ha='center', va='center', fontsize=9, fontweight='bold', color=color)

plt.colorbar(im2, ax=axes2[1], fraction=0.046, pad=0.04, label='Distance (Å)')

# Column 3: Total H-bonds
data_hbond = mut_df[['Total_HBonds']].values
im3 = axes2[2].imshow(data_hbond, aspect='auto', cmap='Greens', vmin=0, vmax=max(6, data_hbond.max()))
axes2[2].set_yticks([])
axes2[2].set_xticks([0])
axes2[2].set_xticklabels(['Total\nH-bonds'], fontsize=10, fontweight='bold')
axes2[2].set_title('H-bonds\n(Total)', fontsize=12, fontweight='bold', pad=10)

for i in range(len(data_hbond)):
    value = int(data_hbond[i,0])
    color = 'white' if value > 3 else 'black'
    axes2[2].text(0, i, f'{value}', 
                ha='center', va='center', fontsize=9, fontweight='bold', color=color)

# Column 4: Active Site H-bonds
data_active_hbond = mut_df[['Active_Site_HBonds']].values
im4 = axes2[3].imshow(data_active_hbond, aspect='auto', cmap='Oranges', 
                      vmin=0, vmax=max(6, data_active_hbond.max()))
axes2[3].set_yticks([])
axes2[3].set_xticks([0])
axes2[3].set_xticklabels(['Active\nSite'], fontsize=10, fontweight='bold')
axes2[3].set_title('H-bonds to\nActive Site', fontsize=12, fontweight='bold', pad=10)

for i in range(len(data_active_hbond)):
    value = int(data_active_hbond[i,0])
    color = 'white' if value > 3 else 'black'
    axes2[3].text(0, i, f'{value}', 
                ha='center', va='center', fontsize=9, fontweight='bold', color=color)

plt.suptitle('Mutation Analysis: Distances and H-bonding Network', 
             fontsize=14, fontweight='bold', y=0.98)

plt.tight_layout()
plt.savefig('mutation_comprehensive_heatmap.png', dpi=300, bbox_inches='tight')
print("✓ Saved: mutation_comprehensive_heatmap.png")

# ==============================================================================
# CREATE H-BOND NETWORK HEATMAP
# ==============================================================================
if not hbond_df.empty:
    print("\n" + "="*70)
    print("CREATING H-BOND NETWORK HEATMAP")
    print("="*70)
    
    fig3, ax = plt.subplots(figsize=(12, 10))
    
    # Create H-bond matrix
    all_residues = sorted(set(hbond_df['Donor_Residue'].astype(str).tolist() + 
                             hbond_df['Acceptor_Residue'].astype(str).tolist()), 
                         key=lambda x: int(x))
    
    n = len(all_residues)
    hbond_matrix = np.zeros((n, n))
    
    # Map residue to index
    res_to_idx = {res: i for i, res in enumerate(all_residues)}
    
    # Fill matrix with H-bond distances
    for _, row in hbond_df.iterrows():
        donor_idx = res_to_idx[str(row['Donor_Residue'])]
        acceptor_idx = res_to_idx[str(row['Acceptor_Residue'])]
        dist = row['Distance_Angstrom']
        
        # Symmetric matrix
        hbond_matrix[donor_idx, acceptor_idx] = dist
        hbond_matrix[acceptor_idx, donor_idx] = dist
    
    # Get residue names
    res_labels = []
    for res in all_residues:
        if res in ['49', '83', '196']:
            # Active site - mark with *
            donor_match = hbond_df[hbond_df['Donor_Residue'].astype(str) == res]
            acceptor_match = hbond_df[hbond_df['Acceptor_Residue'].astype(str) == res]
            
            if len(donor_match) > 0:
                resname = donor_match['Donor_ResName'].iloc[0]
            elif len(acceptor_match) > 0:
                resname = acceptor_match['Acceptor_ResName'].iloc[0]
            else:
                resname = "UNK"
            res_labels.append(f"{resname}{res}*")
        else:
            donor_match = hbond_df[hbond_df['Donor_Residue'].astype(str) == res]
            acceptor_match = hbond_df[hbond_df['Acceptor_Residue'].astype(str) == res]
            
            if len(donor_match) > 0:
                resname = donor_match['Donor_ResName'].iloc[0]
            elif len(acceptor_match) > 0:
                resname = acceptor_match['Acceptor_ResName'].iloc[0]
            else:
                resname = "UNK"
            res_labels.append(f"{resname}{res}")
    
    # Mask zeros
    masked_matrix = np.ma.masked_where(hbond_matrix == 0, hbond_matrix)
    
    # Plot
    im = ax.imshow(masked_matrix, cmap='RdYlGn_r', vmin=2.0, vmax=3.5, aspect='auto')
    
    # Set ticks
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(res_labels, rotation=45, ha='right', fontsize=10)
    ax.set_yticklabels(res_labels, fontsize=10)
    
    # Add values and boxes
    for i in range(n):
        for j in range(n):
            if hbond_matrix[i, j] > 0:
                ax.text(j, i, f'{hbond_matrix[i,j]:.2f}', 
                       ha='center', va='center', fontsize=8, fontweight='bold')
                # Add box around active site interactions
                if (all_residues[i] in ['49', '83', '196'] or 
                    all_residues[j] in ['49', '83', '196']):
                    rect = Rectangle((j-0.5, i-0.5), 1, 1, 
                                   fill=False, edgecolor='red', linewidth=2)
                    ax.add_patch(rect)
    
    plt.colorbar(im, ax=ax, label='H-bond Distance (Å)')
    ax.set_title('Hydrogen Bond Network\n(* = Active Site Residue, Red Box = Active Site Interaction)', 
                fontsize=14, fontweight='bold', pad=15)
    
    plt.tight_layout()
    plt.savefig('hbond_network_heatmap.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: hbond_network_heatmap.png")

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================
print("\n" + "="*70)
print("SUMMARY STATISTICS")
print("="*70)

print("\nMutations closest to active site:")
closest = mut_df.nsmallest(5, 'Min_Distance')[['Mutation_ResName', 'Mutation_Residue', 'Min_Distance', 'Total_HBonds', 'Active_Site_HBonds']]
print(closest.to_string(index=False))

print("\nMutations farthest from active site:")
farthest = mut_df.nlargest(5, 'Min_Distance')[['Mutation_ResName', 'Mutation_Residue', 'Min_Distance', 'Total_HBonds', 'Active_Site_HBonds']]
print(farthest.to_string(index=False))

if not hbond_df.empty:
    print("\nMutations with most H-bonds:")
    top_hbonds = mut_df.nlargest(5, 'Total_HBonds')[['Mutation_ResName', 'Mutation_Residue', 'Total_HBonds', 'Active_Site_HBonds', 'Min_Distance']]
    print(top_hbonds.to_string(index=False))
    
    print("\nMutations with most active site H-bonds:")
    top_active_hbonds = mut_df.nlargest(5, 'Active_Site_HBonds')[['Mutation_ResName', 'Mutation_Residue', 'Active_Site_HBonds', 'Total_HBonds', 'Min_Distance']]
    print(top_active_hbonds.to_string(index=False))

print("\n" + "="*70)
print("ANALYSIS COMPLETE")
print("="*70)
print("\nGenerated files:")
print("  1. mutation_distances_heatmap.png")
print("  2. mutation_comprehensive_heatmap.png")
if not hbond_df.empty:
    print("  3. hbond_network_heatmap.png")
print("="*70 + "\n")

plt.show()