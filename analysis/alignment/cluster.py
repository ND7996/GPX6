"""
GPX6 Evolutionary Distance Analysis with Color-Coded Legend
Clean plot with comprehensive legend showing levels by color intensity
"""

import glob
import numpy as np
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sklearn.manifold import MDS
from scipy import stats
import re

print("="*70)
print("GPX6 Evolutionary Pathway Analysis with Color Legend")
print("="*70)

###########################################
# 1. Load all sequences
###########################################

print("\n[1/6] Loading sequences...")

sequences = {}
seq_data = []

def extract_level(name):
    """Extract level number from sequence name"""
    match = re.search(r'level(\d+)', name, re.IGNORECASE)
    if match:
        return int(match.group(1))
    return None

# Load MOUSE sequences
mouse_files = glob.glob("D:/PhD_Thesis/GPX6/analysis/alignment/MOUSE/*.fasta")
for fasta_file in mouse_files:
    for record in SeqIO.parse(fasta_file, "fasta"):
        name = f"MOUSE_{record.id}"
        level = extract_level(fasta_file)
        sequences[name] = str(record.seq)
        seq_data.append({'name': name, 'type': 'MOUSE', 'seq': str(record.seq), 'level': level})

# Load HUMAN sequences  
human_files = glob.glob("D:/PhD_Thesis/GPX6/analysis/alignment/HUMAN/*.fasta")
for fasta_file in human_files:
    for record in SeqIO.parse(fasta_file, "fasta"):
        name = f"HUMAN_{record.id}"
        level = extract_level(fasta_file)
        sequences[name] = str(record.seq)
        seq_data.append({'name': name, 'type': 'HUMAN', 'seq': str(record.seq), 'level': level})

# Load ancestral node_25
for record in SeqIO.parse("D:/PhD_Thesis/GPX6/analysis/alignment/node_25.fasta", "fasta"):
    name = f"ANCESTOR_{record.id}"
    sequences[name] = str(record.seq)
    seq_data.append({'name': name, 'type': 'ANCESTOR', 'seq': str(record.seq), 'level': None})

print(f"  ✓ Loaded {len(sequences)} sequences")

###########################################
# 2-4. Calculate distances and MDS
###########################################

print("\n[2/6] Calculating pairwise distances...")

seq_names = list(sequences.keys())
n_seqs = len(seq_names)
distance_matrix = np.zeros((n_seqs, n_seqs))

for i in range(n_seqs):
    for j in range(i+1, n_seqs):
        seq1 = sequences[seq_names[i]]
        seq2 = sequences[seq_names[j]]
        min_len = min(len(seq1), len(seq2))
        matches = sum(1 for k in range(min_len) if seq1[k] == seq2[k])
        identity = matches / min_len if min_len > 0 else 0
        distance = 1 - identity
        distance_matrix[i, j] = distance
        distance_matrix[j, i] = distance

print(f"  ✓ Calculated {n_seqs}x{n_seqs} distance matrix")

print("\n[3/6] Performing dimensionality reduction...")
mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42, n_init=4)
coords_2d = mds.fit_transform(distance_matrix)
print("  ✓ Projected sequences into 2D space")

print("\n[4/6] Calculating distances from ancestor...")
ancestor_idx = [i for i, name in enumerate(seq_names) if 'ANCESTOR' in name][0]
ancestor_distances = distance_matrix[ancestor_idx, :]

df_data = pd.DataFrame({
    'sequence': seq_names,
    'distance_from_ancestor': ancestor_distances,
    'type': [s['type'] for s in seq_data],
    'level': [s['level'] for s in seq_data],
    'x': coords_2d[:, 0],
    'y': coords_2d[:, 1]
})

ancestor_x = df_data.loc[df_data['type'] == 'ANCESTOR', 'x'].values[0]
ancestor_y = df_data.loc[df_data['type'] == 'ANCESTOR', 'y'].values[0]

df_save = df_data[df_data['type'] != 'ANCESTOR'][['sequence', 'distance_from_ancestor', 'type', 'level']].copy()
df_save.to_csv("distances_from_ancestor_with_levels.csv", index=False)
print("  ✓ Saved: distances_from_ancestor_with_levels.csv")

###########################################
# 5. Statistical tests
###########################################

print("\n[5/6] Performing statistical tests...")

df_no_ancestor = df_data[df_data['type'] != 'ANCESTOR'].copy()
human_distances = df_no_ancestor[df_no_ancestor['type'] == 'HUMAN']['distance_from_ancestor']
mouse_distances = df_no_ancestor[df_no_ancestor['type'] == 'MOUSE']['distance_from_ancestor']

if len(human_distances) > 1 and len(mouse_distances) > 1:
    u_stat, u_pval = stats.mannwhitneyu(human_distances, mouse_distances, alternative='two-sided')
    print(f"  Human-Mouse separation: p = {u_pval:.4e}")

###########################################
# 6. Create color mapping by distance
###########################################

print("\n[6/6] Creating visualization with color-coded legend...")

# Calculate color for each level based on distance
from matplotlib import colors as mcolors

# Get human and mouse data sorted by level
human_data = df_no_ancestor[df_no_ancestor['type'] == 'HUMAN'].sort_values('level')
mouse_data = df_no_ancestor[df_no_ancestor['type'] == 'MOUSE'].sort_values('level')

# Create color mapping (darker = closer, lighter = farther)
human_colors = {}
mouse_colors = {}

if len(human_data) > 0:
    h_min = human_data['distance_from_ancestor'].min()
    h_max = human_data['distance_from_ancestor'].max()
    for _, row in human_data.iterrows():
        # Normalize distance (0 = closest, 1 = farthest)
        if h_max > h_min:
            norm_dist = (row['distance_from_ancestor'] - h_min) / (h_max - h_min)
        else:
            norm_dist = 0
        # Dark red (close) to light red (far)
        human_colors[row['level']] = mcolors.to_hex((0.5 + 0.5*norm_dist, 0, 0))

if len(mouse_data) > 0:
    m_min = mouse_data['distance_from_ancestor'].min()
    m_max = mouse_data['distance_from_ancestor'].max()
    for _, row in mouse_data.iterrows():
        if m_max > m_min:
            norm_dist = (row['distance_from_ancestor'] - m_min) / (m_max - m_min)
        else:
            norm_dist = 0
        # Dark teal (close) to light teal (far)
        mouse_colors[row['level']] = mcolors.to_hex((0, 0.3 + 0.6*norm_dist, 0.3 + 0.6*norm_dist))

# Create main figure
fig, ax = plt.subplots(figsize=(20, 12))

# Draw pathways
for idx, row in df_data.iterrows():
    if row['type'] != 'ANCESTOR':
        base_color = '#FF6B6B' if row['type'] == 'HUMAN' else '#4ECDC4'
        ax.plot([row['x'], ancestor_x], [row['y'], ancestor_y], 
                color=base_color, alpha=0.15, linewidth=1.5, zorder=1)

# Plot sequences with individual colors
for idx, row in df_data.iterrows():
    if row['type'] == 'HUMAN':
        color = human_colors.get(row['level'], '#FF6B6B')
        ax.scatter(row['x'], row['y'], c=color, s=200, alpha=0.85,
                  edgecolors='black', linewidth=1.5, zorder=3)
    elif row['type'] == 'MOUSE':
        color = mouse_colors.get(row['level'], '#4ECDC4')
        ax.scatter(row['x'], row['y'], c=color, s=200, alpha=0.85,
                  edgecolors='black', linewidth=1.5, zorder=3)

# Plot ancestor as a large circle
ax.scatter(ancestor_x, ancestor_y, c='#FFD700', s=500, 
           edgecolors='black', linewidth=3, zorder=5)

# Create comprehensive legend
legend_elements = []

# Add ancestor
legend_elements.append(mpatches.Patch(facecolor='#FFD700', edgecolor='black', 
                                      label='Ancestor', linewidth=1.5))
legend_elements.append(mpatches.Patch(facecolor='white', edgecolor='white', label=''))  # spacer

# Add human levels header
legend_elements.append(mpatches.Patch(facecolor='white', edgecolor='white', 
                                      label='━━━ Human Levels ━━━'))
for _, row in human_data.iterrows():
    level = int(row['level'])
    color = human_colors[level]
    dist = row['distance_from_ancestor']
    label = f"Level {level:2d} (d={dist:.4f})"
    legend_elements.append(mpatches.Patch(facecolor=color, edgecolor='black', 
                                         label=label, linewidth=1))

legend_elements.append(mpatches.Patch(facecolor='white', edgecolor='white', label=''))  # spacer

# Add mouse levels header
legend_elements.append(mpatches.Patch(facecolor='white', edgecolor='white', 
                                      label='━━━ Mouse Levels ━━━'))
for _, row in mouse_data.iterrows():
    level = int(row['level'])
    color = mouse_colors[level]
    dist = row['distance_from_ancestor']
    label = f"Level {level:2d} (d={dist:.4f})"
    legend_elements.append(mpatches.Patch(facecolor=color, edgecolor='black', 
                                         label=label, linewidth=1))

# Add legend with larger font
ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1.02, 0.5),
          fontsize=10, frameon=True, fancybox=True, shadow=True, 
          ncol=1, title='Levels & Distances\n(darker = closer to ancestor)',
          title_fontsize=11)

ax.set_xlabel('MDS Dimension 1', fontsize=14, fontweight='bold')
ax.set_ylabel('MDS Dimension 2', fontsize=14, fontweight='bold')
ax.set_title('Evolutionary Pathways to Ancestral Sequence\n(Color intensity indicates distance from ancestor)', 
             fontsize=16, fontweight='bold', pad=20)
ax.grid(True, alpha=0.3, linewidth=0.8)

plt.tight_layout()
plt.savefig("evolutionary_pathways_color_legend.png", dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved: evolutionary_pathways_color_legend.png")

# Create a simpler version with gradient bars
fig, ax = plt.subplots(figsize=(20, 12))

# Draw pathways
for idx, row in df_data.iterrows():
    if row['type'] != 'ANCESTOR':
        base_color = '#FF6B6B' if row['type'] == 'HUMAN' else '#4ECDC4'
        ax.plot([row['x'], ancestor_x], [row['y'], ancestor_y], 
                color=base_color, alpha=0.15, linewidth=1.5, zorder=1)

# Plot sequences
for idx, row in df_data.iterrows():
    if row['type'] == 'HUMAN':
        color = human_colors.get(row['level'], '#FF6B6B')
        ax.scatter(row['x'], row['y'], c=color, s=200, alpha=0.85,
                  edgecolors='black', linewidth=1.5, zorder=3)
    elif row['type'] == 'MOUSE':
        color = mouse_colors.get(row['level'], '#4ECDC4')
        ax.scatter(row['x'], row['y'], c=color, s=200, alpha=0.85,
                  edgecolors='black', linewidth=1.5, zorder=3)

ax.scatter(ancestor_x, ancestor_y, c='#FFD700', s=500,
           edgecolors='black', linewidth=3, zorder=5)

# Simplified legend with just gradient explanation
legend_elements = [
    mpatches.Patch(facecolor='#FFD700', edgecolor='black', 
                  label='Ancestor (node_25)', linewidth=1.5),
    mpatches.Patch(facecolor='none', edgecolor='none', label=''),
    mpatches.Patch(facecolor='#800000', edgecolor='black', 
                  label='Human - Dark Red (close to ancestor)', linewidth=1),
    mpatches.Patch(facecolor='#FF6B6B', edgecolor='black', 
                  label='Human - Light Red (far from ancestor)', linewidth=1),
    mpatches.Patch(facecolor='none', edgecolor='none', label=''),
    mpatches.Patch(facecolor='#006666', edgecolor='black', 
                  label='Mouse - Dark Teal (close to ancestor)', linewidth=1),
    mpatches.Patch(facecolor='#4ECDC4', edgecolor='black', 
                  label='Mouse - Light Teal (far from ancestor)', linewidth=1),
]

ax.legend(handles=legend_elements, loc='upper left', fontsize=12, 
          frameon=True, fancybox=True, shadow=True)

ax.set_xlabel('MDS Dimension 1', fontsize=14, fontweight='bold')
ax.set_ylabel('MDS Dimension 2', fontsize=14, fontweight='bold')
ax.set_title('Evolutionary Pathways to Ancestral Sequence\n(Color intensity indicates distance from ancestor)', 
             fontsize=16, fontweight='bold', pad=20)
ax.grid(True, alpha=0.3, linewidth=0.8)

plt.tight_layout()
plt.savefig("evolutionary_pathways_simple_legend.png", dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved: evolutionary_pathways_simple_legend.png")

###########################################
# Summary
###########################################

print("\n" + "="*70)
print("✓✓✓ ANALYSIS COMPLETE! ✓✓✓")
print("="*70)

print("\n📁 Output Files:")
print("  → distances_from_ancestor_with_levels.csv")
print("  → evolutionary_pathways_color_legend.png (detailed legend with all levels)")
print("  → evolutionary_pathways_simple_legend.png (simple gradient legend)")

print("\n📊 Key Results:")
print(f"\nTotal sequences: {len(df_no_ancestor)}")

if len(human_distances) > 0:
    closest_h = human_data.nsmallest(1, 'distance_from_ancestor').iloc[0]
    farthest_h = human_data.nlargest(1, 'distance_from_ancestor').iloc[0]
    print(f"\nHuman (n={len(human_distances)}):")
    print(f"  Mean distance: {human_distances.mean():.6f}")
    print(f"  Closest level: {int(closest_h['level'])} (d={closest_h['distance_from_ancestor']:.6f}) - DARKEST red")
    print(f"  Farthest level: {int(farthest_h['level'])} (d={farthest_h['distance_from_ancestor']:.6f}) - LIGHTEST red")

if len(mouse_distances) > 0:
    closest_m = mouse_data.nsmallest(1, 'distance_from_ancestor').iloc[0]
    farthest_m = mouse_data.nlargest(1, 'distance_from_ancestor').iloc[0]
    print(f"\nMouse (n={len(mouse_distances)}):")
    print(f"  Mean distance: {mouse_distances.mean():.6f}")
    print(f"  Closest level: {int(closest_m['level'])} (d={closest_m['distance_from_ancestor']:.6f}) - DARKEST teal")
    print(f"  Farthest level: {int(farthest_m['level'])} (d={farthest_m['distance_from_ancestor']:.6f}) - LIGHTEST teal")

if len(human_distances) > 1 and len(mouse_distances) > 1:
    print(f"\nStatistical Test:")
    print(f"  Human-Mouse separation: p = {u_pval:.4e}")
    print(f"  Result: {'SIGNIFICANT' if u_pval < 0.05 else 'NOT significant'}")

print("\n💡 Color Guide:")
print("  - Darker colors = closer to ancestor (more conserved)")
print("  - Lighter colors = farther from ancestor (more diverged)")
print("="*70)