import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime
import matplotlib.patches as mpatches

# ===================== FILE SETTINGS =====================
file_path = r'D:\PhD_Thesis\GPX6\analysis\alignment\MEGA-result1.csv'
save_dir = os.path.dirname(file_path)
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

# ===================== LOAD DATA =====================
df = pd.read_csv(file_path)
print("CSV loaded:", df.shape)

# ===================== PARSE CSV =====================
human_levels, human_distances = [], []
mouse_levels, mouse_distances = [], []

for i in range(len(df)):
    name = df.iloc[i, 0]
    dist = df.iloc[i, -1]

    if name.lower() == "ancestor" or "node" in name.lower():
        continue

    if "HUMAN" in name.upper():
        for part in name.split("_"):
            if "level" in part.lower():
                level = int(part.lower().replace("level", ""))
                human_levels.append(level)
                human_distances.append(dist)
                break
    elif "MOUSE" in name.upper():
        for part in name.split("_"):
            if "level" in part.lower():
                level = int(part.lower().replace("level", ""))
                mouse_levels.append(level)
                mouse_distances.append(dist)
                break

# ===================== SORT =====================
human_data = sorted(zip(human_levels, human_distances))
mouse_data = sorted(zip(mouse_levels, mouse_distances))

hL = np.array([x[0] for x in human_data])
hD = np.array([x[1] for x in human_data])
mL = np.array([x[0] for x in mouse_data])
mD = np.array([x[1] for x in mouse_data])

# ===================== NORMALIZE TO START FROM ZERO =====================
dmin = min(hD.min(), mD.min())
hD0 = hD - dmin
mD0 = mD - dmin

# Find closest to ancestor
human_min_idx = np.argmin(hD)
mouse_min_idx = np.argmin(mD)

print(f"\nClosest to Ancestor:")
print(f"Human: Level {hL[human_min_idx]} (distance = {hD[human_min_idx]:.6f})")
print(f"Mouse: Level {mL[mouse_min_idx]} (distance = {mD[mouse_min_idx]:.6f})")

# ===================== COLOR GRADIENT FUNCTIONS =====================
def get_human_gradient(distance, species_distances):
    """
    Returns color for human pathway based on distance from ancestor:
    - Closest to ancestor (smallest distance) = Dark green (#006400)
    - Farthest from ancestor (largest distance) = Light green (not white, to stay visible)
    """
    min_dist = np.min(species_distances)
    max_dist = np.max(species_distances)
    
    # Normalize distance to 0-1 range
    if max_dist == min_dist:
        normalized = 0
    else:
        normalized = (distance - min_dist) / (max_dist - min_dist)
    
    # Clamp normalized value
    normalized = np.clip(normalized, 0, 1)
    
    # Dark green (0, 0.392, 0) to Light green (0.7, 0.95, 0.7) - NOT pure white
    r = 0.0 + 0.7 * normalized
    g = 0.392 + 0.558 * normalized  # Goes to 0.95
    b = 0.0 + 0.7 * normalized
    
    return (r, g, b)

def get_mouse_gradient(distance, species_distances):
    """
    Returns color for mouse pathway based on distance from ancestor:
    - Closest to ancestor (smallest distance) = Dark green (#006400)
    - Farthest from ancestor (largest distance) = Light green (not white, to stay visible)
    """
    min_dist = np.min(species_distances)
    max_dist = np.max(species_distances)
    
    # Normalize distance to 0-1 range
    if max_dist == min_dist:
        normalized = 0
    else:
        normalized = (distance - min_dist) / (max_dist - min_dist)
    
    # Clamp normalized value
    normalized = np.clip(normalized, 0, 1)
    
    # Dark green (0, 0.392, 0) to Light green (0.7, 0.95, 0.7) - NOT pure white
    r = 0.0 + 0.7 * normalized
    g = 0.392 + 0.558 * normalized  # Goes to 0.95
    b = 0.0 + 0.7 * normalized
    
    return (r, g, b)

# ===================== FIGURE 1: RADIAL LINES FROM ANCESTOR =====================
print("\nCreating Figure 1: Radial pathways from ancestor...")

fig1, ax1 = plt.subplots(figsize=(15, 8))

# Plot connecting lines from ancestor to all points
ancestor_x, ancestor_y = 0, 0

# Human pathway lines - radial from ancestor
for i in range(len(hL)):
    color = get_human_gradient(hD0[i], hD0)
    ax1.plot([ancestor_x, hL[i]], [ancestor_y, hD0[i]], 
            color=color, alpha=0.8, linewidth=2.5, zorder=1)

# Mouse pathway lines - radial from ancestor
for i in range(len(mL)):
    color = get_mouse_gradient(mD0[i], mD0)
    ax1.plot([ancestor_x, mL[i]], [ancestor_y, mD0[i]], 
            color=color, alpha=0.8, linewidth=2.5, zorder=1)

# Plot human points with gradient
for i, (x, y) in enumerate(zip(hL, hD0)):
    color = get_human_gradient(hD0[i], hD0)
    ax1.scatter(x, y, s=400, c=[color], edgecolors='black', 
               linewidths=2, zorder=3, marker='o')
    ax1.text(x, y, str(hL[i]), ha='center', va='center', 
            fontsize=8, fontweight='bold', color='black', zorder=5)
    # Add distance label below the point
    ax1.text(x, y - 0.0015, f'{hD0[i]:.4f}', ha='center', va='top', 
            fontsize=6, color='darkslategray', zorder=5, style='italic')

# Plot mouse points with gradient
for i, (x, y) in enumerate(zip(mL, mD0)):
    color = get_mouse_gradient(mD0[i], mD0)
    ax1.scatter(x, y, s=400, c=[color], edgecolors='black', 
               linewidths=2, zorder=3, marker='o')
    ax1.text(x, y, str(mL[i]), ha='center', va='center', 
            fontsize=8, fontweight='bold', color='black', zorder=5)
    # Add distance label below the point
    ax1.text(x, y - 0.0015, f'{mD0[i]:.4f}', ha='center', va='top', 
            fontsize=6, color='darkslategray', zorder=5, style='italic')

# Plot ancestor
ax1.scatter(0, 0, s=500, c='#FFD700', edgecolors='black', 
           linewidths=2.5, zorder=4, marker='h')
ax1.text(0, 0, 'A', ha='center', va='center', 
        fontsize=11, fontweight='bold', color='black', zorder=5)

# Legend for Figure 1 - Create gradient colored patches with more samples
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

# Create a custom legend with gradient samples showing the full spectrum
# More samples for better visualization of the gradient (showing 8 colors instead of 6)
gradient_samples = [
    (get_human_gradient(0.0, np.array([0, 1])), 'Closest to Ancestor'),
    (get_human_gradient(0.15, np.array([0, 1])), ''),
    (get_human_gradient(0.3, np.array([0, 1])), ''),
    (get_human_gradient(0.45, np.array([0, 1])), ''),
    (get_human_gradient(0.6, np.array([0, 1])), ''),
    (get_human_gradient(0.75, np.array([0, 1])), ''),
    (get_human_gradient(0.9, np.array([0, 1])), ''),
    (get_human_gradient(1.0, np.array([0, 1])), 'Farthest from Ancestor')
]

gradient_patches = []
for color, label in gradient_samples:
    patch = mpatches.Patch(color=color, label=label if label else None, 
                          edgecolor='black', linewidth=0.8)
    gradient_patches.append(patch)

ancestor_marker = plt.Line2D([0], [0], marker='h', color='w', markerfacecolor='#FFD700', 
                            markersize=12, label='Ancestor (A)', markeredgecolor='black', 
                            markeredgewidth=1.2)

legend1 = ax1.legend(handles=gradient_patches + [ancestor_marker], 
                   loc='center right', fontsize=8, framealpha=0.97, 
                   edgecolor='black', title='Legend',
                   title_fontsize=9, borderpad=0.4, labelspacing=0.4,
                   handlelength=2, handletextpad=0.6, columnspacing=1.0)
legend1.get_title().set_fontweight('bold')
legend1.get_title().set_fontweight('bold')

# Add text labels to identify human and mouse sequences
ax1.text(1, hD0[0] + 0.001, 'Human', fontsize=10, fontweight='bold', 
        color='darkgreen', ha='center', va='bottom')
ax1.text(1, mD0[0] - 0.001, 'Mouse', fontsize=10, fontweight='bold', 
        color='darkgreen', ha='center', va='top')

# Formatting Figure 1
ax1.set_xlabel('Level Number', fontsize=13, fontweight='bold')
ax1.set_ylabel('Distance from Ancestor', fontsize=13, fontweight='bold')
ax1.set_title('GPX6 Evolutionary Pathways from Ancestor\n(Radial View - Single Green Gradient)', 
             fontsize=15, fontweight='bold', pad=15)

all_levels = list(range(0, max(max(hL), max(mL)) + 1))
ax1.set_xticks(all_levels)
ax1.set_xticklabels(all_levels, fontsize=10)
ax1.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
ax1.set_xlim(-0.5, max(max(hL), max(mL)) + 0.5)
y_max = max(max(hD0), max(mD0))
ax1.set_ylim(-0.001, y_max * 1.1)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

plt.tight_layout()

# Save Figure 1
fig1_path = os.path.join(save_dir, f'GPX6_radial_colored_{timestamp}.png')
plt.savefig(fig1_path, dpi=300, bbox_inches='tight')
print(f"Figure 1 saved: {fig1_path}")
plt.close(fig1)

# ===================== FIGURE 2: SEQUENTIAL PATHWAY =====================
print("\nCreating Figure 2: Sequential level-to-level pathways...")

fig2, ax2 = plt.subplots(figsize=(15, 8))

# Human pathway - connect points sequentially
if len(hL) > 1:
    # Create full arrays including ancestor
    hL_full = np.concatenate([[0], hL])
    hD0_full = np.concatenate([[0], hD0])
    
    # Plot line segments with colors
    for i in range(len(hL_full) - 1):
        color = get_human_gradient(hD0_full[i+1], hD0)
        ax2.plot([hL_full[i], hL_full[i+1]], [hD0_full[i], hD0_full[i+1]], 
                color=color, alpha=0.85, linewidth=3, zorder=1)

# Mouse pathway - connect points sequentially
if len(mL) > 1:
    # Create full arrays including ancestor
    mL_full = np.concatenate([[0], mL])
    mD0_full = np.concatenate([[0], mD0])
    
    # Plot line segments with colors
    for i in range(len(mL_full) - 1):
        color = get_mouse_gradient(mD0_full[i+1], mD0)
        ax2.plot([mL_full[i], mL_full[i+1]], [mD0_full[i], mD0_full[i+1]], 
                color=color, alpha=0.85, linewidth=3, zorder=1)

# Plot human points with gradient
for i, (x, y) in enumerate(zip(hL, hD0)):
    color = get_human_gradient(hD0[i], hD0)
    ax2.scatter(x, y, s=400, c=[color], edgecolors='black', 
               linewidths=2, zorder=3, marker='o')
    ax2.text(x, y, str(hL[i]), ha='center', va='center', 
            fontsize=8, fontweight='bold', color='black', zorder=5)
    # Add distance label below the point
    ax2.text(x, y - 0.0015, f'{hD0[i]:.4f}', ha='center', va='top', 
            fontsize=6, color='darkslategray', zorder=5, style='italic')

# Plot mouse points with gradient
for i, (x, y) in enumerate(zip(mL, mD0)):
    color = get_mouse_gradient(mD0[i], mD0)
    ax2.scatter(x, y, s=400, c=[color], edgecolors='black', 
               linewidths=2, zorder=3, marker='o')
    ax2.text(x, y, str(mL[i]), ha='center', va='center', 
            fontsize=8, fontweight='bold', color='black', zorder=5)
    # Add distance label below the point
    ax2.text(x, y - 0.0015, f'{mD0[i]:.4f}', ha='center', va='top', 
            fontsize=6, color='darkslategray', zorder=5, style='italic')

# Plot ancestor
ax2.scatter(0, 0, s=500, c='#FFD700', edgecolors='black', 
           linewidths=2.5, zorder=4, marker='h')
ax2.text(0, 0, 'A', ha='center', va='center', 
        fontsize=11, fontweight='bold', color='black', zorder=5)

# Legend for Figure 2
legend2 = ax2.legend(handles=gradient_patches + [ancestor_marker], 
                   loc='center right', fontsize=8, framealpha=0.97, 
                   edgecolor='black', title='Legend',
                   title_fontsize=9, borderpad=0.4, labelspacing=0.4,
                   handlelength=2, handletextpad=0.6, columnspacing=1.0)
legend2.get_title().set_fontweight('bold')
legend2.get_title().set_fontweight('bold')

# Add text labels to identify human and mouse sequences
ax2.text(1, hD0[0] + 0.001, 'Human', fontsize=10, fontweight='bold', 
        color='darkgreen', ha='center', va='bottom')
ax2.text(1, mD0[0] - 0.001, 'Mouse', fontsize=10, fontweight='bold', 
        color='darkgreen', ha='center', va='top')

# Formatting Figure 2
ax2.set_xlabel('Level Number', fontsize=13, fontweight='bold')
ax2.set_ylabel('Distance from Ancestor', fontsize=13, fontweight='bold')
ax2.set_title('GPX6 Evolutionary Pathways from Ancestor)', 
             fontsize=15, fontweight='bold', pad=15)

ax2.set_xticks(all_levels)
ax2.set_xticklabels(all_levels, fontsize=10)
ax2.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
ax2.set_xlim(-0.5, max(max(hL), max(mL)) + 0.5)
ax2.set_ylim(-0.001, y_max * 1.1)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

plt.tight_layout()

# Save Figure 2
fig2_path = os.path.join(save_dir, f'GPX6_sequential_colored_{timestamp}.png')
plt.savefig(fig2_path, dpi=300, bbox_inches='tight')
print(f"Figure 2 saved: {fig2_path}")
plt.show()

# ===================== DETAILED STATISTICS =====================
print("\n" + "="*70)
print("DETAILED ANALYSIS")
print("="*70)

print(f"\nHuman GPX6:")
closest_h_indices = np.argsort(hD0)[:3]
farthest_h_indices = np.argsort(hD0)[-3:]
print(f"  Closest levels: {hL[closest_h_indices].tolist()}")
print(f"  Farthest levels: {hL[farthest_h_indices].tolist()}")
print(f"  Mean distance: {np.mean(hD):.6f}")

print(f"\nMouse GPX6:")
closest_m_indices = np.argsort(mD0)[:3]
farthest_m_indices = np.argsort(mD0)[-3:]
print(f"  Closest levels: {mL[closest_m_indices].tolist()}")
print(f"  Farthest levels: {mL[farthest_m_indices].tolist()}")
print(f"  Mean distance: {np.mean(mD):.6f}")

print("\n" + "="*70)
print(f"\nBoth figures saved to: {save_dir}")
print(f"  - Figure 1 (Radial): GPX6_radial_colored_{timestamp}.png")
print(f"  - Figure 2 (Sequential): GPX6_sequential_colored_{timestamp}.png")
print("="*70)