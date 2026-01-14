"""
GPX6 Distance Plot - Side-by-Side Layout
Human | Ancestor | Mouse
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects

# Read data
df = pd.read_csv("distances_from_ancestor_with_levels.csv")
df['level'] = df['level'].astype(int)

# Separate
human_df = df[df['type'] == 'HUMAN'].sort_values('level')
mouse_df = df[df['type'] == 'MOUSE'].sort_values('level')

# Find closest
h_closest = human_df.nsmallest(1, 'distance_from_ancestor').iloc[0]
m_closest = mouse_df.nsmallest(1, 'distance_from_ancestor').iloc[0]

print("="*70)
print("Creating Side-by-Side Comparison Plot")
print("="*70)

###########################################
# Side-by-Side Layout: Human | Ancestor | Mouse
###########################################

fig, ax = plt.subplots(figsize=(20, 10))

# Positions
ancestor_x = 0
human_x_offset = -1  # Human on LEFT
mouse_x_offset = 1   # Mouse on RIGHT

# Plot HUMAN (left side) with pathways
for _, row in human_df.iterrows():
    level = row['level']
    dist = row['distance_from_ancestor']
    x_pos = human_x_offset * level
    
    # Line to ancestor
    ax.plot([x_pos, ancestor_x], [dist, 0.945], 
           color='#FF6B6B', alpha=0.15, linewidth=1, zorder=1)
    
    # Point
    ax.scatter(x_pos, dist, c='#FF6B6B', s=200, alpha=0.9,
              edgecolors='black', linewidth=1.5, zorder=5)
    
    # Label
    text = ax.text(x_pos, dist, str(level), 
                  ha='center', va='center', 
                  fontsize=9, fontweight='bold', color='white', zorder=6)
    text.set_path_effects([
        path_effects.Stroke(linewidth=2.5, foreground='black'),
        path_effects.Normal()
    ])

# Plot MOUSE (right side) with pathways
for _, row in mouse_df.iterrows():
    level = row['level']
    dist = row['distance_from_ancestor']
    x_pos = mouse_x_offset * level
    
    # Line to ancestor
    ax.plot([x_pos, ancestor_x], [dist, 0.945], 
           color='#4ECDC4', alpha=0.15, linewidth=1, zorder=1)
    
    # Point
    ax.scatter(x_pos, dist, c='#4ECDC4', s=200, alpha=0.9,
              edgecolors='black', linewidth=1.5, marker='s', zorder=5)
    
    # Label
    text = ax.text(x_pos, dist, str(level), 
                  ha='center', va='center', 
                  fontsize=9, fontweight='bold', color='white', zorder=6)
    text.set_path_effects([
        path_effects.Stroke(linewidth=2.5, foreground='black'),
        path_effects.Normal()
    ])

# Plot ANCESTOR (center)
ax.scatter(ancestor_x, 0.945, c='#FFD700', s=600,
          edgecolors='black', linewidth=3, zorder=10)
ax.text(ancestor_x, 0.945, 'A', 
       ha='center', va='center', fontsize=14, fontweight='bold', 
       color='black', zorder=11)

# Add labels for sections
ax.text(-12, 0.974, 'HUMAN', fontsize=16, fontweight='bold', 
       ha='center', color='#FF6B6B')
ax.text(12, 0.974, 'MOUSE', fontsize=16, fontweight='bold', 
       ha='center', color='#4ECDC4')
ax.text(0, 0.9435, 'Ancestor', fontsize=12, fontweight='bold', 
       ha='center', va='top', color='#B8860B')

# Styling
ax.set_xlabel('← Human Levels | Mouse Levels →', fontsize=15, fontweight='bold', labelpad=12)
ax.set_ylabel('Distance from Ancestor', fontsize=15, fontweight='bold', labelpad=12)
ax.set_title('GPX6 Evolutionary Distance: Species Comparison\n(Ancestor in center; branches show divergence)', 
            fontsize=17, fontweight='bold', pad=25)

# Legend (WITHOUT ANCESTOR - no yellow box)
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#FF6B6B', 
           markersize=11, markeredgecolor='black', markeredgewidth=1.5, 
           label='Human (left)', linestyle='None'),
    Line2D([0], [0], marker='s', color='w', markerfacecolor='#4ECDC4', 
           markersize=11, markeredgecolor='black', markeredgewidth=1.5, 
           label='Mouse (right)', linestyle='None'),
]
ax.legend(handles=legend_elements, fontsize=13, loc='upper right', 
         framealpha=0.98, edgecolor='black', shadow=True)

# Grid
ax.grid(True, alpha=0.2, linewidth=0.8, linestyle=':', color='gray')
ax.set_axisbelow(True)

# Axes
ax.set_xlim(-22, 22)
ax.set_ylim(0.944, 0.976)

# Custom x-axis ticks to show actual level numbers
human_levels = human_df['level'].unique()
mouse_levels = mouse_df['level'].unique()

# Create tick positions and labels
tick_positions = []
tick_labels = []

# Human side (negative x positions, but show positive level numbers)
for level in human_levels:
    tick_positions.append(-level)
    tick_labels.append(str(int(level)))

# Ancestor at 0
tick_positions.append(0)
tick_labels.append('0')

# Mouse side (positive x positions, show positive level numbers)
for level in mouse_levels:
    tick_positions.append(level)
    tick_labels.append(str(int(level)))

ax.set_xticks(tick_positions)
ax.set_xticklabels(tick_labels, fontsize=10)

# Add vertical line at ancestor
ax.axvline(x=0, color='gold', linewidth=2, linestyle='--', alpha=0.5, zorder=0)

# Clean spines
for spine in ax.spines.values():
    spine.set_linewidth(2)
ax.spines['top'].set_visible(False)

plt.tight_layout()
plt.savefig("side_by_side_comparison.png", dpi=300, bbox_inches='tight')
plt.close()
print("✓ Saved: side_by_side_comparison.png")

print("\n" + "="*70)
print("✓✓✓ SIDE-BY-SIDE COMPARISON COMPLETE!")
print("="*70)
print("\nLayout: Human (LEFT) ← Ancestor (CENTER) → Mouse (RIGHT)")
print(f"\n📊 Level Ranges:")
print(f"  Human: Levels 1-19 ({len(human_df)} levels)")
print(f"  Mouse: Levels 1-20 ({len(mouse_df)} levels)")
print(f"\n🏆 Closest levels:")
print(f"  Human: Level {h_closest['level']} (d={h_closest['distance_from_ancestor']:.6f})")
print(f"  Mouse: Level {m_closest['level']} (d={m_closest['distance_from_ancestor']:.6f})")
print("="*70)