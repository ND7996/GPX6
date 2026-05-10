import matplotlib.pyplot as plt
import numpy as np

# Data
mouse_to_human = [
    54, 49, 24, 139, 47, 60, 74, 144, 99, 177,
    178, 104, 4, 52, 87, 102, 107, 142, 181, 48, 143
]

human_to_mouse = [
    87, 173, 47, 143, 60, 104, 142, 139, 181, 48,
    3, 54, 144, 177, 102, 24, 99, 52, 178, 74
]

x = np.arange(len(mouse_to_human))

# Colors
mouse_color = '#FF8C42'
human_color = '#00CED1'
highlight_color = '#FFD400'

fig, ax = plt.subplots(figsize=(13, 3.6))

# ---- MAIN LINES ----
ax.plot(x, mouse_to_human, color=mouse_color, lw=3, zorder=2)
ax.plot(x[:len(human_to_mouse)], human_to_mouse,
        color=human_color, lw=3, zorder=1)

# ---- POINTS ----
for i, val in enumerate(mouse_to_human):
    face = highlight_color if val == 49 else mouse_color
    edge = '#B89B00' if val == 49 else '#D96B1F'
    ax.scatter(x[i], val, s=180, color=face,
               edgecolor=edge, linewidth=2.2, zorder=4)

for i, val in enumerate(human_to_mouse):
    ax.scatter(x[i], val, s=180, color=human_color,
               edgecolor='#00A0A8', linewidth=2.2, zorder=3)

# ---- LABELS (FIXED OVERLAPS) ----
for i, val in enumerate(mouse_to_human):
    offset = 30 if val == 54 else 16 + (i % 3) * 6
    face = highlight_color if val == 49 else '#2D2D2D'
    txt_col = 'black' if val == 49 else 'white'

    ax.text(x[i], val + offset, str(val),
            ha='center', va='bottom',
            fontsize=9, fontweight='bold',
            color=txt_col, zorder=6,
            bbox=dict(boxstyle='round,pad=0.3',
                      facecolor=face,
                      edgecolor='black',
                      linewidth=0.6))

for i, val in enumerate(human_to_mouse):
    ax.text(x[i], val - (18 + (i % 3) * 6), str(val),
            ha='center', va='top',
            fontsize=9, fontweight='bold',
            color='white', zorder=6,
            bbox=dict(boxstyle='round,pad=0.3',
                      facecolor='#2D2D2D',
                      edgecolor='black',
                      linewidth=0.6))

# ---- ARROWS (HYSTERESIS) ----
# Mouse → Human
ax.annotate('', xy=(12.4, mouse_to_human[12]),
            xytext=(11.3, mouse_to_human[11]),
            arrowprops=dict(arrowstyle='->',
                            color=mouse_color,
                            lw=3,
                            mutation_scale=18))

# Human → Mouse
ax.annotate('', xy=(9.6, human_to_mouse[9]),
            xytext=(10.7, human_to_mouse[10]),
            arrowprops=dict(arrowstyle='->',
                            color=human_color,
                            lw=3,
                            mutation_scale=18))

# ---- AXES & SCALES (RESTORED) ----
ax.set_xlim(-1, len(x))
ax.set_ylim(-10, 200)

# Left Y-axis (Mouse)
ax.set_ylabel('Mouse\nResidues', color=mouse_color,
              fontsize=11, fontweight='bold', labelpad=12)
ax.tick_params(axis='y', colors=mouse_color, width=2, length=10)
ax.set_yticks(np.arange(0, 201, 20))

# Right Y-axis (Human)
ax2 = ax.twinx()
ax2.set_ylim(ax.get_ylim())
ax2.set_ylabel('Human\nResidues', color=human_color,
               fontsize=11, fontweight='bold', labelpad=12)
ax2.tick_params(axis='y', colors=human_color, width=2, length=10)
ax2.set_yticks(np.arange(0, 201, 20))

# ---- CLEAN LOOK ----
ax.set_xticks([])
for spine in ax.spines.values():
    spine.set_visible(False)
for spine in ax2.spines.values():
    spine.set_visible(False)

ax.set_facecolor('white')
fig.patch.set_facecolor('white')

plt.tight_layout()
plt.savefig("gpx6_hysteresis_with_scales.png", dpi=300)
plt.show()
