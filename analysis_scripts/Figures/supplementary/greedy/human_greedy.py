"""
Human Homolog Alignment - JCIM Style (10pt Fonts)
Saved to Figures_FINAL
"""

import sys
import os

ACS_PATH = r"D:\PhD_Thesis\analysis\FINAL_PUBLICATION_FIGURES"

if ACS_PATH not in sys.path:
    sys.path.append(ACS_PATH)

from acsfonts import set_jcim_style
set_jcim_style()

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# JCIM font settings (all black)
plt.rcParams.update({
    "font.size": 10,
    "font.family": "Arial",
    "text.color": "black",
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black"
})

FONT_FAMILY = "Arial"
FONT_SIZE = 10
FONT_SIZE_TITLE = 12

def visualize_homolog_alignment(alignment_data, positions, title="Homolog 2"):
    fig, ax = plt.subplots(figsize=(20, 12))
    
    colors = {0: 'white', 1: '#333333', 2: '#FFD700'}  # black and gold
    
    n_rows = len(alignment_data)
    n_cols = len(alignment_data[0])
    box_size = 0.85
    
    # Draw grid cells
    for row_idx, row in enumerate(alignment_data):
        for col_idx, value in enumerate(row):
            color = colors.get(value, 'white')
            rect = mpatches.Rectangle((col_idx, n_rows - row_idx - 1), box_size, box_size,
                                       facecolor=color, edgecolor='black', linewidth=0.5)
            ax.add_patch(rect)

        # Row labels (L0, L1, L2...)
        ax.text(-0.7, n_rows - row_idx - 1 + 0.425, f"L{row_idx}",
                ha='right', va='center', fontsize=FONT_SIZE, 
                fontfamily=FONT_FAMILY, color='black')
    
    # Column labels (position numbers)
    for idx, pos in enumerate(positions):
        ax.text(idx + 0.425, n_rows + 0.3, str(pos), ha='center', va='bottom',
                fontsize=FONT_SIZE, fontfamily=FONT_FAMILY, color='black')
    
    ax.set_xlim(-1.5, n_cols + 0.5)
    ax.set_ylim(-0.5, n_rows + 1.5)
    ax.set_aspect('equal')
    ax.axis('off')
    
    ax.text(n_cols/2, n_rows + 1.2, title, ha='center', va='bottom',
            fontsize=FONT_SIZE_TITLE, fontfamily=FONT_FAMILY, color='black')
    
    plt.tight_layout()
    return fig

if __name__ == "__main__":
    positions = [48, 52, 47, 99, 54, 177, 144, 74, 178, 143, 139, 173, 87, 142, 102, 24, 60, 181, 3, 104]
    pos_to_col = {pos: idx for idx, pos in enumerate(positions)}
    
    yellow_positions = [87, 99, 47, 143, 60, 104, 142, 139, 181, 52,
                        48, 144, 54, 177, 102, 24, 3, 173, 178, 74]
    
    alignment_data = []
    selected_positions = set()
    
    for row_idx in range(20):
        row = [0] * 20
        for pos in selected_positions:
            if pos in pos_to_col:
                row[pos_to_col[pos]] = 1
        new_pos = yellow_positions[row_idx]
        if new_pos in pos_to_col:
            row[pos_to_col[new_pos]] = 2
        selected_positions.add(new_pos)
        alignment_data.append(row)
    
    fig = visualize_homolog_alignment(alignment_data, positions, title="Human to Mouse Mutational Path")
    
    # Save to Figures_FINAL
    output_dir = r"D:\PhD_Thesis\analysis\FINAL_PUBLICATION_FIGURES\Figures_FINAL"
    os.makedirs(output_dir, exist_ok=True)
    fig.savefig(os.path.join(output_dir, "human_homolog_alignment.png"), dpi=600, bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, "human_homolog_alignment.pdf"), bbox_inches='tight')
    print(f"Figure saved to: {output_dir}")
    plt.show()