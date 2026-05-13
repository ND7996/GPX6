"""
Mouse Homolog Alignment - JCIM Style (10pt Fonts)
Saved to Figures_FINAL
"""

import sys
import os

ACS_PATH = r"./analysis_scripts/Scripts_to_generate_figures/Figures"

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

def visualize_homolog_alignment(alignment_data, positions, title="Homolog 1",
                                annotate_pos=None, annotate_row=None,
                                annotate_label=None):
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
        weight = 'bold' if (annotate_pos is not None and pos == annotate_pos) else 'normal'
        color = '#cc2200' if (annotate_pos is not None and pos == annotate_pos) else 'black'
        ax.text(idx + 0.425, n_rows + 0.3, str(pos), ha='center', va='bottom',
                fontsize=FONT_SIZE, fontfamily=FONT_FAMILY, 
                fontweight=weight, color=color)

    # Annotation for position 49
    if annotate_pos is not None and annotate_row is not None and annotate_label is not None:
        col_idx = positions.index(annotate_pos)
        row_y = n_rows - annotate_row - 1
        bx = col_idx + 0.425

        bracket_top = n_rows + 1.45
        bracket_base = n_rows + 0.55
        ax.annotate('', xy=(bx, bracket_base), xytext=(bx, bracket_top),
                    arrowprops=dict(arrowstyle='->', color='#cc2200', lw=1.0))

        ax.text(bx, bracket_top + 0.08, annotate_label, ha='center', va='bottom',
                fontsize=9, fontfamily=FONT_FAMILY, color='black',
                bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='#cc2200', lw=0.5))

        # Highlight cell
        highlight = mpatches.FancyBboxPatch((col_idx - 0.05, row_y - 0.05),
                                             box_size + 0.10, box_size + 0.10,
                                             boxstyle="round,pad=0.05", facecolor='none',
                                             edgecolor='#cc2200', linewidth=1.5, zorder=5)
        ax.add_patch(highlight)

        # Vertical guide line
        ax.plot([bx, bx], [-0.3, n_rows + 0.55], color='#cc2200', linewidth=0.5,
                linestyle='--', alpha=0.35, zorder=0)

        # Horizontal arrow
        arrow_y = row_y + 0.425
        ax.annotate('', xy=(col_idx + box_size + 0.1, arrow_y),
                    xytext=(col_idx + 2.2, arrow_y),
                    arrowprops=dict(arrowstyle='->', color='#cc2200', lw=0.8))
        ax.text(col_idx + 2.35, arrow_y, f"Introduced at L{annotate_row}",
                va='center', fontsize=9, fontfamily=FONT_FAMILY, color='black')

    ax.set_xlim(-1.5, n_cols + 0.5)
    ax.set_ylim(-0.5, n_rows + 2.2)
    ax.set_aspect('equal')
    ax.axis('off')

    ax.text(n_cols / 2, n_rows + 2.0, title, ha='center', va='bottom',
            fontsize=FONT_SIZE_TITLE, fontfamily=FONT_FAMILY, color='black')

    plt.tight_layout()
    return fig

if __name__ == "__main__":
    positions = [48, 52, 47, 99, 54, 177, 144, 74, 178, 143, 139, 87, 142, 102, 24, 60, 181, 104, 49, 4, 107]
    pos_to_col = {pos: idx for idx, pos in enumerate(positions)}

    yellow_positions = [54, 49, 24, 139, 47, 60, 74, 144, 99, 177,
                        178, 104, 4, 52, 87, 102, 107, 142, 181, 48, 143]

    alignment_data = []
    selected_positions = set()

    for row_idx in range(21):
        row = [0] * len(positions)
        for pos in selected_positions:
            if pos in pos_to_col:
                row[pos_to_col[pos]] = 1
        new_pos = yellow_positions[row_idx]
        if new_pos in pos_to_col:
            row[pos_to_col[new_pos]] = 2
        selected_positions.add(new_pos)
        alignment_data.append(row)

    pos49_introduced_at_row = yellow_positions.index(49)

    fig = visualize_homolog_alignment(
        alignment_data, positions, title="Mouse to Human Mutational Path",
        annotate_pos=49, annotate_row=pos49_introduced_at_row,
        annotate_label="mouse Sec-GPX6\ncatalytic-residue-swap variant",
    )

    # Save to Figures_FINAL
    output_dir = r"./analysis_scripts/Scripts_to_generate_figures/Figures/Figures_FINAL"
    os.makedirs(output_dir, exist_ok=True)
    fig.savefig(os.path.join(output_dir, "mouse_homolog_alignment.png"), dpi=600, bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, "mouse_homolog_alignment.pdf"), bbox_inches='tight')
    print(f"Figure saved to: {output_dir}")
    plt.show()

