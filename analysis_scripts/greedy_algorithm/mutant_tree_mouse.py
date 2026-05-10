import sys
import os

ACS_PATH = r"D:\PhD_Thesis\analysis\FINAL_PUBLICATION_FIGURES"

if ACS_PATH not in sys.path:
    sys.path.append(ACS_PATH)

from acsfonts import acs_figure, save_acs_figure

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


def visualize_homolog_alignment(alignment_data, positions, title="Homolog 1",
                                annotate_pos=None, annotate_row=None,
                                annotate_label=None):
    fig, ax = plt.subplots(figsize=(20, 12))

    colors = {
        0: 'white',
        1: 'black',
        2: 'yellow'
    }

    n_rows = len(alignment_data)
    n_cols = len(alignment_data[0])

    box_size = 0.85
    for row_idx, row in enumerate(alignment_data):
        for col_idx, value in enumerate(row):
            color = colors.get(value, 'white')
            rect = mpatches.Rectangle(
                (col_idx, n_rows - row_idx - 1),
                box_size, box_size,
                facecolor=color,
                edgecolor='black',
                linewidth=0.8
            )
            ax.add_patch(rect)

        # Add level label on the left of each row
        ax.text(-0.7, n_rows - row_idx - 1 + 0.425, f"L{row_idx}",
                ha='right', va='center', fontsize=8, fontweight='normal')

    # Column position labels at top
    for idx, pos in enumerate(positions):
        weight = 'bold' if (annotate_pos is not None and pos == annotate_pos) else 'normal'
        color  = '#cc2200' if (annotate_pos is not None and pos == annotate_pos) else 'black'
        ax.text(idx + 0.425, n_rows + 0.3, str(pos),
                ha='center', va='bottom', fontsize=8,
                fontweight=weight, color=color)

    # ── Annotation for position 49 ──────────────────────────────────────────
    if annotate_pos is not None and annotate_row is not None and annotate_label is not None:
        col_idx = positions.index(annotate_pos)
        row_y   = n_rows - annotate_row - 1
        bx      = col_idx + 0.425

        # Vertical arrow above the column header
        bracket_top  = n_rows + 1.45
        bracket_base = n_rows + 0.55
        ax.annotate('', xy=(bx, bracket_base), xytext=(bx, bracket_top),
                    arrowprops=dict(arrowstyle='->', color='#cc2200', lw=1.4))

        # Label above the arrow
        ax.text(bx, bracket_top + 0.08, annotate_label,
                ha='center', va='bottom', fontsize=7.5,
                color='#cc2200', fontstyle='italic',
                bbox=dict(boxstyle='round,pad=0.25', fc='#fff4f2',
                          ec='#cc2200', lw=0.8))

        # Red highlight border around the cell where pos 49 first appears
        highlight = mpatches.FancyBboxPatch(
            (col_idx - 0.05, row_y - 0.05),
            box_size + 0.10, box_size + 0.10,
            boxstyle="round,pad=0.05",
            facecolor='none', edgecolor='#cc2200', linewidth=2.0,
            zorder=5
        )
        ax.add_patch(highlight)

        # Dashed vertical guide line through the whole column
        ax.plot([bx, bx], [-0.3, n_rows + 0.55],
                color='#cc2200', linewidth=0.7,
                linestyle='--', alpha=0.35, zorder=0)

        # Horizontal arrow + label pointing to the highlighted row
        arrow_y = row_y + 0.425
        ax.annotate('',
                    xy=(col_idx + box_size + 0.1, arrow_y),
                    xytext=(col_idx + 2.2, arrow_y),
                    arrowprops=dict(arrowstyle='->', color='#cc2200', lw=1.2))
        ax.text(col_idx + 2.35, arrow_y,
                f"Introduced at L{annotate_row}",
                va='center', fontsize=7.5, color='#cc2200', fontstyle='italic')

    ax.set_xlim(-1.5, n_cols + 0.5)
    ax.set_ylim(-0.5, n_rows + 2.2)
    ax.set_aspect('equal')
    ax.axis('off')

    ax.text(n_cols / 2, n_rows + 2.0, title,
            ha='center', va='bottom', fontsize=14, fontweight='bold')

    plt.tight_layout()
    return fig


if __name__ == "__main__":
    positions = [48, 52, 47, 99, 54, 177, 144, 74, 178, 143, 139, 87, 142, 102, 24, 60, 181, 104, 49, 4, 107]
    pos_to_col = {pos: idx for idx, pos in enumerate(positions)}

    yellow_positions = [
        54, 49, 24, 139, 47, 60, 74, 144, 99, 177,
        178, 104, 4, 52, 87, 102, 107, 142, 181, 48, 143,
    ]

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

    pos49_introduced_at_row = yellow_positions.index(49)  # == 1

    fig = visualize_homolog_alignment(
        alignment_data,
        positions,
        title="Homolog 1",
        annotate_pos=49,
        annotate_row=pos49_introduced_at_row,
        annotate_label="mouse Sec-GPX6\ncatalytic-residue-swap variant of mouse",
    )

    output_dir = r"D:\PhD_Thesis\analysis\FINAL_PUBLICATION_FIGURES\mutant_tree_mouse"
    os.makedirs(output_dir, exist_ok=True)
    fig.savefig(os.path.join(output_dir, "mouse_to_human_mutants.png"), dpi=300, bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, "mouse_to_human_mutants.pdf"), bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, "mouse_to_human_mutants.svg"), bbox_inches='tight')
    print(f"Figure saved to: {output_dir}")
    plt.show()