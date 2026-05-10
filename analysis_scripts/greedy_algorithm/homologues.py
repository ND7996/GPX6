import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

# =================== ACS PUBLICATION STYLE ===================
plt.rcParams.update({
    'font.family':       'sans-serif',
    'font.sans-serif':   ['Liberation Sans', 'Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size':          7,
    'axes.labelsize':     8,
    'xtick.labelsize':    7,
    'ytick.labelsize':    7,
    'pdf.fonttype':       42,
    'ps.fonttype':        42,
})

OUTPUT_DIR = r"D:\PhD_Thesis\GPX6\analysis\greedy_algorithm"
os.makedirs(OUTPUT_DIR, exist_ok=True)


def visualize_homolog_alignment(alignment_data, positions, title="Homolog"):
    n_rows = len(alignment_data)
    n_cols = len(alignment_data[0])

    cell  = 0.30
    fig_w = max(5.0, n_cols * cell + 1.2)
    fig_h = max(3.0, n_rows * cell + 1.0)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    colors = {0: 'white', 1: '#222222', 2: '#FFD700'}

    box_size = 0.82
    for row_idx, row in enumerate(alignment_data):
        for col_idx, value in enumerate(row):
            color = colors.get(value, 'white')
            rect = mpatches.Rectangle(
                (col_idx, n_rows - row_idx - 1),
                box_size, box_size,
                facecolor=color,
                edgecolor='#555555',
                linewidth=0.5
            )
            ax.add_patch(rect)

    # Position labels — rotated 90° for readability
    for idx, pos in enumerate(positions):
        ax.text(idx + 0.41, n_rows + 0.15, str(pos),
                ha='center', va='bottom',
                fontsize=6, fontweight='normal',
                rotation=90)

    # Row (step) labels on the left
    for row_idx in range(n_rows):
        ax.text(-0.55, n_rows - row_idx - 0.59,
                str(row_idx + 1),
                ha='right', va='center', fontsize=6)

    # Axis labels
    ax.text(-1.4, n_rows / 2, 'Mutation step',
            ha='center', va='center', fontsize=7, rotation=90)
    ax.text(n_cols / 2, n_rows + 1.4, 'Residue position',
            ha='center', va='bottom', fontsize=7)

    # Title
    ax.text(n_cols / 2, n_rows + 2.4, title,
            ha='center', va='bottom', fontsize=8, fontweight='bold')

    # Legend (outside right)
    legend_x = n_cols + 0.8
    legend_y = n_rows - 1
    for color, label in [('#222222', 'Previously mutated'),
                          ('#FFD700', 'Current mutation'),
                          ('white',   'Not mutated')]:
        rect = mpatches.Rectangle((legend_x, legend_y), 0.7, 0.7,
                                   facecolor=color, edgecolor='#555555', linewidth=0.5)
        ax.add_patch(rect)
        ax.text(legend_x + 0.9, legend_y + 0.35, label,
                ha='left', va='center', fontsize=6)
        legend_y -= 1.1

    ax.set_xlim(-1.8, n_cols + 6.0)
    ax.set_ylim(-0.8, n_rows + 3.0)
    ax.set_aspect('equal')
    ax.axis('off')

    plt.tight_layout(pad=0.3)
    return fig


def build_alignment(positions, yellow_positions):
    pos_to_col = {pos: idx for idx, pos in enumerate(positions)}
    alignment_data = []
    selected = set()
    for row_idx in range(len(yellow_positions)):
        row = [0] * len(positions)
        for pos in selected:
            if pos in pos_to_col:
                row[pos_to_col[pos]] = 1
        new_pos = yellow_positions[row_idx]
        if new_pos in pos_to_col:
            row[pos_to_col[new_pos]] = 2
        selected.add(new_pos)
        alignment_data.append(row)
    return alignment_data


def save_fig(fig, name):
    for ext in ('png', 'pdf', 'svg'):
        path = os.path.join(OUTPUT_DIR, f"{name}.{ext}")
        kw = dict(bbox_inches='tight', facecolor='white')
        if ext == 'png':
            kw['dpi'] = 600
        fig.savefig(path, **kw)
    print(f"Saved: {name}")


# =================== MOUSE → HUMAN ===================
positions_h1 = [48,52,47,99,54,177,144,74,178,143,139,87,142,102,24,60,181,104,49,4,107]
yellow_h1    = [54,49,24,139,47,60,74,144,99,177,178,104,4,52,87,102,107,142,181,48,143]

fig1 = visualize_homolog_alignment(
    build_alignment(positions_h1, yellow_h1),
    positions_h1,
    #title="Mouse → Human GPX6 mutation pathway"
)
save_fig(fig1, "mouse_to_human_mutants")
plt.show()


# =================== HUMAN → MOUSE ===================
positions_h2 = [48,52,47,99,54,177,144,74,178,143,139,173,87,142,102,24,60,181,3,104]
yellow_h2    = [87,173,47,143,60,104,142,139,181,48,3,54,144,177,102,24,99,52,178,74]

fig2 = visualize_homolog_alignment(
    build_alignment(positions_h2, yellow_h2),
    positions_h2,
    #title="Human → Mouse GPX6 mutation pathway"
)
save_fig(fig2, "human_to_mouse_mutants")
plt.show()