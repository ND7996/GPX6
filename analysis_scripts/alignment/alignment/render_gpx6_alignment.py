#!/usr/bin/env python3
"""Render the human/mouse GPX6 Clustal alignment used in the SI."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle


ROOT = Path(__file__).resolve().parents[1]
ALIGNMENT_FASTA = ROOT / "alignment" / "gpx6_human_mouse_clustalo.fasta"
OUTPUTS = [
    ROOT / "alignment" / "gpx6_human_mouse_clustalo.png",
    ROOT / "Figures" / "FigS11.png",
]

# EVB mutation files use the mature/PDB numbering that starts at Pro25 of the
# UniProt sequence. The Clustal figure uses full UniProt coordinates.
PDB_TO_UNIPROT_OFFSET = 24

# Residues fixed along the calculated paths, in mature/PDB numbering.
MOUSE_TO_HUMAN_PDB = {
    4, 24, 47, 48, 52, 54, 60, 74, 87, 99,
    102, 104, 107, 139, 142, 143, 144, 177, 178, 181,
}
HUMAN_TO_MOUSE_PDB = {
    3, 24, 47, 48, 52, 54, 60, 87, 99, 102,
    104, 139, 142, 143, 144, 173, 177, 178, 181,
}
CATALYTIC_PDB = 49


AA_COLORS = {
    **{aa: "#4b4b4b" for aa in "AVLIMFWYPG"},
    **{aa: "#2b8c5a" for aa in "STNQCU"},
    **{aa: "#2b61c2" for aa in "KRH"},
    **{aa: "#d54032" for aa in "DE"},
}
IDENTICAL_FILL = "#e9e9e9"
SUBSTITUTION_FILL = "#fff0a8"
HUMAN_MARK = "#1f77b4"
MOUSE_MARK = "#2ca25f"
CATALYTIC_MARK = "#b00020"


def read_fasta(path):
    records = []
    header = None
    seq = []
    for line in path.read_text().splitlines():
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq)))
            header = line[1:]
            seq = []
        else:
            seq.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def label_from_header(header):
    if "GPX6_HUMAN" in header:
        return "Human GPX6 (P59796)"
    if "GPX6_MOUSE" in header:
        return "Mouse GPX6 (Q91WR8)"
    return header.split()[0]


def to_uniprot_positions(pdb_positions):
    return {position + PDB_TO_UNIPROT_OFFSET for position in pdb_positions}


def draw_alignment(records):
    if len(records) != 2:
        raise ValueError("Expected exactly two aligned sequences")

    labels = [label_from_header(header) for header, _ in records]
    human, mouse = [seq.replace("-", "") for _, seq in records]
    if len(human) != len(mouse):
        raise ValueError("Human and mouse sequences must have the same length")

    human_marks = to_uniprot_positions(HUMAN_TO_MOUSE_PDB)
    mouse_marks = to_uniprot_positions(MOUSE_TO_HUMAN_PDB)
    catalytic_position = CATALYTIC_PDB + PDB_TO_UNIPROT_OFFSET

    block_size = 60
    cell = 1.0
    row_step = 1.05
    block_step = 3.55
    x0 = 16.2
    top_pad = 0.7
    legend_y = 14.6

    fig, ax = plt.subplots(figsize=(34.5, 16.2), dpi=100)
    ax.set_xlim(0, x0 + block_size + 8.2)
    ax.set_ylim(0, 16.2)
    ax.invert_yaxis()
    ax.axis("off")

    n = len(human)
    for block_start in range(0, n, block_size):
        block_idx = block_start // block_size
        block_end = min(block_start + block_size, n)
        y_h = top_pad + block_idx * block_step + 0.9
        y_m = y_h + row_step

        ax.text(0.3, y_h + 0.54, labels[0], ha="left", va="center", fontsize=26)
        ax.text(0.3, y_m + 0.54, labels[1], ha="left", va="center", fontsize=26)

        tick_positions = [block_start + 1]
        tick_positions.extend(
            p for p in range(((block_start + 10) // 10) * 10, block_end + 1, 10)
            if p >= block_start + 1
        )
        if block_end == n and n not in tick_positions:
            tick_positions.append(n)
        tick_positions = sorted(set(tick_positions))

        for pos in tick_positions:
            col = pos - block_start - 1
            x = x0 + col + 0.5
            ax.text(x, y_h - 0.45, str(pos), ha="center", va="bottom", fontsize=20, color="#555555")
            ax.plot([x, x], [y_h - 0.35, y_h - 0.18], color="#999999", lw=1.4)

        for col, pos in enumerate(range(block_start + 1, block_end + 1)):
            x = x0 + col
            h_aa = human[pos - 1]
            m_aa = mouse[pos - 1]
            fill = IDENTICAL_FILL if h_aa == m_aa else SUBSTITUTION_FILL

            for y, aa, mark_positions, mark_color in (
                (y_h, h_aa, human_marks, HUMAN_MARK),
                (y_m, m_aa, mouse_marks, MOUSE_MARK),
            ):
                ax.add_patch(Rectangle((x, y), cell * 0.92, cell * 0.82, facecolor=fill, edgecolor="white", lw=1.3))
                ax.text(
                    x + 0.46,
                    y + 0.43,
                    aa,
                    ha="center",
                    va="center",
                    fontsize=21,
                    fontweight="bold",
                    color=AA_COLORS.get(aa, "#4b4b4b"),
                    family="DejaVu Sans Mono",
                )
                if pos in mark_positions:
                    ax.add_patch(
                        Rectangle(
                            (x + 0.04, y + 0.04),
                            cell * 0.84,
                            cell * 0.74,
                            facecolor="none",
                            edgecolor=mark_color,
                            lw=2.4,
                        )
                    )

            if pos == catalytic_position:
                ax.add_patch(
                    Rectangle(
                        (x - 0.02, y_h - 0.08),
                        cell * 0.98,
                        row_step + cell * 0.88,
                        facecolor="none",
                        edgecolor=CATALYTIC_MARK,
                        lw=3.0,
                    )
                )
                ax.text(
                    x + 0.46,
                    y_h - 0.65,
                    str(CATALYTIC_PDB),
                    ha="center",
                    va="center",
                    fontsize=20,
                    fontweight="bold",
                    color=CATALYTIC_MARK,
                )

        ax.text(x0 + block_end - block_start + 1.2, y_h + 0.54, str(block_end), ha="left", va="center", fontsize=20, color="#555555")
        ax.text(x0 + block_end - block_start + 1.2, y_m + 0.54, str(block_end), ha="left", va="center", fontsize=20, color="#555555")

    handles = [
        Patch(facecolor=IDENTICAL_FILL, edgecolor="none", label="Identical"),
        Patch(facecolor=SUBSTITUTION_FILL, edgecolor="none", label="Substitution"),
        Patch(facecolor="white", edgecolor=CATALYTIC_MARK, linewidth=2.5, label="Catalytic site 49"),
        Patch(facecolor="white", edgecolor=HUMAN_MARK, linewidth=2.5, label="H->M fixed residues (19)"),
        Patch(facecolor="white", edgecolor=MOUSE_MARK, linewidth=2.5, label="M->H fixed residues (20)"),
    ]
    ax.legend(
        handles=handles,
        loc="lower left",
        bbox_to_anchor=(0.0, -0.02),
        ncol=5,
        frameon=False,
        fontsize=20,
        handlelength=2.2,
        columnspacing=2.2,
    )

    fig.tight_layout(pad=0.2)
    return fig


def main():
    records = read_fasta(ALIGNMENT_FASTA)
    fig = draw_alignment(records)
    for output in OUTPUTS:
        fig.savefig(output, dpi=100, bbox_inches="tight", pad_inches=0.15)
        print(f"Wrote {output}")
    plt.close(fig)


if __name__ == "__main__":
    main()
