"""
Distance-based trilateration from Hamming distances â€” GPX6 protein sequences
Adapted for a single FASTA input with references: MOUSE_GPX6_WT, HUMAN_GPX6_WT, ANCESTOR_node_25
"""

from pathlib import Path
import csv
import math
from pprint import pprint

import matplotlib.pyplot as plt
import numpy as np

# â”€â”€ Input / Output â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
INPUT_FASTA = Path(
    r"./analysis_scripts/alignment/gpx6_human_mouse_input.fasta"
)
OUTPUT_PREFIX = Path("results/gpx6_trilateration")

# Substrings that identify the three reference sequences (case-insensitive)
REFERENCE_IDS = ["MOUSE_GPX6_WT", "HUMAN_GPX6_WT", "ANCESTOR_node_25"]


# â”€â”€ Parsers â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
def parse_fasta(path: Path):
    records = []
    header = None
    chunks = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks).upper()))
                header = line[1:].strip() or f"seq_{len(records) + 1}"
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks).upper()))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


# â”€â”€ Distance metric â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
def levenshtein_distance(seq_a: str, seq_b: str) -> int:
    if len(seq_a) < len(seq_b):
        seq_a, seq_b = seq_b, seq_a
    previous = list(range(len(seq_b) + 1))
    for i, char_a in enumerate(seq_a, start=1):
        current = [i]
        for j, char_b in enumerate(seq_b, start=1):
            current.append(min(
                current[j - 1] + 1,
                previous[j] + 1,
                previous[j - 1] + (char_a != char_b)
            ))
        previous = current
    return previous[-1]


def sequence_distance(seq_a: str, seq_b: str) -> float:
    """Normalised Levenshtein â€” handles unequal-length sequences, returns [0, 1]."""
    scale = max(len(seq_a), len(seq_b))
    if scale == 0:
        return 0.0
    return levenshtein_distance(seq_a, seq_b) / scale


# â”€â”€ Geometry â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
def trilaterate(reference_xy, distances):
    p1 = np.asarray(reference_xy[0], dtype=float)
    p2 = np.asarray(reference_xy[1], dtype=float)
    p3 = np.asarray(reference_xy[2], dtype=float)
    d1, d2, d3 = [float(v) for v in distances]
    matrix = 2.0 * np.array([
        [p2[0] - p1[0], p2[1] - p1[1]],
        [p3[0] - p1[0], p3[1] - p1[1]],
    ])
    rhs = np.array([
        d1**2 - d2**2 + np.dot(p2, p2) - np.dot(p1, p1),
        d1**2 - d3**2 + np.dot(p3, p3) - np.dot(p1, p1),
    ])
    solution, _, _, _ = np.linalg.lstsq(matrix, rhs, rcond=None)
    residuals = np.array([
        abs(np.linalg.norm(solution - pt) - d)
        for pt, d in zip(reference_xy, distances)
    ])
    return solution, residuals


def write_csv(path: Path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


# â”€â”€ Grouping for plot colours â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
def get_group(seq_id: str) -> str:
    sid = seq_id.upper()
    if "MOUSE" in sid:
        return "Mouse GPX6 variants"
    elif "HUMAN" in sid:
        return "Human GPX6 variants"
    else:
        return "Other"

GROUP_COLORS = {
    "Mouse GPX6 variants": "#2166ac",
    "Human GPX6 variants": "#d6604d",
    "Other":               "#4dac26",
}


# â”€â”€ Main â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
def main():
    all_records = parse_fasta(INPUT_FASTA)
    print(f"Loaded {len(all_records)} sequences from {INPUT_FASTA.name}")

    # Split references vs queries
    reference_records = []
    query_records = []
    for seq_id, seq in all_records:
        matched = any(ref_id.lower() in seq_id.lower() for ref_id in REFERENCE_IDS)
        if matched:
            reference_records.append((seq_id, seq))
        else:
            query_records.append((seq_id, seq))

    if len(reference_records) != 3:
        found = [r[0] for r in reference_records]
        raise ValueError(
            f"Expected exactly 3 references matching {REFERENCE_IDS}, "
            f"but found {len(reference_records)}: {found}"
        )

    print(f"References : {[r[0] for r in reference_records]}")
    print(f"Queries    : {len(query_records)} sequences")

    reference_ids       = [r[0] for r in reference_records]
    reference_sequences = [r[1] for r in reference_records]

    # Equilateral triangle anchors
    reference_xy = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.5, math.sqrt(3.0) / 2.0],
    ])

    # Compute distances and trilaterate
    rows = []
    for query_id, query_seq in query_records:
        distances = np.array([
            sequence_distance(query_seq, ref_seq)
            for ref_seq in reference_sequences
        ], dtype=float)
        xy, residuals = trilaterate(reference_xy, distances)
        rows.append({
            "sequence_id"   : query_id,
            "distance_to_A" : float(distances[0]),
            "distance_to_B" : float(distances[1]),
            "distance_to_C" : float(distances[2]),
            "x"             : float(xy[0]),
            "y"             : float(xy[1]),
            "residual_A"    : float(residuals[0]),
            "residual_B"    : float(residuals[1]),
            "residual_C"    : float(residuals[2]),
            "mean_residual" : float(residuals.mean()),
        })

    # Export CSV
    output_csv = Path(f"{OUTPUT_PREFIX}_coordinates.csv")
    write_csv(output_csv, rows, fieldnames=list(rows[0].keys()))
    print(f"\nSaved CSV : {output_csv}")

    # â”€â”€ Plot â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    fig, ax = plt.subplots(figsize=(8, 8))

    # Triangle outline
    triangle = np.vstack((reference_xy, reference_xy[0]))
    ax.plot(triangle[:, 0], triangle[:, 1], color="black", linewidth=1.5, zorder=1)

    # Reference anchor markers
    ax.scatter(reference_xy[:, 0], reference_xy[:, 1],
               s=220, color="#d04f2a", zorder=5, marker="^")

    # Reference labels â€” short names, positioned outside triangle corners
    short_names = [rid.replace("_chainA", "").replace("_chain", "") for rid in reference_ids]
    anchor_label_kwargs = [
        dict(ha="left",   va="top",    x=reference_xy[0][0] - 0.03, y=reference_xy[0][1] - 0.05),
        dict(ha="right",  va="top",    x=reference_xy[1][0] + 0.03, y=reference_xy[1][1] - 0.05),
        dict(ha="center", va="bottom", x=reference_xy[2][0],         y=reference_xy[2][1] + 0.04),
    ]
    for name, kw in zip(short_names, anchor_label_kwargs):
        ax.text(kw["x"], kw["y"], name,
                ha=kw["ha"], va=kw["va"],
                fontsize=9, fontweight="bold", color="#d04f2a")

    # Query points â€” colored by species group, NO per-point text labels
    groups_present = sorted(set(get_group(r["sequence_id"]) for r in rows))
    for group in groups_present:
        grows = [r for r in rows if get_group(r["sequence_id"]) == group]
        ax.scatter(
            [r["x"] for r in grows],
            [r["y"] for r in grows],
            color=GROUP_COLORS[group],
            s=90, edgecolor="white", linewidth=0.6,
            zorder=3, label=group, alpha=0.9
        )

    ax.legend(loc="upper right", framealpha=0.9, fontsize=9,
              title="Sequence group", title_fontsize=9)

    # # ax.set_title(
    # #     "GPX6 â€” Distance-based trilateration\n"
    # #     "(normalised Levenshtein distance, protein sequences)",
    # #     fontsize=12
    # )
    ax.set_aspect("equal")
    ax.set_xlim(-0.3, 1.3)
    ax.set_ylim(-0.15, 1.05)
    ax.axis("off")   # hide ticks â€” triangle is the coordinate system
    fig.tight_layout()

    output_png = Path(f"{OUTPUT_PREFIX}.png")
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.show()
    print(f"Saved plot: {output_png}")
    pprint(rows)


if __name__ == "__main__":
    main()
