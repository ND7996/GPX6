from pathlib import Path
import csv
import math
from pprint import pprint

import matplotlib.pyplot as plt
import numpy as np

# ── Input / Output ────────────────────────────────────────────────────────────
INPUT_FASTA = Path(
    r"D:\PhD_Thesis\analysis\alignment\all_sequences_for_selection.fasta"
)
OUTPUT_PREFIX = Path("results/gpx6_trilateration")

REFERENCE_IDS = ["MOUSE_GPX6_WT", "HUMAN_GPX6_WT", "ANCESTOR_node_25"]


# ── Parsers ───────────────────────────────────────────────────────────────────
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
    scale = max(len(seq_a), len(seq_b))
    if scale == 0:
        return 0.0
    return levenshtein_distance(seq_a, seq_b) / scale


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


def main():
    all_records = parse_fasta(INPUT_FASTA)
    print(f"Loaded {len(all_records)} sequences from {INPUT_FASTA.name}")

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

    # MOUSE=bottom-left, HUMAN=bottom-right, ANCESTOR=top
    # Reference order from REFERENCE_IDS: [MOUSE, HUMAN, ANCESTOR]
    reference_xy = np.array([
        [0.0, 0.0],                      # MOUSE  — bottom-left
        [1.0, 0.0],                      # HUMAN  — bottom-right
        [0.5, math.sqrt(3.0) / 2.0],    # ANCESTOR — top
    ])

    rows = []
    for query_id, query_seq in query_records:
        distances = np.array([
            sequence_distance(query_seq, ref_seq)
            for ref_seq in reference_sequences
        ], dtype=float)
        xy, residuals = trilaterate(reference_xy, distances)
        closest_idx = int(np.argmin(distances))   # closest reference by distance
        rows.append({
            "sequence_id"      : query_id,
            "group"            : get_group(query_id),
            "closest_ref"      : reference_ids[closest_idx],
            "closest_xy"       : reference_xy[closest_idx].tolist(),
            "distance_to_A"    : float(distances[0]),
            "distance_to_B"    : float(distances[1]),
            "distance_to_C"    : float(distances[2]),
            "x"                : float(xy[0]),
            "y"                : float(xy[1]),
            "residual_A"       : float(residuals[0]),
            "residual_B"       : float(residuals[1]),
            "residual_C"       : float(residuals[2]),
            "mean_residual"    : float(residuals.mean()),
        })

    output_csv = Path(f"{OUTPUT_PREFIX}_coordinates.csv")
    write_csv(output_csv,
              [{k:v for k,v in r.items() if k!="closest_xy"} for r in rows],
              fieldnames=[k for k in rows[0].keys() if k!="closest_xy"])
    print(f"\nSaved CSV : {output_csv}")

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect("equal")
    ax.axis("off")

    # Triangle outline
    triangle = np.vstack((reference_xy, reference_xy[0]))
    ax.plot(triangle[:, 0], triangle[:, 1], color="black", linewidth=2, zorder=1)

    # Reference anchor dots
    ax.scatter(reference_xy[:, 0], reference_xy[:, 1],
               s=200, color="#d04f2a", zorder=10, edgecolor="white", linewidth=0.8)

    # Lines from each query point to its closest reference vertex
    # drawn BEFORE the scatter dots so dots sit on top
    for r in rows:
        col = GROUP_COLORS[r["group"]]
        cxy = r["closest_xy"]
        ax.plot([r["x"], cxy[0]], [r["y"], cxy[1]],
                color=col, linewidth=0.5, alpha=0.30, zorder=2)

    # Query scatter — grouped, drawn largest group first so small groups visible
    groups_present = sorted(set(r["group"] for r in rows),
                            key=lambda g: -sum(1 for r in rows if r["group"]==g))
    for group in groups_present:
        grows = [r for r in rows if r["group"] == group]
        ax.scatter(
            [r["x"] for r in grows],
            [r["y"] for r in grows],
            color=GROUP_COLORS[group],
            s=60, edgecolor="white", linewidth=0.5,
            zorder=5, label=f"{group} (n={len(grows)})", alpha=0.90
        )

    # Reference labels outside corners
    short_names = [rid.replace("_chainA","").replace("_chain","") for rid in reference_ids]
    label_cfg = [
        dict(ha="right", va="top",    x=reference_xy[0][0]-0.04, y=reference_xy[0][1]-0.05),
        dict(ha="left",  va="top",    x=reference_xy[1][0]+0.04, y=reference_xy[1][1]-0.05),
        dict(ha="center",va="bottom", x=reference_xy[2][0],       y=reference_xy[2][1]+0.04),
    ]
    for name, kw in zip(short_names, label_cfg):
        ax.text(kw["x"], kw["y"], name,
                ha=kw["ha"], va=kw["va"],
                fontsize=9, color="#222222")

    ax.legend(loc="upper right", framealpha=0.9, fontsize=8,
              title="Sequence group", title_fontsize=8)

    ax.set_xlim(-0.25, 1.25)
    ax.set_ylim(-0.15, 1.05)
    fig.tight_layout()

    output_png = Path(f"{OUTPUT_PREFIX}.png")
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.show()
    print(f"Saved plot: {output_png}")
    pprint(rows[:3])


if __name__ == "__main__":
    main()