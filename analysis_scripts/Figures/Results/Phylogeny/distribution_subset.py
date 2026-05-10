import sys
import os

ACS_PATH = r"D:\PhD_Thesis\analysis\FINAL_PUBLICATION_FIGURES"

if ACS_PATH not in sys.path:
    sys.path.append(ACS_PATH)

from acsfonts import acs_figure, save_acs_figure

from pathlib import Path
import csv
import math

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

# ── Paths ──────────────────────────────────────────────────────────────────────
# All three FASTA files are loaded and merged
FASTA_MOUSE    = Path(r"D:\PhD_Thesis\GPX6\prep_structures\MOUSE\mutant_pdbs\sequences_mouse.fasta")
FASTA_HUMAN    = Path(r"D:\PhD_Thesis\GPX6\prep_structures\HUMAN\mutant_pdbs\sequences_human.fasta")
FASTA_ANCESTOR = Path(r"D:\PhD_Thesis\analysis\alignment\ancestor_25.fasta")

OUTPUT_PREFIX = Path("results/gpx6")

# ── Reference IDs: substrings matched against FASTA headers (case-insensitive) ─
# Order: A = Mouse (BL), B = Human (BR), C = Ancestor (top)
REFERENCE_IDS = ["original_mousecys", "original_humansec", "ANCESTOR_node_25"]

# Triangle vertices: A=Mouse BL, B=Human BR, C=Ancestor top
VERTICES = np.array([
    [0.0, 0.0],                    # A – Mouse
    [1.0, 0.0],                    # B – Human
    [0.5, math.sqrt(3) / 2],       # C – Ancestor
])
REF_LABELS  = ["Mouse\nGPX6 WT", "Human\nGPX6 WT", "Ancestor\n(node_25)"]
REF_COLORS  = ["#2166ac", "#d6604d", "#4a9c6d"]

_LABEL_OFFSETS = [
    dict(dx=-0.10, dy=-0.10, ha="right",  va="top"),    # Mouse – bottom left
    dict(dx= 0.10, dy=-0.10, ha="left",   va="top"),    # Human – bottom right
    dict(dx= 0.00, dy= 0.08, ha="center", va="bottom"), # Ancestor – top
]

SOFTMAX_TEMP = 0.03
JITTER_SCALE = 0.006
RNG          = np.random.default_rng(seed=42)

# ── FASTA parser ───────────────────────────────────────────────────────────────
def parse_fasta(path: Path):
    records, header, chunks = [], None, []
    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks).upper()))
                header = line[1:].strip() or f"seq_{len(records)+1}"
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks).upper()))
    if not records:
        raise ValueError(f"No FASTA records in {path}")
    return records

def load_all_records():
    """Load and merge sequences from all three FASTA files, deduplicating by ID."""
    seen = {}
    for fasta_path in [FASTA_MOUSE, FASTA_HUMAN, FASTA_ANCESTOR]:
        records = parse_fasta(fasta_path)
        print(f"  Loaded {len(records):>3} sequences from {fasta_path.name}")
        for sid, seq in records:
            if sid not in seen:
                seen[sid] = seq
            else:
                print(f"    [WARN] Duplicate ID skipped: '{sid}'")
    return list(seen.items())

# ── Core edit distance ─────────────────────────────────────────────────────────
def levenshtein_distance(a: str, b: str) -> int:
    if len(a) < len(b):
        a, b = b, a
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        curr = [i]
        for j, cb in enumerate(b, 1):
            curr.append(min(curr[j-1]+1, prev[j]+1, prev[j-1]+(ca != cb)))
        prev = curr
    return prev[-1]

# ── Distance metrics ───────────────────────────────────────────────────────────
def dist_p(a: str, b: str, scale: float = None) -> float:
    if len(a) == len(b):
        cols = [(ca, cb) for ca, cb in zip(a, b) if ca != '-' and cb != '-']
    else:
        cols = list(zip(a, b))
    return sum(ca != cb for ca, cb in cols) / len(cols) if cols else 0.0

def dist_jc(a: str, b: str, scale: float = None) -> float:
    p   = min(dist_p(a, b), 0.94)
    arg = 1.0 - (20.0 * p / 19.0)
    return -(19.0 / 20.0) * math.log(arg) if arg > 0 else 3.0

def dist_lev(a: str, b: str, scale: float = None) -> float:
    denom = scale if scale is not None else max(len(a), len(b))
    return 0.0 if denom == 0 else levenshtein_distance(a, b) / denom

def dist_sqrt_lev(a: str, b: str, scale: float = None) -> float:
    return math.sqrt(dist_lev(a, b, scale))

METRICS = [
    ("p-distance",          dist_p),
    ("Jukes-Cantor (JC20)", dist_jc),
    ("Norm. Levenshtein",   dist_lev),
    ("Sqrt Levenshtein",    dist_sqrt_lev),
]

# ── Softmax barycentric ────────────────────────────────────────────────────────
def distances_to_xy(d_A, d_B, d_C, temp=SOFTMAX_TEMP):
    ds  = np.array([d_A, d_B, d_C], dtype=float)
    ds -= ds.min()
    w   = np.exp(-ds / temp)
    w  /= w.sum()
    xy  = w[0]*VERTICES[0] + w[1]*VERTICES[1] + w[2]*VERTICES[2]
    return float(xy[0]), float(xy[1]), w[0], w[1], w[2], int(np.argmin([d_A, d_B, d_C]))

# ── Grouping ───────────────────────────────────────────────────────────────────
# Mouse variants: GPX6_level01_ … GPX6_level20_  (prefixed with "GPX6_")
# Human variants: level01_ … level19_            (no GPX6_ prefix)
GROUP_STYLE = {
    "Mouse GPX6 variants": dict(color="#2166ac", marker="o", zorder=4),
    "Human GPX6 variants": dict(color="#d6604d", marker="o", zorder=5),
    "Other":               dict(color="#4dac26", marker="o", zorder=3),
}

def get_group(seq_id: str) -> str:
    sid = seq_id.lower()
    if "mouse" in sid or sid.startswith("gpx6_level"):
        return "Mouse GPX6 variants"
    if "human" in sid or (sid.startswith("level") and not sid.startswith("gpx6_level")):
        return "Human GPX6 variants"
    return "Other"

# ── Plot helpers ───────────────────────────────────────────────────────────────
def draw_triangle(ax):
    tri = np.vstack([VERTICES, VERTICES[0]])
    ax.plot(tri[:, 0], tri[:, 1], color="#333333", linewidth=1.8, zorder=2)
    for (vx, vy), col, label, off in zip(VERTICES, REF_COLORS, REF_LABELS, _LABEL_OFFSETS):
        ax.scatter(vx, vy, s=300, color=col, zorder=8,
                   edgecolors="white", linewidths=2.0)
        ax.text(vx + off["dx"], vy + off["dy"], label,
                ha=off["ha"], va=off["va"],
                fontsize=9, color=col,
                multialignment="center", zorder=9)

def scatter_groups(ax, rows):
    for group, style in GROUP_STYLE.items():
        pts = [r for r in rows if r["group"] == group]
        if not pts:
            continue
        xs = np.array([r["x"] for r in pts]) + RNG.uniform(-JITTER_SCALE, JITTER_SCALE, len(pts))
        ys = np.array([r["y"] for r in pts]) + RNG.uniform(-JITTER_SCALE, JITTER_SCALE, len(pts))
        ax.scatter(xs, ys, color=style["color"], marker=style["marker"],
                   s=80, edgecolor="white", linewidth=0.7,
                   zorder=style["zorder"], alpha=0.92, label=group)

def finalise_panel(ax, title):
    ax.set_aspect("equal")
    ax.set_xlim(-0.30, 1.30)
    ax.set_ylim(-0.25, 1.15)
    ax.axis("off")

def write_csv(path: Path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    print("Loading FASTA files...")
    all_records = load_all_records()
    print(f"Total unique sequences: {len(all_records)}\n")

    # ── Debug: show all headers ────────────────────────────────────────────────
    print("All sequence IDs:")
    for sid, seq in all_records:
        print(f"  '{sid}'  (len={len(seq)})")
    print()

    # ── Split into references and queries ─────────────────────────────────────
    ref_records, query_records = [], []
    for sid, seq in all_records:
        if any(r.lower() in sid.lower() for r in REFERENCE_IDS):
            ref_records.append((sid, seq))
        else:
            query_records.append((sid, seq))

    if len(ref_records) != 3:
        found   = [r[0] for r in ref_records]
        missing = [r for r in REFERENCE_IDS
                   if not any(r.lower() in sid.lower() for sid, _ in ref_records)]
        raise ValueError(
            f"Expected 3 references, found {len(ref_records)}.\n"
            f"  Found:   {found}\n"
            f"  Missing substrings: {missing}\n"
            f"  Update REFERENCE_IDS to match the headers printed above."
        )

    # Order references to match REFERENCE_IDS (A=Mouse, B=Human, C=Ancestor)
    ordered_refs = []
    for target in REFERENCE_IDS:
        match = next((r for r in ref_records if target.lower() in r[0].lower()), None)
        if match is None:
            raise ValueError(f"Reference '{target}' not found.")
        ordered_refs.append(match)

    ref_ids  = [r[0] for r in ordered_refs]
    ref_seqs = [r[1] for r in ordered_refs]

    max_ref_len = max(len(s) for s in ref_seqs)
    print(f"  A=Mouse:    {ref_ids[0]}  (len={len(ref_seqs[0])})")
    print(f"  B=Human:    {ref_ids[1]}  (len={len(ref_seqs[1])})")
    print(f"  C=Ancestor: {ref_ids[2]}  (len={len(ref_seqs[2])})")
    print(f"  Lev scale (max ref len): {max_ref_len}")
    print(f"  Queries: {len(query_records)} | T={SOFTMAX_TEMP}\n")

    FIG10_NAMES = ["Fig10_1a", "Fig10_1b", "Fig10_1c", "Fig10_1d"]
    all_csv = []

    for idx, (metric_label, metric_fn) in enumerate(METRICS):
        rows = []
        nearest_counts = {0: 0, 1: 0, 2: 0}

        for qid, qseq in query_records:
            shared_scale = max(max_ref_len, len(qseq))
            d_A = metric_fn(qseq, ref_seqs[0], shared_scale)
            d_B = metric_fn(qseq, ref_seqs[1], shared_scale)
            d_C = metric_fn(qseq, ref_seqs[2], shared_scale)
            x, y, w_A, w_B, w_C, nearest = distances_to_xy(d_A, d_B, d_C)
            nearest_counts[nearest] += 1
            rows.append({
                "sequence_id": qid, "metric": metric_label,
                "group":   get_group(qid),
                "nearest": REF_LABELS[nearest].replace("\n", " "),
                "d_Mouse": d_A, "d_Human": d_B, "d_Ancestor": d_C,
                "w_Mouse": w_A, "w_Human": w_B, "w_Ancestor": w_C,
                "x": x, "y": y,
            })

        all_csv.extend(rows)

        print(f"[{metric_label}] nearest counts → "
              f"Mouse:{nearest_counts[0]}  Human:{nearest_counts[1]}  Ancestor:{nearest_counts[2]}")

        # Save each metric as Fig10a–d
        fig_s, ax_s = plt.subplots(figsize=(7, 7))
        draw_triangle(ax_s)
        scatter_groups(ax_s, rows)
        finalise_panel(ax_s, metric_label)
        ax_s.legend(loc="upper right", framealpha=0.9, fontsize=8,
                    title="Sequence group", title_fontsize=8)
        fig_s.tight_layout()

        fig_name = FIG10_NAMES[idx] + ".png"
        out_s = Path(OUTPUT_PREFIX.parent / fig_name)
        out_s.parent.mkdir(parents=True, exist_ok=True)
        fig_s.savefig(out_s, dpi=300, bbox_inches="tight")
        plt.close(fig_s)
        print(f"  Saved: {out_s}\n")

    # Save CSV
    write_csv(Path(f"{OUTPUT_PREFIX}_coordinates.csv"), all_csv,
              ["sequence_id", "metric", "group", "nearest",
               "d_Mouse", "d_Human", "d_Ancestor",
               "w_Mouse", "w_Human", "w_Ancestor", "x", "y"])
    print(f"Saved CSV: {OUTPUT_PREFIX}_coordinates.csv")


if __name__ == "__main__":
    main()