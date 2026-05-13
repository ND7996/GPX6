import sys
import os

from pathlib import Path
import csv
import math

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import numpy as np
from adjustText import adjust_text

# â”€â”€ JCIM figure fonts (Arial 10 pt, black) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
mpl.rcParams.update({
    "font.family":        "sans-serif",
    "font.sans-serif":    ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":          10,
    "axes.labelsize":     10,
    "axes.titlesize":     10,
    "xtick.labelsize":    10,
    "ytick.labelsize":    10,
    "legend.fontsize":    10,
    "legend.title_fontsize": 10,
    "text.color":         "black",
    "axes.labelcolor":    "black",
    "xtick.color":        "black",
    "ytick.color":        "black",
    "pdf.fonttype":       42,   # embeds fonts as TrueType in PDF/EPS
    "ps.fonttype":        42,
})

# â”€â”€ Paths â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
FASTA_MOUSE    = Path(r"./prep_structures\MOUSE\mutant_pdbs\sequences_mouse.fasta")
FASTA_HUMAN    = Path(r"./prep_structures\HUMAN\mutant_pdbs\sequences_human.fasta")
FASTA_ANCESTOR = Path(r"./analysis_scripts/alignment/gpx6_human_mouse_input.fasta")

OUTPUT_PREFIX = Path("results/gpx6")

REFERENCE_IDS = ["original_mousecys", "original_humansec", "ANCESTOR_node_25"]

VERTICES = np.array([
    [0.0, 0.0],
    [1.0, 0.0],
    [0.5, math.sqrt(3) / 2],
])
REF_LABELS = ["Mouse\nGPX6 WT", "Human\nGPX6 WT", "Ancestor\n(node_25)"]
REF_COLORS = ["#2166ac", "#E87722", "#4a9c6d"]   # blue, orange, green

_LABEL_OFFSETS = [
    dict(dx=-0.10, dy=-0.10, ha="right",  va="top"),
    dict(dx= 0.10, dy=-0.10, ha="left",   va="top"),
    dict(dx= 0.00, dy= 0.08, ha="center", va="bottom"),
]

SOFTMAX_TEMP = 0.20   # controls spread along Mouseâ€“Human axis
ARC_HEIGHT   = 0.30   # how high dots bow toward Ancestor (0=flat, ~0.15=pronounced arc)
RNG          = np.random.default_rng(seed=42)

# â”€â”€ Mutation table â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
MUTATION_TABLE = {
    3:   ("K", "N"),  4:   ("S", "R"),  16:  ("V", "I"),  22:  ("N", "L"),
    24:  ("I", "L"),  25:  ("D", "N"),  27:  ("G", "E"),  29:  ("F", "Y"),
    30:  ("V", "I"),  31:  ("N", "Q"),  33:  ("Q", "K"),  35:  ("Y", "F"),
    40:  ("I", "V"),  47:  ("S", "A"),  48:  ("F", "Y"),  49:  ("C", "U"),
    52:  ("T", "A"),  54:  ("T", "Q"),  60:  ("T", "A"),  67:  ("P", "N"),
    69:  ("N", "G"),  71:  ("T", "I"),  74:  ("G", "A"),  87:  ("K", "T"),
    99:  ("R", "C"),  102: ("G", "S"),  104: ("Y", "F"),  107: ("N", "S"),
    119: ("D", "E"),  120: ("N", "K"),  126: ("S", "T"),  137: ("E", "D"),
    139: ("F", "L"),  142: ("P", "S"),  143: ("E", "S"),  144: ("H", "Q"),
    148: ("D", "E"),  173: ("R", "H"),  177: ("H", "Q"),  178: ("T", "A"),
    181: ("R", "S"),  182: ("I", "T"),  184: ("Q", "K"),  188: ("M", "L"),
    192: ("N", "K"),  194: ("T", "F"),  195: ("S", "N"),
}

MOUSE_YELLOW_POSITIONS = [54, 49, 24, 139, 47, 60, 74, 144, 99, 177,
                           178, 104, 4, 52, 87, 102, 107, 142, 181, 48, 143]

HUMAN_YELLOW_POSITIONS = [87, 99, 47, 143, 60, 104, 142, 139, 181,
                           52, 48, 144, 54, 177, 102, 24, 3, 173, 178, 74]

def make_label(pos: int, is_mouse: bool) -> str:
    if pos not in MUTATION_TABLE:
        return str(pos)
    mouse_aa, human_aa = MUTATION_TABLE[pos]
    return f"{mouse_aa}{pos}{human_aa}" if is_mouse else f"{human_aa}{pos}{mouse_aa}"

def build_label_lookup(query_records):
    mouse_rows = [(sid, seq) for sid, seq in query_records
                  if "mouse" in sid.lower() or sid.lower().startswith("gpx6_level")]
    human_rows = [(sid, seq) for sid, seq in query_records
                  if "human" in sid.lower() or
                  (sid.lower().startswith("level") and not sid.lower().startswith("gpx6_level"))]
    label_lookup = {}
    for row_idx, (sid, _) in enumerate(mouse_rows):
        if row_idx < len(MOUSE_YELLOW_POSITIONS):
            label_lookup[sid] = make_label(MOUSE_YELLOW_POSITIONS[row_idx], is_mouse=True)
    for row_idx, (sid, _) in enumerate(human_rows):
        if row_idx < len(HUMAN_YELLOW_POSITIONS):
            label_lookup[sid] = make_label(HUMAN_YELLOW_POSITIONS[row_idx], is_mouse=False)
    return label_lookup

# â”€â”€ FASTA parser â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
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

# â”€â”€ Distance metrics â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
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

# â”€â”€ Softmax barycentric â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
def distances_to_xy(d_A, d_B, d_C, temp=SOFTMAX_TEMP):
    ds  = np.array([d_A, d_B, d_C], dtype=float)
    ds -= ds.min()
    w   = np.exp(-ds / temp)
    w  /= w.sum()
    xy  = w[0]*VERTICES[0] + w[1]*VERTICES[1] + w[2]*VERTICES[2]
    return float(xy[0]), float(xy[1]), w[0], w[1], w[2], int(np.argmin([d_A, d_B, d_C]))

# â”€â”€ Arc jitter â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
def arc_jitter(x, y, rng, arc_height=ARC_HEIGHT, noise_scale=0.008):
    """
    Lift a point toward the Ancestor vertex using a parabolic arc,
    then add a tiny amount of noise for visual separation.

    The arc is parameterised by t = x (position along the Mouseâ€“Human
    baseline, 0=Mouse, 1=Human). The parabolic term 4*t*(1-t) is exactly
    zero at both endpoints and peaks at t=0.5, so dots near either WT stay
    pinned to the baseline while midpoint dots bow upward.
    """
    t = np.clip(x, 0.0, 1.0)

    # Baseline point (linear interpolation between Mouse and Human vertices)
    base = (1.0 - t) * VERTICES[0] + t * VERTICES[1]

    # Direction from baseline point toward the Ancestor vertex
    direction = VERTICES[2] - base
    dist_to_ancestor = np.linalg.norm(direction)
    if dist_to_ancestor > 1e-9:
        direction = direction / dist_to_ancestor

    # Parabolic lift: peaks at t=0.5, zero at t=0 and t=1
    lift = 4.0 * t * (1.0 - t) * arc_height

    # Tiny noise for visual scatter
    noise = rng.normal(0.0, noise_scale, size=2)

    px = base[0] + lift * direction[0] + noise[0]
    py = base[1] + lift * direction[1] + noise[1]

    return float(px), float(py)

# â”€â”€ Grouping â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
GROUP_STYLE = {
    "Mouse GPX6 variants": dict(color="#2166ac", marker="o", zorder=4),
    "Human GPX6 variants": dict(color="#E87722", marker="o", zorder=5),
    "Other":               dict(color="#4dac26", marker="o", zorder=3),
}

GROUP_LABEL_COLOR = {
    "Mouse GPX6 variants": "#2166ac",
    "Human GPX6 variants": "#E87722",
    "Other":               "#4dac26",
}

def get_group(seq_id: str) -> str:
    sid = seq_id.lower()
    if "mouse" in sid or sid.startswith("gpx6_level"):
        return "Mouse GPX6 variants"
    if "human" in sid or (sid.startswith("level") and not sid.startswith("gpx6_level")):
        return "Human GPX6 variants"
    return "Other"

# â”€â”€ Plot helpers â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
def draw_triangle(ax):
    tri = np.vstack([VERTICES, VERTICES[0]])
    ax.plot(tri[:, 0], tri[:, 1], color="#333333", linewidth=1.8, zorder=2)
    for (vx, vy), col, label, off in zip(VERTICES, REF_COLORS, REF_LABELS, _LABEL_OFFSETS):
        ax.scatter(vx, vy, s=300, color=col, zorder=8,
                   edgecolors="white", linewidths=2.0)
        ax.text(vx + off["dx"], vy + off["dy"], label,
                ha=off["ha"], va=off["va"],
                fontsize=10, color="black",          # JCIM: 10 pt, black
                fontfamily="Arial",
                multialignment="center", zorder=9)

def scatter_groups(ax, rows, label_lookup: dict):
    texts  = []
    dot_xs = []
    dot_ys = []

    for group, style in GROUP_STYLE.items():
        pts = [r for r in rows if r["group"] == group]
        if not pts:
            continue

        for r in pts:
            px, py = arc_jitter(r["x"], r["y"], RNG)
            dot_xs.append(px)
            dot_ys.append(py)

            seq_label   = label_lookup.get(r["sequence_id"])
            is_labelled = seq_label is not None

            ax.scatter(
                px, py,
                color=style["color"],
                marker=style["marker"],
                s=110 if is_labelled else 55,
                edgecolor="black" if is_labelled else "white",
                linewidth=1.0 if is_labelled else 0.5,
                zorder=7 if is_labelled else style["zorder"],
                alpha=1.0 if is_labelled else 0.82,
            )

        ax.scatter([], [], color=style["color"], marker=style["marker"],
                   s=55, edgecolor="white", linewidth=0.5,
                   alpha=0.82, label=group)

    if texts:
        text_objs = [t for t, _, _ in texts]
        adjust_text(
            text_objs,arial
            x=np.array(dot_xs),
            y=np.array(dot_ys),
            ax=ax,
            expand=(1.5, 1.8),
            force_text=(0.35, 0.55),
            force_points=(0.45, 0.65),
            avoid_self=True,
            arrowprops=dict(
                arrowstyle="-",
                color="#999999",
                lw=0.55,
                shrinkA=2,
                shrinkB=3,
            ),
            min_arrow_len=5,
        )

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

# â”€â”€ Main â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
def main():
    print("Loading FASTA files...")
    all_records = load_all_records()
    print(f"Total unique sequences: {len(all_records)}\n")

    print("All sequence IDs:")
    for sid, seq in all_records:
        print(f"  '{sid}'  (len={len(seq)})")
    print()

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
            f"Expected 3 references, found {len(ref_records)}./n"
            f"  Found:   {found}\n"
            f"  Missing substrings: {missing}\n"
            f"  Update REFERENCE_IDS to match the headers printed above."
        )

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
    print(f"  Queries: {len(query_records)} | T={SOFTMAX_TEMP} | ARC_HEIGHT={ARC_HEIGHT}\n")

    label_lookup = build_label_lookup(query_records)
    print(f"  Labelled sequences: {len(label_lookup)}")
    for sid, lbl in label_lookup.items():
        print(f"    {sid}  â†’  {lbl}")

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
                "group":       get_group(qid),
                "nearest":     REF_LABELS[nearest].replace("\n", " "),
                "d_Mouse": d_A, "d_Human": d_B, "d_Ancestor": d_C,
                "w_Mouse": w_A, "w_Human": w_B, "w_Ancestor": w_C,
                "x": x, "y": y,
            })

        all_csv.extend(rows)

        print(f"[{metric_label}] nearest counts â†’ "
              f"Mouse:{nearest_counts[0]}  Human:{nearest_counts[1]}  Ancestor:{nearest_counts[2]}")

        fig_s, ax_s = plt.subplots(figsize=(9, 9))
        draw_triangle(ax_s)
        scatter_groups(ax_s, rows, label_lookup)
        finalise_panel(ax_s, metric_label)
        ax_s.legend(loc="upper right", framealpha=0.9,
                    title="Sequence group")
        fig_s.tight_layout()

        fig_name = FIG10_NAMES[idx] + ".png"
        out_s = Path(OUTPUT_PREFIX.parent / fig_name)
        out_s.parent.mkdir(parents=True, exist_ok=True)
        fig_s.savefig(out_s, dpi=300, bbox_inches="tight")
        plt.close(fig_s)
        print(f"  Saved: {out_s}\n")

    write_csv(Path(f"{OUTPUT_PREFIX}_coordinates.csv"), all_csv,
              ["sequence_id", "metric", "group", "nearest",
               "d_Mouse", "d_Human", "d_Ancestor",
               "w_Mouse", "w_Human", "w_Ancestor", "x", "y"])
    print(f"Saved CSV: {OUTPUT_PREFIX}_coordinates.csv")


if __name__ == "__main__":
    main()

