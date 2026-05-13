import sys
import os

# =================== JCIM PUBLICATION STYLE ===================
ACS_PATH = r"./analysis_scripts/Scripts_to_generate_figures/Figures"
if ACS_PATH not in sys.path:
    sys.path.append(ACS_PATH)

from acsfonts import set_jcim_style
set_jcim_style()

# Force all text black
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.color": "black",
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black"
})

from pathlib import Path
import csv
import math
import numpy as np

# â”€â”€ Paths â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
FASTA_HUMAN    = Path(r"./analysis_scripts/MPNN/Human/GPX6MPNN_generated/all_generated_sequences.fasta")
FASTA_MOUSE    = Path(r"./analysis_scripts/MPNN/Mouse/GPX6MPNN_generated/all_generated_sequences.fasta")
OUTPUT_PREFIX = Path("results/gpx6_mpnn_ternary_sharp")

# â”€â”€ Reference sequences â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
ANC_SEQ = "PQKMKMDCNKGVTGTIYEYGALTLNGEEYIQFKQYAGKHVLFVNVATYGLTAQYPELNALQEELKHFGVIVLGFPCNQFGKQEPGKNSEILSGLKYVRPGGGFVPNFQLFEKGDVNGEKEQKVFTFLKNSCPPTSDLLGSSSQLFWEPMKVHDIRWNFEKFLVGPDGVPVMRWFHRAPVSTVKSDILEYLKQF"
HUM_SEQ = "PQNRKVDNKGVTGTIYEYGALTLNGEEYIQFKQFAGKVLFVNVAAYLAAQYPELNALQEELKNFGVIVLAFPCNQFGKQEPGTNSEILLGLKYVCPGSGFVPSFQLFEKGDVNGEKEQKVFTFLKNSPPTSDLLGSSSQLFWEPMKVDIRWNFEKFLVGPDGVPVMWFQAPVSTVKSDILEYLKQFNT"
MOU_SEQ = "PQKSKVDNKGVTGTVYEYGANTIDGGEFVNFQQYAGKILFVNVASFCGLTATYPELNTLQEELKPFNVTVLGFPCNQFGKQEPGKNSEILLGLKYVRPGGGYVPNFQLFEKGDVNGDNEQKVFSFLKNSPPTSELFGSPELFWDPMKVDIRWNFEKFLVGPDGVPVMRWFTPVRIVQSDIMEYLNQTS"

REF_SEQS    = [MOU_SEQ, HUM_SEQ, ANC_SEQ]
MAX_REF_LEN = max(len(MOU_SEQ), len(HUM_SEQ), len(ANC_SEQ))

VERTICES = np.array([
    [0.0, 0.0],
    [1.0, 0.0],
    [0.5, math.sqrt(3) / 2],
])
REF_LABELS = ["Mouse\nGPX6 WT", "Human\nGPX6 WT", "Ancestor\n(node_25)"]
REF_COLORS = ["#2166ac", "#d6604d", "#4a9c6d"]

# MODIFIED: Adjusted offsets to bring labels closer to vertices
_LABEL_OFFSETS = [
    dict(dx=-0.03, dy=-0.03, ha="right",  va="top"),      # Mouse (bottom-left)
    dict(dx= 0.03, dy=-0.03, ha="left",   va="top"),      # Human (bottom-right)
    dict(dx= 0.00, dy= 0.03, ha="center", va="bottom"),   # Ancestor (top)
]

SOFTMAX_TEMP = 0.03
JITTER_SCALE = 0.006
RNG          = np.random.default_rng(seed=42)

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
    max_len = max(len(a), len(b))
    a_pad   = a.ljust(max_len, '-')
    b_pad   = b.ljust(max_len, '-')
    cols    = [(ca, cb) for ca, cb in zip(a_pad, b_pad)
               if not (ca == '-' and cb == '-')]
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
def distances_to_xy(d_A, d_B, d_C):
    ds  = np.array([d_A, d_B, d_C], dtype=float)
    ds -= ds.min()
    w   = np.exp(-ds / SOFTMAX_TEMP)
    w  /= w.sum()
    xy  = w[0]*VERTICES[0] + w[1]*VERTICES[1] + w[2]*VERTICES[2]
    return float(xy[0]), float(xy[1]), w[0], w[1], w[2], int(np.argmin([d_A, d_B, d_C]))

# â”€â”€ Groups â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
GROUP_STYLE = {
    "Human-designed": dict(color="#d6604d", marker="o", zorder=5),
    "Mouse-designed": dict(color="#2166ac", marker="o", zorder=4),
}

def compute_rows(records, group_label, metric_fn):
    rows = []
    for sid, seq in records:
        shared_scale = max(MAX_REF_LEN, len(seq))
        d_A = metric_fn(seq, REF_SEQS[0], shared_scale)
        d_B = metric_fn(seq, REF_SEQS[1], shared_scale)
        d_C = metric_fn(seq, REF_SEQS[2], shared_scale)
        x, y, w_A, w_B, w_C, nearest = distances_to_xy(d_A, d_B, d_C)
        rows.append({
            "sequence_id": sid, "group": group_label,
            "nearest":  REF_LABELS[nearest].replace("\n", " "),
            "d_Mouse":  d_A, "d_Human": d_B, "d_Ancestor": d_C,
            "w_Mouse":  w_A, "w_Human": w_B, "w_Ancestor": w_C,
            "x": x, "y": y,
        })
    return rows

def draw_triangle(ax):
    tri = np.vstack([VERTICES, VERTICES[0]])
    ax.plot(tri[:, 0], tri[:, 1], color="#333333", linewidth=1.8, zorder=2)
    for (vx, vy), col, label, off in zip(VERTICES, REF_COLORS, REF_LABELS, _LABEL_OFFSETS):
        ax.scatter(vx, vy, s=300, color=col, zorder=8,
                   edgecolors="white", linewidths=2.0)
        # Text in BLACK (JCIM requirement) - now positioned closer to vertices
        ax.text(vx + off["dx"], vy + off["dy"], label,
                ha=off["ha"], va=off["va"],
                color="black",  # BLACK text
                multialignment="center", zorder=9)

def scatter_groups(ax, rows_h, rows_m):
    for rows, group in [(rows_m, "Mouse-designed"), (rows_h, "Human-designed")]:
        if not rows:
            continue
        g  = GROUP_STYLE[group]
        xs = np.array([r["x"] for r in rows]) + RNG.uniform(-JITTER_SCALE, JITTER_SCALE, len(rows))
        ys = np.array([r["y"] for r in rows]) + RNG.uniform(-JITTER_SCALE, JITTER_SCALE, len(rows))
        ax.scatter(xs, ys, color=g["color"], marker=g["marker"],
                   s=70, edgecolor="white", linewidth=0.7,
                   zorder=g["zorder"], alpha=0.92, label=group)

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
    human_records = parse_fasta(FASTA_HUMAN)
    mouse_records = parse_fasta(FASTA_MOUSE)
    print(f"Loaded: Human={len(human_records)}, Mouse={len(mouse_records)}")

    FIG_NAMES = ["Fig11a", "Fig11b", "Fig11c", "Fig11d"]
    all_csv = []

    for idx, (metric_label, metric_fn) in enumerate(METRICS):
        rows_h = compute_rows(human_records, "Human-designed", metric_fn)
        rows_m = compute_rows(mouse_records, "Mouse-designed", metric_fn)
        all_csv.extend(rows_h + rows_m)

        fig, ax = plt.subplots(figsize=(7, 7))
        draw_triangle(ax)
        scatter_groups(ax, rows_h, rows_m)
        finalise_panel(ax, metric_label)
        
        # Legend with black text
        leg = ax.legend(loc="lower right", frameon=False, fontsize=8)
        for text in leg.get_texts():
            text.set_color("black")
        
        fig.tight_layout()

        out_s = Path(f"{OUTPUT_PREFIX}_{FIG_NAMES[idx]}.png")
        out_s.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_s, dpi=600, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {out_s}")

    # Save combined CSV
    write_csv(Path(f"{OUTPUT_PREFIX}_coordinates.csv"), all_csv,
              ["sequence_id","group","nearest",
               "d_Mouse","d_Human","d_Ancestor",
               "w_Mouse","w_Human","w_Ancestor","x","y"])
    print(f"Saved CSV: {OUTPUT_PREFIX}_coordinates.csv")

if __name__ == "__main__":
    main()

