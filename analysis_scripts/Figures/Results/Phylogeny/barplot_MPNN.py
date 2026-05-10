# ═══════════════════════════════════════════════════════════════════
# Horizontal Barplot: Distance from Ancestor Reference
# Exports PNG (300 dpi) + SVG (for Figma / Illustrator)
# ═══════════════════════════════════════════════════════════════════

from pathlib import Path
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
import pandas as pd

matplotlib.rcParams['svg.fonttype'] = 'none'   # keeps text as <text> in SVG (editable in Figma)

# ── PATHS ───────────────────────────────────────────────────────────
FASTA_ANCESTOR = Path(r"D:\PhD_Thesis\MPNN\results\proteinmpnn_all\proteinmpnn_outputs\combined_ANCESTOR_labeled.fasta")
FASTA_HUMAN    = Path(r"D:\PhD_Thesis\MPNN\results\proteinmpnn_all\proteinmpnn_outputs\combined_HUMAN_labeled.fasta")
FASTA_MOUSE    = Path(r"D:\PhD_Thesis\MPNN\results\proteinmpnn_all\proteinmpnn_outputs\combined_MOUSE_labeled.fasta")
OUTPUT_DIR     = Path(r"D:\PhD_Thesis\analysis\FINAL_PUBLICATION_FIGURES\Phylogeny")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ── WT sequences ────────────────────────────────────────────────────
ANC_WT = "PQKMKMDCNKGVTGTIYEYGALTLNGEEYIQFKQYAGKHVLFVNVATYGLTAQYPELNALQEELKHFGVIVLGFPCNQFGKQEPGKNSEILSGLKYVRPGGGFVPNFQLFEKGDVNGEKEQKVFTFLKNSCPPTSDLLGSSSQLFWEPMKVHDIRWNFEKFLVGPDGVPVMRWFHRAPVSTVKSDILEYLKQF"
HUM_WT = "PQNRKVDNKGVTGTIYEYGALTLNGEEYIQFKQFAGKVLFVNVAAYLAAQYPELNALQEELKNFGVIVLAFPCNQFGKQEPGTNSEILLGLKYVCPGSGFVPSFQLFEKGDVNGEKEQKVFTFLKNSPPTSDLLGSSSQLFWEPMKVDIRWNFEKFLVGPDGVPVMWFQAPVSTVKSDILEYLKQFNT"
MOU_WT = "PQKSKVDNKGVTGTVYEYGANTIDGGEFVNFQQYAGKILFVNVASFCGLTATYPELNTLQEELKPFNVTVLGFPCNQFGKQEPGKNSEILLGLKYVRPGGGYVPNFQLFEKGDVNGDNEQKVFSFLKNSPPTSELFGSPELFWDPMKVDIRWNFEKFLVGPDGVPVMRWFTPVRIVQSDIMEYLNQTS"

# ── COLORS ──────────────────────────────────────────────────────────
COLOR_OF = {
    "human":    "#E07B39",
    "mouse":    "#3A7EBF",
    "ancestor": "#9B59B6",
}

# ── Read sequences from FASTA ────────────────────────────────────────
def read_sequences(fasta_path, seq_type):
    sequences = {}
    if fasta_path.exists():
        for record in SeqIO.parse(fasta_path, "fasta"):
            sequences[f"{seq_type}_{record.id}"] = str(record.seq)
    return sequences

# ── BLOSUM62 distance (1 − fractional identity on global alignment) ──
def calculate_distance(seq1, seq2):
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score   = -10
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(seq1, seq2)
    if not alignments:
        return 1.0

    best     = alignments[0]
    al1, al2 = best[0], best[1]
    matches  = sum(a == b for a, b in zip(al1, al2) if a != '-' and b != '-')
    total    = sum(1      for a, b in zip(al1, al2) if a != '-' or  b != '-')
    identity = matches / total if total > 0 else 0
    return 1.0 - identity

# ────────────────────────────────────────────────────────────────────
print("="*60)
print("Generating Barplot: Distance from Ancestor Reference")
print("="*60)

# 1. Collect sequences
print("\n1. Collecting sequences...")
all_sequences = {
    "ancestor_Anc_node": ANC_WT,
    "human_Hs_WT":       HUM_WT,
    "mouse_Mm_WT":       MOU_WT,
}
all_sequences.update(read_sequences(FASTA_ANCESTOR, "ancestor"))
all_sequences.update(read_sequences(FASTA_HUMAN,    "human"))
all_sequences.update(read_sequences(FASTA_MOUSE,    "mouse"))
print(f"   Total sequences: {len(all_sequences)}")

# 2. Calculate distances
print("\n2. Calculating distances from ancestor reference...")
records = []
for seq_name, seq in all_sequences.items():
    if seq_name.startswith("ancestor"):
        sp = "ancestor"
    elif seq_name.startswith("human"):
        sp = "human"
    elif seq_name.startswith("mouse"):
        sp = "mouse"
    else:
        sp = "other"

    bare  = seq_name[len(sp)+1:]
    label = bare if bare else seq_name
    dist  = calculate_distance(ANC_WT, seq)
    records.append({"label": label, "species": sp, "distance": dist})

df = pd.DataFrame(records)
df = df.sort_values("distance", ascending=True).reset_index(drop=True)

print(f"   Calculated {len(df)} distances")
print(f"   Min : {df['distance'].min():.4f}")
print(f"   Max : {df['distance'].max():.4f}")
print(f"   Mean: {df['distance'].mean():.4f}")

df.to_csv(OUTPUT_DIR / "distances_from_ancestor.csv", index=False)

# ── Build figure ─────────────────────────────────────────────────────
print("\n3. Generating figure...")

n      = len(df)
bar_h  = 0.55
fig_h  = max(5, n * 0.32 + 1.2)
fig_w  = 7.5

fig, ax = plt.subplots(figsize=(fig_w, fig_h))

y_pos  = np.arange(n)
colors = [COLOR_OF.get(sp, "#AAAAAA") for sp in df["species"]]
x_max  = df["distance"].max()

ax.barh(y_pos, df["distance"], height=bar_h,
        color=colors, edgecolor="none", linewidth=0)

for i, (dist, sp) in enumerate(zip(df["distance"], df["species"])):
    if dist == 0:
        continue
    ax.text(dist + x_max * 0.005, i, f"{dist:.2f}",
            va="center", ha="left", fontsize=7.5, color="#333333")

ax.set_yticks(y_pos)
ax.set_yticklabels(df["label"], fontsize=9)
ax.tick_params(axis="y", length=0)
ax.set_xlabel("Distance to ancestor (BLOSUM62)", fontsize=10)
ax.set_xlim(0, x_max * 1.13)
ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.2f"))
ax.tick_params(axis="x", labelsize=9)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_linewidth(0.8)
ax.xaxis.grid(True, color="#dddddd", linewidth=0.6, zorder=0)
ax.set_axisbelow(True)

legend_patches = [
    mpatches.Patch(facecolor=COLOR_OF["human"],    label="Human"),
    mpatches.Patch(facecolor=COLOR_OF["mouse"],    label="Mouse"),
    mpatches.Patch(facecolor=COLOR_OF["ancestor"], label="Ancestor"),
]
ax.legend(handles=legend_patches, loc="lower right", fontsize=9,
          frameon=True, framealpha=0.9, edgecolor="#cccccc",
          handlelength=1.2, handleheight=1.0)

plt.tight_layout(pad=0.8)

# ── Save PNG (300 dpi, for publication) ──────────────────────────────
png_path = OUTPUT_DIR / "Barplot_Distance_from_Ancestor.png"
fig.savefig(png_path, dpi=300, bbox_inches="tight")
print(f"   ✓ PNG saved : {png_path}")

# ── Save SVG (Figma / Illustrator — tiny file, fully editable) ───────
svg_path = OUTPUT_DIR / "Barplot_Distance_from_Ancestor.svg"
fig.savefig(svg_path, format="svg", bbox_inches="tight")
print(f"   ✓ SVG saved : {svg_path}")

plt.close()

print("\n" + "="*60)
print("DONE")
print("="*60)
print(f"\nOutput directory: {OUTPUT_DIR}")
print("\nFiles generated:")
print("  Barplot_Distance_from_Ancestor.png  ← publication figure (300 dpi)")
print("  Barplot_Distance_from_Ancestor.svg  ← drag into Figma; all text & bars editable")
print("  distances_from_ancestor.csv         ← raw distance table")