"""
GPX6 Sequence Analysis â€” JCIM Publication-Quality Clustered Heatmap
JCIM compliant: Arial 8pt, black text, inward ticks, frameless legend
Labels: M1-M20 (Mouse), H1-H19 (Human), M_WT, H_WT, Ancestor
"""

import sys
import os

# =================== JCIM PUBLICATION STYLE ===================
ACS_PATH = r"./analysis_scripts/Scripts_to_generate_figures/Figures"
if ACS_PATH not in sys.path:
    sys.path.append(ACS_PATH)

from acsfonts import set_jcim_style
set_jcim_style()

# JCIM standard: Arial 8pt for body text, but labels may need to be larger for readability
import matplotlib.pyplot as plt
plt.rcParams.update({
    "font.size": 8,
    "font.family": "Arial",
    "axes.labelsize": 8,
    "xtick.labelsize": 8,  # Increased from 7 to 8
    "ytick.labelsize": 8,  # Increased from 7 to 8
    "legend.fontsize": 8,
    "text.color": "black",
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "xtick.direction": "in",
    "ytick.direction": "in"
})

from pathlib import Path
from itertools import combinations
import math
import re
import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices

import seaborn as sns
import matplotlib.patches as mpatches

from matplotlib.colors import LinearSegmentedColormap
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage

# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# PATHS
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

FASTA_MOUSE    = Path(r"./prep_structures\MOUSE\mutant_pdbs\sequences_mouse.fasta")
FASTA_HUMAN    = Path(r"./prep_structures\HUMAN\mutant_pdbs\sequences_human.fasta")
FASTA_ANCESTOR = Path(r"./analysis_scripts/alignment/gpx6_human_mouse_input.fasta")

OUTPUT_DIR = Path(r"./analysis_scripts/Scripts_to_generate_figures/Figures/GPX6_original")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# COLORS
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

C_HUMAN    = "#E41A1C"  # Red
C_MOUSE    = "#4DAF4A"  # Green
C_ANCESTOR = "#377EB8"  # Blue

COLOR_OF = {
    "human": C_HUMAN,
    "mouse": C_MOUSE,
    "ancestor": C_ANCESTOR
}

def create_cmap():
    return LinearSegmentedColormap.from_list(
        "gpx6_species_cmap",
        [C_ANCESTOR, C_MOUSE, C_HUMAN],
        N=512
    )

heatmap_cmap = create_cmap()

# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# ALIGNMENT
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

ALIGNER = PairwiseAligner()
ALIGNER.mode = "global"
ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")
ALIGNER.open_gap_score = -10
ALIGNER.extend_gap_score = -0.5

def normsim(a, b):
    ab = float(ALIGNER.score(a, b))
    aa = float(ALIGNER.score(a, a))
    bb = float(ALIGNER.score(b, b))
    denom = math.sqrt(aa * bb)
    return ab / denom if denom else 0.0

def species_of(name):
    n = name.upper()
    if "MOUSE" in n or "MUS" in n or "GPX6_LEVEL" in n:
        return "mouse"
    if "HUMAN" in n or "HOMO" in n or "LEVEL" in n:
        return "human"
    return "ancestor"

def extract_level_number(name):
    """Extract level number from sequence name (e.g., 'Level 1' -> 1)"""
    match = re.search(r'level[_ ]?(\d+)', name, re.IGNORECASE)
    if match:
        return int(match.group(1))
    return None

def get_short_name(idx, name, seq, wt_mouse_idx, wt_human_idx, anc_idx):
    """Generate short name: M1-M20, H1-H19, M_WT, H_WT, Ancestor"""
    species = species_of(name)
    
    # Ancestor
    if idx == anc_idx or species == "ancestor":
        return "Ancestor"
    
    # Mouse WT
    if idx == wt_mouse_idx:
        return "M_WT"
    
    # Human WT
    if idx == wt_human_idx:
        return "H_WT"
    
    # Mouse variants (M1, M2, M3...)
    if species == "mouse":
        level = extract_level_number(name)
        if level is not None:
            return f"M{level}"
        # Try to extract number from name
        numbers = re.findall(r'\d+', name)
        if numbers:
            return f"M{numbers[0]}"
        return "M_var"
    
    # Human variants (H1, H2, H3...)
    if species == "human":
        level = extract_level_number(name)
        if level is not None:
            return f"H{level}"
        numbers = re.findall(r'\d+', name)
        if numbers:
            return f"H{numbers[0]}"
        return "H_var"
    
    return name[:8]

# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# LOAD SEQUENCES
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

print("Loading sequences...")
records = []

for f, prefix in [
    (FASTA_MOUSE, "MOUSE"),
    (FASTA_HUMAN, "HUMAN"),
    (FASTA_ANCESTOR, "ANCESTOR")
]:
    if f.exists():
        for r in SeqIO.parse(f, "fasta"):
            clean_id = r.id.replace("|", "_").replace(" ", "_")
            r.id = f"{prefix}_{clean_id}"
            records.append(r)

if not records:
    raise ValueError("No sequences found!")

names = [r.id for r in records]
seqs = [str(r.seq).replace("U", "C") for r in records]
species = [species_of(n) for n in names]

n = len(seqs)
print(f"\nTotal sequences loaded: {n}")

# Find indices for WT and ancestor
wt_mouse_idx = None
wt_human_idx = None
anc_idx = None

for idx, (name, sp) in enumerate(zip(names, species)):
    if sp == "mouse" and ("original" in name.lower() or "wild" in name.lower() or "mousecys" in name.lower() or idx == 0):
        if wt_mouse_idx is None:
            wt_mouse_idx = idx
            print(f"  Mouse WT: {name}")
    if sp == "human" and ("original" in name.lower() or "wild" in name.lower() or ("U49C" in name) or "humansec" in name.lower()):
        if wt_human_idx is None:
            wt_human_idx = idx
            print(f"  Human WT: {name}")
    if sp == "ancestor":
        if anc_idx is None:
            anc_idx = idx
            print(f"  Ancestor: {name}")

# Generate short names
short_names = []
for idx, (nm, seq, sp) in enumerate(zip(names, seqs, species)):
    short = get_short_name(idx, nm, seq, wt_mouse_idx, wt_human_idx, anc_idx)
    short_names.append(short)

# Print short names for verification
print(f"\nShort names (first 10): {short_names[:10]}")
print(f"  Mouse labels: {[s for s in short_names if s.startswith('M')]}")
print(f"  Human labels: {[s for s in short_names if s.startswith('H')]}")
print(f"  WT: M_WT, H_WT")
print(f"  Ancestor: Ancestor")

mouse_count = sum(1 for sp in species if sp == 'mouse')
human_count = sum(1 for sp in species if sp == 'human')
anc_count = sum(1 for sp in species if sp == 'ancestor')
print(f"\n  Mouse: {mouse_count} (M_WT + M1-M{mouse_count-1})")
print(f"  Human: {human_count} (H_WT + H1-H{human_count-1})")
print(f"  Ancestor: {anc_count}")

# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# DISTANCE MATRIX
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

print("\nComputing pairwise distances...")
sim = np.zeros((n, n))

for i, j in combinations(range(n), 2):
    s = normsim(seqs[i], seqs[j])
    sim[i, j] = sim[j, i] = s

np.fill_diagonal(sim, 1.0)
dist = 1 - sim

df_dist = pd.DataFrame(dist, index=short_names, columns=short_names)

row_colors = [COLOR_OF[s] for s in species]
col_colors = row_colors

# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# CLUSTERING
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

condensed = squareform(dist)
linkage_matrix = linkage(condensed, method="average")

# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# HEATMAP
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

# Larger figure size for 40+ sequences with readable labels
# Increase figure size to accommodate readable fonts
fig_size = max(14, min(18, n / 2.5))  # Made larger for better label visibility
print(f"Figure size: {fig_size} x {fig_size} inches for {n} sequences")

cg = sns.clustermap(
    df_dist,
    row_linkage=linkage_matrix,
    col_linkage=linkage_matrix,
    cmap=heatmap_cmap,
    figsize=(fig_size, fig_size),
    row_colors=row_colors,
    col_colors=row_colors,
    linewidths=0.3,
    cbar_kws={
        "label": "Genetic Distance",
        "ticks": [0, 0.25, 0.5, 0.75, 1],
        "shrink": 0.6,
        "aspect": 30
    },
    cbar_pos=(0.92, 0.35, 0.012, 0.35),
    dendrogram_ratio=0.12,
    colors_ratio=0.02
)

# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# LAYOUT FIX
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

cg.gs.update(left=0.1, right=0.88)
if cg.cax:
    cg.cax.set_position([0.91, 0.35, 0.010, 0.35])

# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# LABEL FORMATTING (READABLE SIZE - 9-10pt for labels)
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

# Use larger font for labels - readability is crucial
label_fs = 9 if n > 40 else 10  # Minimum 10pt for readability

cg.ax_heatmap.tick_params(axis='y', pad=8, colors='black', labelsize=label_fs)
cg.ax_heatmap.tick_params(axis='x', pad=8, colors='black', labelsize=label_fs)

# Rotate x labels if needed, but keep them readable
rotation_angle = 90 if n > 30 else 45
cg.ax_heatmap.set_xticklabels(
    cg.ax_heatmap.get_xticklabels(),
    rotation=rotation_angle,
    fontsize=label_fs,
    color="black",
    ha='center' if rotation_angle == 90 else 'right'
)

cg.ax_heatmap.set_yticklabels(
    cg.ax_heatmap.get_yticklabels(),
    rotation=0,
    fontsize=label_fs,
    color="black"
)

# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# CLEAN COLOR STRIPS
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

cg.ax_row_colors.set_xticks([])
cg.ax_col_colors.set_yticks([])
cg.ax_row_colors.set_frame_on(False)
cg.ax_col_colors.set_frame_on(False)

# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# LEGEND (JCIM STYLE - FRAMELESS, 8-10pt)
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

legend_handles = [
    mpatches.Patch(color=C_HUMAN, label="Human (H)"),
    mpatches.Patch(color=C_MOUSE, label="Mouse (M)"),
    mpatches.Patch(color=C_ANCESTOR, label="Ancestor")
]

legend = cg.fig.legend(
    handles=legend_handles,
    loc="upper right",
    bbox_to_anchor=(1.02, 0.98),
    frameon=False,
    fontsize=9,  # Increased from 8 to 9 for better readability
    handlelength=1.5,
    handleheight=1.0,
    borderaxespad=0.3,
)

for text in legend.get_texts():
    text.set_color("black")

# Colorbar styling
if cg.cax:
    cg.cax.yaxis.label.set_color("black")
    cg.cax.yaxis.label.set_fontsize(9)  # Increased from 8 to 9
    cg.cax.tick_params(colors='black', labelsize=8)

# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# SAVE OUTPUT
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

cg.savefig(OUTPUT_DIR / "GPX6_heatmap.pdf", dpi=600, bbox_inches="tight")
cg.savefig(OUTPUT_DIR / "GPX6_heatmap.png", dpi=600, bbox_inches="tight")
cg.savefig(OUTPUT_DIR / "GPX6_heatmap.tiff", dpi=600, bbox_inches="tight")

plt.close()

print(f"\nâœ… DONE â€” JCIM compliant heatmap with readable labels")
print(f"   Total sequences: {n}")
print(f"   Mouse: M_WT, M1-M{mouse_count-1}")
print(f"   Human: H_WT, H1-H{human_count-1}")
print(f"   Ancestor: Ancestor")
print(f"   Label size: {label_fs}pt (optimized for readability)")
print(f"   Legend size: 10pt")
print(f"   Font: Arial")

