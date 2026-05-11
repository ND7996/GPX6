"""
GPX6 Sequence Analysis — original sequences (Human / Mouse / Ancestor)
ACS publication style:
  - No titles inside figure panels
  - 8 pt font body, 7 pt tick labels
  - Single-column (3.3 in) or double-column (7 in) widths
  - Horizontal barplots coloured by species
  - Sequences sorted by distance to ancestor within each group
"""
import sys
import os

ACS_PATH = r"D:\PhD_Thesis\analysis\FINAL_PUBLICATION_FIGURES"

if ACS_PATH not in sys.path:
    sys.path.append(ACS_PATH)

from acsfonts import acs_figure, save_acs_figure
from pathlib import Path
from itertools import combinations
import csv
import math

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import networkx as nx
import numpy as np

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices


# ═══════════════════════════════════════════════════════════════════
# PATHS
# ═══════════════════════════════════════════════════════════════════

INPUT_FASTA = Path(
    r"D:\PhD_Thesis\analysis\alignment\all_sequences_for_selection.fasta"
)

OUTPUT_DIR  = Path(r"D:\PhD_Thesis\analysis\FINAL_PUBLICATION_FIGURES\GPX6_original")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ═══════════════════════════════════════════════════════════════════
# ACS STYLE
# ═══════════════════════════════════════════════════════════════════

C_HUMAN    = "#E07B39"
C_MOUSE    = "#3A7EBF"
C_ANCESTOR = "#9B59B6"

plt.rcParams.update({
    "font.family"       : "Arial",
    "font.size"         : 8,
    "axes.labelsize"    : 8,
    "axes.titlesize"    : 8,
    "xtick.labelsize"   : 7,
    "ytick.labelsize"   : 7,
    "legend.fontsize"   : 7,
    "legend.frameon"    : True,
    "legend.edgecolor"  : "#cccccc",
    "legend.framealpha" : 0.9,
    "axes.linewidth"    : 0.6,
    "xtick.major.width" : 0.6,
    "ytick.major.width" : 0.6,
    "axes.spines.top"   : False,
    "axes.spines.right" : False,
    "figure.dpi"        : 300,
})


# ═══════════════════════════════════════════════════════════════════
# ALIGNMENT
# ═══════════════════════════════════════════════════════════════════

ALIGNER = PairwiseAligner()
ALIGNER.mode             = "global"
ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")
ALIGNER.open_gap_score   = -10
ALIGNER.extend_gap_score = -0.5


# ═══════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════

def load_records(path):
    recs = list(SeqIO.parse(path, "fasta"))
    if not recs:
        raise ValueError(f"No sequences found in {path}")
    return recs

def normsim(a, b):
    ab = float(ALIGNER.score(a, b))
    aa = float(ALIGNER.score(a, a))
    bb = float(ALIGNER.score(b, b))
    d  = math.sqrt(aa * bb)
    return max(0.0, min(1.0, ab / d)) if d else 0.0

def species_of(name):
    nm = name.upper()
    if "MOUSE" in nm:
        return "mouse"
    elif "HUMAN" in nm:
        return "human"
    else:
        return "ancestor"

def save_csv(path, rows, fields):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(rows)

COLOR_OF = {"human": C_HUMAN, "mouse": C_MOUSE, "ancestor": C_ANCESTOR}


# ═══════════════════════════════════════════════════════════════════
# LOAD
# ═══════════════════════════════════════════════════════════════════

print("Loading sequences …")
records = load_records(INPUT_FASTA)

names   = []
seqs    = []
species = []

for r in records:
    seq = str(r.seq).replace("U", "C")
    names.append(r.id)
    seqs.append(seq)
    species.append(species_of(r.id))

n = len(records)

# Find ancestor reference
anc_idx = next((i for i, s in enumerate(species) if s == "ancestor"), None)
anc_seq = seqs[anc_idx]


# ═══════════════════════════════════════════════════════════════════
# PAIRWISE DISTANCES
# ═══════════════════════════════════════════════════════════════════

print("Computing pairwise distances …")

sim_mat   = np.zeros((n, n))
pair_rows = []

for (i, j) in combinations(range(n), 2):
    s = normsim(seqs[i], seqs[j])
    sim_mat[i, j] = sim_mat[j, i] = s
    pair_rows.append({"seqA": names[i], "seqB": names[j],
                      "similarity": s, "distance": 1 - s})

np.fill_diagonal(sim_mat, 1.0)

dist_mat_full = 1 - sim_mat


# ═══════════════════════════════════════════════════════════════════
# HIERARCHICAL CLUSTERING
# ═══════════════════════════════════════════════════════════════════

print("Performing hierarchical clustering …")

import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform

dist_cond = squareform(dist_mat_full)

Z = sch.linkage(dist_cond, method='ward')


# ═══════════════════════════════════════════════════════════════════
# HEATMAP (HIERARCHICAL CLUSTERED)
# ═══════════════════════════════════════════════════════════════════

print("Generating clustered heatmap …")

import seaborn as sns
import pandas as pd
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

short_names = []
for nm in names:
    if "HUMAN" in nm.upper():
        s = nm.replace("HUMAN_GPX6_", "Hs_").replace("_chainA", "")
    elif "MOUSE" in nm.upper():
        s = nm.replace("MOUSE_GPX6_", "Mm_").replace("_chainA", "")
    else:
        s = "Anc_node"
    short_names.append(s)

dist_df = pd.DataFrame(dist_mat_full, index=names, columns=names)

short_df = dist_df.copy()
short_df.index = short_names
short_df.columns = short_names

row_colors = [COLOR_OF[sp] for sp in species]

fig_height = max(8.5, len(names) * 0.19 + 1.5)

cg = sns.clustermap(
    short_df,
    row_linkage=Z,
    col_linkage=Z,
    cmap="RdYlBu_r",
    vmin=0,
    vmax=1,
    figsize=(7.2, fig_height),
    row_cluster=True,
    col_cluster=True,
    dendrogram_ratio=(0.18, 0.02),
    cbar_pos=(0.92, 0.25, 0.015, 0.55),
    xticklabels=False,
    yticklabels=True,
    row_colors=row_colors,
    colors_ratio=0.025,
    linewidths=0,
    tree_kws={"linewidth": 1.1}
)

reordered = cg.dendrogram_row.reordered_ind
reordered_species = [species[i] for i in reordered]

for label, sp in zip(cg.ax_heatmap.yaxis.get_ticklabels(), reordered_species):
    label.set_color(COLOR_OF[sp])

cg.ax_heatmap.tick_params(axis="y", labelsize=6)

legend_elements = [
    mpatches.Patch(facecolor=C_HUMAN, label="Human"),
    mpatches.Patch(facecolor=C_MOUSE, label="Mouse"),
    mpatches.Patch(facecolor=C_ANCESTOR, label="Ancestor")
]

cg.fig.legend(handles=legend_elements,
              loc="upper right",
              bbox_to_anchor=(0.98, 0.98),
              fontsize=7)

cg.ax_cbar.set_ylabel("Normalised distance\n(BLOSUM62)", fontsize=8)

plt.tight_layout()

plt.savefig(
    OUTPUT_DIR / "gpx6_heatmap_clustered.png",
    dpi=600,
    bbox_inches="tight"
)

plt.close()

print("Heatmap saved.")


# ═══════════════════════════════════════════════════════════════════
# SEQUENCE SIMILARITY NETWORK
# ═══════════════════════════════════════════════════════════════════

print("Building SSN …")

SIMILARITY_THRESHOLD = 0.75

G = nx.Graph()

for nm, sp in zip(names, species):
    G.add_node(nm, species=sp)

for row in pair_rows:
    if row["similarity"] >= SIMILARITY_THRESHOLD:
        G.add_edge(row["seqA"], row["seqB"], weight=row["similarity"])

pos = nx.kamada_kawai_layout(G)

mouse_nodes    = [nd for nd,d in G.nodes(data=True) if d["species"]=="mouse"]
human_nodes    = [nd for nd,d in G.nodes(data=True) if d["species"]=="human"]
ancestor_nodes = [nd for nd,d in G.nodes(data=True) if d["species"]=="ancestor"]

fig, ax = plt.subplots(figsize=(7,7))

nx.draw_networkx_edges(G,pos,alpha=0.3,edge_color="#aaaaaa")

nx.draw_networkx_nodes(G,pos,nodelist=mouse_nodes,node_color=C_MOUSE,node_size=220)
nx.draw_networkx_nodes(G,pos,nodelist=human_nodes,node_color=C_HUMAN,node_size=220)
nx.draw_networkx_nodes(G,pos,nodelist=ancestor_nodes,node_color=C_ANCESTOR,node_size=320)

nx.draw_networkx_labels(G,pos,font_size=5)

ax.axis("off")

plt.savefig(
    OUTPUT_DIR / "gpx6_ssn_network.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("SSN saved.")

print("\nDone.")