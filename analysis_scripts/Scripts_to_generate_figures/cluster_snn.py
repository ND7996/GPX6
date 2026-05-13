"""
GPX6 Sequence Analysis â€” original sequences (Human / Mouse / Ancestor)
ACS publication style.

CHANGES:
  - acsfonts used only for font consistency (no bold anywhere)
  - Improved label placement to avoid overlaps
  - All fontweight="bold" removed throughout
"""
import sys
import os

ACS_PATH = r"./analysis_scripts/Scripts_to_generate_figures/Figures"
if ACS_PATH not in sys.path:
    sys.path.append(ACS_PATH)

from acsfonts import acs_figure, save_acs_figure
from pathlib import Path
from itertools import combinations
import csv
import math

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon as MplPolygon
import networkx as nx
import numpy as np
from scipy.spatial import ConvexHull

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices


# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# PATHS
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

INPUT_FASTA = Path(
    r"./analysis_scripts/alignment/gpx6_human_mouse_input.fasta"
)
OUTPUT_DIR = Path(
    r"./analysis_scripts/Scripts_to_generate_figures/Figures/GPX6_original"
)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# STYLE
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

C_HUMAN    = "#E07B39"
C_MOUSE    = "#3A7EBF"
C_ANCESTOR = "#9B59B6"
COLOR_OF   = {"human": C_HUMAN, "mouse": C_MOUSE, "ancestor": C_ANCESTOR}

HULL_OF = {
    "human":    mcolors.to_rgba(C_HUMAN,    alpha=0.10),
    "mouse":    mcolors.to_rgba(C_MOUSE,    alpha=0.10),
    "ancestor": mcolors.to_rgba(C_ANCESTOR, alpha=0.25),
}

plt.rcParams.update({
    "font.family"       : "Arial",
    "font.size"         : 8,
    "axes.labelsize"    : 8,
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


# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# ALIGNMENT
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

ALIGNER = PairwiseAligner()
ALIGNER.mode                = "global"
ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")
ALIGNER.open_gap_score      = -10
ALIGNER.extend_gap_score    = -0.5


# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# HELPERS
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

def normsim(a, b):
    ab = float(ALIGNER.score(a, b))
    aa = float(ALIGNER.score(a, a))
    bb = float(ALIGNER.score(b, b))
    d  = math.sqrt(aa * bb)
    return max(0.0, min(1.0, ab / d)) if d else 0.0

def species_of(name):
    nm = name.upper()
    if "MOUSE" in nm: return "mouse"
    if "HUMAN" in nm: return "human"
    return "ancestor"

def short_name(nm):
    if "HUMAN" in nm.upper():
        return nm.replace("HUMAN_GPX6_", "Hs_").replace("_chainA", "")
    if "MOUSE" in nm.upper():
        return nm.replace("MOUSE_GPX6_", "Mm_").replace("_chainA", "")
    return "Anc_node"

def is_wt(name):
    return "WT" in name.upper() or "ANCESTOR" in name.upper() or "NODE" in name.upper()

def draw_hull(ax, points, rgba, pad=0.06):
    pts = np.array(points)
    if len(pts) < 3:
        for p in pts:
            ax.add_patch(plt.Circle(p, pad * 2.5, color=rgba, zorder=0))
        return
    try:
        hull  = ConvexHull(pts)
        hverts = pts[hull.vertices]
        centre = pts.mean(axis=0)
        norms  = np.linalg.norm(hverts - centre, axis=1, keepdims=True)
        expanded = centre + (1 + pad / np.maximum(norms, 1e-9)) * (hverts - centre)
        ax.add_patch(MplPolygon(
            expanded, closed=True,
            facecolor=rgba,
            edgecolor=(*rgba[:3], 0.50),
            linewidth=1.0, zorder=0,
        ))
    except Exception:
        pass


# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# LOAD SEQUENCES
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

print("Loading sequences â€¦")
records = list(SeqIO.parse(INPUT_FASTA, "fasta"))
if not records:
    raise ValueError(f"No sequences in {INPUT_FASTA}")

names, seqs, species = [], [], []
for r in records:
    seq = str(r.seq).replace("U", "C")
    names.append(r.id)
    seqs.append(seq)
    species.append(species_of(r.id))

n = len(names)
print(f"  {n} sequences")

anc_idx = next(i for i, s in enumerate(species) if s == "ancestor")
anc_seq = seqs[anc_idx]

short_names     = [short_name(nm) for nm in names]
name_to_short   = dict(zip(names, short_names))
name_to_species = dict(zip(names, species))


# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# PAIRWISE DISTANCES
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

print("Computing pairwise distances â€¦")

sim_mat   = np.zeros((n, n))
pair_rows = []

for i, j in combinations(range(n), 2):
    s = normsim(seqs[i], seqs[j])
    sim_mat[i, j] = sim_mat[j, i] = s
    pair_rows.append({
        "seqA": names[i], "seqB": names[j],
        "similarity": s,  "distance": 1 - s,
    })

np.fill_diagonal(sim_mat, 1.0)
dist_mat_full = 1.0 - sim_mat

d_min = dist_mat_full[dist_mat_full > 0].min()
d_max = dist_mat_full.max()
print(f"  Distance range: {d_min:.4f} â€“ {d_max:.4f}")


# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# HIERARCHICAL CLUSTERING
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

print("Hierarchical clustering â€¦")

import scipy.cluster.hierarchy as sch
from scipy.cluster.hierarchy import fcluster
from scipy.spatial.distance import squareform
import seaborn as sns
import pandas as pd
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

dist_cond      = squareform(dist_mat_full)
Z              = sch.linkage(dist_cond, method="ward")
CLUSTER_CUTOFF = 0.15
cluster_labels = fcluster(Z, CLUSTER_CUTOFF, criterion="distance")
cluster_of     = dict(zip(names, cluster_labels))
n_clusters     = len(set(cluster_labels))
print(f"  {n_clusters} clusters at cutoff {CLUSTER_CUTOFF}")


# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# HEATMAP â€” publication quality
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

print("Generating heatmap â€¦")

warnings.filterwarnings("ignore", category=UserWarning)

from scipy.cluster.hierarchy import dendrogram as _dend
_dinfo      = _dend(Z, no_plot=True)
row_order   = _dinfo["leaves"]

short_reord = [short_names[i] for i in row_order]
orig_reord  = [names[i]       for i in row_order]
sp_reord    = [species[i]     for i in row_order]
cl_reord    = [cluster_of[names[i]] for i in row_order]
mat_reord   = dist_mat_full[np.ix_(row_order, row_order)]

def _blocks(labels):
    blocks, i = [], 0
    while i < len(labels):
        j = i
        while j < len(labels) and labels[j] == labels[i]:
            j += 1
        blocks.append((i, j, labels[i]))
        i = j
    return blocks

blocks = _blocks(cl_reord)

CELL       = 0.20
side       = n * CELL
dend_col_h = 1.2
dend_row_w = 1.4
rlabel_w   = 2.2
clabel_h   = 2.0
cbar_gap   = 0.20
cbar_w_ax  = 0.18
cbar_h_ax  = side * 0.40
clnum_w    = 0.40

fig_w = clnum_w + dend_row_w + side + rlabel_w + cbar_gap + cbar_w_ax + 0.5
fig_h = dend_col_h + side + clabel_h + 0.3

fig = plt.figure(figsize=(fig_w, fig_h))

def _ax(l, b, w, h):
    return fig.add_axes([l/fig_w, b/fig_h, w/fig_w, h/fig_h])

b_heat = clabel_h
l_heat = clnum_w + dend_row_w

ax_col_dend = _ax(l_heat, b_heat + side,          side,       dend_col_h)
ax_row_dend = _ax(clnum_w, b_heat,                dend_row_w, side)
ax_heat     = _ax(l_heat,  b_heat,                side,       side)
ax_cbar     = _ax(l_heat + side + rlabel_w + cbar_gap,
                  b_heat + side * 0.30,
                  cbar_w_ax, cbar_h_ax)

sch.dendrogram(Z, orientation="top", ax=ax_col_dend, no_labels=True,
               color_threshold=0, link_color_func=lambda k: "black")
ax_col_dend.set_xlim(0, n * 10)
ax_col_dend.axis("off")

sch.dendrogram(Z, orientation="left", ax=ax_row_dend, no_labels=True,
               color_threshold=0, link_color_func=lambda k: "black")
ax_row_dend.set_ylim(0, n * 10)
ax_row_dend.invert_xaxis()
ax_row_dend.axis("off")

im = ax_heat.imshow(mat_reord, aspect="auto", origin="upper",
                    cmap="RdBu_r", vmin=0.0, vmax=d_max,
                    interpolation="nearest")

ax_heat.spines[:].set_visible(False)

ax_heat.set_yticks(np.arange(n))
ax_heat.set_yticklabels(short_reord, fontsize=5.5, fontfamily="Arial")
ax_heat.yaxis.set_label_position("right")
ax_heat.yaxis.tick_right()
ax_heat.tick_params(axis="y", length=0, pad=3)
for tick, sp in zip(ax_heat.yaxis.get_ticklabels(), sp_reord):
    tick.set_color(COLOR_OF[sp])

ax_heat.set_xticks(np.arange(n))
ax_heat.set_xticklabels(short_reord, fontsize=5.5, fontfamily="Arial",
                         rotation=90, ha="center", va="top")
ax_heat.xaxis.set_label_position("bottom")
ax_heat.xaxis.tick_bottom()
ax_heat.tick_params(axis="x", length=0, pad=2)
for tick, sp in zip(ax_heat.xaxis.get_ticklabels(), sp_reord):
    tick.set_color(COLOR_OF[sp])

ax_row_dend.set_xlim(ax_row_dend.get_xlim())
dend_xmax = ax_row_dend.get_xlim()[0]
dend_xmin = ax_row_dend.get_xlim()[1]

for start, end, cl_id in blocks:
    y_top = (n - start) * 10
    y_bot = (n - end)   * 10
    for y_line in [y_top, y_bot]:
        ax_row_dend.axhline(y_line, color="black", linewidth=0.8,
                            linestyle=(0, (4, 3)), zorder=5)
    fig.text(
        (clnum_w * 0.5) / fig_w,
        (b_heat + (n - (start + end) / 2) * CELL) / fig_h,
        str(cl_id), ha="center", va="center",
        fontsize=7, color="0.3",
    )

cb = fig.colorbar(im, cax=ax_cbar)
cb.set_label("Distance\n(BLOSUM62)", fontsize=7, labelpad=4)
cb.ax.tick_params(labelsize=6)
cb.set_ticks([0, d_max / 2, d_max])
cb.set_ticklabels(["0\n(identical)", f"{d_max/2:.2f}", f"{d_max:.2f}"])
cb.outline.set_linewidth(0.4)

leg = [mpatches.Patch(facecolor=COLOR_OF[sp], label=sp.capitalize())
       for sp in ("human", "mouse", "ancestor")]
fig.legend(handles=leg, title="Species", title_fontsize=7, fontsize=7,
           loc="upper right", bbox_to_anchor=(0.995, 0.995),
           framealpha=0.92, edgecolor="#cccccc")

fig.savefig(OUTPUT_DIR / "gpx6_heatmap_clustered.png",
            dpi=600, bbox_inches="tight")
plt.close()
print("Heatmap saved.")


# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# SSN â€” FIXED ANCHOR NODES AT TRIANGLE VERTICES
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

print("Building SSN â€¦")

from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse

LOW_T  = 0.80
HIGH_T = 0.92

G_all = nx.Graph()
for nm, sp in zip(names, species):
    G_all.add_node(nm, species=sp, short=short_name(nm), wt=is_wt(nm))

for row in pair_rows:
    sim = row["similarity"]
    if sim >= LOW_T:
        G_all.add_edge(
            row["seqA"], row["seqB"],
            weight = sim,
            strong = sim >= HIGH_T,
            intra  = name_to_species[row["seqA"]] == name_to_species[row["seqB"]],
        )

def find_wt(sp_filter):
    for nm, sp in zip(names, species):
        if sp == sp_filter and is_wt(nm):
            return nm
    for nm, sp in zip(names, species):
        if sp == sp_filter:
            return nm
    raise ValueError(f"No {sp_filter} WT node found")

wt_human = find_wt("human")
wt_mouse = find_wt("mouse")
wt_anc   = names[anc_idx]

print(f"  Anchors: human={wt_human}  mouse={wt_mouse}  ancestor={wt_anc}")

# SWAPPED: Mouse on LEFT (0,0), Human on RIGHT (10,0)
ANCHOR_POS = {
    wt_mouse : np.array([ 0.0,  0.0]),   # LEFT - Mouse
    wt_human : np.array([10.0,  0.0]),   # RIGHT - Human
    wt_anc   : np.array([ 5.0,  8.5]),   # TOP - Ancestor
}

SPECIES_CENTRE = {
    "mouse"   : ANCHOR_POS[wt_mouse],
    "human"   : ANCHOR_POS[wt_human],
    "ancestor": ANCHOR_POS[wt_anc],
}

rng = np.random.default_rng(42)
init_pos = {}

for nm in G_all.nodes():
    if nm in ANCHOR_POS:
        init_pos[nm] = ANCHOR_POS[nm].copy()
    else:
        sp  = name_to_species[nm]
        ctr = SPECIES_CENTRE[sp]
        angle = rng.uniform(0, 2 * np.pi)
        r     = rng.uniform(0.2, 1.0)
        init_pos[nm] = ctr + r * np.array([np.cos(angle), np.sin(angle)])

pos = nx.spring_layout(
    G_all,
    pos        = init_pos,
    fixed      = list(ANCHOR_POS.keys()),
    k          = 1.0,
    iterations = 600,
    seed       = 42,
    weight     = "weight",
)

mouse_nodes    = [nd for nd, d in G_all.nodes(data=True) if d["species"] == "mouse"]
human_nodes    = [nd for nd, d in G_all.nodes(data=True) if d["species"] == "human"]
ancestor_nodes = [nd for nd, d in G_all.nodes(data=True) if d["species"] == "ancestor"]

def nsizes(lst):
    return [520 if G_all.nodes[nd]["wt"] else 180 for nd in lst]

def elist(t):
    return [(u, v) for u, v, _ in t]

def ewidths(t, base, scale):
    return [base + scale * (d["weight"] - LOW_T) / (1 - LOW_T)
            for _, _, d in t]


# â”€â”€ Draw â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
fig, ax = plt.subplots(figsize=(14, 11))
ax.set_facecolor("white")

# SWAPPED: Triangle order now Mouse, Human, Ancestor
tri_pts = np.array([pos[wt_mouse], pos[wt_human], pos[wt_anc]])
tri_loop = np.vstack([tri_pts, tri_pts[0]])
ax.plot(tri_loop[:, 0], tri_loop[:, 1],
        color="black", linewidth=1.0, linestyle=(0, (8, 6)),
        alpha=0.20, zorder=0)

NODE_R = 0.06

for sp, nlist, col in [
        ("mouse",    mouse_nodes,    C_MOUSE),
        ("human",    human_nodes,    C_HUMAN),
        ("ancestor", ancestor_nodes, C_ANCESTOR)]:

    pts = np.array([pos[nd] for nd in nlist if nd in pos])
    if len(pts) == 0:
        continue

    xmin, ymin = pts.min(axis=0)
    xmax, ymax = pts.max(axis=0)
    PAD = max((xmax - xmin) * 0.22, (ymax - ymin) * 0.22, 0.15)
    cx  = (xmin + xmax) / 2
    cy  = (ymin + ymax) / 2
    w   = xmax - xmin + PAD * 2
    h   = ymax - ymin + PAD * 2

    ellipse = Ellipse(
        (cx, cy), width=w, height=h,
        facecolor = mcolors.to_rgba(col, alpha=0.08),
        edgecolor = col,
        linewidth = 2.2,
        linestyle = (0, (8, 4)),
        zorder    = 0,
    )
    ax.add_patch(ellipse)

    ax.text(cx, cy + h / 2 + PAD * 0.5, sp.capitalize(),
            ha="center", va="bottom", fontsize=12,
            color=col, zorder=8)

INTER_COLORS = {
    frozenset({"human", "mouse"})   : "#9B6B4E",
    frozenset({"human", "ancestor"}): "#C47AC0",
    frozenset({"mouse", "ancestor"}): "#5E8FC0",
}

all_inter = [(u, v, d) for u, v, d in G_all.edges(data=True)
             if not d["intra"]]
for u, v, d in all_inter:
    sp_pair = frozenset({name_to_species[u], name_to_species[v]})
    col_e   = INTER_COLORS.get(sp_pair, "#888888")
    lw      = 1.8 if d["strong"] else 0.8
    alpha   = 0.65 if d["strong"] else 0.28
    ls      = "solid" if d["strong"] else (0, (5, 3))
    nx.draw_networkx_edges(G_all, pos, edgelist=[(u, v)],
                           width=lw, alpha=alpha, edge_color=col_e,
                           style=ls, ax=ax)

weak_intra = [(u, v, d) for u, v, d in G_all.edges(data=True)
              if not d["strong"] and d["intra"]]
if weak_intra:
    for sp, col in [("human", C_HUMAN), ("mouse", C_MOUSE),
                    ("ancestor", C_ANCESTOR)]:
        sub = [(u, v, d) for u, v, d in weak_intra
               if name_to_species[u] == sp]
        if sub:
            nx.draw_networkx_edges(G_all, pos, edgelist=elist(sub),
                                   width=ewidths(sub, 0.5, 0.8),
                                   alpha=0.32, edge_color=col, ax=ax)

strong_intra = [(u, v, d) for u, v, d in G_all.edges(data=True)
                if d["strong"] and d["intra"]]
for sp, col in [("human", C_HUMAN), ("mouse", C_MOUSE),
                ("ancestor", C_ANCESTOR)]:
    sub = [(u, v, d) for u, v, d in strong_intra
           if name_to_species[u] == sp]
    if sub:
        nx.draw_networkx_edges(G_all, pos, edgelist=elist(sub),
                               width=ewidths(sub, 1.5, 2.5),
                               alpha=0.85, edge_color=col, ax=ax)

for nlist, col in [(mouse_nodes, C_MOUSE), (human_nodes, C_HUMAN),
                   (ancestor_nodes, C_ANCESTOR)]:
    nx.draw_networkx_nodes(G_all, pos, nodelist=nlist,
                           node_color=col, node_size=nsizes(nlist),
                           edgecolors="white", linewidths=1.4, ax=ax)

# â”€â”€ Node labels â€” overlap-aware placement â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
# Build list of (node, x, y, label, fontsize) sorted so WT nodes
# are placed first (they anchor the layout).
all_xy   = np.array(list(pos.values()))
xy_range = all_xy.max(axis=0) - all_xy.min(axis=0)

label_entries = []
for nd in G_all.nodes():
    x, y = pos[nd]
    lbl  = G_all.nodes[nd]["short"]
    fs   = 7.0 if G_all.nodes[nd]["wt"] else 5.5
    wt   = G_all.nodes[nd]["wt"]
    label_entries.append((nd, x, y, lbl, fs, wt))

# Sort: WT nodes first so they grab the cleanest positions
label_entries.sort(key=lambda e: (0 if e[5] else 1, e[0]))

# Candidate offset directions (angle steps) â€” try these in order
_offsets = [
    ( 0.0,  1.0),   # above
    ( 0.0, -1.0),   # below
    ( 1.0,  0.5),   # right-up
    (-1.0,  0.5),   # left-up
    ( 1.0, -0.5),   # right-down
    (-1.0, -0.5),   # left-down
    ( 1.0,  0.0),   # right
    (-1.0,  0.0),   # left
]

# Approximate character dimensions in data units
# We use a small fraction of the layout range as the text bounding box
CHAR_W = xy_range[0] * 0.018   # per character
LINE_H = xy_range[1] * 0.040   # line height

placed_boxes = []   # list of (x0, y0, x1, y1) already-placed label boxes

def _box(cx, cy, label, fs):
    """Return approximate bounding box centred at (cx, cy)."""
    w = len(label) * CHAR_W * (fs / 6.0)
    h = LINE_H * (fs / 6.0)
    return (cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2)

def _overlaps(b1, b2, margin=0.005):
    mx = xy_range[0] * margin
    my = xy_range[1] * margin
    return not (b1[2] + mx < b2[0] or b2[2] + mx < b1[0] or
                b1[3] + my < b2[1] or b2[3] + my < b1[1])

def _any_overlap(box):
    return any(_overlaps(box, pb) for pb in placed_boxes)

BASE_OFF = xy_range[1] * 0.038   # base offset distance

for nd, nx_, ny_, lbl, fs, wt in label_entries:
    best_xy   = None
    best_dist = -1.0

    for dx_unit, dy_unit in _offsets:
        scale = BASE_OFF * (1.4 if wt else 1.0)
        tx = nx_ + dx_unit * scale
        ty = ny_ + dy_unit * scale
        box = _box(tx, ty, lbl, fs)
        if not _any_overlap(box):
            best_xy   = (tx, ty)
            best_dist = scale
            break

    # Fallback: nudge further out until clear (up to 4Ã— base offset)
    if best_xy is None:
        for mult in [1.6, 2.2, 3.0, 4.0]:
            for dx_unit, dy_unit in _offsets:
                scale = BASE_OFF * mult
                tx = nx_ + dx_unit * scale
                ty = ny_ + dy_unit * scale
                box = _box(tx, ty, lbl, fs)
                if not _any_overlap(box):
                    best_xy = (tx, ty)
                    break
            if best_xy:
                break

    if best_xy is None:
        # Last resort: place directly above with no collision check
        best_xy = (nx_, ny_ + BASE_OFF * 1.2)

    tx, ty = best_xy
    final_box = _box(tx, ty, lbl, fs)
    placed_boxes.append(final_box)

    ax.text(
        tx, ty, lbl,
        fontsize   = fs,
        fontweight = "normal",
        ha         = "center",
        va         = "center",
        color      = "0.10",
        zorder     = 12,
        clip_on    = False,
        bbox       = dict(boxstyle="round,pad=0.1", facecolor="white",
                          edgecolor="none", alpha=0.72),
    )

# â”€â”€ Legend â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
legend_elements = [
    mpatches.Patch(facecolor=C_HUMAN,    edgecolor="white", label="Human"),
    mpatches.Patch(facecolor=C_MOUSE,    edgecolor="white", label="Mouse"),
    mpatches.Patch(facecolor=C_ANCESTOR, edgecolor="white", label="Ancestor"),
    Line2D([0], [0], color="#555", lw=2.0,
           label=f"Strong intra (sim â‰¥ {HIGH_T})"),
    Line2D([0], [0], color="#aaa", lw=0.8,
           label=f"Weak intra (sim â‰¥ {LOW_T})"),
    Line2D([0], [0], color="#9B6B4E", lw=1.5, label="Humanâ€“Mouse inter"),
    Line2D([0], [0], color="#C47AC0", lw=1.5, label="Humanâ€“Ancestor inter"),
    Line2D([0], [0], color="#5E8FC0", lw=1.5, label="Mouseâ€“Ancestor inter"),
    Line2D([0], [0], color="black",   lw=1.0, linestyle=(0, (8, 6)),
           alpha=0.4, label="Anchor triangle"),
]
ax.legend(handles=legend_elements, title="Legend", title_fontsize=8,
          fontsize=7, loc="lower left", framealpha=0.92,
          edgecolor="#cccccc", ncol=1)

ax.axis("off")
ax.autoscale_view()
ax.margins(0.16)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "gpx6_ssn_network.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("SSN saved.")


# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# DISTANCE TO ANCESTOR â€” barplot
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

print("Distance barplot â€¦")

anc_dist = []
for i, seq in enumerate(seqs):
    s = normsim(seq, anc_seq)
    anc_dist.append({
        "sequence"            : names[i],
        "short"               : short_names[i],
        "species"             : species[i],
        "distance_to_ancestor": 1 - s,
    })

by_sp = {
    sp: sorted(
        [x for x in anc_dist if x["species"] == sp],
        key=lambda x: x["distance_to_ancestor"],
    )
    for sp in ("human", "mouse", "ancestor")
}

rows_lbl, rows_val, rows_col = [], [], []
for sp in ("human", "mouse", "ancestor"):
    for x in by_sp[sp]:
        rows_lbl.append(x["short"])
        rows_val.append(x["distance_to_ancestor"])
        rows_col.append(COLOR_OF[sp])

fig_h = max(4, len(rows_lbl) * 0.22)
fig, ax = plt.subplots(figsize=(3.5, fig_h))

ypos = np.arange(len(rows_lbl))
ax.barh(ypos, rows_val, color=rows_col, edgecolor="white")
ax.set_yticks(ypos)
ax.set_yticklabels(rows_lbl, fontsize=6)
ax.set_xlabel("Distance to ancestor (BLOSUM62)")
ax.invert_yaxis()

legend_elements = [
    mpatches.Patch(facecolor=C_HUMAN,    label="Human"),
    mpatches.Patch(facecolor=C_MOUSE,    label="Mouse"),
    mpatches.Patch(facecolor=C_ANCESTOR, label="Ancestor"),
]
ax.legend(handles=legend_elements, fontsize=6, loc="lower right")

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "distance_to_ancestor_barplot.png",
            dpi=400, bbox_inches="tight")
plt.close()


# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# SUMMARY BARPLOT
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

print("Summary plot â€¦")

summary = []
for sp in ("human", "mouse", "ancestor"):
    vals = [x["distance_to_ancestor"] for x in anc_dist if x["species"] == sp]
    summary.append((sp, np.mean(vals), np.std(vals)))

fig, ax = plt.subplots(figsize=(3.3, 2.5))
labels, means, stds, cols = [], [], [], []
for sp, m, s in summary:
    labels.append(sp); means.append(m); stds.append(s); cols.append(COLOR_OF[sp])

ypos = np.arange(len(labels))
ax.barh(ypos, means, xerr=stds, color=cols)
ax.set_yticks(ypos)
ax.set_yticklabels(labels)
ax.set_xlabel("Mean distance to ancestor")
ax.invert_yaxis()

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "distance_summary.png",
            dpi=400, bbox_inches="tight")
plt.close()


# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
# TERNARY / BARYCENTRIC AFFINITY PLOT
# â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

print("Ternary affinity plot â€¦")

import csv as _csv
from matplotlib.patches import FancyArrowPatch
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

def find_anchor(kw, sp_filter):
    for i, (nm, sp) in enumerate(zip(names, species)):
        if sp == sp_filter and kw in nm.upper():
            return i
    for i, sp in enumerate(species):
        if sp == sp_filter:
            return i
    raise ValueError(f"No {sp_filter} sequence found")

idx_A = find_anchor("WT", "human")
idx_B = find_anchor("WT", "mouse")
idx_C = anc_idx

anchor_idxs  = {idx_A, idx_B, idx_C}
anchor_seqs  = [seqs[idx_A], seqs[idx_B], seqs[idx_C]]
anchor_names = [short_name(names[idx_A]),
                short_name(names[idx_B]),
                "Ancestor"]

V = np.array([
    [0.0,                   0.0],
    [1.0,                   0.0],
    [0.5, math.sqrt(3.0) / 2.0],
])

def barycentric_xy(sims):
    sims  = np.asarray(sims, dtype=float)
    total = sims.sum()
    if total <= 0:
        return V.mean(axis=0)
    lam = sims / total
    return lam @ V

tri_rows = []
for i, (nm, seq, sp) in enumerate(zip(names, seqs, species)):
    sims = np.array([normsim(seq, ref) for ref in anchor_seqs])
    xy   = barycentric_xy(sims)
    tri_rows.append({
        "name"         : nm,
        "short"        : short_name(nm),
        "species"      : sp,
        "is_anchor"    : i in anchor_idxs,
        "sim_human_wt" : float(sims[0]),
        "sim_mouse_wt" : float(sims[1]),
        "sim_ancestor" : float(sims[2]),
        "aff_human_wt" : float(sims[0] / sims.sum()),
        "aff_mouse_wt" : float(sims[1] / sims.sum()),
        "aff_ancestor" : float(sims[2] / sims.sum()),
        "x"            : float(xy[0]),
        "y"            : float(xy[1]),
    })

tri_csv = OUTPUT_DIR / "gpx6_ternary_coordinates.csv"
with tri_csv.open("w", newline="") as f:
    w = _csv.DictWriter(f, fieldnames=list(tri_rows[0].keys()))
    w.writeheader(); w.writerows(tri_rows)

query_rows  = [r for r in tri_rows if not r["is_anchor"]]
anchor_rows = [r for r in tri_rows if r["is_anchor"]]

fig, ax = plt.subplots(figsize=(8, 7.5))

tri_loop = np.vstack([V, V[0]])
ax.plot(tri_loop[:, 0], tri_loop[:, 1], color="black",
        linewidth=1.8, zorder=1)

for frac in (1/3, 2/3):
    for a, b, c in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
        p1 = frac * V[a] + (1 - frac) * V[b]
        p2 = frac * V[a] + (1 - frac) * V[c]
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]],
                color="0.82", linewidth=0.5, linestyle="--", zorder=0)

offsets  = [(-0.09, -0.07), (0.09, -0.07), (0.0, 0.06)]
v_colors = [C_HUMAN, C_MOUSE, C_ANCESTOR]
for label, v, off, col in zip(anchor_names, V, offsets, v_colors):
    ax.text(v[0] + off[0], v[1] + off[1], label,
            ha="center", va="center", fontsize=11,
            color=col, zorder=8)

for v, col in zip(V, v_colors):
    ax.scatter(*v, s=280, color=col, marker="D",
               edgecolors="white", linewidths=1.2, zorder=7)

for sp, col in [("human", C_HUMAN), ("mouse", C_MOUSE),
                ("ancestor", C_ANCESTOR)]:
    sub = [r for r in query_rows if r["species"] == sp]
    if not sub:
        continue
    xs    = np.array([r["x"] for r in sub])
    ys    = np.array([r["y"] for r in sub])
    sizes = [140 if "WT" in r["name"].upper() else 70 for r in sub]
    ax.scatter(xs, ys, color=col, s=sizes,
               edgecolors="white", linewidths=0.6,
               zorder=5, label=sp.capitalize(), alpha=0.90)

# â”€â”€ Ternary labels â€” overlap-aware placement â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
tern_range = np.array([1.44, 1.23])   # approx data range of the ternary plot
TCHAR_W  = tern_range[0] * 0.018
TLINE_H  = tern_range[1] * 0.040

t_placed_boxes = []

def _tbox(cx, cy, label, fs):
    w = len(label) * TCHAR_W * (fs / 5.5)
    h = TLINE_H * (fs / 5.5)
    return (cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2)

def _tany_overlap(box):
    return any(_overlaps(box, pb) for pb in t_placed_boxes)

T_BASE_OFF = tern_range[1] * 0.032
centroid   = V.mean(axis=0)

# Sort: WT first
query_sorted = sorted(query_rows,
                      key=lambda r: (0 if "WT" in r["name"].upper() else 1,
                                     r["name"]))

for r in query_sorted:
    rx, ry = r["x"], r["y"]
    lbl    = r["short"]
    fs     = 4.8

    # Default direction: away from centroid
    default_dx = rx - centroid[0]
    default_dy = ry - centroid[1]
    norm = math.sqrt(default_dx**2 + default_dy**2) or 1.0
    primary = [(default_dx / norm, default_dy / norm)]
    all_dirs = primary + [d for d in _offsets
                          if (d[0], d[1]) != (round(default_dx / norm, 1),
                                              round(default_dy / norm, 1))]

    best_xy = None
    for dx_unit, dy_unit in all_dirs:
        for mult in [1.0, 1.5, 2.2, 3.0]:
            scale = T_BASE_OFF * mult
            tx = rx + dx_unit * scale
            ty = ry + dy_unit * scale
            box = _tbox(tx, ty, lbl, fs)
            if not _tany_overlap(box):
                best_xy = (tx, ty)
                break
        if best_xy:
            break

    if best_xy is None:
        best_xy = (rx, ry + T_BASE_OFF * 1.2)

    tx, ty = best_xy
    t_placed_boxes.append(_tbox(tx, ty, lbl, fs))

    ax.text(tx, ty, lbl,
            fontsize=fs, ha="center", va="center",
            color="0.25", zorder=6)

edge_label_kwargs = dict(fontsize=8, color="0.4", ha="center", va="center")
ax.text(0.25, -0.10, "â† more Human WT â†’",  rotation=0,   **edge_label_kwargs)
ax.text(0.75, -0.10, "â† more Mouse WT â†’",  rotation=0,   **edge_label_kwargs)
ax.text(0.17,  0.44, "â† more Ancestor",    rotation=60,  **edge_label_kwargs)
ax.text(0.83,  0.44, "more Ancestor â†’",    rotation=-60, **edge_label_kwargs)

legend_elements = [
    mpatches.Patch(facecolor=C_HUMAN,    label="Human variants"),
    mpatches.Patch(facecolor=C_MOUSE,    label="Mouse variants"),
    mpatches.Patch(facecolor=C_ANCESTOR, label="Ancestor"),
    plt.scatter([], [], s=140, c="grey", edgecolors="white",
                linewidths=0.6, label="WT (larger marker)"),
]
ax.legend(handles=legend_elements, fontsize=7, loc="upper left",
          framealpha=0.92, edgecolor="#cccccc",
          bbox_to_anchor=(0.0, 1.0))

ax.set_aspect("equal")
ax.set_xlim(-0.22, 1.22)
ax.set_ylim(-0.18, 1.05)
ax.axis("off")

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "gpx6_ternary.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("Ternary plot saved.")

print("\nDone.")

