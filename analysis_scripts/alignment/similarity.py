import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

# =========================================================
# USER PATHS (EDIT ONLY IF YOU MOVE FILES)
# =========================================================

BASE_DIR = r"D:\PhD_Thesis\GPX6\analysis\alignment"

DATA_FILE = os.path.join(BASE_DIR, "MEGA-result1.csv")
OUTPUT_FILE = os.path.join(
    BASE_DIR, "GPX6_blossom_distance_heatmap_with_ancestor.png"
)

# =========================================================
# SAFETY CHECK
# =========================================================

print("Looking for distance matrix at:")
print(DATA_FILE)

if not os.path.exists(DATA_FILE):
    raise FileNotFoundError(
        f"\nERROR: Distance matrix not found.\n"
        f"Expected location:\n{DATA_FILE}\n\n"
        f"Fix the filename or path and rerun."
    )

# =========================================================
# LOAD DISTANCE MATRIX
# =========================================================

df = pd.read_csv(DATA_FILE, index_col=0)

# Replace missing values (especially diagonal)
df = df.fillna(0.0)

# Enforce symmetry (important for clustering)
df = (df + df.T) / 2

# =========================================================
# HIERARCHICAL CLUSTERING
# =========================================================

# Convert square matrix → condensed distance matrix
condensed_dist = squareform(df.values)

# Ward linkage (good for blossom-like structure)
Z = linkage(condensed_dist, method="ward")

# =========================================================
# PLOT CLUSTERED HEATMAP
# =========================================================

sns.set(style="white")

g = sns.clustermap(
    df,
    row_linkage=Z,
    col_linkage=Z,
    cmap="magma_r",
    linewidths=0.3,
    figsize=(20, 20),
    cbar_kws={"label": "Pairwise structural distance"},
)

# =========================================================
# HIGHLIGHT ANCESTRAL NODE
# =========================================================

ANCESTOR_LABEL = "node_25"

# Highlight tick labels
for label in g.ax_heatmap.get_yticklabels():
    if ANCESTOR_LABEL in label.get_text():
        label.set_color("cyan")
        label.set_fontweight("bold")

for label in g.ax_heatmap.get_xticklabels():
    if ANCESTOR_LABEL in label.get_text():
        label.set_color("cyan")
        label.set_fontweight("bold")

# Draw box around ancestor cell (diagonal)
labels = [lbl.get_text() for lbl in g.ax_heatmap.get_yticklabels()]

if ANCESTOR_LABEL in labels:
    idx = labels.index(ANCESTOR_LABEL)
    g.ax_heatmap.add_patch(
        plt.Rectangle(
            (idx, idx), 1, 1,
            fill=False, edgecolor="cyan", linewidth=3
        )
    )

# =========================================================
# TITLE
# =========================================================

g.fig.suptitle(
    "Blossom Distance Matrix of GPX6 Variants\n"
    "(Human, Mouse, and Ancestral node_25)",
    fontsize=18,
    y=1.03,
)

# =========================================================
# SAVE FIGURE
# =========================================================

plt.savefig(OUTPUT_FILE, dpi=600, bbox_inches="tight")
plt.close()

print("\nSUCCESS ✅")
print(f"Heatmap saved at:\n{OUTPUT_FILE}")
