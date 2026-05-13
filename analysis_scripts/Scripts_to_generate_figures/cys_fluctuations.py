"""
plot_U49C_pathway.py

Reads: human_HBONDS.csv
Plots U49C mutation across all available levels
Shows Î”G* fluctuation with mutation level

JCIM compliant
- ALL text forced black
- Arial 8 pt
- inward ticks (bottom and left only)
- Professional styling
- Horizontal layout
"""

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.text as mtext

# ==========================================================
# JCIM STYLE
# ==========================================================
ACS_PATH = r"./analysis_scripts/Scripts_to_generate_figures/Figures"
if ACS_PATH not in sys.path:
    sys.path.append(ACS_PATH)

from acsfonts import set_jcim_style
set_jcim_style()

# ==========================================================
# GLOBAL FONT / TEXT SETTINGS (ALL TEXT BLACK)
# ==========================================================
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "text.color": "black",
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "axes.edgecolor": "black"
})

# ==========================================================
# FILES
# ==========================================================
INPUT_FILE = r"./analysis_scripts\human_HBONDS.csv"
OUTPUT_FILE = r"./analysis_scripts/Scripts_to_generate_figures/Figures/U49C_fluctuation.png"

os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)

# ==========================================================
# COLORS
# ==========================================================
C_U49C = "#0B5CAD"

# ==========================================================
# LOAD DATA
# ==========================================================
df = pd.read_csv(INPUT_FILE)
df = df[df["mutation"] == "U49C"].copy()

df["level_num"] = df["level"].str.replace("Level", "", regex=False).astype(int)
df = df.sort_values("level_num")

# ==========================================================
# FIGURE (WIDE AND SHORT - HORIZONTAL)
# ==========================================================
fig, ax = plt.subplots(figsize=(6.5, 2.8))

# ----------------------------------------------------------
# Main U49C pathway
# ----------------------------------------------------------
ax.errorbar(
    df["level_num"],
    df["dg_star"],
    yerr=df["dg_star_error"],
    color=C_U49C,
    linewidth=1.5,
    marker="o",
    markersize=4.5,
    markerfacecolor="white",
    markeredgewidth=1.0,
    capsize=2.0,
    capthick=0.8,
    elinewidth=0.8,
    zorder=5
)

# ----------------------------------------------------------
# AXES FORMATTING
# ----------------------------------------------------------
ax.set_xlabel("Mutation level", fontsize=8)
ax.set_ylabel(r"$\Delta G^{\ddagger}$ (kcal mol$^{-1}$)", fontsize=8)

# Set x-ticks to show ALL levels (0 through 20)
all_levels = list(range(df["level_num"].min(), df["level_num"].max() + 1))
ax.set_xticks(all_levels)
ax.set_xticklabels(all_levels, fontsize=6)  # Rotate for readability

# Y-axis ticks - bottom and left ONLY (no top, no right)
ax.tick_params(
    axis="both",
    direction="in",
    length=3.5,
    width=0.6,
    colors="black",
    top=False,
    right=False
)

# Set limits with padding
ax.set_xlim(df["level_num"].min() - 0.5, df["level_num"].max() + 0.5)

# Light grid
ax.grid(True, axis="y", linestyle="-", linewidth=0.4, alpha=0.25, color="lightgray")
ax.set_axisbelow(True)

# Spines - keep only bottom and left
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["bottom"].set_linewidth(0.5)
ax.spines["left"].set_linewidth(0.5)

# ----------------------------------------------------------
# FORCE ALL TEXT BLACK
# ----------------------------------------------------------
ax.xaxis.label.set_color("black")
ax.yaxis.label.set_color("black")
ax.title.set_color("black")

for t in ax.get_xticklabels():
    t.set_color("black")

for t in ax.get_yticklabels():
    t.set_color("black")

for obj in fig.findobj(match=mtext.Text):
    obj.set_color("black")

plt.tight_layout(pad=0.2)

# ==========================================================
# SAVE
# ==========================================================
plt.savefig(
    OUTPUT_FILE,
    dpi=600,
    bbox_inches="tight",
    facecolor="white"
)

print(f"Saved: {OUTPUT_FILE}")

