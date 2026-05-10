"""
gpx6_all_pathways_fixed.py
FIXES APPLIED:
  - Panels A/B gap eliminated (hspace=0.25, tight margins)
  - All fonts increased to 10 pt (was 8 pt)
  - Axis labels, tick labels, legend, annotations all larger
  - Figure margins tightened (left/right/top/bottom)
  - Mutation labels larger (9 pt, was 7 pt)
  - Legend font 9 pt (was 7 pt)
  - Removed acsfonts dependency — pure matplotlib/Arial
  - All text forced black
  - Inward ticks retained
  - Colored reference lines + legend retained
"""

import sys
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.text as mtext

# ==========================================================
# JCIM STYLE — no external dependency
# ==========================================================
plt.rcParams.update({
    "font.family":      "Arial",
    "font.size":        10,
    "axes.labelsize":   10,
    "axes.titlesize":   10,
    "xtick.labelsize":  10,
    "ytick.labelsize":  10,
    "legend.fontsize":  9,
    "text.color":       "black",
    "axes.labelcolor":  "black",
    "xtick.color":      "black",
    "ytick.color":      "black",
    "axes.edgecolor":   "black",
    "figure.dpi":       600,
    "savefig.dpi":      600,
    "pdf.fonttype":     42,
    "ps.fonttype":      42,
})

# ==========================================================
# FILES — update these paths for your machine
# ==========================================================
HUMAN_FILE = r"D:\PhD_Thesis\analysis\dipole\human_HBONDS.csv"
MOUSE_FILE = r"D:\PhD_Thesis\analysis\dipole\mouse_HBONDS.csv"

HUMAN_REFS = ["humansec", "U49C"]
MOUSE_REFS = ["mousecys", "C49U"]

OUTPUT_FILE = r"D:\PhD_Thesis\analysis\NetworkX\gpx6_all_pathways.png"
os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)

# ==========================================================
# COLORS
# ==========================================================
C_GREEDY = "#0B5CAD"

REF_STYLES = {
    "hsec": {"ls": "--",        "lw": 1.6, "name": "Human Sec WT",        "color": "#1f77b4"},
    "hcys": {"ls": "-.",        "lw": 1.6, "name": "Human Cys WT",        "color": "#ff7f0e"},
    "msec": {"ls": ":",         "lw": 1.6, "name": "Mouse Sec WT (C49U)", "color": "#2ca02c"},
    "mcys": {"ls": (0, (5, 1)), "lw": 1.6, "name": "Mouse Cys WT",        "color": "#d62728"},
}

# ==========================================================
# LOAD DATA
# ==========================================================
def load(path, refs):
    df  = pd.read_csv(path)
    df2 = df[~df["mutation"].isin(refs)].copy()
    df2["level_num"] = df2["level"].str.replace("Level", "").astype(int)
    levels = sorted(df2["level_num"].unique())

    mut_paths = {}
    for mut in df2["mutation"].unique():
        rows = df2[df2["mutation"] == mut].sort_values("level_num")
        mut_paths[mut] = {
            "levels": rows["level_num"].values,
            "dg":     rows["dg_star"].values,
            "err":    rows["dg_star_error"].values,
        }

    greedy_levels, greedy_muts, greedy_dg, greedy_err = [], [], [], []
    for lv in levels:
        lv_data = df2[df2["level_num"] == lv]
        best    = lv_data.loc[lv_data["dg_star"].idxmin()]
        greedy_levels.append(lv)
        greedy_muts.append(best["mutation"])
        greedy_dg.append(best["dg_star"])
        greedy_err.append(best["dg_star_error"])

    return {
        "levels":        np.array(levels),
        "mut_paths":     mut_paths,
        "greedy_levels": np.array(greedy_levels),
        "greedy_muts":   greedy_muts,
        "greedy_dg":     np.array(greedy_dg),
        "greedy_err":    np.array(greedy_err),
        "all_dg":        df2["dg_star"].values,
    }

# ==========================================================
# READ FILES
# ==========================================================
h = load(HUMAN_FILE, HUMAN_REFS)
m = load(MOUSE_FILE, MOUSE_REFS)

h_raw = pd.read_csv(HUMAN_FILE)
m_raw = pd.read_csv(MOUSE_FILE)

def get_ref(df, name):
    return float(df.loc[df["mutation"] == name, "dg_star"].values[0])

hsec_dg = get_ref(h_raw, "humansec")
hcys_dg = get_ref(h_raw, "U49C")
mcys_dg = get_ref(m_raw, "mousecys")
msec_dg = get_ref(m_raw, "C49U")

REF_ORDER = [
    ("hsec", hsec_dg),
    ("hcys", hcys_dg),
    ("msec", msec_dg),
    ("mcys", mcys_dg),
]

# ==========================================================
# GLOBAL SCALE
# ==========================================================
global_vmin = min(min(h["all_dg"]), min(m["all_dg"]))
global_vmax = max(max(h["all_dg"]), max(m["all_dg"]))

# ==========================================================
# FIGURE — tighter layout, panels closer together
# ==========================================================
fig = plt.figure(figsize=(12, 9))

gs = fig.add_gridspec(
    2, 1,
    hspace  = 0.25,   # ← was 0.58; much tighter gap between A and B
    left    = 0.09,
    right   = 0.80,
    top     = 0.95,
    bottom  = 0.07,
)

ax_top = fig.add_subplot(gs[0])
ax_bot = fig.add_subplot(gs[1])

PANELS = [
    (ax_top, h, "a"),
    (ax_bot, m, "b"),
]

# ==========================================================
# PLOT
# ==========================================================
for ax, d, panel in PANELS:

    levels = d["levels"]
    cmap   = plt.colormaps["Blues"]
    norm   = mcolors.Normalize(vmin=global_vmin - 2, vmax=global_vmax + 2)

    # Background trajectories
    for mut, path in d["mut_paths"].items():
        mean_dg = path["dg"].mean()
        ax.plot(
            path["levels"], path["dg"],
            lw=0.8, alpha=0.30,
            color=cmap(norm(mean_dg)), zorder=1,
        )

    viability_dg = mcys_dg if ax == ax_top else hcys_dg

    # Reference lines
    ref_handles = []
    for key, val in REF_ORDER:
        s    = REF_STYLES[key]
        line = ax.axhline(
            val, color=s["color"], lw=s["lw"], ls=s["ls"],
            alpha=0.85, label=f"{s['name']} ({val:.2f})",
        )
        ref_handles.append(line)

    # Viability zone
    ax.axhspan(viability_dg, global_vmax + 10,
               color="lightgray", alpha=0.25)

    # Greedy path
    ax.errorbar(
        d["greedy_levels"], d["greedy_dg"],
        yerr=d["greedy_err"],
        color=C_GREEDY, lw=2.3,
        marker="o", markersize=5.2,
        markerfacecolor="white", markeredgewidth=1.4,
        capsize=3.0, zorder=5,
    )

    # Mutation labels — larger font
    ytop = viability_dg + 1.6
    ybot = global_vmin  - 1.2

    for i, (lv, mut, dg) in enumerate(
        zip(d["greedy_levels"], d["greedy_muts"], d["greedy_dg"])
    ):
        ytxt = ytop if i % 2 == 0 else ybot
        va   = "bottom" if i % 2 == 0 else "top"
        ax.annotate(
            mut,
            xy=(lv, dg), xytext=(lv, ytxt),
            fontsize=9, color="black",      # ← was 7
            ha="center", va=va,
            arrowprops=dict(arrowstyle="-", color="black", lw=0.6),
        )

    # Panel label
    ax.text(
        -0.06, 1.02, panel,
        transform=ax.transAxes,
        fontsize=10, fontweight="bold", color="black",   # ← was 12, now bold
    )

    # Axes formatting
    ax.set_xlim(-0.5, levels[-1] + 0.5)
    ax.set_ylim(global_vmin - 2.2, global_vmax + 3)
    ax.set_xticks(levels)
    ax.set_xticklabels(levels, fontsize=10)

    ax.set_xlabel("Mutation level", fontsize=10)
    ax.set_ylabel(r"Activation free energy $\Delta G^*$ (kcal mol$^{-1}$)", fontsize=10)

    ax.tick_params(
        axis="both", direction="in",
        length=4, width=0.8, colors="black",
    )

    # Legend
    leg = ax.legend(
        handles=ref_handles,
        loc="upper left",
        fontsize=9,          # ← was 7
        frameon=False,
    )
    for text in leg.get_texts():
        text.set_color("black")

# ==========================================================
# FORCE ALL TEXT BLACK
# ==========================================================
for ax in fig.axes:
    ax.xaxis.label.set_color("black")
    ax.yaxis.label.set_color("black")
    ax.title.set_color("black")
    for t in ax.get_xticklabels() + ax.get_yticklabels():
        t.set_color("black")
    for txt in ax.texts:
        txt.set_color("black")

for obj in fig.findobj(match=mtext.Text):
    obj.set_color("black")

# ==========================================================
# SAVE
# ==========================================================
for ext in ("png", "pdf"):
    out = OUTPUT_FILE.replace(".png", f".{ext}")
    plt.savefig(out, dpi=600, bbox_inches="tight",
                pad_inches=0.02, facecolor="white")
    print("Saved:", out)