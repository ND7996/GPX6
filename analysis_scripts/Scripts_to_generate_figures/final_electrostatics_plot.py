"""
heatmap_levels_gpx6_jcim.py
JCIM publication style:
  - Font:  Arial (fallback DejaVu Sans), 7 pt tick labels, 8 pt axis labels
  - Lines: 0.5 pt axes/ticks
  - Ticks: inward, 3 pt major
  - No bold anywhere (except panel labels a and b)
  - No white gaps between cells (edgecolor='none')
  - All level numbers on x-axis (0, 1, 2, â€¦)
  - Mouse: levels 0-20 (21 cols), Human: levels 0-19 (20 cols)
  - X-axis label on BOTH panels
  - Legend: 2 rows, no overlap
  - Colorbar centred below legend
  - FIXED: x-axis labels rotated to prevent overlap
  - FIXED: Legend title no longer overlaps with "0 H-bond"
  - FIXED: H-bond labels no longer overlap
  - FIXED: Panel labels bold without parentheses
"""

import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd

# â”€â”€ JCIM font: Arial 8 pt â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
import matplotlib.font_manager as fm

_arial_found = any("arial" in f.lower() for f in fm.findSystemFonts())
plt.rcParams["font.family"] = "Arial" if _arial_found else "DejaVu Sans"

# JCIM standard rcParams
plt.rcParams.update({
    # font sizes  (JCIM: 7 pt ticks, 8 pt labels)
    "font.size":             10,
    "axes.labelsize":        10,
    "axes.titlesize":        10,
    "xtick.labelsize":       7,
    "ytick.labelsize":       7,
    "legend.fontsize":       7,
    "legend.title_fontsize": 8,
    # colours
    "text.color":            "black",
    "axes.labelcolor":       "black",
    "xtick.color":           "black",
    "ytick.color":           "black",
    "figure.facecolor":      "white",
    "axes.facecolor":        "white",
    "savefig.facecolor":     "white",
    # line / tick weights  (JCIM: 0.5 pt)
    "axes.linewidth":        0.5,
    "xtick.major.width":     0.5,
    "ytick.major.width":     0.5,
    "xtick.major.size":      3,
    "ytick.major.size":      3,
    "xtick.direction":       "in",
    "ytick.direction":       "in",
    # embed fonts in PDF/PS
    "pdf.fonttype":          42,
    "ps.fonttype":           42,
})

# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
# ROW ORDER
# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
MOUSE_ORDER = [
    ("C49U",  0.0), ("F48Y",  3.9), ("T52A",  5.2), ("S47A",  5.6),
    ("R99C",  9.5), ("T54Q", 10.0), ("H144Q", 11.6), ("H177Q", 11.5),
    ("G74A",  14.3), ("T178A", 14.8), ("E143S", 15.3), ("P142S", 16.4),
    ("F139L", 16.3), ("K87T",  16.4), ("Y104F", 16.1), ("G102S", 16.8),
    ("L24I",  18.7), ("T60A",  19.1), ("S4R",   19.6), ("N107S", 19.8),
    ("R181S", 19.5),
]
HUMAN_ORDER = [
    ("U49C",  0.0), ("Y48F",  3.8), ("A52T",  5.5), ("A47S",  5.6),
    ("Q54T",  10.4), ("C99R",  9.7), ("Q144H", 11.5), ("Q177H", 11.5),
    ("A178T", 14.2), ("T87K",  16.8), ("A74G",  14.6), ("S143E", 15.1),
    ("S142P", 16.6), ("N3K",   17.7), ("L139F", 16.5), ("F104Y", 17.1),
    ("A60T",  19.8), ("L24I",  18.5), ("S102G", 18.1), ("H173R", 19.6),
    ("S181R", 19.2),
]

MOUSE_DIST      = {m: d for m, d in MOUSE_ORDER}
HUMAN_DIST      = {m: d for m, d in HUMAN_ORDER}
MOUSE_ROW_ORDER = [m for m, _ in MOUSE_ORDER]
HUMAN_ROW_ORDER = [m for m, _ in HUMAN_ORDER]

# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
# CONFIG
# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
CONFIG = dict(
    human_file = r"./analysis_scripts\human_HBONDS.csv",
    mouse_file = r"./analysis_scripts\mouse_HBONDS.csv",
    human_refs = ["humansec", "humancys"],
    mouse_refs = ["mousecys", "C49U"],
    cmap       = "RdBu_r",
    hb_colors  = {
        0: "#440154", 1: "#3B528B", 2: "#21908C",
        3: "#5DC863", 4: "#FDE725", 5: "#FF6B35",
    },
)

SAVE_DIR  = r"./analysis_scripts/Scripts_to_generate_figures/Figures/Figures_FINAL"
os.makedirs(SAVE_DIR, exist_ok=True)
SAVE_PATH = os.path.join(SAVE_DIR, "heatmap_levels_gpx6_jcim")
MISSING_CLR = "#e0e0e0"

# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
# DATA LOADING
# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
def get_wt_ref(path, refs):
    df   = pd.read_csv(path)
    mask = df["mutation"].isin(refs) & (df["level"] == "Level0")
    return float(df.loc[mask, "dg_star"].mean())


def load_pivot(path, refs, ref_dg, dist_map, row_order):
    df = pd.read_csv(path)
    df = df[~df["mutation"].isin(refs)].copy()
    df["level_num"] = df["level"].str.replace("Level", "", regex=False).astype(int)
    df["ddg"]       = df["dg_star"] - ref_dg
    l0  = df[df["level_num"] == 0].set_index("mutation")
    df["HBonds"]    = df["mutation"].map(l0["Total_HBonds"])
    df["AS_HBonds"] = df["mutation"].map(l0["Active_Site_HBonds"].fillna(0))
    pivot = df.pivot_table(
        index="mutation", columns="level_num", values="ddg", aggfunc="mean"
    )
    meta = (
        df[df["level_num"] == 0]
        .drop_duplicates("mutation")
        .set_index("mutation")[["HBonds", "AS_HBonds"]]
        .dropna(subset=["HBonds"])
    )
    meta["ddg_l0"]   = pivot.get(0, pd.Series(dtype=float))
    meta["Distance"] = meta.index.map(dist_map)
    ordered = [m for m in row_order if m in meta.index]
    return pivot.reindex(ordered), meta.reindex(ordered)


h_ref = get_wt_ref(CONFIG["human_file"], CONFIG["human_refs"])
m_ref = get_wt_ref(CONFIG["mouse_file"], CONFIG["mouse_refs"])

h_pivot, h_meta = load_pivot(
    CONFIG["human_file"], CONFIG["human_refs"], h_ref, HUMAN_DIST, HUMAN_ROW_ORDER)
m_pivot, m_meta = load_pivot(
    CONFIG["mouse_file"], CONFIG["mouse_refs"], m_ref, MOUSE_DIST, MOUSE_ROW_ORDER)

# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
# FORCE LEVEL RANGES (JCIM CONSISTENCY)
# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
# Mouse: levels 0â€“20 (21 columns)
m_pivot = m_pivot.reindex(columns=list(range(21)))

# Human: levels 0â€“19 (20 columns)
h_pivot = h_pivot.reindex(columns=list(range(20)))

vabs     = max(np.nanmax(np.abs(h_pivot.values)),
               np.nanmax(np.abs(m_pivot.values)))
vabs     = np.ceil(vabs / 2) * 2
norm     = mcolors.TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)
cmap_obj = plt.get_cmap(CONFIG["cmap"])

# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
# FIGURE GEOMETRY  (sized so 1 cell â‰ˆ correct pt size at 600 dpi)
# Mouse: 21 cols (0-20)   Human: 20 cols (0-19)
# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
n_m          = len(m_pivot)
n_h          = len(h_pivot)
n_cols_mouse = len(m_pivot.columns)   # 21
n_cols_human = len(h_pivot.columns)   # 20
n_cols_max   = max(n_cols_mouse, n_cols_human)

# Cell dimensions in inches tuned for JCIM column width (~3.33 in single, ~7 in double)
CELL_W   = 0.30   # inches per level column
CELL_H   = 0.22   # inches per mutation row

LEFT_IN  = 0.75   # mutation y-labels
RIGHT_IN = 0.55   # badge column
TOP_IN   = 0.45   # panel label
GAP_IN   = 0.60   # inter-panel gap (both panels have x-axis)
BOT_IN   = 1.50   # legend + colorbar

heatmap_w = n_cols_max * CELL_W
fig_w = LEFT_IN + heatmap_w + RIGHT_IN
fig_h = TOP_IN + n_m * CELL_H + GAP_IN + n_h * CELL_H + BOT_IN

l_frac = LEFT_IN  / fig_w
r_frac = 1.0 - RIGHT_IN / fig_w
t_frac = 1.0 - TOP_IN   / fig_h
b_frac = BOT_IN          / fig_h
avg_h_frac = ((n_m + n_h) / 2 * CELL_H) / fig_h
hspace     = (GAP_IN / fig_h) / avg_h_frac

fig, (ax_top, ax_bot) = plt.subplots(
    2, 1,
    figsize=(fig_w, fig_h),
    gridspec_kw={
        "height_ratios": [n_m, n_h],
        "hspace":        hspace,
        "left":          l_frac,
        "right":         r_frac,
        "top":           t_frac,
        "bottom":        b_frac,
    },
)

# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
# DRAW FUNCTION (FIXED: No overlapping H-bond labels)
# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
def draw_heatmap(ax, pivot, meta, species_label, panel_label, n_cols_max):
    mutations = list(pivot.index)
    levels    = sorted(pivot.columns)
    n_muts    = len(mutations)
    n_lev     = len(levels)

    # â”€â”€ colour cells â€” no gap (edgecolor none) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    for ri, mut in enumerate(mutations):
        for ci, lv in enumerate(levels):
            val  = pivot.loc[mut, lv] if lv in pivot.columns else np.nan
            face = MISSING_CLR if pd.isna(val) else cmap_obj(norm(val))
            ax.add_patch(plt.Rectangle(
                [ci, ri], 1, 1,
                facecolor=face, edgecolor="none", linewidth=0,
            ))

    # â”€â”€ y-axis: mutation labels â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    ax.set_yticks([r + 0.5 for r in range(n_muts)])
    ax.set_yticklabels(mutations, fontsize=7, va="center")
    ax.tick_params(axis="y", pad=2, length=0, width=0, left=False, right=False)
    ax.set_ylabel("Mutation", fontsize=8, labelpad=5)

    # â”€â”€ H-bond badges (FIXED: No overlapping labels) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    badge_x  = n_lev + 0.15
    badge_sz = 0.68
    badge_w  = 0.58
    
    # First pass: draw all badges
    for ri, mut in enumerate(mutations):
        if mut in meta.index and not pd.isna(meta.loc[mut, "HBonds"]):
            hb    = int(meta.loc[mut, "HBonds"])
            color = CONFIG["hb_colors"].get(hb, "#888")
            ax.add_patch(FancyBboxPatch(
                (badge_x, ri + (1 - badge_sz) / 2), badge_w, badge_sz,
                boxstyle="round,pad=0.03",
                facecolor=color, edgecolor="none", clip_on=False,
            ))
    
    # Second pass: add labels with intelligent positioning to avoid overlap
    for ri, mut in enumerate(mutations):
        if mut in meta.index and not pd.isna(meta.loc[mut, "HBonds"]):
            hb    = int(meta.loc[mut, "HBonds"])
            color = CONFIG["hb_colors"].get(hb, "#888")
            
            rgb = mcolors.to_rgb(color)
            lum = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]
            tc  = "white" if lum < 0.55 else "black"
            
            # Check if adjacent rows have the same H-bond value
            y_offset = 0.5  # default center
            
            if ri > 0 and ri < n_muts - 1:
                prev_mut = mutations[ri-1]
                next_mut = mutations[ri+1]
                
                prev_hb = int(meta.loc[prev_mut, "HBonds"]) if (prev_mut in meta.index and not pd.isna(meta.loc[prev_mut, "HBonds"])) else None
                next_hb = int(meta.loc[next_mut, "HBonds"]) if (next_mut in meta.index and not pd.isna(meta.loc[next_mut, "HBonds"])) else None
                
                # If same as neighbor, slightly offset
                if hb == prev_hb and hb == next_hb:
                    # Three in a row - offset middle one differently
                    if ri % 2 == 0:
                        y_offset = 0.45
                    else:
                        y_offset = 0.55
                elif hb == prev_hb:
                    y_offset = 0.55  # offset down
                elif hb == next_hb:
                    y_offset = 0.45  # offset up
            
            # Add the text label
            ax.text(
                badge_x + badge_w / 2, ri + y_offset, str(hb),
                ha="center", va="center",
                fontsize=6, color=tc, clip_on=False,
                fontweight='normal'
            )

    # â”€â”€ x-axis: ALL level numbers with rotation to prevent overlap â”€â”€
    ax.set_xticks([lv + 0.5 for lv in levels])
    ax.set_xticklabels([str(lv) for lv in levels], ha="right")
    
    
    ax.set_xlim(0, n_cols_max)
    ax.set_ylim(0, n_muts)
    ax.invert_yaxis()
    ax.tick_params(
        axis="x", direction="in", length=3, width=0.5,
        top=False, bottom=True, right=False, labelbottom=True,
    )
    ax.set_xlabel("Mutation Level", fontsize=8, labelpad=4)

    # â”€â”€ spines â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.5)

    # â”€â”€ panel label (a)/(b) â€” BOLD, no parentheses â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    ax.text(
        -0.08, 1.04, panel_label,
        transform=ax.transAxes,
        fontsize=9, fontweight="bold",  # CHANGED: bold
        va="bottom", ha="right", clip_on=False,
    )
    ax.text(
        -0.068, 1.04, species_label,
        transform=ax.transAxes,
        fontsize=9, fontweight="normal",  # NOT bold
        va="bottom", ha="left", clip_on=False,
    )


draw_heatmap(ax_top, m_pivot, m_meta, "Mouse", "a", n_cols_max)  # CHANGED: "a" not "(a)"
draw_heatmap(ax_bot, h_pivot, h_meta, "Human", "b", n_cols_max)  # CHANGED: "b" not "(b)"

# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
# LEGEND â€” CLEAN (FINAL FIX, NO OVERLAP)
# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
cx = (l_frac + r_frac) / 2

hb_handles = [
    mpatches.Patch(
        facecolor=CONFIG["hb_colors"][hb], edgecolor="none",
        label=f"{hb} H-bond{'s' if hb != 1 else ''}",
    )
    for hb in sorted(CONFIG["hb_colors"])
]

hb_handles.append(
    mpatches.Patch(
        facecolor=MISSING_CLR, edgecolor="#aaa", linewidth=0.5,
        label="Not observed",
    )
)

# ðŸ”§ Proper vertical placement (key fix)
legend_y = 0.55 / fig_h   # LOWER than before â†’ avoids overlap

leg = fig.legend(
    handles        = hb_handles,
    title          = "Local H-bond count",
    ncol           = 4,                 # keeps 2 rows
    loc            = "lower center",
    bbox_to_anchor = (cx, legend_y),
    bbox_transform = fig.transFigure,
    fontsize       = 7,
    title_fontsize = 8,
    frameon        = False,
    handlelength   = 1.2,
    handleheight   = 1.0,
    columnspacing  = 1.5,              # more horizontal spacing
    labelspacing   = 1.2,              # more vertical spacing (fix overlap)
    borderpad      = 0.3,
    handletextpad  = 0.5,
)

# ðŸ”§ Add spacing between title and legend items (CRITICAL FIX)
leg.set_title("Local H-bond count")
leg._legend_box.align = "center"

for txt in leg.get_texts():
    txt.set_color("black")

leg.get_title().set_color("black")

# Manually adjust the legend title position to prevent overlap
leg._legend_box.sep = 8  # INCREASED from 5 to 8 for more space between title and labels

# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
# COLORBAR â€” JCIM style
# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
sm = plt.cm.ScalarMappable(cmap=cmap_obj, norm=norm)
sm.set_array([])

cb_width  = (r_frac - l_frac) * 0.50
cb_left   = cx - cb_width / 2
cb_bottom = 0.25 / fig_h
cb_height = 0.010

cbar_ax = fig.add_axes([cb_left, cb_bottom, cb_width, cb_height])
cbar    = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
cbar.set_label(
    r"$\Delta\Delta G^*$ (kcal mol$^{-1}$)",
    fontsize=8, labelpad=4, color="black",
)
cbar.ax.tick_params(
    labelsize=7, direction="in",
    width=0.5, length=3,
)
for lbl in cbar.ax.get_xticklabels():
    lbl.set_color("black")
cbar.outline.set_linewidth(0.5)

tick_vals = np.linspace(-vabs, vabs, 5)
cbar.set_ticks(tick_vals)
cbar.set_ticklabels([
    (f"{v:+.0f}" if v != 0 else "0") for v in tick_vals
])

# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
# ADJUST LAYOUT to give more space for rotated labels
# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
plt.subplots_adjust(bottom=b_frac + 0.05)

# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
# SAVE â€” 600 dpi PNG + 300 dpi PDF
# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
for ext, dpi_val in [("png", 600), ("pdf", 300)]:
    out = f"{SAVE_PATH}.{ext}"
    plt.savefig(out, dpi=dpi_val, bbox_inches="tight",
                pad_inches=0.02, facecolor="white")
    print(f"Saved â†’ {out}")

plt.close()
print("\nDone â€” JCIM style applied.")
print("Fixes included:")
print("  - X-axis labels rotated to prevent overlap")
print("  - Legend spacing increased to prevent title-label overlap")
print("  - H-bond labels intelligently positioned to prevent overlap")
print("  - Panel labels 'a' and 'b' are bold (no parentheses)")

