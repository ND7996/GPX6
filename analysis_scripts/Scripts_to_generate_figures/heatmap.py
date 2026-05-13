"""
GPX6 Mutant Heatmap â€” JCIM publication version
Final improvements:
1. One vertical distance strip only (left side)
2. Exact distance values shown beside strip
3. Removed redundant bottom distance colorbar
4. Shared Î”G* / Î”G0 scales across species
5. JCIM double-column layout
6. Clean typography + tight spacing
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.font_manager as fm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import TwoSlopeNorm
import seaborn as sns


# =============================================================================
# FONT
# =============================================================================
_available = {f.name for f in fm.fontManager.ttflist}
FONT = "Arial" if "Arial" in _available else "Liberation Sans"


# =============================================================================
# MATPLOTLIB STYLE (JCIM)
# =============================================================================
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": [FONT, "DejaVu Sans"],
    "font.size": 6,
    "axes.labelsize": 7,
    "axes.titlesize": 8,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.4,
    "ytick.major.width": 0.4,
    "xtick.major.size": 2,
    "ytick.major.size": 2,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


# =============================================================================
# FIGURE PARAMETERS
# =============================================================================
FIG_W = 9.5
FIG_H = 3.2
WSPACE = 0.15  # Reduced from 0.45 to minimize whitespace

FS_BASE = 5
FS_LABEL = 6
FS_TITLE = 7

LW_CELL = 0.15
LW_SPINE = 0.5
LW_TICK = 0.4


# =============================================================================
# COLORS
# =============================================================================
CMAP = sns.diverging_palette(240, 10, s=75, l=50, as_cmap=True)

_cool = CMAP(np.linspace(0.0, 0.42, 256))
DIST_CMAP = mcolors.LinearSegmentedColormap.from_list("dist", _cool)


# =============================================================================
# FILES
# =============================================================================
MOUSE_CSV = r"./analysis_scripts\mouse_HBONDS.csv"
HUMAN_CSV = r"./analysis_scripts\human_HBONDS.csv"

OUT_DIR = r"./analysis_scripts/Scripts_to_generate_figures/Figures/Results\Mutational_map"
os.makedirs(OUT_DIR, exist_ok=True)


# =============================================================================
# DISTANCE OVERRIDES
# =============================================================================
MOUSE_DIST_OVERRIDE = {
    "C49U":0.00,"F48Y":3.88,"T52A":5.20,"S47A":5.59,"R99C":9.50,
    "T54Q":10.00,"H144Q":11.60,"H177Q":11.51,"G74A":14.29,"T178A":14.81,
    "E143S":15.29,"P142S":16.44,"F139L":16.28,"K87T":16.44,"Y104F":16.10,
    "G102S":16.77,"L24I":18.69,"T60A":19.05,"S4R":19.64,"N107S":19.80,
    "R181S":19.54
}

HUMAN_DIST_OVERRIDE = {
    "U49C":0.00,"Y48F":3.84,"A52T":5.49,"A47S":5.62,"Q54T":10.42,
    "C99R":9.72,"Q144H":11.47,"Q177H":11.46,"A178T":14.21,"A74G":14.58,
    "S143E":15.12,"S142P":16.59,"L139F":16.54,"T87K":16.84,"F104Y":17.08,
    "N3K":17.66,"S102G":18.15,"L24I":18.51,"H173R":19.61,"S181R":19.23,
    "A60T":19.79
}


# =============================================================================
# ORDER
# =============================================================================
MOUSE_ORDER = [
"C49U","F48Y","T52A","S47A","R99C","T54Q","H144Q","H177Q","G74A",
"T178A","E143S","P142S","F139L","K87T","Y104F","G102S","L24I",
"T60A","S4R","N107S","R181S"
]

HUMAN_ORDER = [
"U49C","Y48F","A52T","A47S","Q54T","C99R","Q144H","Q177H","A178T",
"T87K","A74G","S143E","S142P","N3K","L139F","F104Y","A60T",
"L24I","S102G","H173R","S181R"
]


# =============================================================================
# FUNCTIONS
# =============================================================================
def load_data(path, wt_name, dist_override):

    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()

    wt = df[(df.level=="Level0") & (df.mutation==wt_name)].iloc[0]

    wt_star = float(wt.dg_star)
    wt_dg0  = float(wt.dg0)

    df = df[df.mutation != wt_name].copy()

    agg = df.groupby(["mutation","level"],as_index=False).agg(
        dg_star=("dg_star","mean"),
        dg0=("dg0","mean"),
        dist=("dist_to_Sec/Cys49","first")
    )

    star = agg.pivot(index="mutation",columns="level",values="dg_star")
    dg0  = agg.pivot(index="mutation",columns="level",values="dg0")

    dist = agg.groupby("mutation")["dist"].first()

    for k,v in dist_override.items():
        if k in dist.index:
            dist[k]=v

    levels = sorted(star.columns,key=lambda x:int(x.replace("Level","")))
    return star[levels], dg0[levels], wt_star, wt_dg0, dist


def reorder(obj, order):
    return obj.reindex(list(reversed(order)))


def global_norm(p1,p2,wt1,wt2):
    vals = np.concatenate([p1.values.ravel(),p2.values.ravel()])
    vals = vals[~np.isnan(vals)]
    vc = np.mean([wt1,wt2])
    half = max(abs(vals.min()-vc),abs(vals.max()-vc))
    return TwoSlopeNorm(vmin=vc-half,vcenter=vc,vmax=vc+half)


def draw_panel(ax,data,norm,title,label_y=True):

    nrows,ncols = data.shape

    mesh = ax.pcolormesh(
        data.values,
        cmap=CMAP,
        norm=norm,
        edgecolors="white",
        linewidths=LW_CELL,
        rasterized=True
    )

    ax.set_xlim(0,ncols)
    ax.set_ylim(0,nrows)

    ax.set_xticks(np.arange(ncols)+0.5)
    ax.set_xticklabels(
        [x.replace("Level","L") for x in data.columns],
        ha="right"
    )

    ax.set_yticks(np.arange(nrows)+0.5)
    
    # Improved y-axis labels: keep original mutation names without rotation
    ax.set_yticklabels(data.index, ha="right", fontsize=FS_BASE)
    
    # Ensure all y-labels are visible
    ax.tick_params(axis='y', which='major', pad=1, left=True, right=False)

    ax.set_xlabel("Level")
    if label_y:
        ax.set_ylabel("Mutation", labelpad=8)

    ax.set_title(title,pad=3)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    div = make_axes_locatable(ax)
    cax = div.append_axes("right",size="3%",pad=0.04)
    cb = plt.colorbar(mesh,cax=cax)
    cb.ax.tick_params(labelsize=5)
    
    # Add label for colorbar
    if "Î”G*" in title:
        cb.set_label("kcal/mol", fontsize=5, labelpad=1)
    else:
        cb.set_label("kcal/mol", fontsize=5, labelpad=1)

    return nrows


def draw_distance_strip(fig, ax, dist_vals, vmax):

    fig.canvas.draw()
    pos = ax.get_position()

    strip_w = 0.015
    gap = 0.045  # Slightly reduced gap

    axd = fig.add_axes([
        pos.x0-gap-strip_w,
        pos.y0,
        strip_w,
        pos.height
    ])

    norm = mcolors.Normalize(vmin=0,vmax=vmax)

    axd.pcolormesh(
        dist_vals.reshape(-1,1),
        cmap=DIST_CMAP,
        norm=norm
    )

    n = len(dist_vals)

    axd.set_xlim(0,1)
    axd.set_ylim(0,n)

    axd.set_xticks([])

    axd.set_yticks(np.arange(n)+0.5)
    axd.set_yticklabels(
        [f"{x:.1f}" for x in dist_vals],
        fontsize=5
    )

    axd.tick_params(left=False,right=False,pad=1)

    for s in axd.spines.values():
        s.set_visible(False)

    axd.set_title("dist./n(Ã…)",fontsize=5,pad=2)


def make_figure(star,dg0,norm_star,norm_dg0,dist,outfile):

    fig, (ax1, ax2) = plt.subplots(
        1, 2,
        figsize=(FIG_W, FIG_H),
        gridspec_kw={"wspace": WSPACE, "width_ratios": [1, 1]}
    )

    # Adjust left margin to accommodate distance strip
    fig.subplots_adjust(left=0.14, right=0.96)

    draw_panel(ax1, star, norm_star, r"$\Delta G^*$", True)
    draw_panel(ax2, dg0, norm_dg0, r"$\Delta G^0$", False)

    # Remove y-label from right panel to reduce redundancy
    ax2.set_ylabel("")

    vmax = np.nanmax(dist.values)

    vals = dist.reindex(star.index).values
    draw_distance_strip(fig, ax1, vals, vmax)

    for ext in [".png", ".pdf"]:
        fig.savefig(outfile+ext, bbox_inches="tight", pad_inches=0.02)

    plt.close(fig)


# =============================================================================
# MAIN
# =============================================================================
if __name__ == "__main__":

    m_star,m_dg0,m_wts,m_wtd,m_dist = load_data(
        MOUSE_CSV,"mousecys",MOUSE_DIST_OVERRIDE
    )

    h_star,h_dg0,h_wts,h_wtd,h_dist = load_data(
        HUMAN_CSV,"humansec",HUMAN_DIST_OVERRIDE
    )

    m_star = reorder(m_star,MOUSE_ORDER)
    m_dg0  = reorder(m_dg0,MOUSE_ORDER)
    m_dist = reorder(m_dist,MOUSE_ORDER)

    h_star = reorder(h_star,HUMAN_ORDER)
    h_dg0  = reorder(h_dg0,HUMAN_ORDER)
    h_dist = reorder(h_dist,HUMAN_ORDER)

    norm_star = global_norm(m_star,h_star,m_wts,h_wts)
    norm_dg0  = global_norm(m_dg0,h_dg0,m_wtd,h_wtd)

    make_figure(
        m_star,m_dg0,norm_star,norm_dg0,m_dist,
        os.path.join(OUT_DIR,"mutant_heatmap_mouse")
    )

    make_figure(
        h_star,h_dg0,norm_star,norm_dg0,h_dist,
        os.path.join(OUT_DIR,"mutant_heatmap_human")
    )

    print("Done.")

