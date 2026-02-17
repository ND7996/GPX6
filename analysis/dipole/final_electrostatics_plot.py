"""
heatmap_levels_gpx6.py
======================
Two heatmaps stacked vertically: Human (top) | Mouse (bottom)
Rows    = mutations, sorted by H-bond count then Level-0 ΔΔG*
Columns = mutation level (0 → max)
Colour  = ΔΔG* (diverging RdBu_r, white = 0)
Annotations: H-bond count badge on left, active-site marker on row label
Grey cells = mutation not observed at that level
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from scipy import stats

CONFIG = dict(
    human_file = r"D:\PhD_Thesis\GPX6\analysis\dipole\human_mutants_merged.csv",
    mouse_file = r"D:\PhD_Thesis\GPX6\analysis\dipole\mouse_mutants_merged.csv",
    human_refs = ["humansec", "humancys"],
    mouse_refs = ["mousecys", "C49U"],
    output     = r"D:\PhD_Thesis\GPX6\analysis\dipole\fig_heatmap_levels",
    cmap       = "RdBu_r",
    c_human    = "#1A6FAF",
    c_mouse    = "#D4860A",
    # H-bond badge colours (matches plasma scale concept)
    hb_colors  = {0:"#440154", 1:"#3B528B", 2:"#21908C",
                  3:"#5DC863", 4:"#FDE725", 5:"#FF6B35"},
)

plt.rcParams.update({
    "figure.facecolor" : "white",
    "axes.facecolor"   : "white",
    "font.family"      : "DejaVu Sans",
    "font.size"        : 8.5,
    "figure.dpi"       : 300,
})

# ── LOAD ──────────────────────────────────────────────────────────────────────
def get_wt_ref(path, refs):
    df = pd.read_csv(path)
    return df.loc[df["mutation"].isin(refs) & (df["level"] == "Level0"),
                  "dg_star"].mean()

def load_pivot(path, refs, ref_dg):
    df = pd.read_csv(path)
    df = df[~df["mutation"].isin(refs)].copy()
    df["level_num"] = df["level"].str.replace("Level","").astype(int)
    df["ddg"]       = df["dg_star"] - ref_dg

    # Propagate HBond / active-site info from Level0
    l0 = df[df["level_num"] == 0].set_index("mutation")
    df["HBonds"]    = df["mutation"].map(l0["Total_HBonds"])
    df["AS_HBonds"] = df["mutation"].map(
                          l0["Active_Site_HBonds"].fillna(0))

    # Build pivot: rows=mutation, cols=level
    pivot = df.pivot_table(index="mutation", columns="level_num",
                           values="ddg", aggfunc="mean")

    # Mutation metadata (sorted by HBonds, then Level-0 ΔΔG*)
    meta = (df[df["level_num"] == 0]
            .drop_duplicates("mutation")
            .set_index("mutation")[["HBonds","AS_HBonds"]]
            .dropna(subset=["HBonds"]))
    meta["ddg_l0"] = pivot.get(0, pd.Series(dtype=float))
    meta = meta.sort_values(["HBonds","ddg_l0"])

    # Reindex pivot rows to match sorted metadata
    pivot = pivot.reindex(meta.index)
    return pivot, meta

h_ref = get_wt_ref(CONFIG["human_file"], CONFIG["human_refs"])
m_ref = get_wt_ref(CONFIG["mouse_file"], CONFIG["mouse_refs"])
h_pivot, h_meta = load_pivot(CONFIG["human_file"], CONFIG["human_refs"], h_ref)
m_pivot, m_meta = load_pivot(CONFIG["mouse_file"], CONFIG["mouse_refs"], m_ref)

# Unified colour scale across both species
vabs = max(np.nanmax(np.abs(h_pivot.values)),
           np.nanmax(np.abs(m_pivot.values)))
norm = mcolors.TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)
cmap = plt.get_cmap(CONFIG["cmap"])
missing_color = "#e0e0e0"

# ── DRAW ONE HEATMAP ──────────────────────────────────────────────────────────
def draw_heatmap(ax, pivot, meta, species, sp_color):
    mutations = list(pivot.index)
    levels    = sorted(pivot.columns)
    n_muts    = len(mutations)
    n_levels  = len(levels)

    # Cell rectangles
    for ri, mut in enumerate(mutations):
        for ci, lv in enumerate(levels):
            val = pivot.loc[mut, lv] if lv in pivot.columns else np.nan
            if pd.isna(val):
                face = missing_color
            else:
                face = cmap(norm(val))
            rect = plt.Rectangle([ci, ri], 1, 1,
                                  facecolor=face,
                                  edgecolor="white",
                                  linewidth=0.5, zorder=1)
            ax.add_patch(rect)

            # Annotate cell value only at Level 0 for readability
            if lv == 0 and not pd.isna(val):
                lum = 0.299*face[0] + 0.587*face[1] + 0.114*face[2]
                tc  = "white" if lum < 0.45 else "#1a1a1a"
                ax.text(ci + 0.5, ri + 0.5, f"{val:+.1f}",
                        ha="center", va="center",
                        fontsize=6.5, color=tc, zorder=2,
                        fontweight="bold")

    # Y-axis: mutation labels with H-bond badge + active-site star
    ax.set_yticks([r + 0.5 for r in range(n_muts)])
    ylabels = []
    for mut in mutations:
        hb  = int(meta.loc[mut, "HBonds"]) if mut in meta.index else "?"
        ast = " ★" if (mut in meta.index and
                       meta.loc[mut, "AS_HBonds"] > 0) else ""
        ylabels.append(f"{mut}{ast}")
    ax.set_yticklabels(ylabels, fontsize=8)

    # Colour the y-tick labels by H-bond count
    for ri, (tick, mut) in enumerate(zip(ax.get_yticklabels(), mutations)):
        if mut in meta.index:
            hb = int(meta.loc[mut, "HBonds"])
            tick.set_color(CONFIG["hb_colors"].get(hb, "#333333"))
            if mut in meta.index and meta.loc[mut, "AS_HBonds"] > 0:
                tick.set_fontweight("bold")

    # H-bond count badge strip on the right
    badge_x = n_levels + 0.15
    for ri, mut in enumerate(mutations):
        if mut not in meta.index: continue
        hb    = int(meta.loc[mut, "HBonds"])
        color = CONFIG["hb_colors"].get(hb, "#888888")
        badge = FancyBboxPatch(
            (badge_x, ri + 0.15), 0.65, 0.70,
            boxstyle="round,pad=0.05",
            facecolor=color, edgecolor="white",
            linewidth=0.4, zorder=3,
            clip_on=False,
        )
        ax.add_patch(badge)
        lum = 0.299 * plt.cm.colors.to_rgb(color)[0] \
              if hasattr(plt.cm, 'colors') else 0.5
        face_rgb = mcolors.to_rgb(color)
        lum = 0.299*face_rgb[0] + 0.587*face_rgb[1] + 0.114*face_rgb[2]
        tc  = "white" if lum < 0.55 else "#1a1a1a"
        ax.text(badge_x + 0.325, ri + 0.5, str(hb),
                ha="center", va="center",
                fontsize=6.5, color=tc, fontweight="bold",
                zorder=4, clip_on=False)

    # X-axis: level numbers, every other tick
    shown_levels = [l for l in levels if l % 2 == 0]
    shown_pos    = [levels.index(l) + 0.5 for l in shown_levels]
    ax.set_xticks(shown_pos)
    ax.set_xticklabels([str(l) for l in shown_levels], fontsize=7.5)

    ax.set_xlim(0, n_levels)
    ax.set_ylim(0, n_muts)
    ax.invert_yaxis()

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)

    # Species label bar
    ax.text(-0.5, -0.8, species,
            transform=ax.transData, fontsize=11,
            fontweight="bold", color=sp_color,
            va="bottom", ha="left")

    # Thin separator line below title

# ── FIGURE ────────────────────────────────────────────────────────────────────
n_h = len(h_pivot)
n_m = len(m_pivot)
n_lv_h = len(h_pivot.columns)
n_lv_m = len(m_pivot.columns)

cell_w = 0.38
cell_h = 0.34
badge_w = 0.5

fig_w = max(n_lv_h, n_lv_m) * cell_w + 3.5
fig_h = (n_h + n_m) * cell_h + 3.2

fig, (ax_h, ax_m) = plt.subplots(
    2, 1, figsize=(fig_w, fig_h),
    gridspec_kw={"height_ratios": [n_h, n_m], "hspace": 0.55},
)

draw_heatmap(ax_h, h_pivot, h_meta, "Human", CONFIG["c_human"])
draw_heatmap(ax_m, m_pivot, m_meta, "Mouse", CONFIG["c_mouse"])

# Shared x-label
for ax in [ax_h, ax_m]:
    ax.set_xlabel("Mutation Level", fontsize=9, fontweight="bold", labelpad=6)

ax_h.set_ylabel("Mutation", fontsize=9, fontweight="bold", labelpad=6)
ax_m.set_ylabel("Mutation", fontsize=9, fontweight="bold", labelpad=6)

# ── COLORBAR ──────────────────────────────────────────────────────────────────
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar_ax = fig.add_axes([1.01, 0.25, 0.018, 0.50])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label("ΔΔG*  (kcal mol⁻¹)", fontsize=8.5, labelpad=8)
cbar.ax.tick_params(labelsize=7.5)
cbar.outline.set_linewidth(0.4)
cbar.ax.axhline(norm(0), color="#333333", lw=1.2, ls="--", alpha=0.6)

# ── H-BOND LEGEND ─────────────────────────────────────────────────────────────
hb_patches = [
    mpatches.Patch(facecolor=CONFIG["hb_colors"][hb],
                   edgecolor="white", linewidth=0.5,
                   label=f"{hb} H-bond{'s' if hb!=1 else ''}")
    for hb in sorted(CONFIG["hb_colors"])
]
hb_patches += [
    mpatches.Patch(facecolor=missing_color, edgecolor="#cccccc",
                   linewidth=0.5, label="Not observed at level"),
]
fig.legend(
    handles=hb_patches,
    title="H-bond count\n(badge + label colour)",
    title_fontsize=8,
    loc="lower right",
    bbox_to_anchor=(1.14, 0.02),
    fontsize=7.5, frameon=True,
    facecolor="white", edgecolor="#cccccc",
    handlelength=1.2, labelspacing=0.35,
)

# ── ANNOTATIONS ───────────────────────────────────────────────────────────────
fig.text(1.01, 0.78,
         "★ = active-site\n    H-bond mutant\n\nBold label =\n    active-site",
         fontsize=7.5, color="#333333", va="top")
fig.text(0.5, 0.995,
         "Level-0 ΔΔG* annotated in first column  ·  "
         "Grey = mutation absent at that level  ·  "
         "Rows sorted by H-bond count then Level-0 ΔΔG*",
         ha="center", fontsize=8, color="#555555")

# ── TITLE ─────────────────────────────────────────────────────────────────────
fig.suptitle(
    "ΔΔG* Across All Mutation Levels — Human vs Mouse GPX6",
    fontsize=12, fontweight="bold", y=1.012,
)

plt.savefig(f"{CONFIG['output']}.png", dpi=300, bbox_inches="tight",
            facecolor="white")
plt.savefig(f"{CONFIG['output']}.pdf", bbox_inches="tight", facecolor="white")
print(f"Saved: {CONFIG['output']}.png + .pdf")