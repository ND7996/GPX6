"""
gpx6_all_pathways.py
====================
- NO colorbars
- Unified y-axis scale across both panels
- Legend placed BELOW each panel (no overlap with data)
- All 4 reference lines use FIXED colour + linestyle in BOTH panels:
    Human Sec WT  → red    dashed   (--)
    Human Cys WT  → orange dash-dot (-.)
    Mouse Sec WT  → green  dotted   (: )
    Mouse Cys WT  → purple dense dashes
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D

# ── CONFIG ────────────────────────────────────────────────────────────────────
HUMAN_FILE  = r"D:\PhD_Thesis\GPX6\analysis\dipole\human_HBONDS.csv"
MOUSE_FILE  = r"D:\PhD_Thesis\GPX6\analysis\dipole\mouse_HBONDS.csv"
HUMAN_REFS  = ["humansec", "humancys"]
MOUSE_REFS  = ["mousecys", "C49U"]
OUTPUT_FILE = r"D:\PhD_Thesis\GPX6\analysis\NetworkX\gpx6_all_pathways.png"

os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)

plt.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor":   "white",
    "font.family":      "DejaVu Sans",
    "font.size":        9,
    "figure.dpi":       300,
})

# ── COLOURS / STYLES — fixed per species, identical in both panels ─────────────
C_GREEDY = "#1A6FAF"

REF_STYLES = {
    "hsec": {"col": "#E63946", "ls": "--",        "lw": 2.0, "label": "Human Sec WT"},
    "hcys": {"col": "#FF9F1C", "ls": "-.",         "lw": 2.0, "label": "Human Cys WT"},
    "msec": {"col": "#2DC653", "ls": ":",          "lw": 2.0, "label": "Mouse Sec WT (C49U)"},
    "mcys": {"col": "#7B2D8B", "ls": (0, (5, 1)), "lw": 2.5, "label": "Mouse Cys WT"},
}

# ── DATA ──────────────────────────────────────────────────────────────────────
def load(path, refs):
    df = pd.read_csv(path)
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
        best = lv_data.loc[lv_data["dg_star"].idxmin()]
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

h = load(HUMAN_FILE, HUMAN_REFS)
m = load(MOUSE_FILE, MOUSE_REFS)

h_raw = pd.read_csv(HUMAN_FILE)
m_raw = pd.read_csv(MOUSE_FILE)

def get_ref(df, name):
    return float(df.loc[df["mutation"] == name, "dg_star"].values[0])

hsec_dg = get_ref(h_raw, "humansec")
hcys_dg = get_ref(h_raw, "humancys")
mcys_dg = get_ref(m_raw, "mousecys")
msec_dg = get_ref(m_raw, "C49U")

# Draw order — same every panel
REF_ORDER = [("hsec", hsec_dg), ("hcys", hcys_dg),
             ("msec", msec_dg), ("mcys", mcys_dg)]

# ── GLOBAL Y SCALE ────────────────────────────────────────────────────────────
global_vmin = min(min(h["all_dg"]), min(m["all_dg"]))
global_vmax = max(max(h["all_dg"]), max(m["all_dg"]))

# ── STAGGER HELPER ─────────────────────────────────────────────────────────────
def stagger_labels(items, min_gap=0.85):
    sorted_items = sorted(items, key=lambda x: x[0])
    label_y = [r[0] for r in sorted_items]

    for i in range(1, len(label_y)):
        if label_y[i] - label_y[i-1] < min_gap:
            label_y[i] = label_y[i-1] + min_gap

    for i in range(len(label_y)-2, -1, -1):
        if label_y[i+1] - label_y[i] < min_gap:
            label_y[i] = label_y[i+1] - min_gap

    return {r[1]: (label_y[j], r[0]) for j, r in enumerate(sorted_items)}

# ── FIGURE ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 1, figsize=(26, 18),
                          gridspec_kw={"hspace": 0.60})
fig.subplots_adjust(left=0.055, right=0.76, top=0.94, bottom=0.10)

for ax, d, species, direction, panel_label in [
    (axes[0], h, "Human GPX6", "Sec → Mouse-like", "A"),
    (axes[1], m, "Mouse GPX6", "Cys → Human-like", "B"),
]:
    levels          = d["levels"]
    greedy_muts_set = set(d["greedy_muts"])
    cmap            = plt.colormaps["Blues"]
    norm            = mcolors.Normalize(vmin=global_vmin - 2, vmax=global_vmax + 2)

    # ── Mutation trajectories ─────────────────────────────────────────────
    for mut, path in d["mut_paths"].items():
        mean_dg = path["dg"].mean()
        lw      = 1.2 if mut in greedy_muts_set else 0.9
        alpha   = 0.55 if mut in greedy_muts_set else 0.35
        zorder  = 3    if mut in greedy_muts_set else 2
        ax.plot(path["levels"], path["dg"],
                color=cmap(norm(mean_dg)), lw=lw, alpha=alpha, zorder=zorder)

    # ── Viability limit (for shading only — no text label) ───────────────
    viability_dg = mcys_dg if panel_label == "A" else hcys_dg

    # ── Draw all 4 reference lines with FIXED styles ──────────────────────
    for key, dg_val in REF_ORDER:
        s = REF_STYLES[key]
        ax.axhline(dg_val, color=s["col"], lw=s["lw"],
                   ls=s["ls"], alpha=0.95, zorder=6)

    # ── Staggered inline labels (right side, no viability suffix) ────────
    label_map  = stagger_labels([(dg, k) for k, dg in REF_ORDER], min_gap=0.85)
    x_line_end = levels[-1] + 0.4
    x_text     = levels[-1] + 0.65

    for key, dg_val in REF_ORDER:
        s  = REF_STYLES[key]
        ly, ty = label_map[key]
        ax.annotate(
            f"{s['label']}\n{dg_val:.2f} kcal mol⁻¹",
            xy=(x_line_end, ty), xytext=(x_text, ly),
            fontsize=8.5, color=s["col"], fontweight="bold",
            va="center", ha="left", clip_on=False, zorder=20,
            arrowprops=dict(arrowstyle="-", color=s["col"],
                            lw=1.0, alpha=0.85, shrinkA=0, shrinkB=0),
            annotation_clip=False,
        )

    # ── Non-viable shading (no text label on it) ──────────────────────────
    ax.axhspan(viability_dg, global_vmax + 8, alpha=0.06, color="#d62728", zorder=0)

    # ── Greedy path ───────────────────────────────────────────────────────
    ax.errorbar(d["greedy_levels"], d["greedy_dg"],
                yerr=d["greedy_err"],
                color=C_GREEDY, lw=3.5, zorder=10,
                marker="o", markersize=7,
                markerfacecolor="white", markeredgecolor=C_GREEDY,
                markeredgewidth=2.5, elinewidth=1.5,
                capsize=4, capthick=1.5, ecolor=C_GREEDY)

    # ── Mutation labels (alternating above / below) ───────────────────────
    y_all_min     = global_vmin - 1.5
    label_above_y = [viability_dg + 1.5, viability_dg + 3.0]
    label_below_y = [y_all_min - 0.3,    y_all_min - 1.8]

    for idx, (lv, mut, dg_val) in enumerate(
            zip(d["greedy_levels"], d["greedy_muts"], d["greedy_dg"])):
        if idx % 2 == 0:
            y_tx = label_above_y[(idx // 2) % 2]; va = "bottom"
        else:
            y_tx = label_below_y[(idx // 2) % 2]; va = "top"
        ax.annotate(mut, xy=(lv, dg_val), xytext=(lv, y_tx),
            fontsize=7.5, ha="center", va=va,
            color=C_GREEDY, fontweight="bold",
            arrowprops=dict(arrowstyle="-", color=C_GREEDY,
                            lw=0.9, alpha=0.6, shrinkA=0, shrinkB=2),
            zorder=15, annotation_clip=False)

    # ── Legend below axes ─────────────────────────────────────────────────
    legend_handles = [
        Line2D([0], [0], color=C_GREEDY, lw=3.0, marker="o",
               markerfacecolor="white", markeredgecolor=C_GREEDY,
               markeredgewidth=2, markersize=6,
               label="Greedy path (lowest dG* at each level)"),
    ]
    for key, dg_val in REF_ORDER:
        s = REF_STYLES[key]
        legend_handles.append(
            Line2D([0], [0], color=s["col"], lw=2.2, ls=s["ls"],
                   label=f"{s['label']}  ({dg_val:.2f} kcal mol⁻¹)")
        )
    ax.legend(
        handles=legend_handles,
        fontsize=8.5, framealpha=0.95,
        loc="lower center",
        bbox_to_anchor=(0.38, -0.30),
        handlelength=2.5, borderpad=0.8, ncol=3,
        title=f"{'Human' if panel_label=='A' else 'Mouse'} GPX6 — Reference lines",
        title_fontsize=9.0,
    )

    # ── Axes ──────────────────────────────────────────────────────────────
    y_bot    = label_below_y[1] - 0.8
    y_top_ax = label_above_y[1] + 1.2
    ax.set_ylim(y_bot, y_top_ax)
    ax.set_xlim(-0.5, levels[-1] + 0.5)
    ax.set_xticks(levels)
    ax.set_xticklabels([str(int(l)) for l in levels], fontsize=9)
    ax.set_xlabel("Mutation level", fontsize=12, fontweight="bold")
    ax.set_ylabel("Activation free energy  ΔG* (kcal mol⁻¹)",
                  fontsize=12, fontweight="bold")
    ax.set_title(
        f"{panel_label}   {species}    ({direction})",
        fontweight="bold", color=C_GREEDY, fontsize=13, pad=10, loc="left",
    )
    for spine in ax.spines.values():
        spine.set_linewidth(0.6)
    ax.tick_params(axis="y", labelsize=9.5)

# ── Unify y-axis ──────────────────────────────────────────────────────────────
y_lims     = [ax.get_ylim() for ax in axes]
shared_bot = min(y[0] for y in y_lims)
shared_top = max(y[1] for y in y_lims)
for ax in axes:
    ax.set_ylim(shared_bot, shared_top)

# ── Main title ────────────────────────────────────────────────────────────────
fig.suptitle(
    "All EVB-measured mutational pathways connecting Human and Mouse GPX6\n"
    "Each line = one mutation tracked across all levels in its epistatic background  ·  "
    "Bold line = greedy path (minimum ΔG* selected at each level)",
    fontsize=11.5, fontweight="bold", y=0.98,
)

plt.savefig(OUTPUT_FILE, dpi=300, bbox_inches="tight", facecolor="white")
print(f"Saved: {OUTPUT_FILE}")