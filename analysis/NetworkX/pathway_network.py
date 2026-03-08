"""
gpx6_all_pathways.py
====================
Shows ALL measured dG* pathways from the CSV.
At each level, every remaining mutation was measured by EVB in that background.
Each mutation traces a path across levels — these are the real possible pathways.
The greedy path (minimum dG* at each level) is highlighted for comparison.
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

C_HU = "#1A6FAF"
C_MO = "#D4860A"

# ── DATA ──────────────────────────────────────────────────────────────────────
def load(path, refs):
    df = pd.read_csv(path)
    wt_sec = df.loc[df["mutation"] == refs[0], "dg_star"].values[0]
    wt_sec_err = df.loc[df["mutation"] == refs[0], "dg_star_error"].values[0]
    wt_cys = df.loc[df["mutation"] == refs[1], "dg_star"].values[0]
    wt_cys_err = df.loc[df["mutation"] == refs[1], "dg_star_error"].values[0]

    df2 = df[~df["mutation"].isin(refs)].copy()
    df2["level_num"] = df2["level"].str.replace("Level", "").astype(int)

    levels = sorted(df2["level_num"].unique())

    # For each mutation: collect (level, dg_star, dg_star_error) across all levels
    mutations = df2["mutation"].unique()
    mut_paths = {}
    for mut in mutations:
        rows = df2[df2["mutation"] == mut].sort_values("level_num")
        mut_paths[mut] = {
            "levels": rows["level_num"].values,
            "dg":     rows["dg_star"].values,
            "err":    rows["dg_star_error"].values,
        }

    # Greedy path: at each level pick mutation with lowest dg_star
    greedy_levels, greedy_muts, greedy_dg, greedy_err = [], [], [], []
    for lv in levels:
        lv_data = df2[df2["level_num"] == lv]
        best = lv_data.loc[lv_data["dg_star"].idxmin()]
        greedy_levels.append(lv)
        greedy_muts.append(best["mutation"])
        greedy_dg.append(best["dg_star"])
        greedy_err.append(best["dg_star_error"])

    return {
        "wt_sec": wt_sec, "wt_sec_err": wt_sec_err,
        "wt_cys": wt_cys, "wt_cys_err": wt_cys_err,
        "levels": np.array(levels),
        "mut_paths": mut_paths,
        "greedy_levels": np.array(greedy_levels),
        "greedy_muts":   greedy_muts,
        "greedy_dg":     np.array(greedy_dg),
        "greedy_err":    np.array(greedy_err),
        "all_dg": df2["dg_star"].values,
    }

h = load(HUMAN_FILE, HUMAN_REFS)
m = load(MOUSE_FILE, MOUSE_REFS)

# ── FIGURE: two tall wide panels stacked vertically ──────────────────────────
fig, axes = plt.subplots(2, 1, figsize=(22, 16),
                          gridspec_kw={"hspace": 0.52})

for ax, d, sp_color, species, direction, panel_label in [
    (axes[0], h, C_HU, "Human GPX6", "Sec → Mouse-like", "A"),
    (axes[1], m, C_MO, "Mouse GPX6", "Cys → Human-like", "B"),
]:
    levels      = d["levels"]
    wt_sec      = d["wt_sec"]
    wt_cys      = d["wt_cys"]
    greedy_muts_set = set(d["greedy_muts"])

    # Colour each non-greedy mutation line by its mean dG* across levels
    # (yellow = low barrier throughout, red = high barrier throughout)
    all_mean_dg = [d["mut_paths"][m]["dg"].mean()
                   for m in d["mut_paths"] if m not in greedy_muts_set
                   or True]  # include all for colour range
    vmin = min(d["all_dg"])
    vmax = max(d["all_dg"])
    cmap = plt.colormaps["RdYlBu_r"]
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    # ── Draw each mutation's trajectory ──────────────────────────────────
    greedy_mut_sequence = d["greedy_muts"]  # ordered list

    for mut, path in d["mut_paths"].items():
        lv   = path["levels"]
        dg   = path["dg"]
        mean_dg = dg.mean()

        is_greedy_chosen = mut in greedy_muts_set

        if not is_greedy_chosen:
            # Non-greedy: thin coloured line
            ax.plot(lv, dg,
                    color=cmap(norm(mean_dg)),
                    lw=1.0, alpha=0.45, zorder=2)
        else:
            # Will draw greedy on top separately
            pass

    # ── Draw individual mutation paths that were greedy-selected ─────────
    # (these are the same mutations, just highlighted at their selected level)
    for mut, path in d["mut_paths"].items():
        if mut in greedy_muts_set:
            lv  = path["levels"]
            dg  = path["dg"]
            mean_dg = dg.mean()
            ax.plot(lv, dg,
                    color=cmap(norm(mean_dg)),
                    lw=1.2, alpha=0.55, zorder=3)

    # ── Reference lines ───────────────────────────────────────────────────
    ax.axhline(wt_sec, color=sp_color, lw=1.8, ls="--", alpha=0.85, zorder=6,
               label=f"Sec WT  {wt_sec:.2f} kcal mol⁻¹")
    ax.axhline(wt_cys, color="#333333", lw=2.0, ls="-.", alpha=0.9, zorder=6,
               label=f"Cys WT (viability limit)  {wt_cys:.2f} kcal mol⁻¹")

    # ── Non-viable shading ────────────────────────────────────────────────
    y_top = vmax + 5
    ax.axhspan(wt_cys, y_top, alpha=0.05, color="#d62728", zorder=0)
    ax.text(levels[-1] * 0.92, wt_cys + 0.35,
            "Non-viable  (dG* > Cys WT)",
            color="#c0392b", fontsize=9, alpha=0.85,
            ha="right", style="italic", zorder=7)

    # ── Greedy path bold on top with error bars ───────────────────────────
    ax.errorbar(d["greedy_levels"], d["greedy_dg"],
                yerr=d["greedy_err"],
                color=sp_color, lw=3.5, zorder=10,
                marker="o", markersize=7,
                markerfacecolor="white",
                markeredgecolor=sp_color,
                markeredgewidth=2.5,
                elinewidth=1.5, capsize=4, capthick=1.5,
                ecolor=sp_color,
                label="Greedy path (lowest dG* at each level)")

    # ── ALL mutation labels on greedy path — alternating above/below ──────
    y_all_min = vmin - 1.5
    y_all_max = vmax + 5

    # Two label rows above (above non-viable zone) and two below data
    label_above_y = [wt_cys + 1.4, wt_cys + 2.7]
    label_below_y = [y_all_min - 0.2, y_all_min - 1.5]

    for idx, (lv, mut, dg_val) in enumerate(
            zip(d["greedy_levels"], d["greedy_muts"], d["greedy_dg"])):
        if idx % 2 == 0:
            row  = (idx // 2) % 2
            y_tx = label_above_y[row]
            va   = "bottom"
        else:
            row  = (idx // 2) % 2
            y_tx = label_below_y[row]
            va   = "top"

        ax.annotate(
            mut,
            xy=(lv, dg_val),
            xytext=(lv, y_tx),
            fontsize=7.5,
            ha="center", va=va,
            color=sp_color,
            fontweight="bold",
            arrowprops=dict(
                arrowstyle="-",
                color=sp_color,
                lw=0.9, alpha=0.6,
                shrinkA=0, shrinkB=2,
            ),
            zorder=15,
            annotation_clip=False,
        )

    # ── Axes formatting ───────────────────────────────────────────────────
    y_bot = label_below_y[1] - 0.8
    y_top_ax = label_above_y[1] + 1.0
    ax.set_ylim(y_bot, y_top_ax)
    ax.set_xlim(-0.5, levels[-1] + 0.5)
    ax.set_xticks(levels)
    ax.set_xticklabels([str(int(l)) for l in levels], fontsize=9)
    ax.set_xlabel("Mutation level", fontsize=12, fontweight="bold")
    ax.set_ylabel("Activation free energy  ΔG* (kcal mol⁻¹)",
                  fontsize=12, fontweight="bold")
    ax.set_title(
        f"{panel_label}   {species}    ({direction})",
        fontweight="bold", color=sp_color,
        fontsize=13, pad=10, loc="left",
    )
    ax.legend(fontsize=9, framealpha=0.9,
              loc="upper right", handlelength=2.0,
              borderpad=0.6)
    for spine in ax.spines.values():
        spine.set_linewidth(0.6)
    ax.tick_params(axis="y", labelsize=9.5)

    # ── Shared colorbar per panel ─────────────────────────────────────────
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.55, pad=0.01, aspect=25)
    cbar.set_label("Mean ΔG* of mutation across levels\n(kcal mol⁻¹)",
                   fontsize=8.5, labelpad=6)
    cbar.ax.tick_params(labelsize=8)
    cbar.outline.set_linewidth(0.4)

# ── Main title ─────────────────────────────────────────────────────────────
fig.suptitle(
    "All EVB-measured mutational pathways connecting Human and Mouse GPX6\n"
    "Each line = one mutation tracked across all levels in its epistatic background  ·  "
    "Bold line = greedy path (minimum ΔG* selected at each level)",
    fontsize=11.5, fontweight="bold", y=1.01,
)

plt.savefig(OUTPUT_FILE, dpi=300, bbox_inches="tight", facecolor="white")
print(f"Saved: {OUTPUT_FILE}")