import sys
import os

ACS_PATH = r"./analysis_scripts/Scripts_to_generate_figures/Figures"

if ACS_PATH not in sys.path:
    sys.path.append(ACS_PATH)

# ===== JCIM formatting =====
from acsfonts import jcim_figure, save_jcim_figure

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ==========================================================
# OUTPUT DIRECTORY
# ==========================================================
OUTDIR = Path(r"./analysis_scripts/Scripts_to_generate_figures/Figures/parameter_plots")
OUTDIR.mkdir(parents=True, exist_ok=True)

jcim_figure()

# ==========================================================
# FLAT-BOTTOM HARMONIC POTENTIAL
# ==========================================================
def flat_bottom_harmonic_potential(d, d_min, d_max, k):
    if d_min <= d <= d_max:
        return 0.0
    d_ref = d_min if d < d_min else d_max
    return 0.5 * k * (d - d_ref) ** 2


# ==========================================================
# TOLERANCE SETTINGS
# ==========================================================
tolerance_factor = 0.10
min_tolerance = 0.20


def compute_bounds(distances):
    bounds = {}
    for pair, d in distances.items():
        delta = max(d * tolerance_factor, min_tolerance)
        bounds[pair] = (round(d - delta, 2), round(d + delta, 2))
    return bounds


# ==========================================================
# K ASSIGNMENT RULES
# ==========================================================
def assign_k(pair_name):

    pair = pair_name.upper()

    # weak labile S-H
    if "S-H" in pair or "SG-HG1" in pair:
        return 0.5

    # moderate proton donor bonds
    if "O-H" in pair or "SE-H" in pair:
        return 5.0

    # strong heavy atom donor-acceptor distances
    return 10.0


# ==========================================================
# PLOTTING
# ==========================================================
def plot_and_save(restraints, stem, figsize=(4.6, 2.6)):
    fig, ax = plt.subplots(figsize=figsize)

    for name, d_min, d_max, k in restraints:
        x = np.linspace(d_min - 1.0, d_max + 1.0, 250)
        y = [flat_bottom_harmonic_potential(v, d_min, d_max, k) for v in x]

        # legend now ONLY shows pair name
        ax.plot(x, y, lw=1.0, label=name)

        ax.axvline(d_min, ls="--", lw=0.6, alpha=0.35, color="gray")
        ax.axvline(d_max, ls="--", lw=0.6, alpha=0.35, color="gray")

    ax.set_xlabel(r"Distance ($\AA$)")
    ax.set_ylabel(r"Potential Energy (kcal mol$^{-1}$)")
    ax.tick_params(labelsize=7)
    ax.grid(alpha=0.25)

    ax.legend(
        fontsize=5.0,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        framealpha=0.8
    )

    plt.tight_layout()
    plt.subplots_adjust(right=0.62)

    save_jcim_figure(fig, OUTDIR / f"{stem}.png")
    save_jcim_figure(fig, OUTDIR / f"{stem}.pdf")

    plt.show()


# ==========================================================
# FIGURE 1 : Sulfur full labels
# ==========================================================
distances1 = {
    "49.SG-196.O1": 3.40,
    "83.OE1-196.O1": 3.20,
    "49.SG-49.HG1": 1.50,
    "196.O1-196.H1": 1.00,
    "83.OE1-196.H1": 3.70,
    "196.O1-49.HG1": 2.80,
}

bounds1 = compute_bounds(distances1)

r1 = []
for pair in distances1:
    dmin, dmax = bounds1[pair]
    r1.append((pair, dmin, dmax, assign_k(pair)))

plot_and_save(r1, "bounds_fig1_sulfur_full")


# ==========================================================
# FIGURE 2 : Sulfur short labels
# ==========================================================
distances2 = {
    "S-O": 1.70,
    "O-O": 1.48,
    "S-H": 1.34,
    "O-H": 0.96,
}

bounds2 = compute_bounds(distances2)

r2 = []
for pair in distances2:
    dmin, dmax = bounds2[pair]
    r2.append((pair, dmin, dmax, assign_k(pair)))

plot_and_save(r2, "bounds_fig2_sulfur_short")


# ==========================================================
# FIGURE 3 : Selenium
# ==========================================================
distances3 = {
    "Se-O": 1.50,
    "O-O": 1.48,
    "Se-H": 1.47,
    "O-H": 0.96,
}

bounds3 = compute_bounds(distances3)

r3 = []
for pair in distances3:
    dmin, dmax = bounds3[pair]
    r3.append((pair, dmin, dmax, assign_k(pair)))

plot_and_save(r3, "bounds_fig3_selenium")


print(f"\nAll figures saved to: {OUTDIR}")
