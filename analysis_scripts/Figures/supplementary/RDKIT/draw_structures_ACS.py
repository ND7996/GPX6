import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe

# ── JCIM Style (Uniform Font Settings) ─────────────────────────────────────────
plt.rcParams.update({
    "font.family":       "Arial",
    "font.size":         9,
    "axes.labelsize":    9,
    "axes.titlesize":    9,
    "xtick.labelsize":   9,
    "ytick.labelsize":   9,
    "legend.fontsize":   9,
    "text.color":        "black",
    "axes.labelcolor":   "black",
    "xtick.color":       "black",
    "ytick.color":       "black",
    "figure.dpi":        600,
    "savefig.dpi":       600,
    "savefig.facecolor": "white",
    "axes.facecolor":    "white",
    "figure.facecolor":  "white",
})

# ── Parameters ────────────────────────────────────────────────────────────────
ACS_DPI       = 600
FIG_W         = 3.4   # Width for single panel
FIG_H         = 3.5   # Height for single panel
BOND_LW       = 1.2
BOND_SEP      = 0.055
ATOM_FS       = 9.0
ATOM_FS_SMALL = 7.5
FONT          = "Arial"
PANEL_LABEL_FS = 12
TITLE_FS      = 9

COLORS = {
    "nitrogen": "#1a6faf",
    "oxygen":   "#c0392b",
    "sulfur":   "#b7950b",
    "selenium": "#6c3483",
    "carbon":   "#2d2d2d",
    "hydrogen": "#7f8c8d",
}

# ── Drawing Helpers ──────────────────────────────────────────────────────────
def bond(ax, p1, p2, style="single", lw=None):
    """Draw a chemical bond between two points."""
    lw = lw or BOND_LW
    p1, p2 = np.array(p1, float), np.array(p2, float)
    if style == "single":
        ax.add_line(Line2D([p1[0], p2[0]], [p1[1], p2[1]], 
                          color="black", lw=lw, solid_capstyle="round", zorder=1))
    elif style == "double":
        d = p2 - p1
        L = np.linalg.norm(d)
        if L > 0:
            n = np.array([-d[1], d[0]]) / L * BOND_SEP
            for sgn in (-1, 1):
                off = sgn * n
                ax.add_line(Line2D([p1[0]+off[0], p2[0]+off[0]], 
                                  [p1[1]+off[1], p2[1]+off[1]], 
                                  color="black", lw=lw, solid_capstyle="round", zorder=1))

def label_atom(ax, xy, text, element="carbon", size="normal"):
    """Add atom label with white stroke outline for readability."""
    color = COLORS.get(element, "black")
    fs = ATOM_FS if size == "normal" else ATOM_FS_SMALL
    # White stroke (outline) effect
    ax.text(xy[0], xy[1], text, fontsize=fs, fontfamily=FONT, 
            ha="center", va="center", color="white", zorder=3,
            path_effects=[pe.withStroke(linewidth=4.0, foreground="white")])
    # Colored text on top
    ax.text(xy[0], xy[1], text, fontsize=fs, fontfamily=FONT, 
            ha="center", va="center", color=color, zorder=4)

def setup_ax(ax, xlim, ylim):
    """Set up axis with equal aspect and remove ticks."""
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

def add_title(ax, title):
    """Add centered title below the chemical structure."""
    ax.text(0.5, -0.15, title, transform=ax.transAxes,
            fontsize=TITLE_FS, fontfamily=FONT, ha="center", va="top", zorder=10)

# ── Panel Draw Functions ────────────────────────────────────────────────────
def draw_A(ax):
    """Cysteine (Cys) - with sulfur"""
    setup_ax(ax, (-1.8, 4.6), (-1.8, 2.8))
    
    # Atom positions
    N = np.array([0, 0])
    CA = np.array([1.2, 0])
    CB = np.array([1.95, 0.9])
    SG = np.array([3.3, 0.9])
    HN = np.array([-0.65, 0.5])
    HA = np.array([1.05, -0.7])
    HB1 = np.array([1.2, 1.65])
    HB2 = np.array([2.5, 1.65])
    HG1 = np.array([3.5, 0.2])
    
    # Bonds
    for p in [(HN, N), (N, CA), (CA, HA), (CA, CB), 
              (CB, HB1), (CB, HB2), (CB, SG), (SG, HG1)]:
        bond(ax, p[0], p[1])
    
    # Labels
    label_atom(ax, N, "N", "nitrogen")
    label_atom(ax, CA, "CA")
    label_atom(ax, CB, "CB")
    label_atom(ax, SG, "SG", "sulfur")
    label_atom(ax, HN, "HN", "hydrogen", "small")
    label_atom(ax, HA, "HA", "hydrogen", "small")
    label_atom(ax, HB1, "HB1", "hydrogen", "small")
    label_atom(ax, HB2, "HB2", "hydrogen", "small")
    label_atom(ax, HG1, "HG1", "hydrogen", "small")
    
    add_title(ax, "Cysteine (Cys)")

def draw_B(ax):
    """Selenocysteine (Sec) - with selenium"""
    setup_ax(ax, (-1.8, 4.8), (-1.8, 2.8))
    
    # Atom positions
    N = np.array([0, 0])
    CA = np.array([1.2, 0])
    CB = np.array([1.95, 0.9])
    SE = np.array([3.45, 0.9])
    HN = np.array([-0.65, 0.5])
    HA = np.array([1.05, -0.7])
    HB1 = np.array([1.2, 1.65])
    HB2 = np.array([2.5, 1.65])
    HG1 = np.array([3.65, 0.2])
    
    # Bonds
    for p in [(HN, N), (N, CA), (CA, HA), (CA, CB), 
              (CB, HB1), (CB, HB2), (CB, SE), (SE, HG1)]:
        bond(ax, p[0], p[1])
    
    # Labels
    label_atom(ax, N, "N", "nitrogen")
    label_atom(ax, CA, "CA")
    label_atom(ax, CB, "CB")
    label_atom(ax, SE, "SE", "selenium")
    label_atom(ax, HN, "HN", "hydrogen", "small")
    label_atom(ax, HA, "HA", "hydrogen", "small")
    label_atom(ax, HB1, "HB1", "hydrogen", "small")
    label_atom(ax, HB2, "HB2", "hydrogen", "small")
    label_atom(ax, HG1, "HG1", "hydrogen", "small")
    
    add_title(ax, "Selenocysteine (Sec)")

def draw_C(ax):
    """Hydrogen Peroxide (H₂O₂)"""
    setup_ax(ax, (-1.8, 4.0), (-1.2, 2.0))
    
    # Atom positions
    O1 = np.array([0.7, 0.4])
    O2 = np.array([1.8, 0.4])
    H1 = np.array([-0.05, 1.1])
    H2 = np.array([2.55, -0.3])
    
    # Bonds
    bond(ax, H1, O1)
    bond(ax, O1, O2)
    bond(ax, O2, H2)
    
    # Labels
    label_atom(ax, O1, "O1", "oxygen")
    label_atom(ax, O2, "O2", "oxygen")
    label_atom(ax, H1, "H1", "hydrogen", "small")
    label_atom(ax, H2, "H2", "hydrogen", "small")
    
    add_title(ax, "Hydrogen Peroxide")

def draw_D(ax):
    """Glutamine (Gln)"""
    setup_ax(ax, (-2.5, 2.5), (-3.0, 3.0))
    
    # Atom positions
    CA = np.array([0, 2])
    CB = np.array([0, 0.8])
    CG = np.array([0, -0.4])
    CD = np.array([0, -1.6])
    OE1 = np.array([1.2, -1.6])
    NE2 = np.array([-1.2, -1.6])
    HA = np.array([-0.6, 2.4])
    HB1 = np.array([-0.6, 0.4])
    HB2 = np.array([0.6, 0.4])
    HG1 = np.array([-0.6, -0.8])
    HG2 = np.array([0.6, -0.8])
    HE21 = np.array([-1.9, -1.0])
    HE22 = np.array([-1.9, -2.2])
    
    # Bonds
    bond(ax, CA, CB)
    bond(ax, CB, CG)
    bond(ax, CG, CD)
    bond(ax, CD, OE1, style="double")
    bond(ax, CD, NE2)
    bond(ax, CA, HA)
    bond(ax, CB, HB1)
    bond(ax, CB, HB2)
    bond(ax, CG, HG1)
    bond(ax, CG, HG2)
    bond(ax, NE2, HE21)
    bond(ax, NE2, HE22)
    
    # Labels
    label_atom(ax, CA, "CA")
    label_atom(ax, CB, "CB")
    label_atom(ax, CG, "CG")
    label_atom(ax, CD, "CD")
    label_atom(ax, OE1, "OE1", "oxygen")
    label_atom(ax, NE2, "NE2", "nitrogen")
    label_atom(ax, HA, "HA", "hydrogen", "small")
    label_atom(ax, HB1, "HB1", "hydrogen", "small")
    label_atom(ax, HB2, "HB2", "hydrogen", "small")
    label_atom(ax, HG1, "HG1", "hydrogen", "small")
    label_atom(ax, HG2, "HG2", "hydrogen", "small")
    label_atom(ax, HE21, "HE21", "hydrogen", "small")
    label_atom(ax, HE22, "HE22", "hydrogen", "small")
    
    add_title(ax, "Glutamine (Gln)")

# ── Main Execution: Save 4 Separate Figures ─────────────────────────────────
def main():
    # Create output directory
    outdir = "FigS2_output"
    os.makedirs(outdir, exist_ok=True)
    
    # Figure A: Cysteine
    fig_A, ax_A = plt.subplots(figsize=(FIG_W, FIG_H))
    draw_A(ax_A)
    plt.tight_layout(pad=0.5)
    path_A = os.path.join(outdir, "FigS2_A_Cysteine.png")
    fig_A.savefig(path_A, dpi=ACS_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig_A)
    print(f"Saved: {path_A}")
    
    # Figure B: Selenocysteine
    fig_B, ax_B = plt.subplots(figsize=(FIG_W, FIG_H))
    draw_B(ax_B)
    plt.tight_layout(pad=0.5)
    path_B = os.path.join(outdir, "FigS2_B_Selenocysteine.png")
    fig_B.savefig(path_B, dpi=ACS_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig_B)
    print(f"Saved: {path_B}")
    
    # Figure C: Hydrogen Peroxide
    fig_C, ax_C = plt.subplots(figsize=(FIG_W, FIG_H))
    draw_C(ax_C)
    plt.tight_layout(pad=0.5)
    path_C = os.path.join(outdir, "FigS2_C_HydrogenPeroxide.png")
    fig_C.savefig(path_C, dpi=ACS_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig_C)
    print(f"Saved: {path_C}")
    
    # Figure D: Glutamine
    fig_D, ax_D = plt.subplots(figsize=(FIG_W, FIG_H))
    draw_D(ax_D)
    plt.tight_layout(pad=0.5)
    path_D = os.path.join(outdir, "FigS2_D_Glutamine.png")
    fig_D.savefig(path_D, dpi=ACS_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig_D)
    print(f"Saved: {path_D}")
    
    print(f"\nAll 4 figures saved to: {outdir}/")

if __name__ == "__main__":
    main()