"""
acs_fonts.py  —  JCIM publication style
Place this file in your figures folder and import it in every plotting script:

    from acs_fonts import set_jcim_style, jcim_figure, jcim_subplots, save_jcim_figure
"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def set_jcim_style():
    """
    Apply JCIM-compliant rcParams globally.
    Call ONCE at the top of every script, before creating any figure.
    Do NOT set rcParams anywhere else.
    """
    plt.rcParams.update({
        # ── Fonts ──────────────────────────────────────────────────────────
        "font.family":           "sans-serif",
        "font.sans-serif":       ["Arial", "Helvetica", "Liberation Sans"],
        "font.size":             8,
        "axes.labelsize":        8,
        "axes.titlesize":        8,
        "xtick.labelsize":       8,
        "ytick.labelsize":       8,
        "legend.fontsize":       8,
        "text.usetex":           False,

        # ── Ticks ──────────────────────────────────────────────────────────
        "xtick.direction":       "in",
        "ytick.direction":       "in",
        "xtick.major.width":     0.75,
        "ytick.major.width":     0.75,
        "xtick.major.size":      4,
        "ytick.major.size":      4,
        "xtick.minor.size":      2,
        "ytick.minor.size":      2,
        "xtick.minor.visible":   True,
        "ytick.minor.visible":   True,

        # ── Lines & axes ───────────────────────────────────────────────────
        "axes.linewidth":        0.75,
        "axes.edgecolor":        "black",
        "lines.linewidth":       1.0,
        "lines.markersize":      4,
        "lines.markeredgewidth": 0.5,

        # ── Legend ─────────────────────────────────────────────────────────
        "legend.frameon":        False,
        "legend.handlelength":   1.5,
        "legend.handletextpad":  0.8,

        # ── Spines / grid ──────────────────────────────────────────────────
        "axes.spines.top":       False,
        "axes.spines.right":     False,
        "axes.grid":             False,
        "grid.alpha":            0.3,

        # ── Output quality ─────────────────────────────────────────────────
        "figure.dpi":            300,
        "savefig.dpi":           600,
        "savefig.bbox":          "tight",
        "savefig.pad_inches":    0.05,

        # ── Font embedding (required for PDF submission) ────────────────────
        "pdf.fonttype":          42,
        "ps.fonttype":           42,

        # ── Colour cycle ───────────────────────────────────────────────────
        "axes.prop_cycle": plt.cycler(color=[
            "#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd",
            "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
        ]),

        # ── Image defaults ─────────────────────────────────────────────────
        "image.cmap":            "viridis",
        "image.interpolation":   "none",
        "image.aspect":          "auto",
    })


# ── Convenience constructors ────────────────────────────────────────────────

def jcim_figure(width=3.25, height=2.5):
    """
    Single-panel JCIM figure.
    Default: 3.25 in wide (single column).
    """
    set_jcim_style()
    fig, ax = plt.subplots(figsize=(width, height))
    return fig, ax


def jcim_subplots(nrows=1, ncols=1, width=3.25, height=2.5, **kwargs):
    """
    Multi-panel JCIM figure.
    Total figure size = (width * ncols, height * nrows).
    """
    set_jcim_style()
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(width * ncols, height * nrows),
        **kwargs,
    )
    return fig, axes


def jcim_gridspec(fig, nrows, ncols, width_ratios=None, **kwargs):
    """
    Return a GridSpec attached to *fig* with JCIM defaults.
    Pass width_ratios / height_ratios / wspace / hspace as needed.
    """
    return gridspec.GridSpec(
        nrows, ncols,
        width_ratios=width_ratios,
        **kwargs,
    )


# ── Save helper ─────────────────────────────────────────────────────────────

def save_jcim_figure(fig, filename, formats=None, fmt=None):
    """
    Save *fig* to one or more formats.

    Parameters
    ----------
    fig      : matplotlib Figure
    filename : str or Path  — output path WITHOUT extension
    formats  : list of str  — default ['png', 'pdf']
    """
    if formats is None and fmt is not None:
        formats = fmt if isinstance(fmt, list) else [fmt]
    if formats is None:
        formats = ["png", "pdf"]

    _params = {
        "png":  {"dpi": 600, "bbox_inches": "tight"},
        "pdf":  {"dpi": 600, "bbox_inches": "tight"},
        "eps":  {"bbox_inches": "tight"},
        "tiff": {"dpi": 600, "bbox_inches": "tight"},
        "svg":  {"bbox_inches": "tight"},
    }

    for fmt in formats:
        if fmt in _params:
            out = f"{filename}.{fmt}"
            fig.savefig(out, **_params[fmt])
            print(f"Saved: {out}")
        else:
            print(f"Warning: unsupported format '{fmt}'")