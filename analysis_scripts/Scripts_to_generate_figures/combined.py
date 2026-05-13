<<<<<<< HEAD
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from pathlib import Path
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# =================== ACS PUBLICATION STYLE ===================
ACS_PATH = r"./analysis_scripts/Scripts_to_generate_figures/Figures"
if ACS_PATH not in sys.path:
    sys.path.append(ACS_PATH)

from acsfonts import set_acs_style, acs_figure, save_acs_figure

set_acs_style()
sns.set_style("whitegrid", {'axes.linewidth': 0.8, 'grid.linewidth': 0.4})
set_acs_style()

# =================== OUTPUT FOLDER ===================
outdir = Path(r"./analysis_scripts/Scripts_to_generate_figures/Figures/Results\LFER")
outdir.mkdir(parents=True, exist_ok=True)

# =================== SEM HELPER ===================
def sem(series):
    n = series.dropna().shape[0]
    if n < 2:
        return np.nan
    return series.std(ddof=1) / np.sqrt(n)

# =================== AUTO-FIND DISTANCE COLUMN ===================
def find_dist_col(df):
    possible = [
        'dist_to_Sec/Cys49', 'distance_to_Cys49', 'dist_to_Cys49',
        'distance_Cys49', 'dist_C49', 'Distance_to_Cys49',
        'distance_to_Sec49', 'dist_to_Sec49', 'distance_Sec49',
    ]
    for col in possible:
        if col in df.columns:
            return col
    return None

# =================== DATA LOADING & Î”Î”G CALCULATION ===================
def load_and_compute(csv_path, ref_mutation):
    df = pd.read_csv(csv_path)
    for c in ['dg_star', 'dg_star_error', 'dg0', 'dg0_error']:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    df['level_num'] = df['level'].str.replace('Level', '', regex=False).astype(int)

    ref = df[df['mutation'].str.lower() == ref_mutation.lower()].set_index('level')
    if len(ref) == 0:
        raise ValueError(f"No '{ref_mutation}' reference found in {csv_path}!")
    fallback = ref.iloc[0]

    all_data = []
    for level in df['level'].unique():
        sub = df[df['level'] == level].copy()
        r = ref.loc[level] if level in ref.index else fallback
        sub['Î”Î”Gâ€¡'] = sub['dg_star'] - r['dg_star']
        sub['Î”Î”Gâ°']  = sub['dg0']     - r['dg0']
        sub['Î”Î”Gâ€¡_err'] = np.sqrt(sub['dg_star_error']**2 + r['dg_star_error']**2)
        sub['Î”Î”Gâ°_err']  = np.sqrt(sub['dg0_error']**2    + r['dg0_error']**2)
        all_data.append(sub)

    dfp = pd.concat(all_data)
    dfp = dfp[dfp['mutation'].str.lower() != ref_mutation.lower()].copy()

    dist_col = find_dist_col(dfp)
    if dist_col:
        dfp[dist_col] = pd.to_numeric(dfp[dist_col], errors='coerce')
    print(f"[{csv_path.name}] Distance column: '{dist_col}'")

    means_per_mut = dfp.groupby('mutation').agg(
        Î”Î”G0_mean           = ('Î”Î”Gâ°',     'mean'),
        Î”Î”G0_sem            = ('Î”Î”Gâ°',      sem),
        Î”Î”Gddagger_mean     = ('Î”Î”Gâ€¡',     'mean'),
        Î”Î”Gddagger_sem      = ('Î”Î”Gâ€¡',      sem),
        Î”Î”Gddagger_err_mean = ('Î”Î”Gâ€¡_err', 'mean'),
    ).reset_index()

    if dist_col:
        means_per_mut[dist_col] = dfp.groupby('mutation')[dist_col].mean().values

    mask_reg = dfp[['Î”Î”Gâ°', 'Î”Î”Gâ€¡']].notna().all(axis=1)
    x_all = dfp.loc[mask_reg, 'Î”Î”Gâ°'].values
    y_all = dfp.loc[mask_reg, 'Î”Î”Gâ€¡'].values
    slope, intercept, r_value, p_value, std_err_reg = stats.linregress(x_all, y_all)

    print(f"  LFER: Î±={slope:.3f}Â±{std_err_reg:.3f}, RÂ²={r_value**2:.3f}, p={p_value:.3e}, n={len(x_all)}")

    return dict(
        dfp=dfp, means_per_mut=means_per_mut, dist_col=dist_col,
        x_all=x_all, y_all=y_all,
        slope=slope, intercept=intercept, r_value=r_value,
        p_value=p_value, std_err=std_err_reg, n_pts=len(x_all)
    )

# =================== LOAD BOTH DATASETS ===================
print("\n--- Loading human ---")
H = load_and_compute(
    Path(r"./analysis_scripts\human_HBONDS.csv"),
    ref_mutation='humansec'
)
# also exclude humancys from human
H['dfp'] = H['dfp'][H['dfp']['mutation'].str.lower() != 'humancys'].copy()

print("\n--- Loading mouse ---")
M = load_and_compute(
    Path(r"./analysis_scripts\mouse_HBONDS.csv"),
    ref_mutation='mousecys'
)
# also exclude C49U from mouse if present
M['dfp'] = M['dfp'][M['dfp']['mutation'].str.lower() != 'c49u'].copy()

# =================== SHARED AXIS LIMITS ===================
# Collect all Î”Î”Gâ° and Î”Î”Gâ€¡ from both datasets
all_x = np.concatenate([H['x_all'], M['x_all']])
all_y = np.concatenate([H['y_all'], M['y_all']])

PAD = 1.0
shared_xlim = (all_x.min() - PAD, all_x.max() + PAD)
shared_ylim = (all_y.min() - PAD, all_y.max() + PAD)

# Shared fit line range
x_fit = np.linspace(shared_xlim[0], shared_xlim[1], 400)

# Shared distance colormap range (if both have dist_col)
dist_H = H['dist_col']
dist_M = M['dist_col']

all_dist_vals = []
if dist_H:
    all_dist_vals += H['dfp'][dist_H].dropna().tolist()
if dist_M:
    all_dist_vals += M['dfp'][dist_M].dropna().tolist()

if all_dist_vals:
    shared_vmin = min(all_dist_vals)
    shared_vmax = max(all_dist_vals)
    shared_norm = plt.Normalize(vmin=shared_vmin, vmax=shared_vmax)
    cmap_v = plt.cm.viridis
else:
    shared_norm = None
    cmap_v = None

print(f"\nShared Î”Î”Gâ° range: {shared_xlim[0]:.1f} to {shared_xlim[1]:.1f}")
print(f"Shared Î”Î”Gâ€¡ range: {shared_ylim[0]:.1f} to {shared_ylim[1]:.1f}")
if all_dist_vals:
    print(f"Shared distance range: {shared_vmin:.1f} to {shared_vmax:.1f} Ã…")

# =================== PLOTTING FUNCTION ===================
def plot_all_figures(data, label, ref_label, outdir, shared_xlim, shared_ylim,
                     x_fit, shared_norm, cmap_v):

    dfp           = data['dfp']
    means_per_mut = data['means_per_mut']
    dist_col      = data['dist_col']
    slope         = data['slope']
    intercept     = data['intercept']
    r_value       = data['r_value']
    p_value       = data['p_value']
    std_err       = data['std_err']
    n_pts         = data['n_pts']

    y_fit = slope * x_fit + intercept

    # Valid data for scatter
    if dist_col and shared_norm is not None:
        mask_valid = dfp[['Î”Î”Gâ°', 'Î”Î”Gâ€¡', dist_col]].notna().all(axis=1)
        dfp_valid  = dfp[mask_valid].copy()
        mutation_distances = dfp_valid.groupby('mutation')[dist_col].mean().sort_values()
    else:
        dfp_valid = dfp[dfp[['Î”Î”Gâ°', 'Î”Î”Gâ€¡']].notna().all(axis=1)].copy()
        mutation_distances = None

    # ------------------------------------------------------------------
    # FIGURE 1 â€” LFER
    # ------------------------------------------------------------------
    fig = plt.figure(figsize=(7, 3.5))
    gs  = fig.add_gridspec(
        1, 2,
        left=0.10, right=0.98, bottom=0.14, top=0.95,
        width_ratios=[4.5, 1], wspace=0.08
    )
    ax_main   = fig.add_subplot(gs[0])
    ax_legend = fig.add_subplot(gs[1])

    # Raw scatter
    ax_main.scatter(
        dfp_valid['Î”Î”Gâ°'], dfp_valid['Î”Î”Gâ€¡'],
        color='gray', edgecolors='none',
        s=6, alpha=0.30, zorder=2
    )

    # Regression line
    ax_main.plot(x_fit, y_fit, color='#d62728', lw=1.5, zorder=3, alpha=0.85)

    # Per-mutation means
    if dist_col and dist_col in means_per_mut.columns and shared_norm is not None:
        sc = ax_main.scatter(
            means_per_mut['Î”Î”G0_mean'], means_per_mut['Î”Î”Gddagger_mean'],
            c=means_per_mut[dist_col], cmap='viridis', norm=shared_norm,
            s=40, edgecolors='black', linewidth=0.6, alpha=0.95, zorder=6
        )
    else:
        sc = ax_main.scatter(
            means_per_mut['Î”Î”G0_mean'], means_per_mut['Î”Î”Gddagger_mean'],
            color='#4c72b0',
            s=40, edgecolors='black', linewidth=0.6, alpha=0.95, zorder=6
        )

    # Error bars
    ax_main.errorbar(
        means_per_mut['Î”Î”G0_mean'], means_per_mut['Î”Î”Gddagger_mean'],
        xerr=means_per_mut['Î”Î”G0_sem'],
        yerr=means_per_mut['Î”Î”Gddagger_err_mean'],
        fmt='none', elinewidth=0.8, capsize=2.5, capthick=0.8,
        color='black', alpha=0.65, zorder=5
    )

    # RÂ² annotation
    ax_main.text(
        0.04, 0.96,
        f'$R^2$ = {r_value**2:.2f}',
        transform=ax_main.transAxes,
        fontsize=7, fontweight='normal', va='top', ha='left',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                  edgecolor='black', linewidth=0.8, alpha=0.95)
    )

    ax_main.set_xlabel(r'$\Delta\Delta G^{\circ}$ (kcal mol$^{-1}$)', fontsize=8)
    ax_main.set_ylabel(r'$\Delta\Delta G^{\ddagger}$ (kcal mol$^{-1}$)', fontsize=8)
    ax_main.set_xlim(shared_xlim)
    ax_main.set_ylim(shared_ylim)

    if dist_col and shared_norm is not None:
        cbar = plt.colorbar(sc, ax=ax_main, pad=0.02, aspect=25, shrink=0.85)
        cbar.set_label('Distance to residue 49 (Ã…)', rotation=270, labelpad=12, fontsize=7)
        cbar.ax.tick_params(labelsize=6)

    # Legend panel
    ax_legend.axis('off')
    ax_legend.set_xlim(0, 1)
    ax_legend.set_ylim(0, 1)
    ax_legend.text(0.5, 0.99, 'Mutations', fontsize=7, ha='center', va='top')

    if dist_col and mutation_distances is not None and shared_norm is not None:
        y_pos = 0.93
        line_height = min(0.028, 0.88 / max(len(mutation_distances), 1))
        for mutation, distance in mutation_distances.items():
            color = cmap_v(shared_norm(distance))
            ax_legend.add_patch(
                plt.Rectangle((0.03, y_pos - 0.010), 0.12, 0.018,
                               facecolor=color, edgecolor='black', linewidth=0.5)
            )
            ax_legend.text(0.18, y_pos, mutation, fontsize=6, va='center', ha='left')
            ax_legend.text(0.98, y_pos, f'{distance:.1f} Ã…', fontsize=6,
                           va='center', ha='right', color='#444444')
            y_pos -= line_height
            if y_pos < 0.02:
                break

    legend_elements = [
        plt.scatter([], [], color='gray', edgecolors='none', s=6, alpha=0.40,
                    label='Per-level data'),
        plt.scatter([], [], color='white', edgecolors='black', s=40,
                    linewidth=0.6, label='Mutation mean $\\pm$ SEM'),
        plt.Line2D([0], [0], color='#d62728', lw=1.5, label='Linear fit'),
    ]
    ax_main.legend(handles=legend_elements, loc='lower right', fontsize=6,
                   framealpha=0.9, edgecolor='black', frameon=True,
                   handlelength=1.5, borderpad=0.5)

    save_acs_figure(fig, outdir / f"1_LFER_{label}")
    plt.show()
    print(f"Saved Figure 1 ({label})")

    # ------------------------------------------------------------------
    # FIGURE 2 â€” MEAN EFFECT PER EVOLUTIONARY LEVEL
    # ------------------------------------------------------------------
    fig, ax = acs_figure(width=3.3, height=2.5)

    sns.stripplot(
        data=dfp, x='level_num', y='Î”Î”Gâ€¡',
        color='#4c72b0', size=3.5, alpha=0.45,
        jitter=True, ax=ax, zorder=2
    )
    sns.pointplot(
        data=dfp, x='level_num', y='Î”Î”Gâ€¡',
        color='black', markers='D', linestyles='',
        errorbar='se', capsize=0.15,
        markersize=5, linewidth=1.0,
        ax=ax, zorder=4
    )
    ax.axhline(0, color='gray', lw=0.6, linestyle='--', alpha=0.6)
    ax.set_ylim(shared_ylim)
    ax.set_ylabel(r'$\Delta\Delta G^{\ddagger}$ (kcal mol$^{-1}$)', fontsize=8)
    ax.set_xlabel('Evolutionary level', fontsize=8)
    ax.tick_params(axis='both', labelsize=7)

    strip_patch  = mpatches.Patch(color='#4c72b0', alpha=0.5, label='Individual mutations')
    point_handle = plt.Line2D([0], [0], marker='D', color='black', lw=0,
                               markersize=4, label='Mean $\\pm$ SEM')
    ax.legend(handles=[strip_patch, point_handle], fontsize=6,
              loc='upper left', framealpha=0.9, edgecolor='black')

    fig.tight_layout()
    save_acs_figure(fig, outdir / f"2_Mean_by_level_{label}")
    plt.show()
    print(f"Saved Figure 2 ({label})")

    # ------------------------------------------------------------------
    # FIGURE 3 â€” DISTANCE VS Î”Î”Gâ€¡
    # ------------------------------------------------------------------
    if dist_col and shared_norm is not None:
        fig, ax = acs_figure(width=3.3, height=2.8)

        ax.scatter(
            dfp_valid[dist_col], dfp_valid['Î”Î”Gâ€¡'],
            color='gray', edgecolors='none', s=5, alpha=0.25, zorder=1
        )
        sc2 = ax.scatter(
            means_per_mut[dist_col], means_per_mut['Î”Î”Gddagger_mean'],
            c=means_per_mut[dist_col], cmap='viridis', norm=shared_norm,
            s=35, edgecolors='black', linewidth=0.5, alpha=0.92, zorder=4
        )
        ax.errorbar(
            means_per_mut[dist_col], means_per_mut['Î”Î”Gddagger_mean'],
            yerr=means_per_mut['Î”Î”Gddagger_err_mean'],
            fmt='none', elinewidth=0.8, capsize=2.5, capthick=0.8,
            color='black', alpha=0.60, zorder=3
        )

        try:
            from adjustText import adjust_text
            texts = [
                ax.text(row[dist_col], row['Î”Î”Gddagger_mean'], row['mutation'],
                        fontsize=5.5, ha='center')
                for _, row in means_per_mut.iterrows()
            ]
            adjust_text(texts, arrowprops=dict(arrowstyle='->', color='gray', lw=0.4))
        except ImportError:
            for _, row in means_per_mut.iterrows():
                ax.text(row[dist_col], row['Î”Î”Gddagger_mean'], ' ' + row['mutation'],
                        fontsize=5.5, ha='left', va='center',
                        bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                                  edgecolor='none', alpha=0.70))

        cbar2 = plt.colorbar(sc2, ax=ax, pad=0.02)
        cbar2.set_label('Distance to residue 49 (Ã…)', rotation=270, labelpad=12, fontsize=7)
        cbar2.ax.tick_params(labelsize=6)

        ax.set_ylim(shared_ylim)
        ax.set_xlabel('Distance to residue 49 (Ã…)', fontsize=8)
        ax.set_ylabel(r'$\Delta\Delta G^{\ddagger}$ (kcal mol$^{-1}$)', fontsize=8)
        ax.tick_params(axis='both', labelsize=7)

        fig.tight_layout()
        save_acs_figure(fig, outdir / f"3_Distance_vs_effect_{label}")
        plt.show()
        print(f"Saved Figure 3 ({label})")

    # ------------------------------------------------------------------
    # FIGURE 4 â€” HEATMAP
    # ------------------------------------------------------------------
    pivot = dfp.pivot_table(
        values='Î”Î”Gâ€¡', index='level_num', columns='mutation', aggfunc='mean'
    )

    # Shared colormap limits for heatmap
    heatmap_vmin = min(H['dfp']['Î”Î”Gâ€¡'].min(), M['dfp']['Î”Î”Gâ€¡'].min())
    heatmap_vmax = max(H['dfp']['Î”Î”Gâ€¡'].max(), M['dfp']['Î”Î”Gâ€¡'].max())
    heatmap_abs  = max(abs(heatmap_vmin), abs(heatmap_vmax))

    fig, ax = plt.subplots(figsize=(7, 3.0))
    set_acs_style()

    sns.heatmap(
        pivot, cmap='RdBu_r',
        vmin=-heatmap_abs, vmax=heatmap_abs,
        center=0, linewidths=0.3,
        cbar_kws={'label': r'$\Delta\Delta G^{\ddagger}$ (kcal mol$^{-1}$)', 'shrink': 0.8},
        ax=ax
    )
    ax.set_xlabel('Mutation', fontsize=8)
    ax.set_ylabel('Evolutionary level', fontsize=8)
    ax.tick_params(axis='both', labelsize=6)
    ax.collections[0].colorbar.ax.tick_params(labelsize=6)
    ax.collections[0].colorbar.set_label(
        r'$\Delta\Delta G^{\ddagger}$ (kcal mol$^{-1}$)', fontsize=7
    )

    fig.tight_layout()
    save_acs_figure(fig, outdir / f"4_Heatmap_{label}")
    plt.show()
    print(f"Saved Figure 4 ({label})")

    # ------------------------------------------------------------------
    # FIGURE 5 â€” PER-MUTATION SUMMARY BAR
    # ------------------------------------------------------------------
    mut_summary = dfp.groupby('mutation').agg(
        Î”Î”G0_mean       = ('Î”Î”Gâ°', 'mean'),
        Î”Î”Gddagger_mean = ('Î”Î”Gâ€¡', 'mean'),
        Î”Î”Gddagger_sem  = ('Î”Î”Gâ€¡',  sem),
        n               = ('Î”Î”Gâ€¡', 'count'),
    ).reset_index()
    mut_summary = mut_summary.sort_values('Î”Î”Gddagger_mean')
    mut_summary['ci95'] = (
        stats.t.ppf(0.975, df=(mut_summary['n'] - 1).clip(lower=1))
        * mut_summary['Î”Î”Gddagger_sem']
    )

    fig, ax = acs_figure(width=3.3, height=max(2.5, 0.22 * len(mut_summary)))

    colors = ['#d62728' if v > 0 else '#1f77b4' for v in mut_summary['Î”Î”Gddagger_mean']]
    ax.barh(
        mut_summary['mutation'], mut_summary['Î”Î”Gddagger_mean'],
        xerr=mut_summary['ci95'],
        color=colors, edgecolor='black', linewidth=0.5,
        error_kw=dict(elinewidth=0.8, capsize=2.5, capthick=0.8, ecolor='black'),
        alpha=0.80, height=0.65
    )
    ax.axvline(0, color='black', lw=0.8)
    ax.set_xlim(shared_ylim)   # same scale as Î”Î”Gâ€¡ axis
    ax.set_xlabel(r'$\Delta\Delta G^{\ddagger}$ (kcal mol$^{-1}$)', fontsize=8)
    ax.tick_params(axis='both', labelsize=6)

    pos_patch = mpatches.Patch(color='#d62728', alpha=0.8, label=r'$\Delta\Delta G^{\ddagger}$ > 0')
    neg_patch = mpatches.Patch(color='#1f77b4', alpha=0.8, label=r'$\Delta\Delta G^{\ddagger}$ < 0')
    ax.legend(handles=[pos_patch, neg_patch], fontsize=6,
              loc='lower right', framealpha=0.9, edgecolor='black')

    fig.tight_layout()
    save_acs_figure(fig, outdir / f"5_Mutation_summary_{label}")
    plt.show()
    print(f"Saved Figure 5 ({label})")

    # ------------------------------------------------------------------
    # SUMMARY STATS
    # ------------------------------------------------------------------
    print(f"\n{'='*60}")
    print(f"SUMMARY STATISTICS â€” {label.upper()}")
    print(f"{'='*60}")
    print(f"LFER slope Î±  = {slope:.3f} Â± {std_err:.3f}")
    print(f"RÂ²            = {r_value**2:.4f}")
    print(f"p-value       = {p_value:.3e}")
    print(f"n (data pts)  = {n_pts}")
    print(mut_summary[['mutation', 'Î”Î”G0_mean', 'Î”Î”Gddagger_mean',
                        'Î”Î”Gddagger_sem', 'ci95', 'n']].to_string(index=False))

# =================== RUN BOTH ===================
print("\n" + "="*60)
print("PLOTTING HUMAN")
print("="*60)
plot_all_figures(H, label='human', ref_label='humansec',
                 outdir=outdir,
                 shared_xlim=shared_xlim, shared_ylim=shared_ylim,
                 x_fit=x_fit, shared_norm=shared_norm, cmap_v=cmap_v)

print("\n" + "="*60)
print("PLOTTING MOUSE")
print("="*60)
plot_all_figures(M, label='mouse', ref_label='mousecys',
                 outdir=outdir,
                 shared_xlim=shared_xlim, shared_ylim=shared_ylim,
                 x_fit=x_fit, shared_norm=shared_norm, cmap_v=cmap_v)

=======
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from pathlib import Path
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# =================== ACS PUBLICATION STYLE ===================
ACS_PATH = r"./analysis_scripts/Scripts_to_generate_figures/Figures"
if ACS_PATH not in sys.path:
    sys.path.append(ACS_PATH)

from acsfonts import set_acs_style, acs_figure, save_acs_figure

set_acs_style()
sns.set_style("whitegrid", {'axes.linewidth': 0.8, 'grid.linewidth': 0.4})
set_acs_style()

# =================== OUTPUT FOLDER ===================
outdir = Path(r"./analysis_scripts/Scripts_to_generate_figures/Figures/Results\LFER")
outdir.mkdir(parents=True, exist_ok=True)

# =================== SEM HELPER ===================
def sem(series):
    n = series.dropna().shape[0]
    if n < 2:
        return np.nan
    return series.std(ddof=1) / np.sqrt(n)

# =================== AUTO-FIND DISTANCE COLUMN ===================
def find_dist_col(df):
    possible = [
        'dist_to_Sec/Cys49', 'distance_to_Cys49', 'dist_to_Cys49',
        'distance_Cys49', 'dist_C49', 'Distance_to_Cys49',
        'distance_to_Sec49', 'dist_to_Sec49', 'distance_Sec49',
    ]
    for col in possible:
        if col in df.columns:
            return col
    return None

# =================== DATA LOADING & Î”Î”G CALCULATION ===================
def load_and_compute(csv_path, ref_mutation):
    df = pd.read_csv(csv_path)
    for c in ['dg_star', 'dg_star_error', 'dg0', 'dg0_error']:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    df['level_num'] = df['level'].str.replace('Level', '', regex=False).astype(int)

    ref = df[df['mutation'].str.lower() == ref_mutation.lower()].set_index('level')
    if len(ref) == 0:
        raise ValueError(f"No '{ref_mutation}' reference found in {csv_path}!")
    fallback = ref.iloc[0]

    all_data = []
    for level in df['level'].unique():
        sub = df[df['level'] == level].copy()
        r = ref.loc[level] if level in ref.index else fallback
        sub['Î”Î”Gâ€¡'] = sub['dg_star'] - r['dg_star']
        sub['Î”Î”Gâ°']  = sub['dg0']     - r['dg0']
        sub['Î”Î”Gâ€¡_err'] = np.sqrt(sub['dg_star_error']**2 + r['dg_star_error']**2)
        sub['Î”Î”Gâ°_err']  = np.sqrt(sub['dg0_error']**2    + r['dg0_error']**2)
        all_data.append(sub)

    dfp = pd.concat(all_data)
    dfp = dfp[dfp['mutation'].str.lower() != ref_mutation.lower()].copy()

    dist_col = find_dist_col(dfp)
    if dist_col:
        dfp[dist_col] = pd.to_numeric(dfp[dist_col], errors='coerce')
    print(f"[{csv_path.name}] Distance column: '{dist_col}'")

    means_per_mut = dfp.groupby('mutation').agg(
        Î”Î”G0_mean           = ('Î”Î”Gâ°',     'mean'),
        Î”Î”G0_sem            = ('Î”Î”Gâ°',      sem),
        Î”Î”Gddagger_mean     = ('Î”Î”Gâ€¡',     'mean'),
        Î”Î”Gddagger_sem      = ('Î”Î”Gâ€¡',      sem),
        Î”Î”Gddagger_err_mean = ('Î”Î”Gâ€¡_err', 'mean'),
    ).reset_index()

    if dist_col:
        means_per_mut[dist_col] = dfp.groupby('mutation')[dist_col].mean().values

    mask_reg = dfp[['Î”Î”Gâ°', 'Î”Î”Gâ€¡']].notna().all(axis=1)
    x_all = dfp.loc[mask_reg, 'Î”Î”Gâ°'].values
    y_all = dfp.loc[mask_reg, 'Î”Î”Gâ€¡'].values
    slope, intercept, r_value, p_value, std_err_reg = stats.linregress(x_all, y_all)

    print(f"  LFER: Î±={slope:.3f}Â±{std_err_reg:.3f}, RÂ²={r_value**2:.3f}, p={p_value:.3e}, n={len(x_all)}")

    return dict(
        dfp=dfp, means_per_mut=means_per_mut, dist_col=dist_col,
        x_all=x_all, y_all=y_all,
        slope=slope, intercept=intercept, r_value=r_value,
        p_value=p_value, std_err=std_err_reg, n_pts=len(x_all)
    )

# =================== LOAD BOTH DATASETS ===================
print("\n--- Loading human ---")
H = load_and_compute(
    Path(r"./analysis_scripts\human_HBONDS.csv"),
    ref_mutation='humansec'
)
# also exclude humancys from human
H['dfp'] = H['dfp'][H['dfp']['mutation'].str.lower() != 'humancys'].copy()

print("\n--- Loading mouse ---")
M = load_and_compute(
    Path(r"./analysis_scripts\mouse_HBONDS.csv"),
    ref_mutation='mousecys'
)
# also exclude C49U from mouse if present
M['dfp'] = M['dfp'][M['dfp']['mutation'].str.lower() != 'c49u'].copy()

# =================== SHARED AXIS LIMITS ===================
# Collect all Î”Î”Gâ° and Î”Î”Gâ€¡ from both datasets
all_x = np.concatenate([H['x_all'], M['x_all']])
all_y = np.concatenate([H['y_all'], M['y_all']])

PAD = 1.0
shared_xlim = (all_x.min() - PAD, all_x.max() + PAD)
shared_ylim = (all_y.min() - PAD, all_y.max() + PAD)

# Shared fit line range
x_fit = np.linspace(shared_xlim[0], shared_xlim[1], 400)

# Shared distance colormap range (if both have dist_col)
dist_H = H['dist_col']
dist_M = M['dist_col']

all_dist_vals = []
if dist_H:
    all_dist_vals += H['dfp'][dist_H].dropna().tolist()
if dist_M:
    all_dist_vals += M['dfp'][dist_M].dropna().tolist()

if all_dist_vals:
    shared_vmin = min(all_dist_vals)
    shared_vmax = max(all_dist_vals)
    shared_norm = plt.Normalize(vmin=shared_vmin, vmax=shared_vmax)
    cmap_v = plt.cm.viridis
else:
    shared_norm = None
    cmap_v = None

print(f"\nShared Î”Î”Gâ° range: {shared_xlim[0]:.1f} to {shared_xlim[1]:.1f}")
print(f"Shared Î”Î”Gâ€¡ range: {shared_ylim[0]:.1f} to {shared_ylim[1]:.1f}")
if all_dist_vals:
    print(f"Shared distance range: {shared_vmin:.1f} to {shared_vmax:.1f} Ã…")

# =================== PLOTTING FUNCTION ===================
def plot_all_figures(data, label, ref_label, outdir, shared_xlim, shared_ylim,
                     x_fit, shared_norm, cmap_v):

    dfp           = data['dfp']
    means_per_mut = data['means_per_mut']
    dist_col      = data['dist_col']
    slope         = data['slope']
    intercept     = data['intercept']
    r_value       = data['r_value']
    p_value       = data['p_value']
    std_err       = data['std_err']
    n_pts         = data['n_pts']

    y_fit = slope * x_fit + intercept

    # Valid data for scatter
    if dist_col and shared_norm is not None:
        mask_valid = dfp[['Î”Î”Gâ°', 'Î”Î”Gâ€¡', dist_col]].notna().all(axis=1)
        dfp_valid  = dfp[mask_valid].copy()
        mutation_distances = dfp_valid.groupby('mutation')[dist_col].mean().sort_values()
    else:
        dfp_valid = dfp[dfp[['Î”Î”Gâ°', 'Î”Î”Gâ€¡']].notna().all(axis=1)].copy()
        mutation_distances = None

    # ------------------------------------------------------------------
    # FIGURE 1 â€” LFER
    # ------------------------------------------------------------------
    fig = plt.figure(figsize=(7, 3.5))
    gs  = fig.add_gridspec(
        1, 2,
        left=0.10, right=0.98, bottom=0.14, top=0.95,
        width_ratios=[4.5, 1], wspace=0.08
    )
    ax_main   = fig.add_subplot(gs[0])
    ax_legend = fig.add_subplot(gs[1])

    # Raw scatter
    ax_main.scatter(
        dfp_valid['Î”Î”Gâ°'], dfp_valid['Î”Î”Gâ€¡'],
        color='gray', edgecolors='none',
        s=6, alpha=0.30, zorder=2
    )

    # Regression line
    ax_main.plot(x_fit, y_fit, color='#d62728', lw=1.5, zorder=3, alpha=0.85)

    # Per-mutation means
    if dist_col and dist_col in means_per_mut.columns and shared_norm is not None:
        sc = ax_main.scatter(
            means_per_mut['Î”Î”G0_mean'], means_per_mut['Î”Î”Gddagger_mean'],
            c=means_per_mut[dist_col], cmap='viridis', norm=shared_norm,
            s=40, edgecolors='black', linewidth=0.6, alpha=0.95, zorder=6
        )
    else:
        sc = ax_main.scatter(
            means_per_mut['Î”Î”G0_mean'], means_per_mut['Î”Î”Gddagger_mean'],
            color='#4c72b0',
            s=40, edgecolors='black', linewidth=0.6, alpha=0.95, zorder=6
        )

    # Error bars
    ax_main.errorbar(
        means_per_mut['Î”Î”G0_mean'], means_per_mut['Î”Î”Gddagger_mean'],
        xerr=means_per_mut['Î”Î”G0_sem'],
        yerr=means_per_mut['Î”Î”Gddagger_err_mean'],
        fmt='none', elinewidth=0.8, capsize=2.5, capthick=0.8,
        color='black', alpha=0.65, zorder=5
    )

    # RÂ² annotation
    ax_main.text(
        0.04, 0.96,
        f'$R^2$ = {r_value**2:.2f}',
        transform=ax_main.transAxes,
        fontsize=7, fontweight='normal', va='top', ha='left',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                  edgecolor='black', linewidth=0.8, alpha=0.95)
    )

    ax_main.set_xlabel(r'$\Delta\Delta G^{\circ}$ (kcal mol$^{-1}$)', fontsize=8)
    ax_main.set_ylabel(r'$\Delta\Delta G^{\ddagger}$ (kcal mol$^{-1}$)', fontsize=8)
    ax_main.set_xlim(shared_xlim)
    ax_main.set_ylim(shared_ylim)

    if dist_col and shared_norm is not None:
        cbar = plt.colorbar(sc, ax=ax_main, pad=0.02, aspect=25, shrink=0.85)
        cbar.set_label('Distance to residue 49 (Ã…)', rotation=270, labelpad=12, fontsize=7)
        cbar.ax.tick_params(labelsize=6)

    # Legend panel
    ax_legend.axis('off')
    ax_legend.set_xlim(0, 1)
    ax_legend.set_ylim(0, 1)
    ax_legend.text(0.5, 0.99, 'Mutations', fontsize=7, ha='center', va='top')

    if dist_col and mutation_distances is not None and shared_norm is not None:
        y_pos = 0.93
        line_height = min(0.028, 0.88 / max(len(mutation_distances), 1))
        for mutation, distance in mutation_distances.items():
            color = cmap_v(shared_norm(distance))
            ax_legend.add_patch(
                plt.Rectangle((0.03, y_pos - 0.010), 0.12, 0.018,
                               facecolor=color, edgecolor='black', linewidth=0.5)
            )
            ax_legend.text(0.18, y_pos, mutation, fontsize=6, va='center', ha='left')
            ax_legend.text(0.98, y_pos, f'{distance:.1f} Ã…', fontsize=6,
                           va='center', ha='right', color='#444444')
            y_pos -= line_height
            if y_pos < 0.02:
                break

    legend_elements = [
        plt.scatter([], [], color='gray', edgecolors='none', s=6, alpha=0.40,
                    label='Per-level data'),
        plt.scatter([], [], color='white', edgecolors='black', s=40,
                    linewidth=0.6, label='Mutation mean $\\pm$ SEM'),
        plt.Line2D([0], [0], color='#d62728', lw=1.5, label='Linear fit'),
    ]
    ax_main.legend(handles=legend_elements, loc='lower right', fontsize=6,
                   framealpha=0.9, edgecolor='black', frameon=True,
                   handlelength=1.5, borderpad=0.5)

    save_acs_figure(fig, outdir / f"1_LFER_{label}")
    plt.show()
    print(f"Saved Figure 1 ({label})")

    # ------------------------------------------------------------------
    # FIGURE 2 â€” MEAN EFFECT PER EVOLUTIONARY LEVEL
    # ------------------------------------------------------------------
    fig, ax = acs_figure(width=3.3, height=2.5)

    sns.stripplot(
        data=dfp, x='level_num', y='Î”Î”Gâ€¡',
        color='#4c72b0', size=3.5, alpha=0.45,
        jitter=True, ax=ax, zorder=2
    )
    sns.pointplot(
        data=dfp, x='level_num', y='Î”Î”Gâ€¡',
        color='black', markers='D', linestyles='',
        errorbar='se', capsize=0.15,
        markersize=5, linewidth=1.0,
        ax=ax, zorder=4
    )
    ax.axhline(0, color='gray', lw=0.6, linestyle='--', alpha=0.6)
    ax.set_ylim(shared_ylim)
    ax.set_ylabel(r'$\Delta\Delta G^{\ddagger}$ (kcal mol$^{-1}$)', fontsize=8)
    ax.set_xlabel('Evolutionary level', fontsize=8)
    ax.tick_params(axis='both', labelsize=7)

    strip_patch  = mpatches.Patch(color='#4c72b0', alpha=0.5, label='Individual mutations')
    point_handle = plt.Line2D([0], [0], marker='D', color='black', lw=0,
                               markersize=4, label='Mean $\\pm$ SEM')
    ax.legend(handles=[strip_patch, point_handle], fontsize=6,
              loc='upper left', framealpha=0.9, edgecolor='black')

    fig.tight_layout()
    save_acs_figure(fig, outdir / f"2_Mean_by_level_{label}")
    plt.show()
    print(f"Saved Figure 2 ({label})")

    # ------------------------------------------------------------------
    # FIGURE 3 â€” DISTANCE VS Î”Î”Gâ€¡
    # ------------------------------------------------------------------
    if dist_col and shared_norm is not None:
        fig, ax = acs_figure(width=3.3, height=2.8)

        ax.scatter(
            dfp_valid[dist_col], dfp_valid['Î”Î”Gâ€¡'],
            color='gray', edgecolors='none', s=5, alpha=0.25, zorder=1
        )
        sc2 = ax.scatter(
            means_per_mut[dist_col], means_per_mut['Î”Î”Gddagger_mean'],
            c=means_per_mut[dist_col], cmap='viridis', norm=shared_norm,
            s=35, edgecolors='black', linewidth=0.5, alpha=0.92, zorder=4
        )
        ax.errorbar(
            means_per_mut[dist_col], means_per_mut['Î”Î”Gddagger_mean'],
            yerr=means_per_mut['Î”Î”Gddagger_err_mean'],
            fmt='none', elinewidth=0.8, capsize=2.5, capthick=0.8,
            color='black', alpha=0.60, zorder=3
        )

        try:
            from adjustText import adjust_text
            texts = [
                ax.text(row[dist_col], row['Î”Î”Gddagger_mean'], row['mutation'],
                        fontsize=5.5, ha='center')
                for _, row in means_per_mut.iterrows()
            ]
            adjust_text(texts, arrowprops=dict(arrowstyle='->', color='gray', lw=0.4))
        except ImportError:
            for _, row in means_per_mut.iterrows():
                ax.text(row[dist_col], row['Î”Î”Gddagger_mean'], ' ' + row['mutation'],
                        fontsize=5.5, ha='left', va='center',
                        bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                                  edgecolor='none', alpha=0.70))

        cbar2 = plt.colorbar(sc2, ax=ax, pad=0.02)
        cbar2.set_label('Distance to residue 49 (Ã…)', rotation=270, labelpad=12, fontsize=7)
        cbar2.ax.tick_params(labelsize=6)

        ax.set_ylim(shared_ylim)
        ax.set_xlabel('Distance to residue 49 (Ã…)', fontsize=8)
        ax.set_ylabel(r'$\Delta\Delta G^{\ddagger}$ (kcal mol$^{-1}$)', fontsize=8)
        ax.tick_params(axis='both', labelsize=7)

        fig.tight_layout()
        save_acs_figure(fig, outdir / f"3_Distance_vs_effect_{label}")
        plt.show()
        print(f"Saved Figure 3 ({label})")

    # ------------------------------------------------------------------
    # FIGURE 4 â€” HEATMAP
    # ------------------------------------------------------------------
    pivot = dfp.pivot_table(
        values='Î”Î”Gâ€¡', index='level_num', columns='mutation', aggfunc='mean'
    )

    # Shared colormap limits for heatmap
    heatmap_vmin = min(H['dfp']['Î”Î”Gâ€¡'].min(), M['dfp']['Î”Î”Gâ€¡'].min())
    heatmap_vmax = max(H['dfp']['Î”Î”Gâ€¡'].max(), M['dfp']['Î”Î”Gâ€¡'].max())
    heatmap_abs  = max(abs(heatmap_vmin), abs(heatmap_vmax))

    fig, ax = plt.subplots(figsize=(7, 3.0))
    set_acs_style()

    sns.heatmap(
        pivot, cmap='RdBu_r',
        vmin=-heatmap_abs, vmax=heatmap_abs,
        center=0, linewidths=0.3,
        cbar_kws={'label': r'$\Delta\Delta G^{\ddagger}$ (kcal mol$^{-1}$)', 'shrink': 0.8},
        ax=ax
    )
    ax.set_xlabel('Mutation', fontsize=8)
    ax.set_ylabel('Evolutionary level', fontsize=8)
    ax.tick_params(axis='both', labelsize=6)
    ax.collections[0].colorbar.ax.tick_params(labelsize=6)
    ax.collections[0].colorbar.set_label(
        r'$\Delta\Delta G^{\ddagger}$ (kcal mol$^{-1}$)', fontsize=7
    )

    fig.tight_layout()
    save_acs_figure(fig, outdir / f"4_Heatmap_{label}")
    plt.show()
    print(f"Saved Figure 4 ({label})")

    # ------------------------------------------------------------------
    # FIGURE 5 â€” PER-MUTATION SUMMARY BAR
    # ------------------------------------------------------------------
    mut_summary = dfp.groupby('mutation').agg(
        Î”Î”G0_mean       = ('Î”Î”Gâ°', 'mean'),
        Î”Î”Gddagger_mean = ('Î”Î”Gâ€¡', 'mean'),
        Î”Î”Gddagger_sem  = ('Î”Î”Gâ€¡',  sem),
        n               = ('Î”Î”Gâ€¡', 'count'),
    ).reset_index()
    mut_summary = mut_summary.sort_values('Î”Î”Gddagger_mean')
    mut_summary['ci95'] = (
        stats.t.ppf(0.975, df=(mut_summary['n'] - 1).clip(lower=1))
        * mut_summary['Î”Î”Gddagger_sem']
    )

    fig, ax = acs_figure(width=3.3, height=max(2.5, 0.22 * len(mut_summary)))

    colors = ['#d62728' if v > 0 else '#1f77b4' for v in mut_summary['Î”Î”Gddagger_mean']]
    ax.barh(
        mut_summary['mutation'], mut_summary['Î”Î”Gddagger_mean'],
        xerr=mut_summary['ci95'],
        color=colors, edgecolor='black', linewidth=0.5,
        error_kw=dict(elinewidth=0.8, capsize=2.5, capthick=0.8, ecolor='black'),
        alpha=0.80, height=0.65
    )
    ax.axvline(0, color='black', lw=0.8)
    ax.set_xlim(shared_ylim)   # same scale as Î”Î”Gâ€¡ axis
    ax.set_xlabel(r'$\Delta\Delta G^{\ddagger}$ (kcal mol$^{-1}$)', fontsize=8)
    ax.tick_params(axis='both', labelsize=6)

    pos_patch = mpatches.Patch(color='#d62728', alpha=0.8, label=r'$\Delta\Delta G^{\ddagger}$ > 0')
    neg_patch = mpatches.Patch(color='#1f77b4', alpha=0.8, label=r'$\Delta\Delta G^{\ddagger}$ < 0')
    ax.legend(handles=[pos_patch, neg_patch], fontsize=6,
              loc='lower right', framealpha=0.9, edgecolor='black')

    fig.tight_layout()
    save_acs_figure(fig, outdir / f"5_Mutation_summary_{label}")
    plt.show()
    print(f"Saved Figure 5 ({label})")

    # ------------------------------------------------------------------
    # SUMMARY STATS
    # ------------------------------------------------------------------
    print(f"\n{'='*60}")
    print(f"SUMMARY STATISTICS â€” {label.upper()}")
    print(f"{'='*60}")
    print(f"LFER slope Î±  = {slope:.3f} Â± {std_err:.3f}")
    print(f"RÂ²            = {r_value**2:.4f}")
    print(f"p-value       = {p_value:.3e}")
    print(f"n (data pts)  = {n_pts}")
    print(mut_summary[['mutation', 'Î”Î”G0_mean', 'Î”Î”Gddagger_mean',
                        'Î”Î”Gddagger_sem', 'ci95', 'n']].to_string(index=False))

# =================== RUN BOTH ===================
print("\n" + "="*60)
print("PLOTTING HUMAN")
print("="*60)
plot_all_figures(H, label='human', ref_label='humansec',
                 outdir=outdir,
                 shared_xlim=shared_xlim, shared_ylim=shared_ylim,
                 x_fit=x_fit, shared_norm=shared_norm, cmap_v=cmap_v)

print("\n" + "="*60)
print("PLOTTING MOUSE")
print("="*60)
plot_all_figures(M, label='mouse', ref_label='mousecys',
                 outdir=outdir,
                 shared_xlim=shared_xlim, shared_ylim=shared_ylim,
                 x_fit=x_fit, shared_norm=shared_norm, cmap_v=cmap_v)

>>>>>>> 52a0899df1d7ce42275b260475b162758c100d86
print(f"\nAll figures saved to: {outdir}")

