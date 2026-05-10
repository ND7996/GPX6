import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# =================== JCIM PUBLICATION STYLE ===================
ACS_PATH = r"D:\PhD_Thesis\analysis\FINAL_PUBLICATION_FIGURES"
if ACS_PATH not in sys.path:
    sys.path.append(ACS_PATH)

from acsfonts import set_jcim_style
set_jcim_style()

# Force all text black
plt.rcParams.update({
    "text.color": "black",
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black"
})

# -------------------------------
# Load data
# -------------------------------
df = pd.read_csv(r'D:\PhD_Thesis\analysis\dipole\human_HBONDS.csv')

levels = sorted(df['level'].unique(), key=lambda x: int(x.replace('Level', '')))
levels = levels[:20]

TOL = 0.01

# -------------------------------
# Global Y range
# -------------------------------
all_vals = list(df['dg_star']) + list(df['dg0'])
y_min = min(all_vals) - 5
y_max = max(all_vals) + 5

# -------------------------------
# Figure layout
# -------------------------------
nrows = 5
ncols = 4

fig, axes = plt.subplots(
    nrows,
    ncols,
    figsize=(20,18),
    dpi=300
)

axes = axes.flatten()

# -------------------------------
# Plot levels
# -------------------------------
for idx, level in enumerate(levels):

    ax = axes[idx]

    level_data = df[df['level'] == level].copy()

    wt_row = level_data[level_data['mutation'] == 'U49C']
    wt_dg_star = wt_row['dg_star'].values[0] if len(wt_row) > 0 else None

    mutations = []
    dg_star_vals = []
    dg_star_errs = []
    dg0_vals = []
    dg0_errs = []
    colors_dg_star = []

    for _, row in level_data.iterrows():

        mut = row['mutation']
        dg_star = row['dg_star']

        mutations.append(mut)
        dg_star_vals.append(dg_star)
        dg_star_errs.append(row['dg_star_error'])
        dg0_vals.append(row['dg0'])
        dg0_errs.append(row['dg0_error'])

        if mut == 'U49C':
            colors_dg_star.append('red')

        elif wt_dg_star is not None:

            if dg_star < wt_dg_star - TOL:
                colors_dg_star.append('green')

            elif dg_star > wt_dg_star + TOL:
                colors_dg_star.append('blue')

            else:
                colors_dg_star.append('green')

        else:
            colors_dg_star.append('blue')

    x_pos = np.arange(len(mutations))
    width = 0.35

    # ΔG* bars
    ax.bar(
        x_pos - width/2,
        dg_star_vals,
        width,
        yerr=dg_star_errs,
        capsize=3,
        color=colors_dg_star,
        alpha=0.9,
        error_kw={'elinewidth':1,'capthick':1}
    )

    # ΔG0 bars
    ax.bar(
        x_pos + width/2,
        dg0_vals,
        width,
        yerr=dg0_errs,
        capsize=3,
        color='grey',
        alpha=0.6,
        error_kw={'elinewidth':1,'capthick':1}
    )

    # WT reference line
    if wt_dg_star is not None:
        ax.axhline(
            y=wt_dg_star,
            color='red',
            linestyle='--',
            linewidth=1.2,
            alpha=0.7
        )

    # Zero line
    ax.axhline(
        y=0,
        color='black',
        linestyle='-',
        linewidth=0.8,
        alpha=0.3
    )

    ax.set_xticks(x_pos)
    ax.set_xticklabels(
        mutations,
        rotation=65,
        ha='right'
    )

    ax.set_title(level)

    ax.set_ylim(y_min, y_max)

    ax.grid(True, alpha=0.25, linestyle='-', axis='y')
    ax.set_axisbelow(True)

    if idx % ncols == 0:
        ax.set_ylabel('Free Energy (kcal/mol)')

# -------------------------------
# Legend
# -------------------------------
legend_elements = [
    plt.Rectangle((0,0),1,1,fc='red',alpha=0.9,label='WT (U49C)'),
    plt.Rectangle((0,0),1,1,fc='green',alpha=0.9,label='ΔG* < WT'),
    plt.Rectangle((0,0),1,1,fc='blue',alpha=0.9,label='ΔG* > WT'),
    plt.Rectangle((0,0),1,1,fc='grey',alpha=0.6,label='ΔG0')
]

fig.legend(
    handles=legend_elements,
    loc='lower center',
    ncol=4,
    frameon=False
)

# -------------------------------
# Layout
# -------------------------------
fig.subplots_adjust(
    left=0.06,
    right=0.98,
    top=0.96,
    bottom=0.12,
    hspace=0.55,
    wspace=0.35
)

# -------------------------------
# Save
# -------------------------------
save_path = r'D:\PhD_Thesis\analysis\free_energy\mutation_human_all_levels.png'

plt.savefig(save_path, dpi=600, bbox_inches='tight', facecolor='white')

plt.show()

print(f"Saved: {save_path}")