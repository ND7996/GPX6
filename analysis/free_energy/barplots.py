import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the data
df = pd.read_csv(r'D:\PhD_Thesis\GPX6\analysis\human_mutants_FINAL_with_distances.csv')

# Get unique levels and sort them properly
levels = df['level'].unique()
levels = sorted(levels, key=lambda x: int(x.replace('Level', '')))

# Create figure with 6x4 subplots (publication-quality DPI + spacing)
fig, axes = plt.subplots(6, 4, figsize=(48, 72), dpi=100)
axes = axes.flatten()

# Process each level
for idx, level in enumerate(levels):
    if idx >= 24:
        break
    
    ax = axes[idx]
    level_data = df[df['level'] == level].copy()
    
    # Wild-type info
    wt_data = level_data[level_data['mutation'] == 'humansec']
    if len(wt_data) > 0:
        wt_dg_star = wt_data['dg_star'].values[0]
        wt_dg0 = wt_data['dg0'].values[0]
    else:
        wt_dg_star = None
        wt_dg0 = None
    
    # Prepare data
    mutations = []
    dg_star_vals = []
    dg_star_errs = []
    dg0_vals = []
    dg0_errs = []
    colors_dg_star = []
    
    for _, row in level_data.iterrows():
        mutations.append(row['mutation'])
        dg_star_vals.append(row['dg_star'])
        dg_star_errs.append(row['dg_star_error'])
        dg0_vals.append(row['dg0'])
        dg0_errs.append(row['dg0_error'])
        
        if row['mutation'] == 'humansec':
            colors_dg_star.append('red')
        elif wt_dg_star is not None:
            if row['dg_star'] > wt_dg_star:
                colors_dg_star.append('blue')
            else:
                colors_dg_star.append('green')
        else:
            colors_dg_star.append('red')
    
    x_pos = np.arange(len(mutations))
    width = 0.35
    
    # Bars
    ax.bar(x_pos - width/2, dg_star_vals, width,
           yerr=dg_star_errs, capsize=3, color=colors_dg_star,
           alpha=0.85, error_kw={'elinewidth': 1, 'capthick': 1})
    
    ax.bar(x_pos + width/2, dg0_vals, width,
           yerr=dg0_errs, capsize=3, color='grey',
           alpha=0.6, error_kw={'elinewidth': 1, 'capthick': 1})
    
    # WT marker line
    if wt_dg_star is not None:
        ax.axhline(y=wt_dg_star, color='red', linestyle='--',
                   linewidth=1.2, alpha=0.6)
    
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.8, alpha=0.3)
    
    # Labels & titles
    ax.set_xticks(x_pos)
    ax.set_xticklabels(mutations, rotation=90, ha='right', fontsize=12)
    ax.set_ylabel('Free Energy (kcal/mol)', fontsize=14, fontweight='bold')
    
    ax.set_title(f'{level}_barplot\nFree Energy Changes for Mouse → Human Mutations',
                 fontsize=16, fontweight='bold', pad=20)
    
    ax.grid(True, alpha=0.25, linestyle='-', axis='y')
    ax.set_axisbelow(True)
    
    y_min = min(min(dg0_vals), min(dg_star_vals)) - 10
    y_max = max(max(dg0_vals), max(dg_star_vals)) + 10
    ax.set_ylim(y_min, y_max)
    
    # Legend (only first subplot)
    if idx == 0:
        legend_elements = [
            plt.Rectangle((0, 0), 1, 1, fc='red', alpha=0.85, label='Threshold (humansec) ΔG* > 17.68 kcal/mol'),
            plt.Rectangle((0, 0), 1, 1, fc='blue', alpha=0.85, label='ΔG* > WT'),
            plt.Rectangle((0, 0), 1, 1, fc='green', alpha=0.85, label='ΔG* < WT'),
            plt.Rectangle((0, 0), 1, 1, fc='grey', alpha=0.6, label='Mean ΔG0')
        ]
        
        ax.legend(handles=legend_elements,
                  loc='upper center',
                  fontsize=14,
                  framealpha=0.95,
                  bbox_to_anchor=(1.10, -0.35))  # Place in white space
    
    ax.set_xlabel('Mutation', fontsize=14, fontweight='bold')
    ax.tick_params(axis='both', which='major', labelsize=12)

# Hide empty subplots
for idx in range(len(levels), 24):
    axes[idx].axis('off')

# Publication-quality spacing
fig.subplots_adjust(
    left=0.06,
    right=0.98,
    top=0.96,
    bottom=0.05,
    hspace=0.75,
    wspace=0.45
)

# Save (high DPI)
plt.savefig(
    r'D:\PhD_Thesis\GPX6\analysis\figures\mutation_plots_human.png',
    dpi=300,
    bbox_inches='tight',
    facecolor='white'
)

plt.show()

print("Plot saved successfully!")
print(f"Total levels plotted: {len(levels)}")
print(f"Levels found: {levels}")
