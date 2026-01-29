import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

def visualize_homolog_alignment(alignment_data, positions, title="Homolog 2"):
    """
    Visualize sequence alignment with colored boxes.
    Yellow boxes show the new additions at each step.
    
    Parameters:
    -----------
    alignment_data : list of lists
        Each inner list represents a sequence row with values:
        0 = white/empty, 1 = black/filled, 2 = yellow/new addition
    positions : list
        Position numbers to display at the top
    title : str
        Title for the plot
    """
    # Create figure and axis
    fig, ax = plt.subplots(figsize=(20, 12))
    
    # Define colors
    colors = {
        0: 'white',    # Empty box
        1: 'black',    # Previously filled box
        2: 'yellow'    # New addition at this step
    }
    
    n_rows = len(alignment_data)
    n_cols = len(alignment_data[0])
    
    # Draw boxes
    box_size = 0.85
    for row_idx, row in enumerate(alignment_data):
        for col_idx, value in enumerate(row):
            color = colors.get(value, 'white')
            
            # Draw filled box
            rect = mpatches.Rectangle(
                (col_idx, n_rows - row_idx - 1),
                box_size, box_size,
                facecolor=color,
                edgecolor='black',
                linewidth=0.8
            )
            ax.add_patch(rect)
    
    # Add position labels at top
    for idx, pos in enumerate(positions):
        ax.text(idx + 0.425, n_rows + 0.3, str(pos),
                ha='center', va='bottom', fontsize=8, fontweight='normal')
    
    # Set axis properties
    ax.set_xlim(-0.5, n_cols + 0.5)
    ax.set_ylim(-0.5, n_rows + 1.5)
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Add title
    ax.text(n_cols/2, n_rows + 1.2, title,
            ha='center', va='bottom', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    return fig


if __name__ == "__main__":
    # Position numbers from the image
    positions = [48, 52, 47, 99, 54, 177, 144, 74, 178, 143, 139, 173, 87, 142, 102, 24, 60, 181, 3, 104]
    
    # Create a mapping of position to column index
    pos_to_col = {pos: idx for idx, pos in enumerate(positions)}
    
    # Target positions for each row (newly selected in yellow)
    yellow_positions = [
        87,   # row 0
        173,  # row 1
        47,   # row 2
        143,  # row 3
        60,   # row 4
        104,  # row 5
        142,  # row 6
        139,  # row 7
        181,  # row 8
        48,   # row 9
        3,    # row 10
        54,   # row 11
        144,  # row 12
        177,  # row 13
        102,  # row 14
        24,   # row 15
        99,   # row 16
        52,   # row 17
        178,  # row 18
        74,   # row 19
    ]
    
    # Build alignment data with cumulative effect
    alignment_data = []
    selected_positions = set()  # Track all positions selected so far
    
    for row_idx in range(20):
        row = [0] * 20  # Start with all white
        
        # Mark all previously selected positions as black
        for pos in selected_positions:
            if pos in pos_to_col:
                col_idx = pos_to_col[pos]
                row[col_idx] = 1  # Black
        
        # Add the new yellow box for this row
        new_pos = yellow_positions[row_idx]
        if new_pos in pos_to_col:
            col_idx = pos_to_col[new_pos]
            row[col_idx] = 2  # Yellow
        
        # Add this position to the selected set for future rows
        selected_positions.add(new_pos)
        
        alignment_data.append(row)
    
    fig = visualize_homolog_alignment(alignment_data, positions)
    
    # Save the figure
    output_dir = r"D:\PhD_Thesis\Article-GPX6-EVB\Figures"
    
    # Create directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Save in multiple formats
    output_path_png = os.path.join(output_dir, "human_to_mouse_mutants.png")
    output_path_pdf = os.path.join(output_dir, "human_to_mouse_mutants.pdf")
    output_path_svg = os.path.join(output_dir, "human_to_mouse_mutants.svg")
    
    fig.savefig(output_path_png, dpi=300, bbox_inches='tight')
    fig.savefig(output_path_pdf, bbox_inches='tight')
    fig.savefig(output_path_svg, bbox_inches='tight')
    
    print(f"Figure saved to:")
    print(f"  - {output_path_png}")
    print(f"  - {output_path_pdf}")
    print(f"  - {output_path_svg}")
    
    plt.show()