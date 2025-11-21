import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

def visualize_homolog_alignment(alignment_data, positions, title="Homolog 1"):
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
        0: 'white',      # Empty box
        1: 'black',      # Previously filled box
        2: 'yellow'      # New addition at this step
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

# Data extracted from the image - corrected based on cumulative additions
if __name__ == "__main__":
    # Position numbers from the image
    positions = [48, 52, 47, 99, 54, 177, 144, 74, 178, 143, 139, 104, 87, 142, 
                 102, 24, 60, 181, 4, 107]
    
    # Alignment data matching the image (0=white, 1=black, 2=yellow)
    # Yellow shows NEW additions at each row
    alignment_data = [
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  # Row 0: empty
        [0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  # Row 1: +T54Q
        [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0],  # Row 2: +T54Q + mousesec
        [0,0,0,0,1,0,0,0,0,0,2,0,0,0,0,1,0,0,0,0],  # Row 3: +I24L
        [0,0,2,0,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0],
        [0,0,1,0,1,0,0,0,0,0,1,0,0,0,0,1,2,0,0,0],  # Row 5: +S47A
        [0,0,1,0,1,0,0,2,0,0,1,0,0,0,0,1,1,2,0,0],  # Row 6: +T60A
        [0,0,0,1,1,1,0,2,0,0,0,0,0,0,1,0,1,1,0,0],  # Row 7: +G74A
        [0,0,0,1,1,1,0,1,0,0,0,0,0,0,1,0,1,1,2,0],  # Row 8: +H144Q
        [0,0,0,1,1,1,0,1,0,0,0,0,0,0,1,0,1,1,1,2],  # Row 9: +R89C
        [0,0,0,1,1,1,0,1,0,0,0,0,2,0,1,0,1,1,1,1],  # Row 10: +H177Q
        [0,0,0,1,1,1,0,1,0,2,0,0,1,0,1,0,1,1,1,1],  # Row 11: +T178A
        [0,0,0,1,1,1,0,1,0,1,0,0,1,0,1,0,1,1,1,2],  # Row 12: +Y104F
        [0,0,0,1,1,1,0,1,0,1,0,0,1,0,1,2,1,1,1,1],  # Row 13: +S4R
        [0,0,0,1,1,1,0,1,0,1,0,0,1,0,1,1,1,1,1,2],  # Row 14: +T52A
        [0,0,0,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1],  # Row 15: +K87T
        [0,0,0,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1],  # Row 16: +G102S
        [0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2],  # Row 17: +N107S
        [0,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1],  # Row 18: +P142S
        [0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1],  # Row 19: +R181S
        [2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],  # Row 20: +F48Y
    ]
    
    fig = visualize_homolog_alignment(alignment_data, positions)
    plt.show()