import os
import csv
import numpy as np
import pandas as pd
import matplotlib
# Force matplotlib to use a non-GUI backend to avoid Qt issues
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LinearSegmentedColormap

def visualize_all_dipoles(summary_file, output_file, plot_type="3d"):
    """
    Create a comprehensive visualization of all dipoles
    
    Parameters:
    summary_file - Path to the dipole summary CSV
    output_file - Path to save the visualization
    plot_type - Type of plot: "3d", "2d", "rosette", or "all"
    """
    print(f"Creating combined visualization of all dipoles from {summary_file}")
    
    # Load summary data
    if not os.path.exists(summary_file):
        print(f"ERROR: Summary file not found: {summary_file}")
        return
    
    df = pd.read_csv(summary_file)
    
    if len(df) == 0:
        print("No dipole data to visualize!")
        return
    
    # Determine the plot type to generate
    if plot_type == "all":
        # Create a figure with multiple subplots
        fig = plt.figure(figsize=(20, 16))
        
        # 3D plot
        ax1 = fig.add_subplot(221, projection='3d')
        create_3d_dipole_plot(ax1, df)
        
        # 2D plots for principal planes
        ax2 = fig.add_subplot(222)
        create_2d_dipole_plot(ax2, df, plane="xy")
        
        ax3 = fig.add_subplot(223)
        create_2d_dipole_plot(ax3, df, plane="xz")
        
        ax4 = fig.add_subplot(224)
        create_2d_dipole_plot(ax4, df, plane="yz")
        
        plt.tight_layout()
        
    elif plot_type == "3d":
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')
        create_3d_dipole_plot(ax, df)
        
    elif plot_type == "2d":
        # Create three 2D plots for the principal planes
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        create_2d_dipole_plot(axes[0], df, plane="xy")
        create_2d_dipole_plot(axes[1], df, plane="xz")
        create_2d_dipole_plot(axes[2], df, plane="yz")
        plt.tight_layout()
        
    elif plot_type == "rosette":
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='polar')
        create_dipole_rosette(ax, df)
    
    # Save the plot
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Combined visualization saved to {output_file}")
    plt.close(fig)  # Explicitly close the figure to free memory

def create_3d_dipole_plot(ax, df):
    """Create a 3D plot of all dipole vectors"""
    # Create origin point (0,0,0)
    origin = np.zeros((len(df), 3))
    
    # Create colormap based on magnitudes
    magnitudes = df['dipole_magnitude'].values
    norm = plt.Normalize(magnitudes.min(), magnitudes.max())
    cmap = plt.get_cmap('viridis')
    colors = cmap(norm(magnitudes))
    
    # Plot vectors
    for i, (_, row) in enumerate(df.iterrows()):
        x, y, z = row['dipole_x'], row['dipole_y'], row['dipole_z']
        replica = row['replica']
        magnitude = row['dipole_magnitude']
        
        # Plot the dipole vector
        ax.quiver(0, 0, 0, x, y, z, color=colors[i], alpha=0.8, 
                 linewidth=2, arrow_length_ratio=0.1)
        
        # Add text label for replica number at tip of arrow
        ax.text(x*1.1, y*1.1, z*1.1, replica, color=colors[i], fontsize=8)
    
    # Calculate the maximum magnitude for setting axis limits
    max_mag = df['dipole_magnitude'].max() * 1.2
    
    # Set axis limits
    ax.set_xlim([-max_mag, max_mag])
    ax.set_ylim([-max_mag, max_mag])
    ax.set_zlim([-max_mag, max_mag])
    
    # Add labels and title
    ax.set_xlabel('X', fontsize=12)
    ax.set_ylabel('Y', fontsize=12)
    ax.set_zlabel('Z', fontsize=12)
    ax.set_title('3D Visualization of All Dipole Vectors', fontsize=14)
    
    # Add a colorbar
    scalar_map = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar_map.set_array([])
    cbar = plt.colorbar(scalar_map, ax=ax, label='Dipole Magnitude')
    
    # Add a reference sphere
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    radius = max_mag * 0.5  # Sphere radius
    x_sphere = radius * np.cos(u) * np.sin(v)
    y_sphere = radius * np.sin(u) * np.sin(v)
    z_sphere = radius * np.cos(v)
    ax.plot_wireframe(x_sphere, y_sphere, z_sphere, color='gray', alpha=0.2)

def create_2d_dipole_plot(ax, df, plane="xy"):
    """Create a 2D plot of dipole vectors in the specified plane"""
    # Map plane to axes indices
    plane_map = {
        "xy": (0, 1, "X", "Y", "XY Plane"),
        "xz": (0, 2, "X", "Z", "XZ Plane"),
        "yz": (1, 2, "Y", "Z", "YZ Plane")
    }
    
    idx1, idx2, label1, label2, title = plane_map.get(plane, (0, 1, "X", "Y", "XY Plane"))
    
    # Extract components based on plane
    components = []
    for _, row in df.iterrows():
        if plane == "xy":
            components.append((row['dipole_x'], row['dipole_y']))
        elif plane == "xz":
            components.append((row['dipole_x'], row['dipole_z']))
        else:  # yz
            components.append((row['dipole_y'], row['dipole_z']))
    
    components = np.array(components)
    
    # Create colormap based on magnitudes
    magnitudes = df['dipole_magnitude'].values
    norm = plt.Normalize(magnitudes.min(), magnitudes.max())
    cmap = plt.get_cmap('viridis')
    colors = cmap(norm(magnitudes))
    
    # Plot vectors
    for i, (_, row) in enumerate(df.iterrows()):
        x, y = components[i]
        replica = row['replica']
        
        # Plot the vector
        ax.quiver(0, 0, x, y, color=colors[i], scale=1, scale_units='xy', 
                 angles='xy', width=0.005, alpha=0.8)
        
        # Add text label for replica number
        ax.text(x*1.05, y*1.05, replica, color=colors[i], fontsize=8)
    
    # Calculate the maximum value for setting axis limits
    max_val = max(abs(components.min()), abs(components.max())) * 1.2
    
    # Set axis limits
    ax.set_xlim([-max_val, max_val])
    ax.set_ylim([-max_val, max_val])
    
    # Add a circle for reference
    circle = plt.Circle((0, 0), max_val*0.5, fill=False, color='gray', linestyle='--', alpha=0.5)
    ax.add_artist(circle)
    
    # Add gridlines and axes
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.axhline(y=0, color='k', linestyle='-', alpha=0.3)
    ax.axvline(x=0, color='k', linestyle='-', alpha=0.3)
    
    # Add labels and title
    ax.set_xlabel(label1, fontsize=12)
    ax.set_ylabel(label2, fontsize=12)
    ax.set_title(f'Dipole Vectors: {title}', fontsize=14)
    
    # Make axes equal for proper vector representation
    ax.set_aspect('equal')

def create_dipole_rosette(ax, df):
    """Create a polar plot (rosette) of dipole directions"""
    # Calculate theta (azimuthal angle) and r (magnitude)
    theta = []
    r = []
    
    for _, row in df.iterrows():
        x, y = row['dipole_x'], row['dipole_y']
        magnitude = row['dipole_magnitude']
        
        # Calculate azimuthal angle
        angle = np.arctan2(y, x)
        if angle < 0:
            angle += 2 * np.pi  # Convert to [0, 2π]
            
        theta.append(angle)
        r.append(magnitude)
    
    # Create colormap based on magnitude
    norm = plt.Normalize(min(r), max(r))
    cmap = plt.get_cmap('viridis')
    colors = cmap(norm(r))
    
    # Plot points
    ax.scatter(theta, r, c=colors, s=80, alpha=0.7)
    
    # Add replica labels
    for i, (angle, rad) in enumerate(zip(theta, r)):
        replica = df.iloc[i]['replica']
        ax.text(angle, rad*1.05, replica, fontsize=8, 
               ha='center', va='center')
    
    # Add connecting lines to origin
    for angle, rad in zip(theta, r):
        ax.plot([angle, angle], [0, rad], color='gray', alpha=0.3)
    
    # Set up the plot
    ax.set_theta_zero_location('E')  # 0 degrees points east
    ax.set_theta_direction(-1)  # clockwise
    ax.set_title('Dipole Direction Rosette (XY Plane)', fontsize=14)
    ax.set_rticks([])  # Hide radial ticks
    
    # Add cardinal direction labels
    cardinal_dirs = ['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE']
    angles = np.linspace(0, 2*np.pi, len(cardinal_dirs), endpoint=False)
    ax.set_xticks(angles)
    ax.set_xticklabels(cardinal_dirs)
    
    # Add a colorbar
    scalar_map = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar_map.set_array([])
    plt.colorbar(scalar_map, ax=ax, label='Dipole Magnitude')

def create_heatmap_visualization(summary_file, output_file):
    """Create a heatmap visualization of dipole components"""
    print(f"Creating heatmap visualization from {summary_file}")
    
    # Load summary data
    if not os.path.exists(summary_file):
        print(f"ERROR: Summary file not found: {summary_file}")
        return
    
    df = pd.read_csv(summary_file)
    
    if len(df) == 0:
        print("No dipole data to visualize!")
        return
    
    # Sort by replica number (assuming replica is numeric or can be converted)
    # First check if replica is numeric
    try:
        df['replica_num'] = df['replica'].astype(int)
    except:
        # If not numeric, convert to string and zero pad
        df['replica_num'] = df['replica'].astype(str).str.zfill(3)
    
    df = df.sort_values('replica_num')
    
    # Extract components for the heatmap
    components = df[['dipole_x', 'dipole_y', 'dipole_z', 'dipole_magnitude']].values
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Create heatmap
    im = ax.imshow(components, cmap='coolwarm', aspect='auto')
    
    # Add labels
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df['replica'])
    ax.set_xticks(range(4))
    ax.set_xticklabels(['X', 'Y', 'Z', 'Magnitude'])
    
    # Add colorbar
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('Component Value', rotation=-90, va="bottom")
    
    # Add title
    ax.set_title('Dipole Component Heatmap by Replica', fontsize=14)
    
    # Add grid lines
    ax.set_xticks(np.arange(-.5, 4, 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(df), 1), minor=True)
    ax.grid(which="minor", color="black", linestyle='-', linewidth=1)
    
    # Add values in each cell
    for i in range(len(df)):
        for j in range(4):
            text = ax.text(j, i, f"{components[i, j]:.2f}",
                          ha="center", va="center", color="black")
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Heatmap visualization saved to {output_file}")
    plt.close(fig)

def create_statistical_visualization(summary_file, output_file):
    """Create statistical visualizations of dipole data"""
    print(f"Creating statistical visualization from {summary_file}")
    
    # Load summary data
    if not os.path.exists(summary_file):
        print(f"ERROR: Summary file not found: {summary_file}")
        return
    
    df = pd.read_csv(summary_file)
    
    if len(df) == 0:
        print("No dipole data to visualize!")
        return
    
    # Create figure with multiple subplots
    fig = plt.figure(figsize=(18, 12))
    
    # 1. Histogram of magnitudes
    ax1 = fig.add_subplot(221)
    ax1.hist(df['dipole_magnitude'], bins=10, color='skyblue', edgecolor='black')
    ax1.set_xlabel('Dipole Magnitude')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Distribution of Dipole Magnitudes')
    ax1.grid(True, linestyle='--', alpha=0.7)
    
    # Add mean and median lines
    mean_mag = df['dipole_magnitude'].mean()
    median_mag = df['dipole_magnitude'].median()
    ax1.axvline(mean_mag, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_mag:.2f}')
    ax1.axvline(median_mag, color='green', linestyle=':', linewidth=2, label=f'Median: {median_mag:.2f}')
    ax1.legend()
    
    # 2. Box plots of x, y, z components
    ax2 = fig.add_subplot(222)
    component_data = [df['dipole_x'], df['dipole_y'], df['dipole_z']]
    ax2.boxplot(component_data, labels=['X', 'Y', 'Z'])
    ax2.set_ylabel('Component Value')
    ax2.set_title('Distribution of Dipole Components')
    ax2.grid(True, linestyle='--', alpha=0.7)
    
    # 3. Scatter plot of X vs Y components
    ax3 = fig.add_subplot(223)
    scatter = ax3.scatter(df['dipole_x'], df['dipole_y'], c=df['dipole_z'], 
                         cmap='viridis', s=df['dipole_magnitude']*20, alpha=0.7)
    ax3.set_xlabel('X Component')
    ax3.set_ylabel('Y Component')
    ax3.set_title('X vs Y Components (Size=Magnitude, Color=Z)')
    ax3.grid(True, linestyle='--', alpha=0.7)
    ax3.axhline(y=0, color='k', linestyle='-', alpha=0.3)
    ax3.axvline(x=0, color='k', linestyle='-', alpha=0.3)
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax3)
    cbar.set_label('Z Component')
    
    # 4. Component contributions
    ax4 = fig.add_subplot(224)
    
    # Calculate absolute component contributions
    df['abs_x'] = np.abs(df['dipole_x'])
    df['abs_y'] = np.abs(df['dipole_y'])
    df['abs_z'] = np.abs(df['dipole_z'])
    total_abs = df['abs_x'] + df['abs_y'] + df['abs_z']
    df['x_contrib'] = df['abs_x'] / total_abs * 100
    df['y_contrib'] = df['abs_y'] / total_abs * 100
    df['z_contrib'] = df['abs_z'] / total_abs * 100
    
    # Average contributions
    avg_contribs = [df['x_contrib'].mean(), df['y_contrib'].mean(), df['z_contrib'].mean()]
    
    ax4.bar(['X', 'Y', 'Z'], avg_contribs, color=['red', 'green', 'blue'], alpha=0.7)
    ax4.set_ylabel('Average Contribution (%)')
    ax4.set_title('Average Component Contributions to Dipole')
    ax4.grid(True, linestyle='--', alpha=0.7)
    
    # Add percentage labels
    for i, v in enumerate(avg_contribs):
        ax4.text(i, v + 1, f"{v:.1f}%", ha='center')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Statistical visualization saved to {output_file}")
    plt.close(fig)

def create_directional_analysis(summary_file, output_file):
    """Create a visualization focusing on dipole directions"""
    print(f"Creating directional analysis from {summary_file}")
    
    # Load summary data
    if not os.path.exists(summary_file):
        print(f"ERROR: Summary file not found: {summary_file}")
        return
    
    df = pd.read_csv(summary_file)
    
    if len(df) == 0:
        print("No dipole data to visualize!")
        return
    
    # Create figure with multiple subplots
    fig = plt.figure(figsize=(18, 12))
    
    # 1. Direction distribution for each component (positive vs negative)
    ax1 = fig.add_subplot(221)
    
    # Count positive and negative values for each component
    pos_neg_counts = {
        'X': [sum(df['dipole_x'] > 0), sum(df['dipole_x'] < 0)],
        'Y': [sum(df['dipole_y'] > 0), sum(df['dipole_y'] < 0)],
        'Z': [sum(df['dipole_z'] > 0), sum(df['dipole_z'] < 0)]
    }
    
    # Plot grouped bar chart
    x = np.arange(3)
    width = 0.35
    
    ax1.bar(x - width/2, [pos_neg_counts['X'][0], pos_neg_counts['Y'][0], pos_neg_counts['Z'][0]], 
           width, label='Positive', color='green')
    ax1.bar(x + width/2, [pos_neg_counts['X'][1], pos_neg_counts['Y'][1], pos_neg_counts['Z'][1]], 
           width, label='Negative', color='red')
    
    ax1.set_xlabel('Component')
    ax1.set_ylabel('Count')
    ax1.set_title('Direction Distribution (Positive vs Negative)')
    ax1.set_xticks(x)
    ax1.set_xticklabels(['X', 'Y', 'Z'])
    ax1.legend()
    ax1.grid(True, linestyle='--', alpha=0.7)
    
    # Add labels on bars
    for i, v in enumerate([pos_neg_counts['X'][0], pos_neg_counts['Y'][0], pos_neg_counts['Z'][0]]):
        ax1.text(i - width/2, v + 0.1, str(v), ha='center')
    for i, v in enumerate([pos_neg_counts['X'][1], pos_neg_counts['Y'][1], pos_neg_counts['Z'][1]]):
        ax1.text(i + width/2, v + 0.1, str(v), ha='center')
    
    # 2. Octant distribution (which 3D octant the dipoles point to)
    ax2 = fig.add_subplot(222)
    
    # Calculate octant for each dipole
    def get_octant(x, y, z):
        return (int(x > 0), int(y > 0), int(z > 0))
    
    octants = df.apply(lambda row: get_octant(row['dipole_x'], row['dipole_y'], row['dipole_z']), axis=1)
    octant_names = {
        (0, 0, 0): 'O1 (-,-,-)',
        (1, 0, 0): 'O2 (+,-,-)',
        (0, 1, 0): 'O3 (-,+,-)',
        (1, 1, 0): 'O4 (+,+,-)',
        (0, 0, 1): 'O5 (-,-,+)',
        (1, 0, 1): 'O6 (+,-,+)',
        (0, 1, 1): 'O7 (-,+,+)',
        (1, 1, 1): 'O8 (+,+,+)'
    }
    
    octant_counts = octants.value_counts().reindex(list(octant_names.keys()), fill_value=0)
    
    # Plot octant distribution
    colors = plt.cm.tab10(np.linspace(0, 1, 8))
    bars = ax2.bar(range(8), octant_counts, color=colors)
    
    ax2.set_xlabel('Octant')
    ax2.set_ylabel('Count')
    ax2.set_title('3D Octant Distribution')
    ax2.set_xticks(range(8))
    ax2.set_xticklabels([octant_names[k] for k in octant_counts.index], rotation=45, ha='right')
    ax2.grid(True, linestyle='--', alpha=0.7)
    
    # Add labels on bars
    for i, v in enumerate(octant_counts):
        ax2.text(i, v + 0.1, str(v), ha='center')
    
    # 3. Angular deviation from mean dipole
    ax3 = fig.add_subplot(223)
    
    # Calculate mean dipole vector
    mean_dipole = df[['dipole_x', 'dipole_y', 'dipole_z']].mean().values
    mean_dipole_mag = np.linalg.norm(mean_dipole)
    
    # Calculate angles between each dipole and mean dipole
    def calc_angle(row):
        v1 = np.array([row['dipole_x'], row['dipole_y'], row['dipole_z']])
        v2 = mean_dipole
        cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        # Clip to avoid numerical errors
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        return np.degrees(np.arccos(cos_angle))
    
    df['angle_from_mean'] = df.apply(calc_angle, axis=1)
    
    # Plot histogram of angles
    ax3.hist(df['angle_from_mean'], bins=10, color='skyblue', edgecolor='black')
    ax3.set_xlabel('Angle from Mean Dipole (degrees)')
    ax3.set_ylabel('Frequency')
    ax3.set_title(f'Angular Deviation from Mean Dipole Vector [{mean_dipole[0]:.2f}, {mean_dipole[1]:.2f}, {mean_dipole[2]:.2f}]')
    ax3.grid(True, linestyle='--', alpha=0.7)
    
    # 4. Direction consistency radar chart
    ax4 = fig.add_subplot(224, projection='polar')
    
    # Calculate spherical coordinates for each dipole
    def cart_to_sph(row):
        x, y, z = row['dipole_x'], row['dipole_y'], row['dipole_z']
        r = np.sqrt(x**2 + y**2 + z**2)
        theta = np.arctan2(y, x)  # azimuthal angle
        if theta < 0:
            theta += 2 * np.pi
        phi = np.arccos(z/r)  # polar angle
        return theta, phi, r
    
    # Group dipoles by azimuthal angle into bins
    n_bins = 8
    bin_width = 2 * np.pi / n_bins
    angle_bins = [[] for _ in range(n_bins)]
    
    for _, row in df.iterrows():
        theta, phi, r = cart_to_sph(row)
        bin_idx = int(theta / bin_width)
        angle_bins[bin_idx].append(r)
    
    # Calculate average magnitude in each angular bin
    bin_avgs = [np.mean(bin_vals) if bin_vals else 0 for bin_vals in angle_bins]
    bin_angles = np.linspace(0, 2*np.pi, n_bins, endpoint=False)
    
    # Create radar chart
    ax4.bar(bin_angles, bin_avgs, width=bin_width*0.9, alpha=0.7, color='purple')
    ax4.set_theta_zero_location('N')
    ax4.set_theta_direction(-1)  # clockwise
    ax4.set_title('Average Dipole Magnitude by Direction')
    
    # Add cardinal directions
    directions = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
    ax4.set_xticks(np.linspace(0, 2*np.pi, len(directions), endpoint=False))
    ax4.set_xticklabels(directions)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Directional analysis saved to {output_file}")
    plt.close(fig)

def generate_comprehensive_dipole_visualizations(input_dir="./fixed_data/clean_data", output_dir="/home/hp/nayanika/github/GPX6/figures"):
    """Generate comprehensive visualizations for dipole data"""
    
    # Check if the summary file exists
    summary_file = os.path.join(input_dir, "dipole_summary.csv")
    if not os.path.exists(summary_file):
        print(f"ERROR: Summary file not found: {summary_file}")
        return
    
    # Create visualization directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate all visualizations
    visualize_all_dipoles(
        summary_file, 
        os.path.join(output_dir, "dipole_3d_visualization.png"),
        plot_type="3d"
    )
    
    visualize_all_dipoles(
        summary_file, 
        os.path.join(output_dir, "dipole_2d_visualization.png"),
        plot_type="2d"
    )
    
    visualize_all_dipoles(
        summary_file, 
        os.path.join(output_dir, "dipole_rosette.png"),
        plot_type="rosette"
    )
    
    visualize_all_dipoles(
        summary_file, 
        os.path.join(output_dir, "dipole_all_visualizations.png"),
        plot_type="all"
    )
    
    create_heatmap_visualization(
        summary_file,
        os.path.join(output_dir, "dipole_heatmap.png")
    )
    
    create_statistical_visualization(
        summary_file,
        os.path.join(output_dir, "dipole_statistics.png")
    )
    
    create_directional_analysis(
        summary_file,
        os.path.join(output_dir, "dipole_directional_analysis.png")
    )
    
    print(f"\nAll dipole visualizations have been saved to {output_dir}")
    print("Summary of created visualizations:")
    print("1. 3D visualization of all dipole vectors")
    print("2. 2D projections in XY, XZ, and YZ planes")
    print("3. Dipole direction rosette (polar plot)")
    print("4. Combined visualization with all plot types")
    print("5. Heatmap of dipole components by replica")
    print("6. Statistical analysis of dipole data")
    print("7. Directional analysis of dipoles")

# Modified to use the new output directory
if __name__ == "__main__":
    # Set the matplotlib backend to Agg at the very start
    matplotlib.use('Agg')
    
    # Base directory for processed data
    base_dir = "./fixed_data"
    
    # Output directory specified by the user
    output_dir = "/home/hp/nayanika/github/GPX6/figures"
    
    # Generate all visualizations with the new output directory
    generate_comprehensive_dipole_visualizations(
        input_dir=os.path.join(base_dir, "clean_data"),
        output_dir=output_dir
    )