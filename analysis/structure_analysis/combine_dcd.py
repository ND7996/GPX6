import os
import MDAnalysis as mda
from MDAnalysis.coordinates.DCD import DCDWriter

def combine_dcd(directory, selection="all"):
    """
    Combine DCD files in the specified directory based on selection criteria.
    
    Parameters:
        directory (str): Path to the folder containing DCD files.
        selection (str): "fep", "equil", or "all" to filter DCD files.
    """
    dcd_files = sorted([f for f in os.listdir(directory) if f.endswith(".dcd")])
    
    if selection == "fep":
        dcd_files = [f for f in dcd_files if f.startswith("fep")]
    elif selection == "equil":
        dcd_files = [f for f in dcd_files if f.startswith("equil")]
    
    if not dcd_files:
        print("No matching DCD files found.")
        return
    
    valid_dcd_files = []
    for dcd in dcd_files:
        dcd_path = os.path.join(directory, dcd)
        try:
            u = mda.Universe(dcd_path)
            if len(u.trajectory) > 0:
                valid_dcd_files.append(dcd_path)
        except Exception as e:
            print(f"Skipping {dcd_path}: {e}")
    
    if not valid_dcd_files:
        print("No valid DCD files found.")
        return
    
    first_dcd = valid_dcd_files[0]
    output_dcd = os.path.join(directory, "combined.dcd")
    u = mda.Universe(first_dcd)
    
    with DCDWriter(output_dcd, n_atoms=u.atoms.n_atoms) as writer:
        for dcd_path in valid_dcd_files:
            print(f"Processing: {dcd_path}")
            u = mda.Universe(first_dcd, dcd_path)
            for ts in u.trajectory:
                writer.write(u.atoms)
    
    print(f"Combined DCD saved as {output_dcd}")

# Example usage
directory = "/home/hp/results/MOUSE_OLD/MOUSE/level0/R181S/traj00"
combine_dcd(directory, selection="fep")  # Change selection to "fep" or "equil" as needed
