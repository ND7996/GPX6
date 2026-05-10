cimport pymol
import pandas as pd
import os
import glob
import time

# Define the mapping of mutant codes to their actual residue numbers
residue_dict = {
    "D148E": 148, "F139L": 139, "G102S": 102, "H177Q": 177, "K87T": 87, 
    "P142S": 142, "R99C": 99, "S4R": 4, "T178A": 178, "T54Q": 54, 
    "Y104F": 104, "E143S": 143, "F48Y": 48, "H144Q": 144, "I24L": 24, 
    "N107S": 107, "R181S": 181, "S47A": 47, "T52A": 52, "T60A": 60
}

# Define charge values for each residue
charge_dict = {
    "D148E": -1, "F139L": 0, "G102S": 0, "H177Q": 0, "K87T": +1, 
    "P142S": 0, "R99C": +1, "S4R": 0, "T178A": 0, "T54Q": 0, 
    "Y104F": 0, "E143S": -1, "F48Y": 0, "H144Q": 0, "I24L": 0, 
    "N107S": 0, "R181S": +1, "S47A": 0, "T52A": 0, "T60A": 0
}

# Define dg_star values
free_energy_data = {
    "F48Y": 17.16, "T52A": 17.83, "S47A": 19.07, "R99C": 17.29, 
    "H177Q": 16.09, "H144Q": 16.60, "T178A": 15.41, "A74G": 17.88, "E143S": 17.52, "F139L": 17.01, 
    "K87T": 20.21, "P142S": 18.30, "G102S": 18.72, "Y104F": 15.75, "L24I": 17.88, "T60A": 19.22, 
    "R181S": 19.38, "S4R": 15.73, "N107S": 21.25
}

# Initialize PyMOL
pymol.finish_launching(["pymol", "-cq"])

# Define the list of mutants
mutants = list(residue_dict.keys())

# Create an empty list to store data
data = []

# Function to find available PDB files
def find_pdb_file(mutant):
    base_dir = f"/home/hp/results/MOUSE/level1/{mutant}/"
    patterns = ["replica000_fep_050_0.000.pdb", "replica001_fep_050_0.000.pdb", "*.pdb"]
    
    for pattern in patterns:
        files = glob.glob(os.path.join(base_dir, pattern))
        if files:
            return files[0]
    return None

# Iterate through the mutants
for mutant in mutants:
    pdb_file = find_pdb_file(mutant)
    if not pdb_file:
        print(f"Warning: No PDB file found for {mutant}")
        continue
    
    print(f"Processing {mutant} using file: {pdb_file}")
    pymol.cmd.delete(mutant)
    pymol.cmd.load(pdb_file, mutant)
    residue_number = residue_dict[mutant]
    
    try:
        dist_value = pymol.cmd.get_distance(f"{mutant} and resi 49 and name CA", 
                                            f"{mutant} and resi {residue_number} and name CA")
        charge = charge_dict[mutant]
        dg_star = free_energy_data.get(mutant, None)
        data.append([49, residue_number, dist_value, charge, dg_star])
        print(f"Successfully calculated distance for {mutant}: {dist_value}Å")
    except pymol.CmdException as e:
        print(f"Error calculating distance for {mutant}: {str(e)}")
    finally:
        pymol.cmd.delete(mutant)
        time.sleep(0.1)

# Create a pandas dataframe
df = pd.DataFrame(data, columns=["Residue 49", "Mutant Residue", "Distance", "Charge", "dg_star"])

df.to_csv("mutant_distances_with_charges_dg_star.csv", index=False)
print("\nFinal Results:")
print(df)

# Quit PyMOL when done
pymol.cmd.quit()
