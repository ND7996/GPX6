
import MDAnalysis as mda
from MDAnalysis.coordinates.DCD import DCDWriter

# List all the .dcd files in the order they should be combined
dcd_files = [
    "equil_000_1.000.dcd", "equil_001_1.000.dcd", "equil_002_1.000.dcd", "equil_003_1.000.dcd", "equil_004_1.000.dcd",
    "fep_000_1.000.dcd", "fep_001_0.980.dcd", "fep_002_0.960.dcd", "fep_003_0.940.dcd", "fep_004_0.920.dcd",
    "fep_005_0.900.dcd", "fep_006_0.880.dcd", "fep_007_0.860.dcd", "fep_008_0.840.dcd", "fep_009_0.820.dcd",
    "fep_010_0.800.dcd", "fep_011_0.780.dcd", "fep_012_0.760.dcd", "fep_013_0.740.dcd", "fep_014_0.720.dcd",
    "fep_015_0.700.dcd", "fep_016_0.680.dcd", "fep_017_0.660.dcd", "fep_018_0.640.dcd", "fep_019_0.620.dcd",
    "fep_020_0.600.dcd", "fep_021_0.580.dcd", "fep_022_0.560.dcd", "fep_023_0.540.dcd", "fep_024_0.520.dcd",
    "fep_025_0.500.dcd", "fep_026_0.480.dcd", "fep_027_0.460.dcd", "fep_028_0.440.dcd", "fep_029_0.420.dcd",
    "fep_030_0.400.dcd", "fep_031_0.380.dcd", "fep_032_0.360.dcd", "fep_033_0.340.dcd", "fep_034_0.320.dcd",
    "fep_035_0.300.dcd", "fep_036_0.280.dcd", "fep_037_0.260.dcd", "fep_038_0.240.dcd", "fep_039_0.220.dcd",
    "fep_040_0.200.dcd", "fep_041_0.180.dcd", "fep_042_0.160.dcd", "fep_043_0.140.dcd", "fep_044_0.120.dcd",
    "fep_045_0.100.dcd", "fep_046_0.080.dcd", "fep_047_0.060.dcd", "fep_048_0.040.dcd", "fep_049_0.020.dcd",
    "fep_050_0.000.dcd"
]

# Reference PDB file
pdb_file = "/home/hp/results/mousecys/K3N/minim/minim.pdb"  # Replace with the path to your .pdb file

# Load the reference structure using the PDB file
u = mda.Universe(pdb_file, dcd_files[0])  # Initialize with the first .dcd file

# Output file name for the combined trajectory
output_file = "combined_trajectory.dcd"

# Write the combined trajectory
with DCDWriter(output_file, u.atoms.n_atoms) as w:
    for dcd in dcd_files:
        print(f"Processing {dcd}...")
        u.load_new(dcd)  # Load the next .dcd file
        for ts in u.trajectory:  # Iterate over all frames in the trajectory
            w.write(u.atoms)

print(f"Combined trajectory saved as {output_file}.")

