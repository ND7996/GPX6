===== PyMOL Dipole Visualization =====

These files allow you to visualize the reactant and product structures
along with their dipole vectors in PyMOL.

=== How to Use ===

Method 1 - Using the batch script:
1. Run the 'start_pymol.sh' script:
   $ ./start_pymol.sh

Method 2 - Manual loading:
1. Start PyMOL
2. Run the master script:
   PyMOL> run /path/to/master_dipole.pml

3. Load a specific replica:
   PyMOL> load_replica(15)  # Replace 15 with your desired replica number

=== Keyboard Shortcuts ===
Once a replica is loaded:
F1: View all atoms
F2: View only centers of mass
F3: View only dipole vector

From the master script:
F4-F12: Load replicas 0-8

=== Additional Notes ===
- Blue spheres represent reactant atoms
- Red spheres represent product atoms
- The larger blue sphere is the reactant center of mass
- The larger red sphere is the product center of mass
- The purple line represents the dipole vector
