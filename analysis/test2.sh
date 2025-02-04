# PREP_STRUCTURES (Generate PDB files)
echo "Running prep_structure_mouse.py..."
python3 /home/hp/nayanika/github/GPX6/prep_structures/prep_structure_mouse.py
echo "PDB files generated!"

# PREP_TOPOLOGY (Generate topology files inside the respective folders)
# Specify the target directory for saving PDB and Topology files
target_dir="/home/hp/nayanika/github/GPX6/prep_structures/MOUSE/48_47_52_99_54_144_177_74_178_143_87_142_104_102_139_24_181_4_60_107"

# Ensure the directory exists
mkdir -p "$target_dir"

# Running prep.sh inside the folder where it's located
echo "Running prep.sh in $target_dir..."
cd "$target_dir"  # Change to the specific folder

# Run the prep.sh script to generate the topology and PDB files
sh "$target_dir/prep.sh"
echo "Topology and PDB files generated in $target_dir!"

# (CHANGE SOME ATOMS)(EXTRA STEP) 
# Run the atomnames script to generate the atom names in PDB files
sh "$target_dir/atomnames.sh"
echo "Atom names changed in $target_dir!" 

cp -r /home/hp/nayanika/github/GPX6/prep_structures/MOUSE/48_47_52_99_54_144_177_74_178_143_87_142_104_102_139_24_181_4_60_107 /home/hp/results/MOUSE
