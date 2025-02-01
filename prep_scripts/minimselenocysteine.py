import os
import subprocess

def main():
    # Directories setup
    pdb_dir = "/home/hp/nayanika/github/GPX6/prep_structures/MOUSE/48_47_52_99_54_144_177_74_178_143_87_142_104_102_139_24_181_4_60_107"
    base_scr_dir = "/home/hp/results/MOUSE/48_47_52_99_54_144_177_74_178_143_87_142_104_102_139_24_181_4_60_107"
    run_dir = "/home/hp/nayanika/github/GPX6/cluster_scripts"
    genrelax_dir = "/home/hp/nayanika/github/GPX6/input"
    topology_dir = "/home/hp/nayanika/github/GPX6/prep_structures/MOUSE/48_47_52_99_54_144_177_74_178_143_87_142_104_102_139_24_181_4_60_107"
    fep_dir = "/home/hp/nayanika/github/GPX6/input/fep"
    
    # Create the base directory if it does not exist
    os.makedirs(base_scr_dir, exist_ok=True)
    
    # Iterate through all files in the pdb_dir
    for pdb_filename in os.listdir(pdb_dir):
        if pdb_filename.endswith(".pdb"):
            pdb_file_path = os.path.join(pdb_dir, pdb_filename)
            system_name = os.path.splitext(pdb_filename)[0]  # Extract system name

            # Define paths for the topology file and output directories
            topology_file = os.path.join(topology_dir, f"{system_name}.top")
            system_dir = os.path.join(base_scr_dir, system_name)

            # Check if both the PDB and topology files exist
            if not os.path.isfile(topology_file):
                print(f"Warning: Topology file not found for {system_name}: {topology_file}")
                continue  # Skip this file if topology is missing

            # Create the system directory
            os.makedirs(system_dir, exist_ok=True)

            # Construct the q_genrelax.py command
            command = [
                "q_genrelax.py",
                os.path.join(genrelax_dir, "selgenrelax.proc"),
                "--top", topology_file,
                "--pdb", pdb_file_path,
                "--fep", os.path.join(fep_dir, "GPX6_wtmousesec.fep"),
                "--outdir", os.path.join(system_dir, "minim"),
                "--rs", os.path.join(run_dir, "run_qdyn_5.sh")
            ]

            # Run the command
            try:
                subprocess.run(command, check=True)
                print(f"Successfully processed: {system_name}")
            except subprocess.CalledProcessError as e:
                print(f"Error processing {system_name}: {e}")

if __name__ == "__main__":
    main()
