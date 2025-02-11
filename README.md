### STEP 1 - Prepare the Structures
This script automates residue mutations for in a given PDB file using PyMOL, applying predefined Mouse-to-Human substitutions (including SEC incorporation) and saving each mutated structure separately. 

### STEP 2 - Running Qprep
Runs qprep5 and automates the preparation of PDB solvated structures for relaxation by applying solvation and boundary conditions 
Extracts the base name (without the .pdb extension) from each PDB file generated in the previous step.
Uses this base_name to create a .inp file for qprep5.
Runs qprep5 with the generated .inp file.
Saves the final solvated structures as ${base_name}_solvated.pdb and the topology as ${base_name}_solvated.top in the same directory

### STEP 3 - Making FEP file
This step executes the makeFEP.py script within the directory:
[Go to FEP Directory](../../input/fep)
The script utilizes the qmap file and the solvated PDB structure to generate the required FEP file
Inputs:
fepmousecys.qmap → for Mouse WT
fephumansec.qmap → for Human WT

### STEP 4 - Make inputs for relaxation
 runs q_genrelax.py, which requires the following inputs:
--genrelax.proc: input for relax (steps,temp,bath coupling,restrains)  
--top: Specifies the topology file ${base_name}_solvated.top, located in either the 
--pdb: Provides the structure file ${base_name}_solvated.pdb, sourced from the 
--fep: Supplies the Free Energy Perturbation (FEP) file from the 
--outdir minim: Defines the output directory for minimized structures, stored in 
--rs run_qdyn_5.sh: Generates the qdyn script inside the minim folder to execute relaxation on each relax.inp* file.

### STEP 5 - Make minimised pdb for fep calculations
The script in folder iterates over System Directories in 
Checks if each subdirectory contains a minim folder
Searches for a topology (.top) file in corresponding minim folder
Takes the last restart file relax_012.re and executes relaxation (rx $relax_file).
Runs qprep5 specifying the .prm, .lib and relax file to write the output minimized structure called minim.pdb in the respective folder

### STEP 6 - Make inputs for FEP 
runs q_genfeps.py, which requires the following inputs:
--genfeps.proc: input for equil and fep (steps,temp,bath coupling,restrains) 
--pdb: Uses the minim.pdb input file from each minim/.
--repeats : Generates number of independent replicas
--frames 51: 51 frames FEP calculation per replica
--fromlambda 1.0: Sets the starting lambda value for FEP.
--prefix replica: Saves results in the replica folder inside the system directory
--rs run_qdyn_5.sh: Generates the qdyn script inside the each replica folders to execute relaxation on each equil.inp* and fep.inp* file.

### STEP 7 - Running qtools for analysis
Runs q_mapper.py over rep* (all replica directories).
Runs q_analysefeps.py over all replica directories (rep*).
Saves the output in test.out and generates a .json file 
Extract Key Statistics
Finds the line containing "Statistics" in test.out, extracts the next 5 lines and saves them in stats_output.txt
Generates a Latex table and csv file
mean dG* value ± standard error 
mean dG0 value ± standard error

