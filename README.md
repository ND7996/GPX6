TITLE

### STEP 1 - Prepare the Structures

This script automates residue mutations for in a given PDB file using PyMOL, applying predefined Mouse-to-Human substitutions (including SEC incorporation) and saving each mutated structure separately. [Mutation Script Directory](../../prep_structures)

Saves the mutated structures in the appropriate directory:  
[Mouse Mutations](../../prep_structures/MOUSE/) or [Human Mutations](../../prep_structures/HUMAN/)

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

The script runs q_genrelax.py, which requires the following inputs:

--top: Specifies the topology file ${base_name}_solvated.top, located in either the [Mouse Mutations](../../prep_structures/MOUSE/) or [Human Mutations](../../prep_structures/HUMAN/)
--pdb: Provides the structure file ${base_name}_solvated.pdb, sourced from the [Mouse Mutations](../../prep_structures/MOUSE/) or [Human Mutations](../../prep_structures/HUMAN/)
--fep: Supplies the Free Energy Perturbation (FEP) file from the [FEP Directory](../../input/fep)
--outdir minim: Defines the output directory for minimized structures, stored in [Results] 
--rs run_qdyn_5.sh: Generates the qdyn script inside the minim folder to execute relaxation on each relax.inp* file.


### STEP 4 - Make minimised pdb for fep calculations

The script in folder (../../prep_structures) iterates over System Directories in [results folder] /home/hp/results
/home/hp/results.

Checks if each subdirectory contains a minim folder
Searches for a topology (.top) file in corresponding minim folder
Takes the last restart file relax_012.re and executes relaxation (rx $relax_file).
Runs qprep5 specifying the .prm, .lib and relax file to write the output minimized structure called minim.pdb in the respective folder

### STEP 5 - Make inputs for FEP 




### STEP 6 - Running qtools for analysis

Runs q_mapper.py over rep* (all replica directories).

q_mapper.py 13.52 1.30 --bins 50 --skip 100 --min 10 --temp 300.0 --dirs rep* --qfep_exec qfep5

Runs q_analysefeps.py over all replica directories (rep*).

q_analysefeps.py rep* > test.out

Saves the output in test.out and generates a .json file 

Extract Key Statistics

grep Statistics test.out -A 5 > stats_output.txt

Finds the line containing "Statistics" in test.out, extracts the next 5 lines and saves them in stats_output.txt

Generates a Latex table of
mean dG* value ± standard error 
mean dG0 value ± standard error


