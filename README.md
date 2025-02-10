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
--outdir minim: Defines the output directory for minimized structures, stored in [Results] /home/hp/results.
--rs run_qdyn_5.sh: Generates the qdyn script inside the minim folder to execute relaxation on each relax.inp* file.


### STEP 4 - Make minimised pdb for fep calculations



### STEP 5 - Make inputs for FEP 



### STEP 6 - Running qtools for analysis







STEP 1 - Prepare the Structures:
Runs prep_structure_mouse.py to generate mutated PDB files with predefined Mouse-to-Human substitutions.
These PDB files are saved in the MOUSE or HUMAN mutation directory.

STEP 2 - Running Qprep:
Moves to the same directory where the mutated PDBs were generated (target_dir).
Runs prep.sh to generate topology (.top) and modified PDB files.
Uses the base name of each PDB file to create qprep5 input files (.inp).
Runs qprep5 to apply solvation and boundary conditions, preparing the structures for relaxation.

Extra Step - Modifying Atom Names:
Runs atomnames.sh to ensure atom names are correctly formatted.


