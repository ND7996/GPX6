# Empirical Valence Bond Simulation on Glutathione Peroxidase 6

## Running Q

### Parameterization

1. **Prepare the GPX6 Structure**  
   Using the GPX6 structure from MD simulation (AlphaFold for human, PDB ID: 7FC2 for mouse), remove the hydrogens using the command:

**Input File**: `/protein_stepwise/GPX6WT/mousecys/1-prep/no_hyd.pdb`

2. **Obtain Parameters for SEC GPX6**  
Acquire parameters for SEC GPX6 in its active site using Maestro (FFLD server).

---

## STEP 1: Preparation

1. **Read the Library File**  
Load the library file from Qtools: **`qoplsaa.lib`**.

2. **Read the Parameters File**  
Load the parameters file along with selenium parameters: **`qoplsaa_all.prm`**.

3. **Prepare a New GPX.lib File**  
Include library entries that are not present in `qoplsaa.lib`, such as:
- Hydrogen peroxide
- Selenol
- Selenenic Acid
- Selenolate ion
- Protonated glutamine

4. **Set the Sphere of Water Molecules**  
Define the boundary sphere around selenium and cysteine:
- Boundary sphere: `49:SG` and `49:SE` (residue number: residue type)
- Radius: `solvate 49:SG 25` and `49:SE 25`
- Grid: `HOH`

---

## Output

1. Generate the topology file.
2. Write the topology file.
3. Write the PDB file.

---

## STEP 2: Making the QMAP File

1. **Create a File for Quantum Atoms**  
Include all quantum atoms participating in the reaction (State 1 - Reactant, State 2 - Product).  
The file should contain the PDB ID and the corresponding LIB ID.

2. **Run the Script to Generate the FEP File**  
Use the topology, PDB after performing Qprep, along with the original `qoplsaa.lib`, `GPX.lib`, and `qoplsaa_all.prm`.

---

### FEP File

Once the FEP file is generated with the script, it will be used for relaxing the system. The FEP file includes:
- Changes in charges, bonds, torsions, impropers, and angles.
- Soft-core potentials.
- Harmonic and Morse potentials for the bonds that are being made and broken.

---

## STEP 3: Running Relaxation

1. **Set Parameters**  
Define steps, temperature, cut-offs, and restraints using: **`genrelax.proc`**.

2. **Run the Script**  
Use the script **`q_genrelax.py`** from `qtools/qscripts-cli`.

3. **Input Files**  
- Topology file from Qprep.
- PDB file from Qprep.
- FEP file.

---

### Generating Inputs

Run the following command to generate inputs:
```bash
q_genrelax.py genrelax.proc --top GPX6cys_mouse.top --pdb GPX6cys_mouse.pdb --fep GPX6_wtmousecys.fep --outdir minim --rest top --rs run_Q5.10_amd.sh
