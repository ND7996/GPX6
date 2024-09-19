
# Empirical Valence Bond Simulation on Glutathione Peroxidase 6 

## Running Q

## Parameterization

1. Using the gpx6 structure from md simulation (alphafold in case of human and pdb id 7fc2 in case of mouse) remove the hydrogens using the command remove hydro in pymol. input file - /protein_stepwise/GPX6WT/mousecys/1-prep/no_hyd.pdb
2. Obtaining the parameters for SEC gpx6 in its active site using maestro (ffld server)
   
## STEP 1 Preparation 

1. Reading the library file from qtools **qoplsaa.lib**
2. Reading the parameters file along with selenium parameters **qoplsaa_all.prm**
3. Preparing a new **GPX.lib** file which includes the library entries of those not present in the qoplsaa.lib
   - Hydrogen peroxide , Selenol, Selenenic Acid, Selenolate ion and protonated glutamine
4. Setting the sphere of water molecules around the selenium and cystiene 
boundary sphere 49:SG and 49:SE (residue number:residue type) and radius solvate 49:SG 25. 49:SE 25. grid HOH

## Output 

1. makes the topology file 
2. writes topology file 
3. writes pdb file 

##  STEP 2 Making the qmap file

1. Create a file with the all quantum atoms that participate in the reaction (state 1 - reactant), (state 2 - product) (The file contains the PDB ID and corresponding LIB ID for it)
2. Run the script to make the fep file, it takes the topology, pdb after performing qprep and the original qoplsaa.lib, GPX.lib, qoplsaa_all.prm 

## Fep file

Once the fep file is generated with the script, we will use this file further for relaxing the system. 

The fep file includes the change in the charges, bonds, torsions, impropers and angles , soft core potentials, Harmonic and Morse potentials for the bonds that are being made and broken

##  STEP 3 Running Relaxation

1. Setting the steps, temperature, cut-offs, restraints -**genrelax.proc**
2. Use the script **q_genrelax.py** from qtools/qscripts-cli
3. Topology file from Qprep
4. Pdb file from Qprep
5. Fep file 

## Generating inputs

q_genrelax.py genrelax.proc --top GPX6cys_mouse.top  --pdb GPX6cys_mouse.pdb --fep GPX6_wtmousecys.fep --outdir minim --rest top --rs run_Q5.10_amd.sh

Run **sbatch run_Q5.10_amd.sh** in the new folder minim , it runs Qdyn6 on every relax input to obtain the dcd, restart and log files 

## Output 

1. Log file - energy summaries and temperature is written to the standard output log file
2. DCD - trajectory file
3. Restart file

##  STEP 4 Writing a new Pdb file after minimization 

1. Reading the library file qoplsaa.lib from qtools
2. Reading GPX.lib     
3. Reading the parameter file qoplsaa_all.prm
4. Reading the topology file from Qprep
5. Reading the last restart file

## Output

The minimised pdb structure that is used to run FEP calculations

##  STEP 5 Running FEP calculation

1. Setting the steps, temperature, cut-offs, restraints -**genfeps.proc** 
2. Use the script **q_genfeps.py** from qtools/qscripts-cli
3. Topology file from Qprep
4. Pdb structure generated from last restart file after minimizing 
5. Fep file 
6. Repeats - as many replicas as you want

## Generating inputs

q_genfeps.py genfeps.proc --pdb minimised.pdb relax_008.inp relax --repeats 5 --frames 51 --fromlambda 1.0 --prefix replicamousecys --rs run_Q5.10_amd.sh

Run **sbatch run_Q5.10_amd.sh** in the replica folders , it runs Qdyn6 on every replica to obtain the dcd, restart, log and en (energy) files 

## Output of FEP

1. Log file - energy summaries and temperature is written to the standard output log file
2. DCD - trajectory file
3. Restart file
4. q_enfiles.list - all the energy files are listed of all the frames 

## STEP 6 Calibration

## Running automapper that gives an Hij and Alpha value which can be used to run mapper 
1. q_automapper.py 5.0 -10.0 10.0 -200.0 --temp 300 --dirs replica*
   
 5.0  ref_dga Reference activation free energy. (modify for you system), -10.0  ref_dg0  Reference reaction free energy. (modify for you system), 10.0  init_hij  Initial guess for Hij (offdiagonal) (modify  for you system), -200.0 init_alpha Initial guess for alpha (state 2 shift) (modify for you system)

## Running mapper to map free energies

2. q_mapper.py 79.106841 -92.883878 --bins 50 --skip 100 --min 10 --temp 300.0 --dirs replica*

79.106841 hij  Hij coupling constant (modify  for you system), -92.883878 alpha state 2 shift (alpha) (modify  for you system)

## Output 

gives the mean, dG, dG0, SD and error from all the replicas 

## Running analyse

3.  q_analysefeps.py replica* 

## Output 

gives a json file **qaf.PlotData.json** compiling different replicas 

##  STEP 7 Plotting free energies 

q_plot.py qaf.PlotData.json 


