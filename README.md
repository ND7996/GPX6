GPX6 Mutation Analysis Workflow
This document outlines the complete workflow for analyzing GPX6 mutations, including structure preparation, relaxation, and Free Energy Perturbation (FEP) calculations.

Overview
The workflow consists of seven main steps:

Structure preparation with PyMOL
Structure solvation using Qprep
FEP file generation
Relaxation input preparation
Structure minimization
FEP calculation setup
Analysis using QTools

Detailed Steps
1. Structure Preparation
Purpose: Automate residue mutations in PDB files using PyMOL, applying Mouse-to-Human substitutions including SEC incorporation.
Location: prep_structures/
Output Directories:

Mouse mutations: prep_structures/MOUSE/
Human mutations: prep_structures/HUMAN/

2. Structure Solvation
Purpose: Prepare PDB structures for relaxation by applying solvation and boundary conditions.
Process:

Uses qprep5 to process each PDB file
Generates input files automatically based on structure names
Creates solvated structures and topology files

Outputs:

Solvated structure: ${base_name}_solvated.pdb
Topology file: ${base_name}_solvated.top

3. FEP File Generation
Purpose: Generate Free Energy Perturbation (FEP) files using makeFEP.py.
Location: input/fep/
Input Maps:

Mouse WT: fepmousecys.qmap
Human WT: fephumansec.qmap

4. Relaxation Input Preparation
Script: /home/hp/nayanika/github/GPX6/prep_scripts/minimcysteine.ipynb
Required Inputs:

genrelax.proc: Relaxation parameters (steps, temperature, bath coupling, restraints)
Topology file: ${base_name}_solvated.top
Structure file: ${base_name}_solvated.pdb
FEP file: from input/fep/

Output:

Directory: /home/hp/results/minim/
Generated script: run_qdyn_5.sh

5. Structure Minimization
Purpose: Generate minimized structures for FEP calculations.
Process:

Iterates through system directories in /home/hp/results/
Processes directories containing minim folders
Uses the last restart file (relax_012.re)
Executes relaxation using qprep5

Output: minim.pdb in respective system folders
6. FEP Calculation Setup
Script: /home/hp/nayanika/github/GPX6/prep_scripts/fepcysteine.ipynb
Parameters:

Input file: minim.pdb from each minim directory
FEP frames: 51
Starting lambda: 1.0
Output prefix: "replica"

Outputs:

Multiple replica directories
Generated script: run_qdyn_5.sh in each replica folder


