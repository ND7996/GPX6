# Free Energy Perturbation (FEP) Workflow Guide

This comprehensive guide outlines the complete automated workflow for performing Free Energy Perturbation calculations using Q software suite, from initial structure preparation through final statistical analysis of GPX6 protein as an example

![Detailed Workflow](https://raw.githubusercontent.com/ND7996/GPX6/main/analysis/figures/detailed_workflow.drawio.png)

## Workflow Overview

This pipeline automates the calculation of binding free energy differences between protein variants through a systematic 7-step process involving structure preparation, solvation, FEP setup, molecular dynamics relaxation, and statistical analysis.

---

## STEP 1 - Prepare the Structures

### Purpose
This script automates residue mutations in a given PDB file using PyMOL, applying predefined Mouse-to-Human substitutions (including SEC incorporation) and saving each mutated structure separately.

### Process Details
- **Tool**: PyMOL automation script
- **Function**: Systematic mutation of specific residues
- **Mutations Applied**:
  - Mouse-to-Human amino acid substitutions
  - Human-to-Mouse amino acid substitutions

- **Output**: Individual PDB files for each mutated variant

### Key Features
- Automated residue substitution based on predefined mapping
- Preservation of protein backbone structure
- Quality control for mutation compatibility
- Standardized naming convention for output files

---

## STEP 2 - Running Qprep

### Purpose
Runs qprep5 and automates the preparation of PDB solvated structures for relaxation by applying solvation and boundary conditions.

### Process Details
1. **File Processing**:
   - Extracts the base name (without the .pdb extension) from each PDB file generated in the previous step
   - Uses this base_name to create a corresponding .inp file for qprep5

2. **Qprep5 Execution**:
   - Runs qprep5 with the generated .inp file
   - Applies solvation parameters and boundary conditions
   - Configures system for molecular dynamics simulation

3. **Output Generation**:
   - Saves the final solvated structures as `${base_name}_solvated.pdb`
   - Generates topology files as `${base_name}_solvated.top`
   - All files stored in the same directory for organization

### Technical Specifications
- **Solvation**: Explicit water model implementation
- **Boundary Conditions**: Periodic boundary conditions applied
- **Force Field**: Q-compatible parameter set

---

## STEP 3 - Making FEP File

### Purpose
This step executes the makeFEP.py script within the directory to generate the required Free Energy Perturbation file.

### Process Details
- **Script**: makeFEP.py execution
- **Function**: The script utilizes the qmap file and the solvated PDB structure to generate the required FEP file
- **Integration**: Seamless connection between solvated structures and FEP calculations

### Input Requirements
| Input File | Purpose | Application |
|------------|---------|-------------|
| `fepmousecys.qmap` | Mouse WT mapping | For Mouse wild-type calculations |
| `fephumansec.qmap` | Human WT mapping | For Human wild-type calculations |

### Output
- FEP file containing lambda-dependent parameters
- Mapping between initial and final states
- Perturbation pathway definition

---

## STEP 4 - Make Inputs for Relaxation

### Purpose
Runs q_genrelax.py to generate comprehensive input files for molecular dynamics relaxation.

### Required Inputs
- **`--genrelax.proc`**: Input parameters for relaxation including:
  - Number of steps
  - Temperature settings
  - Bath coupling parameters
  - Restraint definitions

- **`--top`**: Specifies the topology file `${base_name}_solvated.top`

- **`--pdb`**: Provides the structure file `${base_name}_solvated.pdb`

- **`--fep`**: Supplies the Free Energy Perturbation (FEP) file from Step 3

- **`--outdir minim`**: Defines the output directory for minimized structures

- **`--rs run_qdyn_5.sh`**: Generates the qdyn execution script inside the minim folder

### Process Flow
1. Parameter validation and input file verification
2. Generation of relaxation input files (`relax.inp*`)
3. Creation of execution scripts for automated processing
4. Directory structure setup for organized output management

---

## STEP 5 - Make Minimized PDB for FEP Calculations

### Purpose
The script iterates over system directories to generate minimized PDB structures suitable for FEP calculations.

### Process Details
1. **Directory Scanning**:
   - Iterates over System Directories
   - Checks if each subdirectory contains a minim folder
   - Searches for topology (.top) files in corresponding minim folders

2. **Relaxation Processing**:
   - Identifies the last restart file `relax_012.re`
   - Executes relaxation using the command: `rx $relax_file`

3. **Structure Generation**:
   - Runs qprep5 specifying the .prm, .lib and relax file
   - Writes the output minimized structure as `minim.pdb`
   - Saves in the respective system folder

### Quality Control
- Verification of convergence criteria
- Energy minimization validation
- Structural integrity checks

---

## STEP 6 - Make Inputs for FEP

### Purpose
Runs q_genfeps.py to generate comprehensive FEP calculation input files with multiple replica support.

### Required Inputs
- **`--genfeps.proc`**: Input parameters for equilibration and FEP including:
  - Number of steps
  - Temperature control
  - Bath coupling parameters
  - Restraint specifications

- **`--pdb`**: Uses the `minim.pdb` input file from each minim/ directory

- **`--repeats`**: Generates specified number of independent replicas for statistical robustness

- **`--frames 51`**: Configures 51 frames FEP calculation per replica

- **`--fromlambda 1.0`**: Sets the starting lambda value for FEP calculations

- **`--prefix replica`**: Saves results in replica folders inside system directories

- **`--rs run_qdyn_5.sh`**: Generates qdyn execution scripts inside each replica folder

## STEP 7 - Running Qtools for Analysis

### Purpose
Comprehensive analysis of FEP calculations using Q software analysis tools.

### Analysis Pipeline

#### 7.1 Data Mapping
- **Tool**: q_mapper.py
- **Scope**: Processes all replica directories (rep*)
- **Function**: Maps FEP data for analysis compatibility

#### 7.2 FEP Analysis
- **Tool**: q_analysefeps.py
- **Scope**: Analyzes all replica directories (rep*)
- **Output**: 
  - Primary results saved in `test.out`
  - Structured data in JSON format

#### 7.3 Statistical Extraction
- **Process**: Extract Key Statistics
- **Method**: 
  - Finds the line containing "Statistics" in `test.out`
  - Extracts the next 5 lines of statistical data
  - Saves processed statistics in `stats_output.txt`

#### 7.4 Report Generation
- **LaTeX Table**: Professional formatting for publication
- **CSV File**: Data export for further processing
- **Key Metrics**:
  - Mean ΔG* value ± standard error
  - Mean ΔG₀ value ± standard error

