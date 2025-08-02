![Detailed Workflow](https://raw.githubusercontent.com/ND7996/GPX6/main/analysis/figures/detailed_workflow.drawio.png)

# Free Energy Perturbation (FEP) Workflow Guide

This guide presents a complete, automated workflow for conducting **Free Energy Perturbation (FEP)** calculations using the **Q software suite**. It covers all steps from initial structure preparation through final statistical analysis of free energy differences between protein variants.

---

## 🧭 Workflow Overview

This pipeline automates the calculation of binding free energy differences through a structured 7-step process involving:

1. Structure Preparation
2. System Solvation
3. FEP Setup
4. Molecular Dynamics (MD) Relaxation
5. Structure Minimization
6. FEP Input Generation
7. Statistical Analysis

---

## 🧱 STEP 1 - Structure Preparation

### 🧪 Purpose
Automates point mutations in a PDB structure using PyMOL.

### ⚙️ Details
- **Tool**: PyMOL automation script
- **Function**: Introduces predefined mutations (e.g., mouse-to-human substitutions or SEC incorporation)
- **Output**: Individual PDB files for each variant

### 🔑 Features
- Automated mutation mapping
- Preserves backbone integrity
- Ensures mutation compatibility
- Consistent output naming

---

## 💧 STEP 2 - Solvation with Qprep

### 🧪 Purpose
Prepares solvated systems and topologies for MD relaxation.

### ⚙️ Process
1. Extracts base name from each PDB file
2. Generates `.inp` file and runs `qprep5`
3. Applies solvation and boundary conditions

### 🧾 Output
- `${base_name}_solvated.pdb`
- `${base_name}_solvated.top`

### ⚙️ Specifications
- **Solvation**: Explicit water model
- **Boundary Conditions**: Periodic
- **Force Field**: Q-compatible

---

## 📁 STEP 3 - FEP File Generation

### 🧪 Purpose
Uses `makeFEP.py` to create lambda-dependent FEP input files.

### 📥 Inputs
| File | Role |
|------|------|
| `fepmousecys.qmap` | Mapping for mouse WT |
| `fephumansec.qmap` | Mapping for human WT |

### 📤 Output
- FEP files with initial/final state mappings and perturbation parameters

---

## 🔧 STEP 4 - Relaxation Input Setup

### 🧪 Purpose
Generates input files for energy minimization and equilibration.

### 📥 Required Inputs
- `--genrelax.proc`: Relaxation parameters
- `--top`: Topology file
- `--pdb`: Solvated structure
- `--fep`: FEP file
- `--outdir minim`: Output folder
- `--rs run_qdyn_5.sh`: Qdyn execution script

### 🔁 Steps
- Validates inputs
- Generates `.inp` files
- Creates run scripts
- Organizes output structure

---

## 🧬 STEP 5 - Minimized PDB for FEP

### 🧪 Purpose
Performs energy minimization and prepares PDBs for FEP simulations.

### ⚙️ Process
1. Scans system directories
2. Locates `relax_012.re` restart files
3. Runs relaxation and prepares minimized structures

### 📤 Output
- `minim.pdb` saved in each system folder

### ✅ Quality Checks
- Convergence validation
- Energy profile checks
- Structural integrity

---

## 🔂 STEP 6 - FEP Input Generation

### 🧪 Purpose
Generates complete FEP simulation input sets with replica support.

### 📥 Required Inputs
- `--genfeps.proc`: Parameters for FEP and equilibration
- `--pdb`: `minim.pdb`
- `--repeats`: Number of replicas
- `--frames`: e.g., `51`
- `--fromlambda`: e.g., `1.0`
- `--prefix`: e.g., `replica`
- `--rs`: Execution script

### 📁 Output
FEP-ready folders per replica inside each system directory.

---

## 📊 STEP 7 - Analysis with Qtools

### 🧪 Purpose
Performs statistical evaluation of FEP results.

### 🧬 Pipeline

#### 7.1 Data Mapping
- **Tool**: `q_mapper.py`
- **Function**: Aggregates lambda data across replicas

#### 7.2 FEP Analysis
- **Tool**: `q_analysefeps.py`
- **Output**: `test.out` + JSON structured data

#### 7.3 Statistical Summary
- Extracts block-averaged statistics
- Stores results in `stats_output.txt`

#### 7.4 Report Generation
- Exports:
  - CSV tables for data science workflows
  - LaTeX tables for publications

### 📈 Key Metrics
- Mean ΔG* ± SE
- Mean ΔG₀ ± SE

---

## 📚 Acknowledgments

This pipeline was built to support automated, reproducible FEP calculations in protein mutagenesis projects using the Q software suite.

