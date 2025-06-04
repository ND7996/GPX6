# Protein Preparation and Parameter Generation Workflow

This guide outlines the complete workflow for preparing protein structures and generating force field parameters using Maestro and Schrödinger tools.

## Prerequisites

- Schrödinger Suite (version 2022-1 or compatible)
- Access to Maestro GUI
- Command line access to Schrödinger utilities
- Python environment with q_ffld2q.py script

## Workflow Overview

### 1. Protein Preparation in Maestro

#### Step 1.1: Open Maestro
- Launch Maestro application
- Import your protein structure file

#### Step 1.2: Protein Preparation Wizard
Navigate to: **Protein Preparation → Preparation Workflow**

#### Step 1.3: Preprocessing
- **Action**: Select "Preprocess"
- **Tasks**:
  - Cap termini (if required)
  - Fill missing side chains
  - Remove unwanted molecules/waters (optional)

#### Step 1.4: Structure Validation
- **Action**: Check structures for inconsistencies
- **Manual adjustments**: If needed, use **EDIT → ADD HYDROGENS** to manually add missing hydrogens

#### Step 1.5: 3D Structure Building
- **Tool**: Use **3D Builder**
- **Purpose**: Adds missing atoms (carbon, oxygen, etc.) and assigns proper charges
- **Note**: This step ensures chemical completeness of the structure

#### Step 1.6: Structure Optimization
- **Action**: Run "Optimize" step
- **Functions**:
  - Checks protonation states using PROPKA
  - Adds hydrogen bond flips where required
  - Optimizes side chain conformations

#### Step 1.7: Export Structure
- **Format**: Export as `.mae` file
- **Filename**: `maestroGPX6_wt.mae` (or your preferred naming convention)
- **Location**: File → Export Structures

### 2. Parameter Generation Using ffld_server

#### Step 2.1: Generate Protein Parameters
```bash
/opt/schrodinger2022-1/utilities/ffld_server \
    -imae maestroGPX6_wt.mae \
    -version 14 \
    -print_parameters \
    -out_file GPX_PARAM.log
```

**Parameters explained**:
- `-imae`: Input Maestro file
- `-version 14`: Force field version (OPLS3e)
- `-print_parameters`: Output force field parameters
- `-out_file`: Specify output log file

#### Step 2.2: Generate Small Molecule Parameters (if applicable)
```bash
/opt/schrodinger2022-1/utilities/ffld_server \
    -ipdb h2o2.pdb \
    -print_parameters \
    -version 14 \
    > h2o2_param.ffld11
```

**Note**: Replace `h2o2.pdb` with your small molecule PDB file

### 3. Converting to Q-Chem Format

#### Step 3.1: Use q_ffld2q.py Script
```bash
q_ffld2q.py h2o2_param.ffld11 h2o2.pdb
```

**Required inputs**:
- `h2o2_param.ffld11`: ffld_server output file
- `h2o2.pdb`: Original PDB structure file used to create the ffld_output

**Generated outputs**:
- `.prm` file: Parameter file
- `.lib` file: Library file
- `.prm.chk` file: Checkpoint file


