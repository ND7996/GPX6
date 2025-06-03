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

## File Organization

```
project_directory/
├── input/
│   ├── raw_protein.pdb
│   └── small_molecule.pdb
├── maestro_output/
│   └── maestroGPX6_wt.mae
├── parameters/
│   ├── GPX_PARAM.log
│   ├── h2o2_param.ffld11
│   ├── output.prm
│   ├── output.lib
│   └── output.prm.chk
└── README.md
```

## Important Notes

1. **Force Field Version**: This workflow uses OPLS3e (version 14). Adjust the version parameter if using different force fields.

2. **File Naming**: Maintain consistent naming conventions throughout the workflow to avoid confusion.

3. **Quality Control**: Always visually inspect the prepared structure in Maestro before proceeding to parameter generation.

4. **Path Dependencies**: Ensure Schrödinger installation path is correct (`/opt/schrodinger2022-1/` in this example).

## Troubleshooting

### Common Issues

- **Missing hydrogens**: Use the manual hydrogen addition in Maestro's EDIT menu
- **Protonation states**: Check PROPKA results and manually adjust if necessary
- **File path errors**: Verify Schrödinger installation directory
- **Parameter generation fails**: Ensure input PDB and ffld_output correspond to the same structure

### Validation Steps

1. Check that all atoms have proper connectivity
2. Verify reasonable bond lengths and angles
3. Ensure proper protonation states for biological pH
4. Validate that generated parameters are complete

## Additional Resources

- Schrödinger Documentation: Protein Preparation Wizard
- OPLS Force Field Documentation
- Q-Chem Parameter Generation Guidelines

## Version Information

- Schrödinger Suite: 2022-1
- Force Field: OPLS3e (version 14)
- Last Updated: [Current Date]
