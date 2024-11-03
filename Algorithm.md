## Protocol for Running Empirical Valence Bond Simulations On Several Mutations

### 1. Create **PDB Structures** of Single and Combinations of Mutations

- **Code**: [View PDB Creation Code](prep_structures/prep_structure.ipynb)

- **Available Files**:
  - [Mouse Structures](prep_structures/mousecys)
  - [Human Structures](prep_structures/humansec)

### 2. Create **Topology Files**
Generate topology files for all PDB files individually:

- **Topology File**: [Done with qprep in Q program]
- (https://github.com/ND7996/GPX6/blob/main/prep_structures/mousecys/prep5.inp)

- **Available Files**:
  - [Mouse Topology Files](prep_structures/mousecys)
  - [Human Topology Files](prep_structures/humansec)

### 3. Create Files for **Minimization**
Further processing files for minimization will be handled in this step.

- **Available Files**:
  - [Cysteine](prep_scripts/minimcysteine.ipynb)
  - [Selenocysteine](prep_scripts/minimselenocysteine.ipynb)
  
