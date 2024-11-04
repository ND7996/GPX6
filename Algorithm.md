## Protocol for Running Empirical Valence Bond Simulations On Several Mutations

### 1. Create **PDB Structures** of Single and Combinations of Mutations

- **Code**: [View PDB Creation Code](prep_structures/prep_structure.ipynb)

- **Available Files**:
  - [Mouse Structures](prep_structures/mousecys)
  - [Human Structures](prep_structures/humansec)

### 2. Create **Topology Files**
Generate topology files for PDB files:

- **Topology File**: Done with qprep in Q program
  - [Cysteine](prep_structures/mousecys/prep5.inp)
  - [Selenocysteine](prep_structures/humansec/prep5.inp)
 
- **Available Files**:
  - [Mouse Topology Files](prep_structures/mousecys)
  - [Human Topology Files](prep_structures/humansec)

### 3. Create Input Files for **Minimization**

- **Code**
  - [Cysteine](prep_scripts/minimcysteine.ipynb)
  - [Selenocysteine](prep_scripts/minimselenocysteine.ipynb)
    
### 4. Create PDB Files for **Equilibration and Free Energy calculations**

  - Done with qprep in Q program
  - [Cysteine](prep_structures/mousecys/prep5.inp)
  - [Selenocysteine](prep_structures/humansec/prep5.inp)

### 5. Create Input Files for **FEP calculations**

- **Code**
  - [Cysteine](prep_scripts/fepcysteine.ipynb)
  - [Selenocysteine](prep_scripts/fepselenocysteine.ipynb)
    
### 6. Analysis

  - [Energy Files](analysis/copy_energy.sh)
  - [Trajectory Files](analysis/copy_dcd.sh)
    
- **Code**
- [Mapping Energies](analysis/mapper.sh)

  
    
  
