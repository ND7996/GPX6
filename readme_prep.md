# Protein Preparation and Parameter Generation Workflow

## 1. Open Maestro
Launch Maestro and follow the steps below to prepare your protein structure.

## 2. Protein Preparation Workflow

Navigate to **Protein Preparation** → **Preparation Workflow** and perform the following:

- **Preprocess**  
  - Cap termini  
  - Fill missing side chains  

- **Check Structures**  
  - Add hydrogens if needed:  
    `Edit` → `Add Hydrogens`

- **3D Builder**  
  - Add missing atoms (C, O), charges, etc.

- **Optimize**  
  - Run **PROPKA**  
  - Apply flips (e.g., HIS, ASN, GLN) where necessary

> 💾 **Export the structure** as:  
> `maestroGPX6_wt.mae`

---

## 3. Generate Force Field Parameters

Use `ffld_server` to convert the `.mae` or `.pdb` file into a parameter log.

### For the Maestro structure:
```bash
/opt/schrodinger2022-1/utilities/ffld_server -imae maestroGPX6_wt.mae -version 14 -print_parameters -out_file GPX_PARAM.log
