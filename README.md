# Free Energy Perturbation (FEP) Workflow Guide

This guide presents a complete, automated workflow for conducting **Free Energy Perturbation (FEP)** calculations using the **Q6 software suite**. It covers all steps from initial structure preparation through final statistical analysis of free energy differences between protein variants.

![Detailed Workflow](detailed_workflow.drawio.png)

# This workflow is generated using GPX6 protein as a model

---

## Environment & Dependencies


| Tool | Purpose |
|------|--------|
| [Q6](https://github.com/qusers/Q6) | FEP/MD engine (`qprep6`, `qdyn6`, `qfep6`) |
| [PyMOL](https://pymol.org/) | Structure preparation and mutation |
| [Qtools](https://github.com/mpurg/qtools) | FEP analysis (`q_mapper.py`, `q_analysefeps.py`) |
| Python 3.8 | Programming environment |
| OPLS-AA | Force field |

---

## Installation

Clone the repository and install the required Python packages:

```bash
git clone https://github.com/ND7996/GPX6.git
cd GPX6
pip install -r requirements.txt
```

> **Note:** For reviewer/reproducible use, the Docker workflow below is recommended. It builds Q6 and Qtools inside the container and uses repo-root relative paths throughout the scripts.

---

## Docker Reproducible Environment

This repository includes Docker files for a reviewer-friendly Linux environment:

| File | Purpose |
|------|---------|
| `Dockerfile` | Builds the GPX6 runtime image |
| `docker-compose.yml` | Mounts this repository at `/workspace/GPX6` |
| `.dockerignore` | Keeps large/generated files out of the image build context |
| `requirements-docker.txt` | Python packages for scripts, notebooks, and figures |
| `docker/entrypoint.sh` | Provides the `check` command for tool verification |

Build the image from the repository root:

```bash
docker compose build
```

The build installs Python packages, PyMOL, LaTeX, MPI/SLURM client tools, Qtools, and Q6. Qtools and Q6 are not copied into the repository folder; they are installed inside the image at:

```bash
/opt/qtools
/opt/Q6
```

The following commands are available on `PATH` inside the container:

```bash
q_genrelax.py
q_genfeps.py
q_mapper.py
q_analysefeps.py
qprep5
qdyn5_r8
qfep5
qprep6
qdyn6
qfep6
```

Check the container:

```bash
docker compose run --rm gpx6 check
```

Open a shell inside the container:

```bash
docker compose run --rm gpx6
```

Inside the shell, always run from the repository root:

```bash
cd /workspace/GPX6
```

All `.sh`, `.py`, and `.ipynb` workflow paths have been converted away from machine-specific paths such as `/home/...`, `D:\...`, and `C:\...`. Scripts now use repo-root relative paths such as:

```bash
./parameters
./prep_structures/HUMAN/level0
./prep_structures/MOUSE/level0
./input
./input/fep
./analysis_scripts
./analysis_scripts/Scripts_to_generate_figures/Figures
```

You can change the target species/level without editing files by passing environment variables, for example:

```bash
docker compose run --rm -e GPX6_SPECIES=HUMAN -e GPX6_LEVEL=0 gpx6 python prep_structures/prep_structure.py
docker compose run --rm -e PDB_DIR=./prep_structures/MOUSE/level0 gpx6 bash prep_structures/prep.sh
docker compose run --rm -e BASE_DIR=./prep_structures/HUMAN/level0 gpx6 bash cluster_scripts/run_minim.sh
```

Example commands:

```bash
docker compose run --rm gpx6 python prep_structures/prep_structure.py
docker compose run --rm gpx6 bash prep_structures/prep.sh
docker compose run --rm gpx6 bash prep_structures/preprestart.sh
docker compose run --rm gpx6 python analysis_scripts/Scripts_to_generate_figures/heatmap.py
docker compose run --rm gpx6 bash analysis_scripts/humanWT/statistics_table.sh
docker compose run --rm gpx6 make -C manuscript_files
```

If using a cluster scheduler, the SLURM scripts can still be submitted with `sbatch`. For local Docker runs, the updated helper scripts call run scripts directly with `bash` where appropriate.

---
## Reviewer Reproduction Guide

This section gives a step-by-step command list for a reviewer or new user starting from a clean Linux machine. The Docker image contains Qtools and Q6, so the reviewer does not need to install those separately on the host.

### 1. Install Docker On Linux

```bash
sudo apt update
sudo apt install -y docker.io docker-compose-plugin git
sudo systemctl enable --now docker
sudo usermod -aG docker "$USER"
newgrp docker
```

### 2. Clone The Repository

```bash
git clone https://github.com/ND7996/GPX6.git
cd GPX6
```

If the repository is already downloaded, simply enter the repository root:

```bash
cd GPX6
```

### 3. Build The Docker Image

```bash
docker compose build
```

This step downloads and installs the reproducible software environment. Qtools and Q6 are installed inside the image at `/opt/qtools` and `/opt/Q6`; they are not stored inside the GPX6 repository folder.

### 4. Check The Docker Environment

```bash
docker compose run --rm gpx6 check
```

The check command reports installed Python packages and executable tools. Expected commands include `q_genrelax.py`, `q_genfeps.py`, `q_mapper.py`, `q_analysefeps.py`, `qprep5`, `qdyn5_r8`, and `qfep5`.

### 5. Enter The Container

```bash
docker compose run --rm gpx6
```

Inside the container, work from the repository root:

```bash
cd /workspace/GPX6
```

All commands below assume the current directory is `/workspace/GPX6`.

### 6. Structure Preparation

Run the default structure-preparation script:

```bash
python prep_structures/prep_structure.py
```

Run a specific species and level without editing the script:

```bash
GPX6_SPECIES=HUMAN GPX6_LEVEL=0 python prep_structures/prep_structure.py
GPX6_SPECIES=MOUSE GPX6_LEVEL=0 python prep_structures/prep_structure.py
```

The script uses repo-relative paths such as `./prep_structures/HUMAN/level0.txt` and writes outputs to `./prep_structures/HUMAN/level0` or the matching mouse folder.

### 7. Solvation With Qprep

Run solvation using the default folder:

```bash
bash prep_structures/prep.sh
```

Run solvation for a specific folder:

```bash
PDB_DIR=./prep_structures/HUMAN/level0 bash prep_structures/prep.sh
PDB_DIR=./prep_structures/MOUSE/level0 bash prep_structures/prep.sh
```

`prep.sh` reads force-field files from:

```bash
./parameters/qoplsaa.lib
./parameters/GPX.lib
./parameters/qoplsaa_all2.prm
```

### 8. Generate Minimized PDBs From Restart Files

```bash
BASE_DIR=./prep_structures/HUMAN/level0 bash prep_structures/preprestart.sh
BASE_DIR=./prep_structures/MOUSE/level0 bash prep_structures/preprestart.sh
```

### 9. Run Minimization Helper Scripts

```bash
BASE_DIR=./prep_structures/HUMAN/level0 bash cluster_scripts/run_minim.sh
BASE_DIR=./prep_structures/MOUSE/level0 bash cluster_scripts/run_minim.sh
```

### 10. Run Replica Helper Scripts

```bash
BASE_DIR=./prep_structures/HUMAN/level0 bash cluster_scripts/run_replica.sh
BASE_DIR=./prep_structures/MOUSE/level0 bash cluster_scripts/run_replica.sh
```

For local Docker use these helpers call `bash run_qdyn_5.sh` directly where appropriate. On an HPC cluster, users may adapt the scripts back to `sbatch` submission if desired.

### 11. Execute Preparation Notebooks

```bash
jupyter nbconvert --to notebook --execute --inplace prep_scripts/minim.ipynb
jupyter nbconvert --to notebook --execute --inplace prep_scripts/fep.ipynb
```

### 12. Run Qtools FEP Analysis

Human WT:

```bash
bash analysis_scripts/humanWT/statistics_table.sh
```

Mouse WT:

```bash
bash analysis_scripts/mouseWT/statistics_table.sh
```

These scripts use `q_mapper.py`, `q_analysefeps.py`, and `qfep5` from the Docker image `PATH`.

### 13. Execute Analysis Notebooks

Human analysis notebooks:

```bash
jupyter nbconvert --to notebook --execute --inplace analysis_scripts/humanWT/distances.ipynb
jupyter nbconvert --to notebook --execute --inplace analysis_scripts/humanWT/energy.ipynb
jupyter nbconvert --to notebook --execute --inplace analysis_scripts/humanWT/plotjson.ipynb
jupyter nbconvert --to notebook --execute --inplace analysis_scripts/humanWT/parameters.ipynb
```

Mouse analysis notebooks:

```bash
jupyter nbconvert --to notebook --execute --inplace analysis_scripts/mouseWT/distances.ipynb
jupyter nbconvert --to notebook --execute --inplace analysis_scripts/mouseWT/energy.ipynb
jupyter nbconvert --to notebook --execute --inplace analysis_scripts/mouseWT/plotjson.ipynb
jupyter nbconvert --to notebook --execute --inplace analysis_scripts/mouseWT/parameters.ipynb
jupyter nbconvert --to notebook --execute --inplace analysis_scripts/mouseWT/dcd_pdb.ipynb
jupyter nbconvert --to notebook --execute --inplace analysis_scripts/mouseWT/dcd_pdb2.ipynb
jupyter nbconvert --to notebook --execute --inplace analysis_scripts/mouseWT/extractsequence.ipynb
jupyter nbconvert --to notebook --execute --inplace analysis_scripts/mouseWT/rmsf.ipynb
```

### 14. Generate Figures

Run figure scripts from the repository root:

```bash
python analysis_scripts/Scripts_to_generate_figures/heatmap.py
python analysis_scripts/Scripts_to_generate_figures/free_energy.py
python analysis_scripts/Scripts_to_generate_figures/free_energy_plots.py
python analysis_scripts/Scripts_to_generate_figures/final_electrostatics_plot.py
python analysis_scripts/Scripts_to_generate_figures/lfer_human.py
python analysis_scripts/Scripts_to_generate_figures/lfer_mouse.py
python analysis_scripts/Scripts_to_generate_figures/combined.py
python analysis_scripts/Scripts_to_generate_figures/cluster_heatmap.py
python analysis_scripts/Scripts_to_generate_figures/cluster_snn.py
python analysis_scripts/Scripts_to_generate_figures/phylogeny.py
python analysis_scripts/Scripts_to_generate_figures/pathway_network.py
python analysis_scripts/Scripts_to_generate_figures/distribution.py
python analysis_scripts/Scripts_to_generate_figures/distribution_subset.py
python analysis_scripts/Scripts_to_generate_figures/distribution_with_labels.py
python analysis_scripts/Scripts_to_generate_figures/distribution_subset_with_labels.py
python analysis_scripts/Scripts_to_generate_figures/trilateration.py
python analysis_scripts/Scripts_to_generate_figures/trilateration_with_lines.py
python analysis_scripts/Scripts_to_generate_figures/toc.py
python analysis_scripts/Scripts_to_generate_figures/bounds.py
python analysis_scripts/Scripts_to_generate_figures/cys_fluctuations.py
python analysis_scripts/Scripts_to_generate_figures/human_greedy.py
python analysis_scripts/Scripts_to_generate_figures/mouse_greedy.py
python analysis_scripts/Scripts_to_generate_figures/high_resolution.py
python analysis_scripts/Scripts_to_generate_figures/high_resolution_2.py
python analysis_scripts/Scripts_to_generate_figures/merge.py
```

MPNN-specific figure scripts require the corresponding MPNN FASTA inputs under `./analysis_scripts/MPNN`:

```bash
python analysis_scripts/Scripts_to_generate_figures/barplot_MPNN.py
python analysis_scripts/Scripts_to_generate_figures/distribution_MPNN.py
python analysis_scripts/Scripts_to_generate_figures/phylogeny_MPNN.py
python analysis_scripts/Scripts_to_generate_figures/trilateration_MPNN.py
python analysis_scripts/Scripts_to_generate_figures/trilateration_with_lines_MPNN.py
```

### 15. Build The Manuscript

```bash
make -C manuscript_files
```

Clean manuscript build files:

```bash
make -C manuscript_files clean
```

### 16. Run Commands Without Opening An Interactive Shell

Any command can be run directly through Docker, for example:

```bash
docker compose run --rm gpx6 python prep_structures/prep_structure.py
docker compose run --rm gpx6 bash prep_structures/prep.sh
docker compose run --rm gpx6 python analysis_scripts/Scripts_to_generate_figures/heatmap.py
docker compose run --rm gpx6 make -C manuscript_files
```

### 17. Reproducibility Notes

- Run commands from `/workspace/GPX6` inside Docker.
- The repository is mounted into the container, so generated outputs remain in the working tree.
- Scripts use relative paths such as `./parameters`, `./prep_structures`, `./input`, and `./analysis_scripts`.
- Qtools and Q6 are installed inside Docker at `/opt/qtools` and `/opt/Q6` and are available through `PATH`.
- Large MD outputs may take substantial runtime and storage depending on the number of replicas and systems executed.

---
## Starting Structures

| System | Source | Identifier |
|--------|--------|------------|
| Mouse GPX6 | Experimental (X-ray) | [PDB: 7FC2](https://www.rcsb.org/structure/7FC2) |
| Human GPX6 | AlphaFold model | As described in the manuscript |

Both structures were mutated using PyMOL automation (Step 1) to generate the variant panel used in this study.

---

## Simulation Parameters

Full simulation protocol, parameters, and convergence analysis are described in detail in the associated manuscript.

---

## Workflow Overview

This pipeline automates the calculation of binding free energy differences through a structured 7-step process involving:

1. Structure Preparation
2. System Solvation
3. FEP Setup
4. Molecular Dynamics (MD) Relaxation
5. Structure Minimization
6. FEP Input Generation
7. Statistical Analysis

---

## STEP 1 - Structure Preparation

### Purpose

Automates point mutations in a PDB structure using PyMOL.

### Details

- **Tool**: PyMOL automation script
- **Function**: Introduces predefined mutations (e.g., mouse-to-human substitutions or SEC incorporation)
- **Output**: Individual PDB files for each variant

### Features

- Automated mutation mapping
- Preserves backbone integrity
- Ensures mutation compatibility
- Consistent output naming

---

## STEP 2 - Solvation with Qprep

### Purpose

Prepares solvated systems and topologies for MD relaxation.

### Process

1. Extracts base name from each PDB file
2. Generates `.inp` file and runs `qprep6`
3. Applies solvation and boundary conditions

### Output

- `${base_name}_solvated.pdb`
- `${base_name}_solvated.top`

### Specifications

- **Solvation**: Explicit water model
- **Boundary Conditions**: Spherical
- **Force Field**: OPLS-AA

---

## STEP 3 - FEP File Generation

### Purpose

Uses `makeFEP.py` to create lambda-dependent FEP input files.

### Inputs

| File | Role |
|------|------|
| `fepmousecys.qmap` | Mapping for mouse WT |
| `fephumansec.qmap` | Mapping for human WT |

### Output

- FEP files with initial/final state mappings and perturbation parameters across 51 lambda windows

---

## STEP 4 - Relaxation Input Setup

### Purpose

Generates input files for energy minimization and equilibration.

### Required Inputs

- `--genrelax.proc`: Relaxation parameters
- `--top`: Topology file
- `--pdb`: Solvated structure
- `--fep`: FEP file
- `--outdir minim`: Output folder
- `--rs run_qdyn6.sh`: Qdyn execution script

### Steps

- Validates inputs
- Generates `.inp` files
- Creates run scripts
- Organizes output structure

---

## STEP 5 - Minimized PDB for FEP

### Purpose

Performs energy minimization and prepares PDBs for FEP simulations.

### Process

1. Scans system directories
2. Locates `relax_012.re` restart files
3. Runs relaxation and prepares minimized structures

### Output

- `minim.pdb` saved in each system folder

### Quality Checks

- Convergence validation
- Energy profile checks
- Structural integrity

---

## STEP 6 - FEP Input Generation

### Purpose

Generates complete FEP simulation input sets with replica support.

### Required Inputs

- `--genfeps.proc`: Parameters for FEP and equilibration
- `--pdb`: `minim.pdb`
- `--repeats`: 15 replicas
- `--frames`: 51
- `--fromlambda`: 1.0
- `--prefix`: `replica`
- `--rs`: Execution script

### Output

FEP-ready folders per replica inside each system directory.

---

## STEP 7 - Analysis with Qtools

### Purpose

Performs statistical evaluation of FEP results.

### Pipeline

#### 7.1 Data Mapping

- **Tool**: `q_mapper.py`
- **Function**: Aggregates lambda data across 15 replicas

#### 7.2 FEP Analysis

- **Tool**: `q_analysefeps.py`
- **Output**: `test.out` + JSON structured data

---

This pipeline is built to support automated FEP calculations in protein mutagenesis projects using the Q6 software suite.

---

## Data Availability

All scripts, input files, and parameters are available in this repository.
MPNN sequence design data are available at [github.com/ND7996/MPNN](https://github.com/ND7996/MPNN).

Raw simulation trajectories (`.dcd`, `.en`, `.log`) are currently stored on an external hard drive
connected to the university cluster and are being deposited to
[CSUC (Consorci de Serveis Universitaris de Catalunya)](https://www.csuc.cat/en) for long-term public archival.
A DOI will be added here upon completion. Raw data are available upon request to the corresponding author.


