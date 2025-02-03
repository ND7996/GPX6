#!/bin/bash

# Base directory where minim folders are located
BASE_SCR_DIR="/home/hp/results/MOUSE/48_47_52_99_54_144_177_74_178_143_87_142_104_102_139_24_181_4_60_107"

# Paths to required files
LIB1="/home/hp/nayanika/github/GPX6/parameters/qoplsaa.lib"
LIB2="/home/hp/nayanika/github/GPX6/parameters/GPX.lib"
PRM="/home/hp/nayanika/github/GPX6/parameters/qoplsaa_all2.prm"
RELAX_FILE="relax_012.re"

# Directories for input files
RUN_DIR="/home/hp/nayanika/github/GPX6/cluster_scripts"
GENFEP_DIR="/home/hp/nayanika/github/GPX6/input"
TOPOLOGY_DIR="/home/hp/nayanika/github/GPX6/prep_structures/MOUSE/48_47_52_99_54_144_177_74_178_143_87_142_104_102_139_24_181_4_60_107"
FEP_DIR="/home/hp/nayanika/github/GPX6/input/fep"

# Iterate over each system folder within the base directory
for SYSTEM_NAME in $(ls "$BASE_SCR_DIR"); do
    SYSTEM_DIR="$BASE_SCR_DIR/$SYSTEM_NAME"
    MINIM_DIR="$SYSTEM_DIR/minim"
    
    # Ensure the minim directory exists
    if [ ! -d "$MINIM_DIR" ]; then
        echo "Skipping $SYSTEM_DIR as minim directory does not exist."
        continue
    fi

    # Find the .top file in the minim directory
    TOP_FILE=$(find "$MINIM_DIR" -maxdepth 1 -name "*.top" | head -n 1)
    
    if [ -f "$TOP_FILE" ]; then
        INPUT_FILE="$MINIM_DIR/generate_minim.inp"
        OUTPUT_PDB="$MINIM_DIR/minim.pdb"
        
        # Write the commands to the input file
        cat <<EOF > "$INPUT_FILE"
! Read in library files
readlib $LIB1
readlib $LIB2
! Read in parameter files
readprm $PRM
! Read in structure
readtop $(basename "$TOP_FILE")
rx $RELAX_FILE
writepdb $(basename "$OUTPUT_PDB") y
quit
EOF
        
        # Navigate to the minim directory and run qprep5
        cd "$MINIM_DIR" || exit
        qprep5 < generate_minim.inp  
        
        # Log the status
        if [ -f "$OUTPUT_PDB" ]; then
            echo "Generated $OUTPUT_PDB successfully."
        else
            echo "Failed to generate $OUTPUT_PDB in $MINIM_DIR."
        fi
        
        # Cleanup
        rm -f "$INPUT_FILE"
    else
        echo "No .top file found in $MINIM_DIR."
        continue
    fi
    
    # Define the paths for the required files
    PDB_FILE="$TOPOLOGY_DIR/$SYSTEM_NAME.pdb"
    TOPOLOGY_FILE="$TOPOLOGY_DIR/$SYSTEM_NAME.top"
    RELAX_INPUT_PATH="$MINIM_DIR/relax_012.inp"
    
    # Check if the required files exist
    if [ ! -f "$PDB_FILE" ]; then
        echo "No .pdb file found for $SYSTEM_NAME at $PDB_FILE"
        continue
    fi
    if [ ! -f "$TOPOLOGY_FILE" ]; then
        echo "No topology file found for $SYSTEM_NAME at $TOPOLOGY_FILE"
        continue
    fi
    if [ ! -f "$RELAX_INPUT_PATH" ]; then
        echo "Error: Relax input file '$RELAX_INPUT_PATH' doesn't exist for $SYSTEM_NAME."
        continue
    fi

    # Construct the command to run q_genfeps.py
    CMD=(
        "q_genfeps.py"
        "$GENFEP_DIR/selgenfeps.proc"
        "--pdb" "$PDB_FILE"
        "$RELAX_INPUT_PATH"
        "relax"
        "--repeats" "20"
        "--frames" "51"
        "--fromlambda" "1.0"
        "--prefix" "$SYSTEM_DIR/replica"
        "--rs" "$RUN_DIR/run_qdyn_5.sh"
    )

    # Run the command
    "${CMD[@]}"
    if [ $? -eq 0 ]; then
        echo "Successfully processed $SYSTEM_NAME with q_genfeps.py."
    else
        echo "Error processing $SYSTEM_NAME."
    fi

done

echo "Processing completed for all systems."
