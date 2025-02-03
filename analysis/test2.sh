#!/bin/bash

# Define variables
MOUSE_DIR="/home/hp/results/MOUSE"
SCRIPT_DIR="/home/hp/nayanika/github/GPX6/prep_scripts"
SHELL_SCRIPT="minimselenocysteine.sh"

# Ensure the MOUSE directory exists before copying
if [ -d "/home/hp/nayanika/github/GPX6/prep_structures/MOUSE" ]; then
    echo "Copying MOUSE directory to results..."
    cp -r "/home/hp/nayanika/github/GPX6/prep_structures/MOUSE" "/home/hp/results"
else
    echo "Error: MOUSE directory not found!"
    exit 1
fi

# RUNNING MINIMIZATION
echo "Running $SHELL_SCRIPT..."
cd "$SCRIPT_DIR" || { echo "Error: Cannot change directory to $SCRIPT_DIR"; exit 1; }

# Check if the shell script exists
if [ ! -f "$SHELL_SCRIPT" ]; then
    echo "Error: Shell script $SHELL_SCRIPT not found in $SCRIPT_DIR!"
    exit 1
fi

# Ensure the shell script is executable
chmod +x "$SHELL_SCRIPT"

# Run the shell script
"$SCRIPT_DIR/$SHELL_SCRIPT"
EXIT_STATUS=$?

if [ $EXIT_STATUS -ne 0 ]; then
    echo "Error: Shell script execution failed!"
    exit $EXIT_STATUS
fi

echo "Minimization completed!"

# Ensure python is correctly linked to python3
if ! command -v python &> /dev/null; then
    echo "Python not found! Creating symlink to python3..."
    sudo ln -s /usr/bin/python3 /usr/bin/python
fi

# Verify Python installation
if ! python3 --version &> /dev/null; then
    echo "Python3 is not installed! Installing now..."
    sudo apt update && sudo apt install -y python3
fi
