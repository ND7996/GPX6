#!/bin/bash

# Iterate over replicas from replica000 to replica015
for i in $(seq -f "replica%03g" 0 15); do
    if [[ -d "$i" && -f "$i/run_qdyn_5.sh" ]]; then
        echo "Submitting job in $i..."
        (cd "$i" && sbatch run_qdyn_5.sh)
    else
        echo "Skipping $i (not found or missing run_qdyn_5.sh)"
    fi
done


