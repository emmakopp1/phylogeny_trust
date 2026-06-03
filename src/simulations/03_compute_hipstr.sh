#!/bin/bash
# ------------------------------------------------------------------------------
# Script Name: compute_hipstr.sh
# Description: compute the HIPSTR tree for all ages, simulations.
#              Change SIMULATION_FOLDER to choose another simulation study
# Usage: bash src/simulations/03_compute_hipstr.sh
# ------------------------------------------------------------------------------

TREEANNOTATOR="/Applications/BEAST v10.5.0/bin/treeannotator"
BURNIN=1000
HEIGHTS="mean"

# to be referenced by the user
# path of the simulation folder you want to analyse
SIMULATION_FOLDER="data/simulated-2025-07-28"
#SIMULATION_FOLDER="data/simulated-2025-07-22-12000"
#SIMULATION_FOLDER="data/simulated-2025-07-22-6000"
#SIMULATION_FOLDER="data/simulated-2025-07-22-1500"

# find all .trees files recursively
TREES_FILES=$(find "$SIMULATION_FOLDER" -name "*.trees")

for INPUT_PATH in $TREES_FILES; do
    # extract tree age from filename (number before .trees)
    BASENAME=$(basename "$INPUT_PATH" .trees)
    TREE_AGE="${BASENAME##*-}"

    # build output path: replace ctmc-strict-bd-{age}.trees -> hipstr-{age}.tree
    DIR=$(dirname "$INPUT_PATH")
    OUTPUT_PATH="${DIR}/hipstr-${TREE_AGE}.tree"

    if [ -f "$OUTPUT_PATH" ]; then
        echo "Skipping (already exists): $OUTPUT_PATH"
        continue
    fi

    echo "Processing: $INPUT_PATH -> $OUTPUT_PATH"

    "$TREEANNOTATOR" -type hipstr -burnin "$BURNIN" -heights "$HEIGHTS" "$INPUT_PATH" "$OUTPUT_PATH"
done

echo "Done."
