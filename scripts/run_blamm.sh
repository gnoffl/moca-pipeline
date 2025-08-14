#!/bin/bash
# This script runs the main blamm workflow. It's called by the Python pipeline.

# These arguments are provided by the Python script:
# $1: Path to the blamm executable
# $2: Path to the input motif file (e.g., trimmed.jaspar)
# $3: Path to the sequences manifest file (sequences.mf)
# $4: P-value threshold for the scan

BLAMM_EXECUTABLE="$1"
MOTIFS_FILE="$2"
SEQUENCES_MANIFEST="$3"
P_VALUE="$4"

# Get the directory where blamm is located to run commands from there
BLAMM_DIR=$(dirname "$BLAMM_EXECUTABLE")
cd "$BLAMM_DIR"

echo "Running BLAMM dictionary creation..."
./$(basename "$BLAMM_EXECUTABLE") dict "$SEQUENCES_MANIFEST"

echo "Running BLAMM histogram creation..."
./$(basename "$BLAMM_EXECUTABLE") hist -e "$MOTIFS_FILE" "$SEQUENCES_MANIFEST"

echo "Running BLAMM scan..."
./$(basename "$BLAMM_EXECUTABLE") scan -rc -pt "$P_VALUE" "$MOTIFS_FILE" "$SEQUENCES_MANIFEST"

echo "BLAMM workflow complete."
