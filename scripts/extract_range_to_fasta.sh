#!/bin/bash
# This script extracts gene ranges from a GFF file.

# These arguments are provided by the Python script:
# $1: Path to the input GFF/GTF file
# $2: Flank size (e.g., 1500)
# $3: Path for the output ranges text file

GFF_FILE="$1"
FLANK_SIZE="$2"
OUTPUT_RANGES_FILE="$3"

# Extract gene ranges from the GFF file
while IFS=$'\t' read -r col1 _ col3 col4 col5 _; do
    if [[ $col3 == *"gene"* ]]; then
        start=$((col4 - FLANK_SIZE))
        end=$((col5 + FLANK_SIZE))
        echo -e "$col1:$start-$end"
    fi
done < "$GFF_FILE" > "$OUTPUT_RANGES_FILE"

echo "Gene range extraction complete."
