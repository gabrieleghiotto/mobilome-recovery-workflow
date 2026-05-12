#!/bin/bash
# Step 1: Mine viral and plasmid sequences using geNomad
#
# This script runs geNomad end-to-end on metagenome sequences to recover
# both viral and plasmid genomes simultaneously.

set -e

# Configuration
GENOMAD_DB="${GENOMAD_DB:-/path/to/genomad_db}"
INPUT_DIR="data/metagenomes"
OUTPUT_DIR="data/genomad"
THREADS=${THREADS:-16}

echo "==========================================================================="
echo "STEP 1: Mining viral and plasmid sequences with geNomad"
echo "==========================================================================="

# Check inputs
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

if [ ! -d "$GENOMAD_DB" ]; then
    echo "ERROR: geNomad database not found: $GENOMAD_DB"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Count input files
n_metagenomes=$(find "$INPUT_DIR" -name "*.fna" -o -name "*.fa" | wc -l)
echo "Found $n_metagenomes metagenome(s) to process"

# Process each metagenome
sample_count=0
for metagenome in "$INPUT_DIR"/*.{fna,fa}; do

    # Skip if file doesn't exist (glob expansion issue)
    if [ ! -f "$metagenome" ]; then
        continue
    fi

    sample=$(basename "$metagenome" | sed 's/\.[^.]*$//')
    sample_output="$OUTPUT_DIR/$sample"

    echo ""
    echo "Processing: $metagenome"
    echo "  Output: $sample_output"

    # Run geNomad
    genomad end-to-end \
        --enable-score-calibration \
        --conservative \
        --lenient-taxonomy \
        --full-ictv-lineage \
        --sensitivity 7.0 \
        --cleanup \
        --threads $THREADS \
        "$metagenome" \
        "$sample_output" \
        "$GENOMAD_DB"

    echo "✓ Completed $sample"
    ((sample_count++))
done

echo ""
echo "==========================================================================="
echo "Step 1 complete: Processed $sample_count samples"
echo "==========================================================================="
echo ""
echo "Next steps:"
echo "  1. Check geNomad outputs in $OUTPUT_DIR"
echo "  2. Merge geNomad summaries: python scripts/merge_summaries.py"
echo "  3. Run Step 2: bash 02_quality_filtering.sh"
