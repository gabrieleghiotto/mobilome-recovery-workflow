#!/bin/bash
# Step 3: ANI-based clustering with vclust
#
# vOTUs: ANI ≥ 95%, qcov ≥ 85%, Leiden resolution = 1
# PTUs: gANI ≥ 35%, Leiden resolution = 0.9

set -e

# Configuration
INPUT_DIR="data/filtered"
OUTPUT_DIR="data/clusters"
THREADS=${THREADS:-16}

echo "==========================================================================="
echo "STEP 3: ANI-based clustering with vclust"
echo "==========================================================================="

mkdir -p "$OUTPUT_DIR/vOTUs"
mkdir -p "$OUTPUT_DIR/PTUs"

# =========================================================================
# vOTU CLUSTERING (Viruses)
# =========================================================================
echo ""
echo "Clustering viruses into vOTUs..."

VIRUS_INPUT="$INPUT_DIR/viruses/viruses_1kb.fna"

if [ ! -f "$VIRUS_INPUT" ]; then
    echo "ERROR: Viral input file not found: $VIRUS_INPUT"
    exit 1
fi

# Step 1: Prefilter (k-mer based)
echo "  [1/3] Prefilter - candidate pair selection (min-ident 0.95)"
vclust prefilter \
    -i "$VIRUS_INPUT" \
    -o "$OUTPUT_DIR/vOTUs/filter.txt" \
    --min-ident 0.95

# Step 2: Align (pairwise ANI calculation)
echo "  [2/3] Align - pairwise ANI computation"
vclust align \
    -i "$VIRUS_INPUT" \
    -o "$OUTPUT_DIR/vOTUs/ani.tsv" \
    --filter "$OUTPUT_DIR/vOTUs/filter.txt" \
    --threads $THREADS

# Step 3: Cluster (Leiden algorithm on ANI graph)
echo "  [3/3] Cluster - Leiden community detection"
vclust cluster \
    -i "$OUTPUT_DIR/vOTUs/ani.tsv" \
    -o "$OUTPUT_DIR/vOTUs/clusters.tsv" \
    --ids "$OUTPUT_DIR/vOTUs/ani.ids.tsv" \
    --algorithm leiden \
    --metric ani \
    --ani 0.95 \
    --qcov 0.85 \
    --leiden-resolution 1.0 \
    --out-repr

n_votus=$(tail -n +2 "$OUTPUT_DIR/vOTUs/clusters.tsv" 2>/dev/null | cut -f1 | sort -u | wc -l || echo "?")
echo "✓ vOTU clustering complete: ~$n_votus clusters"

# =========================================================================
# PTU CLUSTERING (Plasmids)
# =========================================================================
echo ""
echo "Clustering plasmids into PTUs..."

PLASMID_INPUT="$INPUT_DIR/plasmids/plasmids_2kb.fna"

if [ ! -f "$PLASMID_INPUT" ]; then
    echo "WARNING: Plasmid input file not found: $PLASMID_INPUT"
    echo "Skipping PTU clustering"
else
    # Step 1: Prefilter (more lenient for plasmids)
    echo "  [1/3] Prefilter - candidate pair selection"
    vclust prefilter \
        -i "$PLASMID_INPUT" \
        -o "$OUTPUT_DIR/PTUs/filter.txt" \
        --min-kmers 20 \
        --min-ident 0.5

    # Step 2: Align
    echo "  [2/3] Align - pairwise ANI computation"
    vclust align \
        -i "$PLASMID_INPUT" \
        -o "$OUTPUT_DIR/PTUs/ani.tsv" \
        --filter "$OUTPUT_DIR/PTUs/filter.txt" \
        --out-ani 0.70 \
        --out-qcov 0.50 \
        --threads $THREADS

    # Step 3: Cluster (using global ANI)
    echo "  [3/3] Cluster - Leiden community detection (gANI)"
    vclust cluster \
        -i "$OUTPUT_DIR/PTUs/ani.tsv" \
        -o "$OUTPUT_DIR/PTUs/clusters.tsv" \
        --ids "$OUTPUT_DIR/PTUs/ani.ids.tsv" \
        --algorithm leiden \
        --metric gani \
        --gani 0.35 \
        --leiden-resolution 0.9 \
        --out-repr

    n_ptus=$(tail -n +2 "$OUTPUT_DIR/PTUs/clusters.tsv" 2>/dev/null | cut -f1 | sort -u | wc -l || echo "?")
    echo "✓ PTU clustering complete: ~$n_ptus clusters"
fi

# =========================================================================
# Extract representatives for next step
# =========================================================================
echo ""
echo "Extracting cluster representatives..."

# vOTU representatives
if [ -f "$OUTPUT_DIR/vOTUs/clusters.tsv" ]; then
    python scripts/extract_representatives.py \
        --clusters "$OUTPUT_DIR/vOTUs/clusters.tsv" \
        --input "$VIRUS_INPUT" \
        --output "$OUTPUT_DIR/vOTUs/representatives.fna"
    echo "✓ Extracted vOTU representatives"
fi

# PTU representatives
if [ -f "$OUTPUT_DIR/PTUs/clusters.tsv" ]; then
    python scripts/extract_representatives.py \
        --clusters "$OUTPUT_DIR/PTUs/clusters.tsv" \
        --input "$PLASMID_INPUT" \
        --output "$OUTPUT_DIR/PTUs/representatives.fna"
    echo "✓ Extracted PTU representatives"
fi

echo ""
echo "==========================================================================="
echo "Step 3 complete: Clustering finished"
echo "==========================================================================="
echo ""
echo "Output files:"
echo "  vOTU clusters:     $OUTPUT_DIR/vOTUs/clusters.tsv"
echo "  vOTU representatives: $OUTPUT_DIR/vOTUs/representatives.fna"
if [ -f "$OUTPUT_DIR/PTUs/clusters.tsv" ]; then
    echo "  PTU clusters:      $OUTPUT_DIR/PTUs/clusters.tsv"
    echo "  PTU representatives: $OUTPUT_DIR/PTUs/representatives.fna"
fi
echo ""
echo "Next step: bash 04_all_vs_all_ani.sh"
