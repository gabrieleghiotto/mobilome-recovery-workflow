#!/bin/bash
# Step 2: Quality filtering of viral and plasmid sequences
#
# Viruses: length ≥ 1 kb + CheckV quality assessment
# Plasmids: length ≥ 2 kb

set -e

# Configuration
INPUT_DIR="data/genomad"
OUTPUT_DIR="data/filtered"
CHECKV_DB="${CHECKV_DB:-/path/to/checkv-db-v1.5}"
THREADS=${THREADS:-16}

echo "==========================================================================="
echo "STEP 2: Quality filtering"
echo "==========================================================================="

mkdir -p "$OUTPUT_DIR/viruses"
mkdir -p "$OUTPUT_DIR/plasmids"

# =========================================================================
# VIRUSES: Length ≥ 1 kb
# =========================================================================
echo ""
echo "Filtering viruses (length ≥ 1 kb)..."

# Find and concatenate all viral sequences
viral_files=$(find "$INPUT_DIR" -name "*_virus.fna" 2>/dev/null || true)
if [ -z "$viral_files" ]; then
    echo "WARNING: No viral sequences found in $INPUT_DIR"
else
    # Concatenate all viral sequences
    cat $viral_files > "$OUTPUT_DIR/viruses/all_viruses.fna"
    echo "  Concatenated viral sequences"

    # Filter by length
    seqkit seq -m 1000 "$OUTPUT_DIR/viruses/all_viruses.fna" \
        > "$OUTPUT_DIR/viruses/viruses_1kb.fna"

    n_before=$(grep -c "^>" "$OUTPUT_DIR/viruses/all_viruses.fna" || echo 0)
    n_after=$(grep -c "^>" "$OUTPUT_DIR/viruses/viruses_1kb.fna" || echo 0)
    echo "  Length filtered: $n_before → $n_after sequences (≥1 kb)"

    # Run CheckV if database exists
    if [ -d "$CHECKV_DB" ]; then
        echo "  Running CheckV..."
        checkv end_to_end \
            "$OUTPUT_DIR/viruses/viruses_1kb.fna" \
            "$OUTPUT_DIR/viruses/checkv_out" \
            -t $THREADS \
            -d "$CHECKV_DB"
        echo "  ✓ CheckV complete"
    else
        echo "  WARNING: CheckV database not found at $CHECKV_DB"
        echo "  Skipping CheckV analysis. Set CHECKV_DB=/path/to/checkv-db-v1.5"
    fi
fi

# =========================================================================
# PLASMIDS: Length ≥ 2 kb
# =========================================================================
echo ""
echo "Filtering plasmids (length ≥ 2 kb)..."

plasmid_files=$(find "$INPUT_DIR" -name "*_plasmid.fna" 2>/dev/null || true)
if [ -z "$plasmid_files" ]; then
    echo "WARNING: No plasmid sequences found in $INPUT_DIR"
else
    # Concatenate all plasmid sequences
    cat $plasmid_files > "$OUTPUT_DIR/plasmids/all_plasmids.fna"
    echo "  Concatenated plasmid sequences"

    # Filter by length
    seqkit seq -m 2000 "$OUTPUT_DIR/plasmids/all_plasmids.fna" \
        > "$OUTPUT_DIR/plasmids/plasmids_2kb.fna"

    n_before=$(grep -c "^>" "$OUTPUT_DIR/plasmids/all_plasmids.fna" || echo 0)
    n_after=$(grep -c "^>" "$OUTPUT_DIR/plasmids/plasmids_2kb.fna" || echo 0)
    echo "  Length filtered: $n_before → $n_after sequences (≥2 kb)"
fi

# =========================================================================
# Merge geNomad summaries
# =========================================================================
echo ""
echo "Merging geNomad summaries..."

# Find and concatenate viral summaries
viral_summaries=$(find "$INPUT_DIR" -name "*_virus_summary.tsv" 2>/dev/null || true)
if [ ! -z "$viral_summaries" ]; then
    # Concatenate with headers only once
    cat $(echo "$viral_summaries" | head -1) > "$OUTPUT_DIR/viruses/all_virus_summaries.tsv"
    for file in $(echo "$viral_summaries" | tail -n +2); do
        tail -n +2 "$file" >> "$OUTPUT_DIR/viruses/all_virus_summaries.tsv"
    done
    echo "  Merged viral summaries: $(wc -l < "$OUTPUT_DIR/viruses/all_virus_summaries.tsv") records"
fi

# Find and concatenate plasmid summaries
plasmid_summaries=$(find "$INPUT_DIR" -name "*_plasmid_summary.tsv" 2>/dev/null || true)
if [ ! -z "$plasmid_summaries" ]; then
    cat $(echo "$plasmid_summaries" | head -1) > "$OUTPUT_DIR/plasmids/all_plasmid_summaries.tsv"
    for file in $(echo "$plasmid_summaries" | tail -n +2); do
        tail -n +2 "$file" >> "$OUTPUT_DIR/plasmids/all_plasmid_summaries.tsv"
    done
    echo "  Merged plasmid summaries: $(wc -l < "$OUTPUT_DIR/plasmids/all_plasmid_summaries.tsv") records"
fi

echo ""
echo "==========================================================================="
echo "Step 2 complete: Quality filtering finished"
echo "==========================================================================="
echo ""
echo "Output files:"
echo "  Viruses:  $OUTPUT_DIR/viruses/viruses_1kb.fna"
echo "  Plasmids: $OUTPUT_DIR/plasmids/plasmids_2kb.fna"
echo ""
echo "Next step: bash 03_vclust_clustering.sh"
