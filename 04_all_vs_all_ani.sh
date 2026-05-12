#!/bin/bash
# Step 4: All-vs-all ANI validation
#
# Performs BLASTn all-vs-all comparison on cluster representatives
# followed by ANI estimation to detect fragmented representatives.

set -e

# Configuration
CLUSTERS_DIR="data/clusters"
THREADS=${THREADS:-16}

echo "==========================================================================="
echo "STEP 4: All-vs-all ANI validation"
echo "==========================================================================="

# =========================================================================
# vOTU ANI validation
# =========================================================================
if [ -f "$CLUSTERS_DIR/vOTUs/representatives.fna" ]; then
    echo ""
    echo "Processing vOTU representatives..."

    REP_FILE="$CLUSTERS_DIR/vOTUs/representatives.fna"
    DB_FILE="$CLUSTERS_DIR/vOTUs/reps_db"
    BLAST_OUT="$CLUSTERS_DIR/vOTUs/blast_results.tsv"
    ANI_OUT="$CLUSTERS_DIR/vOTUs/allVSall_ani.tsv"

    # Create BLAST database
    echo "  Building BLAST database..."
    makeblastdb -in "$REP_FILE" -dbtype nucl -out "$DB_FILE" -title "vOTU representatives"

    # Run all-vs-all BLASTn
    echo "  Running all-vs-all BLASTn..."
    blastn -query "$REP_FILE" -db "$DB_FILE" \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
        -max_target_seqs 10000 \
        -num_threads $THREADS \
        -evalue 1e-5 \
        -out "$BLAST_OUT"

    n_hits=$(wc -l < "$BLAST_OUT")
    echo "  Generated $n_hits BLASTn hits"

    # Calculate ANI from BLAST results
    echo "  Calculating ANI from BLAST results..."
    python scripts/anicalc.py -i "$BLAST_OUT" -o "$ANI_OUT"

    echo "✓ vOTU ANI validation complete"
fi

# =========================================================================
# PTU ANI validation
# =========================================================================
if [ -f "$CLUSTERS_DIR/PTUs/representatives.fna" ]; then
    echo ""
    echo "Processing PTU representatives..."

    REP_FILE="$CLUSTERS_DIR/PTUs/representatives.fna"
    DB_FILE="$CLUSTERS_DIR/PTUs/reps_db"
    BLAST_OUT="$CLUSTERS_DIR/PTUs/blast_results.tsv"
    ANI_OUT="$CLUSTERS_DIR/PTUs/allVSall_ani.tsv"

    # Create BLAST database
    echo "  Building BLAST database..."
    makeblastdb -in "$REP_FILE" -dbtype nucl -out "$DB_FILE" -title "PTU representatives"

    # Run all-vs-all BLASTn
    echo "  Running all-vs-all BLASTn..."
    blastn -query "$REP_FILE" -db "$DB_FILE" \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
        -max_target_seqs 10000 \
        -num_threads $THREADS \
        -evalue 1e-5 \
        -out "$BLAST_OUT"

    n_hits=$(wc -l < "$BLAST_OUT")
    echo "  Generated $n_hits BLASTn hits"

    # Calculate ANI from BLAST results
    echo "  Calculating ANI from BLAST results..."
    python scripts/anicalc.py -i "$BLAST_OUT" -o "$ANI_OUT"

    echo "✓ PTU ANI validation complete"
fi

echo ""
echo "==========================================================================="
echo "Step 4 complete: All-vs-all ANI validation finished"
echo "==========================================================================="
echo ""
echo "Output files:"
echo "  vOTU ANI:  $CLUSTERS_DIR/vOTUs/allVSall_ani.tsv"
echo "  PTU ANI:   $CLUSTERS_DIR/PTUs/allVSall_ani.tsv (if applicable)"
echo ""
echo "Next step: bash 05_refine_representatives.sh"
