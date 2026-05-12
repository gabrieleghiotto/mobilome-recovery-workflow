#!/bin/bash
# Step 5: Refine cluster representatives
#
# Applies two refinement rules:
# 1. Fragment reassignment
# 2. Circular trumps linear

set -e

# Configuration
CLUSTERS_DIR="data/clusters"
FILTERED_DIR="data/filtered"
OUTPUT_DIR="data/final"
THREADS=${THREADS:-16}

echo "==========================================================================="
echo "STEP 5: Refining cluster representatives"
echo "==========================================================================="

mkdir -p "$OUTPUT_DIR"

# =========================================================================
# vOTU refinement
# =========================================================================
if [ -f "$CLUSTERS_DIR/vOTUs/clusters.tsv" ] && \
   [ -f "$CLUSTERS_DIR/vOTUs/allVSall_ani.tsv" ] && \
   [ -f "$FILTERED_DIR/viruses/all_virus_summaries.tsv" ]; then

    echo ""
    echo "Refining vOTU cluster representatives..."

    python scripts/refine_votus.py \
        --clusters "$CLUSTERS_DIR/vOTUs/clusters.tsv" \
        --ani "$CLUSTERS_DIR/vOTUs/allVSall_ani.tsv" \
        --genomad "$FILTERED_DIR/viruses/all_virus_summaries.tsv" \
        --out "$OUTPUT_DIR/votu_clusters_updated_final.tsv" \
        --pid 95.0 \
        --qcov 85.0 \
        --tcov 70.0 \
        --log "$OUTPUT_DIR/votu_refinement.log"

    echo "✓ vOTU refinement complete"

    # Extract final representatives
    echo ""
    echo "Extracting final vOTU representatives..."
    python scripts/extract_representatives.py \
        --clusters "$OUTPUT_DIR/votu_clusters_updated_final.tsv" \
        --input "$FILTERED_DIR/viruses/viruses_1kb.fna" \
        --output "$OUTPUT_DIR/virus_representatives.fna"

    echo "✓ Extracted $(grep -c "^>" "$OUTPUT_DIR/virus_representatives.fna") vOTU representatives"

fi

# =========================================================================
# PTU refinement
# =========================================================================
if [ -f "$CLUSTERS_DIR/PTUs/clusters.tsv" ] && \
   [ -f "$CLUSTERS_DIR/PTUs/allVSall_ani.tsv" ] && \
   [ -f "$FILTERED_DIR/plasmids/all_plasmid_summaries.tsv" ]; then

    echo ""
    echo "Refining PTU cluster representatives..."

    python scripts/refine_ptus.py \
        --clusters "$CLUSTERS_DIR/PTUs/clusters.tsv" \
        --ani "$CLUSTERS_DIR/PTUs/allVSall_ani.tsv" \
        --genomad "$FILTERED_DIR/plasmids/all_plasmid_summaries.tsv" \
        --out "$OUTPUT_DIR/ptu_clusters_updated_final.tsv" \
        --gani 35.0 \
        --qcov 50.0 \
        --tcov 70.0 \
        --log "$OUTPUT_DIR/ptu_refinement.log"

    echo "✓ PTU refinement complete"

    # Extract final representatives
    echo ""
    echo "Extracting final PTU representatives..."
    python scripts/extract_representatives.py \
        --clusters "$OUTPUT_DIR/ptu_clusters_updated_final.tsv" \
        --input "$FILTERED_DIR/plasmids/plasmids_2kb.fna" \
        --output "$OUTPUT_DIR/plasmid_representatives.fna"

    echo "✓ Extracted $(grep -c "^>" "$OUTPUT_DIR/plasmid_representatives.fna") PTU representatives"

fi

# =========================================================================
# Generate summary statistics
# =========================================================================
echo ""
echo "Generating summary statistics..."

# vOTU stats
if [ -f "$OUTPUT_DIR/votu_clusters_updated_final.tsv" ]; then
    n_votus=$(tail -n +2 "$OUTPUT_DIR/votu_clusters_updated_final.tsv" | wc -l)
    total_uvigs=$(tail -n +2 "$OUTPUT_DIR/votu_clusters_updated_final.tsv" | awk -F'\t' '{sum += $3} END {print sum}')

    echo ""
    echo "vOTU Statistics:"
    echo "  Total vOTUs: $n_votus"
    echo "  Total UViGs: $total_uvigs"
    echo "  Average UViGs/vOTU: $(echo "scale=1; $total_uvigs / $n_votus" | bc)"
fi

# PTU stats
if [ -f "$OUTPUT_DIR/ptu_clusters_updated_final.tsv" ]; then
    n_ptus=$(tail -n +2 "$OUTPUT_DIR/ptu_clusters_updated_final.tsv" | wc -l)
    total_psequences=$(tail -n +2 "$OUTPUT_DIR/ptu_clusters_updated_final.tsv" | awk -F'\t' '{sum += $3} END {print sum}')

    echo ""
    echo "PTU Statistics:"
    echo "  Total PTUs: $n_ptus"
    echo "  Total plasmids: $total_psequences"
    echo "  Average plasmids/PTU: $(echo "scale=1; $total_psequences / $n_ptus" | bc)"
fi

echo ""
echo "==========================================================================="
echo "Step 5 complete: Cluster refinement finished"
echo "==========================================================================="
echo ""
echo "Final output files in $OUTPUT_DIR:"
echo "  ✓ votu_clusters_updated_final.tsv"
echo "  ✓ virus_representatives.fna"
if [ -f "$OUTPUT_DIR/ptu_clusters_updated_final.tsv" ]; then
    echo "  ✓ ptu_clusters_updated_final.tsv"
    echo "  ✓ plasmid_representatives.fna"
fi
echo ""
echo "Workflow complete! 🎉"
