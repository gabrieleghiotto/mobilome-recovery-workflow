#!/usr/bin/env python3
"""
Extract cluster representatives from FASTA file based on cluster assignments.

Reads cluster TSV and extracts the representative sequence for each cluster
from the complete FASTA file.
"""

import argparse
from pathlib import Path
from typing import Dict, List
import logging


def setup_logging(verbose: bool = False) -> logging.Logger:
    """Setup logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger(__name__)


def parse_fasta(fasta_file: str) -> Dict[str, str]:
    """Parse FASTA file into dict {seq_id: sequence}."""
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_file) as f:
        for line in f:
            line = line.rstrip()

            if line.startswith('>'):
                # Save previous sequence
                if current_id:
                    sequences[current_id] = ''.join(current_seq)

                # Parse new header
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_id:
            sequences[current_id] = ''.join(current_seq)

    return sequences


def load_clusters(cluster_file: str) -> Dict[str, str]:
    """Load clusters and return {cluster_id: representative_id}."""
    representatives = {}

    with open(cluster_file) as f:
        # Skip header if present
        header = f.readline()

        for line in f:
            parts = line.strip().split('\t')

            if len(parts) >= 2:
                cluster_id = parts[0]
                rep_id = parts[1]
                representatives[cluster_id] = rep_id

    return representatives


def extract_representatives(
    cluster_file: str,
    input_fasta: str,
    output_fasta: str,
    logger: logging.Logger
) -> None:
    """
    Extract representative sequences from FASTA file.

    Parameters
    ----------
    cluster_file : str
        Cluster assignments file
    input_fasta : str
        Complete FASTA file
    output_fasta : str
        Output FASTA with representatives
    logger : logging.Logger
        Logger instance
    """
    # Load data
    logger.info(f"Loading sequences from {input_fasta}...")
    sequences = parse_fasta(input_fasta)
    logger.info(f"  Loaded {len(sequences)} sequences")

    logger.info(f"Loading cluster assignments from {cluster_file}...")
    representatives = load_clusters(cluster_file)
    logger.info(f"  Loaded {len(representatives)} clusters")

    # Extract representatives
    missing = 0
    found = 0
    rep_sequences = {}

    for cluster_id, rep_id in representatives.items():
        if rep_id not in sequences:
            logger.warning(f"Representative {rep_id} not found in FASTA (cluster {cluster_id})")
            missing += 1
            continue

        rep_sequences[rep_id] = sequences[rep_id]
        found += 1

    logger.info(f"Found {found}/{len(representatives)} representatives ({missing} missing)")

    # Write output
    with open(output_fasta, 'w') as f:
        for seq_id, seq in rep_sequences.items():
            f.write(f">{seq_id}\n")
            # Write sequence in 80-char lines
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')

    logger.info(f"Wrote {len(rep_sequences)} representatives to {output_fasta}")


def main():
    """Main execution."""
    parser = argparse.ArgumentParser(
        description='Extract cluster representatives from FASTA'
    )
    parser.add_argument('--clusters', required=True, help='Cluster assignments TSV')
    parser.add_argument('--input', required=True, help='Input FASTA file')
    parser.add_argument('--output', required=True, help='Output FASTA file')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')

    args = parser.parse_args()

    logger = setup_logging(args.verbose)

    # Verify inputs
    for file in [args.clusters, args.input]:
        if not Path(file).exists():
            logger.error(f"File not found: {file}")
            exit(1)

    # Extract
    extract_representatives(args.clusters, args.input, args.output, logger)

    logger.info("Complete!")


if __name__ == "__main__":
    main()
