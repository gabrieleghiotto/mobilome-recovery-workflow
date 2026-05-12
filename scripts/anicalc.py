#!/usr/bin/env python3
"""
ANI (Average Nucleotide Identity) calculator from BLAST results.

This script calculates ANI values from BLAST tabular output.
Adapted from CheckV scripts for use in this workflow.

ANI is calculated as:
  ANI = (number of matching bases) / (BLAST alignment length) * 100
"""

import argparse
import logging
from pathlib import Path
from typing import Dict, List, Tuple
import pandas as pd


def calculate_ani(
    percent_identity: float,
    alignment_length: int,
    query_length: int,
    subject_length: int
) -> Tuple[float, float, float]:
    """
    Calculate ANI and coverage metrics from BLAST alignment.

    Parameters
    ----------
    percent_identity : float
        Percent identity from BLAST (0-100)
    alignment_length : int
        Length of BLAST alignment
    query_length : int
        Total query sequence length
    subject_length : int
        Total subject sequence length

    Returns
    -------
    tuple
        (ani, qcov, scov) - ANI value and coverages
    """
    ani = percent_identity if alignment_length > 0 else 0.0
    qcov = 100 * alignment_length / query_length if query_length > 0 else 0.0
    scov = 100 * alignment_length / subject_length if subject_length > 0 else 0.0

    return ani, qcov, scov


def process_blast_results(
    blast_file: str,
    output_file: str,
    min_ali_len: int = 0,
    min_pid: float = 0.0
) -> None:
    """
    Process BLAST results to calculate ANI values.

    Parameters
    ----------
    blast_file : str
        Input BLAST tabular file (format: qseqid sseqid pident ... qlen slen)
    output_file : str
        Output ANI results file
    min_ali_len : int
        Minimum alignment length filter
    min_pid : float
        Minimum percent identity filter
    """
    # Read BLAST results
    column_names = [
        'query', 'subject', 'pid', 'length',
        'mismatch', 'gapopen', 'qstart', 'qend',
        'sstart', 'send', 'evalue', 'bitscore',
        'qlen', 'slen'
    ]

    try:
        df = pd.read_csv(blast_file, sep='\t', names=column_names, dtype={
            'query': str, 'subject': str, 'pid': float,
            'length': int, 'qlen': int, 'slen': int
        })
    except Exception as e:
        logging.error(f"Error reading BLAST file: {e}")
        return

    logging.info(f"Loaded {len(df)} BLAST results")

    # Apply filters
    df = df[df['length'] >= min_ali_len]
    df = df[df['pid'] >= min_pid]
    logging.info(f"After filtering: {len(df)} results")

    # Calculate ANI metrics
    ani_data = []
    for _, row in df.iterrows():
        ani, qcov, scov = calculate_ani(
            row['pid'], row['length'], row['qlen'], row['slen']
        )

        # Calculate global ANI (gANI) - symmetric measure
        gani = (ani * min(row['qlen'], row['slen'])) / max(row['qlen'], row['slen'])

        ani_data.append({
            'query': row['query'],
            'subject': row['subject'],
            'ani': round(ani, 2),
            'qcov': round(qcov, 2),
            'tcov': round(scov, 2),  # target coverage (tcov)
            'pid': round(row['pid'], 2),
            'gani': round(gani, 2),
            'ali_len': row['length'],
            'evalue': row['evalue']
        })

    # Write output
    result_df = pd.DataFrame(ani_data)
    result_df.to_csv(output_file, sep='\t', index=False)
    logging.info(f"Wrote {len(result_df)} ANI results to {output_file}")


def main():
    """Main execution."""
    parser = argparse.ArgumentParser(
        description='Calculate ANI from BLAST results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python anicalc.py -i blast_results.tsv -o ani_results.tsv
  python anicalc.py -i blast.tsv -o ani.tsv --min-ali-len 500
        """
    )

    parser.add_argument('-i', '--input', required=True, help='BLAST tabular file')
    parser.add_argument('-o', '--output', required=True, help='Output ANI file')
    parser.add_argument('--min-ali-len', type=int, default=0, help='Minimum alignment length')
    parser.add_argument('--min-pid', type=float, default=0.0, help='Minimum percent identity')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')

    args = parser.parse_args()

    # Setup logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=level, format='%(levelname)s: %(message)s')

    # Check input
    if not Path(args.input).exists():
        logging.error(f"File not found: {args.input}")
        exit(1)

    # Process
    logging.info(f"Processing BLAST results from {args.input}")
    process_blast_results(
        args.input, args.output,
        min_ali_len=args.min_ali_len,
        min_pid=args.min_pid
    )

    logging.info("Complete!")


if __name__ == "__main__":
    main()
