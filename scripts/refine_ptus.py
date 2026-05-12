#!/usr/bin/env python3
"""
Refine PTU (plasmid taxonomic unit) cluster representatives.

This script applies two refinement rules to PTU clusters:

1. Fragment reassignment: If a representative sequence is a fragment of a longer
   sequence (gANI ≥ 35%, adjusted thresholds for plasmid diversity), reassign it
   to the longer sequence. Transitive chains are resolved to terminal targets.

2. Circular trumps linear: If a cluster's linear representative contains circular
   members, promote the longest circular member instead.

Author: Gabriele Ghiotto
"""

import argparse
import logging
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple


def setup_logging(output_file: str = None) -> logging.Logger:
    """Setup logging configuration."""
    log_level = logging.INFO
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

    handlers = [logging.StreamHandler()]
    if output_file:
        handlers.append(logging.FileHandler(output_file))

    logging.basicConfig(
        level=log_level,
        format=log_format,
        handlers=handlers
    )
    return logging.getLogger(__name__)


def load_clusters(cluster_file: str) -> Dict[str, List[str]]:
    """Load cluster assignments from vclust output."""
    clusters = {}
    with open(cluster_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                cluster_id = parts[0]
                member_id = parts[1]
                if cluster_id not in clusters:
                    clusters[cluster_id] = []
                clusters[cluster_id].append(member_id)

    return clusters


def load_ani_results(ani_file: str) -> Dict[Tuple[str, str], Dict]:
    """Load ANI results from all-vs-all BLASTn + anicalc.py."""
    ani_results = {}

    try:
        df = pd.read_csv(ani_file, sep='\t', dtype={'query': str, 'subject': str})
    except FileNotFoundError:
        logging.error(f"ANI file not found: {ani_file}")
        return {}
    except Exception as e:
        logging.error(f"Error reading ANI file {ani_file}: {e}")
        return {}

    for _, row in df.iterrows():
        query = str(row['query']).strip()
        subject = str(row['subject']).strip()
        key = (query, subject)

        ani_results[key] = {
            'gani': float(row.get('gani', 0)),  # Global ANI for plasmids
            'ani': float(row.get('ani', 0)),
            'qcov': float(row.get('qcov', 0)),
            'tcov': float(row.get('tcov', 0)),
            'pid': float(row.get('pid', 0))
        }

    return ani_results


def load_genomad_metadata(genomad_file: str) -> Dict[str, Dict]:
    """Load geNomad summary metadata for plasmids."""
    metadata = {}

    try:
        df = pd.read_csv(genomad_file, sep='\t', dtype={'sequence_name': str})
    except FileNotFoundError:
        logging.error(f"geNomad metadata file not found: {genomad_file}")
        return {}

    for _, row in df.iterrows():
        seq_id = str(row['sequence_name']).strip()
        metadata[seq_id] = {
            'length': int(row.get('length', 0)),
            'topology': str(row.get('topology', 'Unknown')),
            '# plasmid genes': int(row.get('# plasmid genes', 0))
        }

    return metadata


def resolve_fragment_chains(
    fragment_map: Dict[str, str],
    logger: logging.Logger
) -> Dict[str, str]:
    """Resolve transitive fragment assignments."""
    resolved = {}

    for fragment, target in fragment_map.items():
        visited = set()
        current = target

        while current in fragment_map and current not in visited:
            visited.add(current)
            current = fragment_map[current]

        if current in visited and current in fragment_map:
            logger.warning(f"Potential cycle in fragment chain starting at {fragment}")
            resolved[fragment] = target
        else:
            resolved[fragment] = current

    return resolved


def reassign_fragment_representatives(
    clusters: Dict[str, List[str]],
    ani_results: Dict[Tuple[str, str], Dict],
    genomad_metadata: Dict[str, Dict],
    logger: logging.Logger,
    gani_threshold: float = 35.0,
    qcov_threshold: float = 50.0,
    tcov_threshold: float = 70.0
) -> Tuple[Dict[str, str], Dict[str, List[str]]]:
    """
    Identify and resolve fragment plasmid representatives.

    Uses global ANI (gANI) which is appropriate for the greater diversity
    in plasmid sequences compared to viruses.

    Parameters
    ----------
    clusters : dict
        Cluster assignments
    ani_results : dict
        ANI/gANI results from all-vs-all comparison
    genomad_metadata : dict
        Sequence metadata
    logger : logging.Logger
        Logger instance
    gani_threshold : float
        Minimum global ANI for fragment reassignment
    qcov_threshold : float
        Minimum query coverage
    tcov_threshold : float
        Maximum target coverage
    """
    fragment_map = {}
    logger.info(f"\nReassigning fragment representatives...")
    logger.info(f"  Thresholds: gANI≥{gani_threshold}, qcov≥{qcov_threshold}, tcov<{tcov_threshold}")

    n_reassigned = 0

    for cluster_id, members in clusters.items():
        if len(members) < 2:
            continue

        rep = members[0]
        rep_length = genomad_metadata.get(rep, {}).get('length', 0)

        best_target = None
        best_metrics = {'gani': 0, 'qcov': 0, 'tcov': 100}

        for member in members[1:]:
            member_length = genomad_metadata.get(member, {}).get('length', 0)

            if member_length <= rep_length:
                continue

            ani_key = (rep, member)
            if ani_key not in ani_results:
                continue

            metrics = ani_results[ani_key]
            gani = metrics.get('gani', 0)
            qcov = metrics.get('qcov', 0)
            tcov = metrics.get('tcov', 100)

            # For plasmids, use gANI threshold
            if (gani >= gani_threshold and
                qcov >= qcov_threshold and
                tcov < tcov_threshold):

                target_length = member_length
                if best_target is None or target_length > genomad_metadata.get(
                    best_target, {}
                ).get('length', 0):
                    best_target = member
                    best_metrics = metrics

        if best_target:
            logger.info(
                f"  Reassigning {rep} (len={rep_length}) to {best_target} "
                f"(len={genomad_metadata.get(best_target, {}).get('length', 0)}) "
                f"[gANI={best_metrics['gani']:.1f}, qcov={best_metrics['qcov']:.1f}]"
            )
            fragment_map[rep] = best_target
            n_reassigned += 1

    if fragment_map:
        fragment_map = resolve_fragment_chains(fragment_map, logger)

    logger.info(f"  Reassigned {n_reassigned} representatives")

    return fragment_map, clusters


def promote_circular_representatives(
    clusters: Dict[str, List[str]],
    genomad_metadata: Dict[str, Dict],
    logger: logging.Logger
) -> Dict[str, str]:
    """Promote circular members as representatives if cluster rep is linear."""
    promotions = {}
    logger.info(f"\nPromoting circular representatives...")

    n_promoted = 0

    for cluster_id, members in clusters.items():
        if len(members) < 2:
            continue

        rep = members[0]
        rep_topology = genomad_metadata.get(rep, {}).get('topology', 'Unknown')

        if rep_topology != 'Linear':
            continue

        circular_members = [
            m for m in members
            if genomad_metadata.get(m, {}).get('topology', '') in ['DTR', 'ITR', 'Circular']
        ]

        if not circular_members:
            continue

        best_circular = max(
            circular_members,
            key=lambda m: genomad_metadata.get(m, {}).get('length', 0)
        )

        best_circular_length = genomad_metadata.get(best_circular, {}).get('length', 0)
        rep_length = genomad_metadata.get(rep, {}).get('length', 0)

        logger.info(
            f"  Promoting {best_circular} (len={best_circular_length}, circular) "
            f"over {rep} (len={rep_length}, linear) in cluster {cluster_id}"
        )

        promotions[cluster_id] = best_circular
        n_promoted += 1

    logger.info(f"  Promoted {n_promoted} circular representatives")

    return promotions


def write_refined_clusters(
    clusters: Dict[str, List[str]],
    fragment_map: Dict[str, str],
    promotions: Dict[str, str],
    output_file: str,
    logger: logging.Logger
) -> None:
    """Write refined cluster assignments."""
    rows = []

    for cluster_id, members in clusters.items():
        if not members:
            continue

        new_rep = members[0]

        if cluster_id in promotions:
            new_rep = promotions[cluster_id]

        if new_rep in fragment_map:
            new_rep = fragment_map[new_rep]

        rows.append({
            'cluster_id': cluster_id,
            'representative_id': new_rep,
            'n_members': len(members),
            'members': ','.join(members)
        })

    df = pd.DataFrame(rows)
    df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"\nWrote refined clusters to {output_file}")


def main():
    """Main execution."""
    parser = argparse.ArgumentParser(
        description='Refine PTU cluster representatives',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python refine_ptus.py \\
    --clusters PTUs/clusters.tsv \\
    --ani PTUs/allVSall_ani.tsv \\
    --genomad PTUs/all_plasmid_summaries.tsv \\
    --out PTUs/ptu_clusters_updated_final.tsv
        """
    )

    parser.add_argument('--clusters', required=True, help='vclust clusters.tsv')
    parser.add_argument('--ani', required=True, help='ANI results TSV from anicalc.py')
    parser.add_argument('--genomad', required=True, help='geNomad summary metadata TSV')
    parser.add_argument('--out', required=True, help='Output refined clusters TSV')
    parser.add_argument('--gani', type=float, default=35.0, help='Min global ANI (default: 35.0)')
    parser.add_argument('--qcov', type=float, default=50.0, help='Min query coverage (default: 50.0)')
    parser.add_argument('--tcov', type=float, default=70.0, help='Max target coverage (default: 70.0)')
    parser.add_argument('--log', help='Log file (optional)')

    args = parser.parse_args()

    # Setup logging
    logger = setup_logging(args.log)
    logger.info("="*70)
    logger.info("PLASMID PTU CLUSTER REFINEMENT")
    logger.info("="*70)

    # Verify input files
    for file in [args.clusters, args.ani, args.genomad]:
        if not Path(file).exists():
            logger.error(f"File not found: {file}")
            exit(1)

    # Load data
    logger.info(f"\nLoading cluster assignments from {args.clusters}...")
    clusters = load_clusters(args.clusters)
    logger.info(f"  Loaded {len(clusters)} clusters")

    logger.info(f"\nLoading ANI results from {args.ani}...")
    ani_results = load_ani_results(args.ani)
    logger.info(f"  Loaded {len(ani_results)} ANI comparisons")

    logger.info(f"\nLoading geNomad metadata from {args.genomad}...")
    genomad_metadata = load_genomad_metadata(args.genomad)
    logger.info(f"  Loaded metadata for {len(genomad_metadata)} sequences")

    # Apply refinement rules
    fragment_map, _ = reassign_fragment_representatives(
        clusters, ani_results, genomad_metadata, logger,
        gani_threshold=args.gani,
        qcov_threshold=args.qcov,
        tcov_threshold=args.tcov
    )

    promotions = promote_circular_representatives(clusters, genomad_metadata, logger)

    # Write output
    write_refined_clusters(
        clusters, fragment_map, promotions, args.out, logger
    )

    logger.info("\n" + "="*70)
    logger.info("REFINEMENT COMPLETE")
    logger.info("="*70)


if __name__ == "__main__":
    main()
