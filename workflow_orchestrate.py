#!/usr/bin/env python3
"""
Workflow orchestration script.

This script coordinates the entire viral and plasmid sequence recovery,
filtering, clustering, and refinement pipeline.

Can be used to run individual steps or the complete workflow.
"""

import argparse
import subprocess
import sys
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional
import yaml
import os


class WorkflowOrchestrator:
    """Orchestrates the complete workflow."""

    def __init__(self, config_file: str = None, output_dir: str = None, verbose: bool = False):
        """
        Initialize orchestrator.

        Parameters
        ----------
        config_file : str
            Path to configuration YAML file
        output_dir : str
            Override output directory
        verbose : bool
            Verbose logging
        """
        self.config = self._load_config(config_file)
        self.output_dir = output_dir or self.config.get('output_base', 'data')
        self.repo_root = Path(__file__).parent
        self.verbose = verbose

        self._setup_logging()

    def _load_config(self, config_file: str = None) -> Dict[str, Any]:
        """Load YAML configuration."""
        if config_file is None:
            config_file = self.repo_root / 'config' / 'workflow_config.yaml'

        try:
            with open(config_file) as f:
                config = yaml.safe_load(f)
            return config
        except FileNotFoundError:
            logging.warning(f"Config file not found: {config_file}, using defaults")
            return {}

    def _setup_logging(self):
        """Setup logging."""
        level = logging.DEBUG if self.verbose else logging.INFO
        log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(level)
        console_handler.setFormatter(logging.Formatter(log_format))

        # File handler
        log_dir = Path(self.output_dir) / 'logs'
        log_dir.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_dir / 'workflow.log')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter(log_format))

        # Root logger
        logging.basicConfig(
            level=level,
            handlers=[console_handler, file_handler],
            format=log_format
        )

        self.logger = logging.getLogger(__name__)

    def _run_bash_script(self, script_name: str, env: Dict = None) -> int:
        """
        Run a bash script.

        Parameters
        ----------
        script_name : str
            Name of bash script (e.g., '01_genomad.sh')
        env : dict
            Environment variables to set

        Returns
        -------
        int
            Return code (0 = success)
        """
        script_path = self.repo_root / script_name

        if not script_path.exists():
            self.logger.error(f"Script not found: {script_path}")
            return 1

        # Setup environment
        exec_env = os.environ.copy()
        if env:
            exec_env.update(env)

        # Add repo root to PATH so scripts can find other scripts
        exec_env['PATH'] = f"{self.repo_root / 'scripts'}{os.pathsep}{exec_env['PATH']}"

        # Run script
        self.logger.info(f"Running {script_name}...")
        result = subprocess.run(
            ['bash', str(script_path)],
            env=exec_env,
            cwd=self.repo_root
        )

        if result.returncode == 0:
            self.logger.info(f"✓ {script_name} completed successfully")
        else:
            self.logger.error(f"✗ {script_name} failed with code {result.returncode}")

        return result.returncode

    def run_genomad(self) -> int:
        """Run Step 1: geNomad mining."""
        self.logger.info("\n" + "="*70)
        self.logger.info("STEP 1: geNomad mining")
        self.logger.info("="*70)

        env = {
            'GENOMAD_DB': self.config.get('genomad_db', '/path/to/genomad_db'),
            'THREADS': str(self.config.get('execution', {}).get('num_threads', 16))
        }

        return self._run_bash_script('01_genomad.sh', env)

    def run_quality_filtering(self) -> int:
        """Run Step 2: Quality filtering."""
        self.logger.info("\n" + "="*70)
        self.logger.info("STEP 2: Quality filtering")
        self.logger.info("="*70)

        env = {
            'CHECKV_DB': self.config.get('checkv_db', '/path/to/checkv_db'),
            'THREADS': str(self.config.get('execution', {}).get('num_threads', 16))
        }

        return self._run_bash_script('02_quality_filtering.sh', env)

    def run_clustering(self) -> int:
        """Run Step 3: vclust clustering."""
        self.logger.info("\n" + "="*70)
        self.logger.info("STEP 3: Clustering with vclust")
        self.logger.info("="*70)

        env = {
            'THREADS': str(self.config.get('execution', {}).get('num_threads', 16))
        }

        return self._run_bash_script('03_vclust_clustering.sh', env)

    def run_ani_validation(self) -> int:
        """Run Step 4: All-vs-all ANI validation."""
        self.logger.info("\n" + "="*70)
        self.logger.info("STEP 4: All-vs-all ANI validation")
        self.logger.info("="*70)

        env = {
            'THREADS': str(self.config.get('execution', {}).get('num_threads', 16))
        }

        return self._run_bash_script('04_all_vs_all_ani.sh', env)

    def run_refinement(self) -> int:
        """Run Step 5: Cluster refinement."""
        self.logger.info("\n" + "="*70)
        self.logger.info("STEP 5: Cluster refinement")
        self.logger.info("="*70)

        env = {
            'THREADS': str(self.config.get('execution', {}).get('num_threads', 16))
        }

        return self._run_bash_script('05_refine_representatives.sh', env)

    def run_all(self) -> int:
        """Run complete pipeline."""
        self.logger.info("\n" + "="*80)
        self.logger.info("VIRAL & PLASMID SEQUENCE RECOVERY AND CLUSTERING WORKFLOW")
        self.logger.info("="*80)

        steps = [
            ('genomad', self.run_genomad),
            ('quality_filtering', self.run_quality_filtering),
            ('clustering', self.run_clustering),
            ('ani_validation', self.run_ani_validation),
            ('refinement', self.run_refinement)
        ]

        results = {}
        for step_name, step_func in steps:
            try:
                return_code = step_func()
                results[step_name] = return_code

                if return_code != 0:
                    self.logger.error(f"\n✗ Workflow failed at step: {step_name}")
                    return 1

            except KeyboardInterrupt:
                self.logger.error("\nWorkflow interrupted by user")
                return 1
            except Exception as e:
                self.logger.error(f"\nError in step {step_name}: {e}", exc_info=True)
                return 1

        # Success
        self.logger.info("\n" + "="*80)
        self.logger.info("✓ WORKFLOW COMPLETED SUCCESSFULLY")
        self.logger.info("="*80)
        self.logger.info(f"\nFinal output in: {self.output_dir}/final/")

        return 0


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Viral and plasmid workflow orchestrator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run complete pipeline
  python workflow_orchestrate.py

  # Run specific steps
  python workflow_orchestrate.py --step 1           # Just geNomad
  python workflow_orchestrate.py --step 3 4 5       # Clustering onwards

  # Custom configuration
  python workflow_orchestrate.py --config custom_config.yaml

  # Dry run (print commands without executing)
  python workflow_orchestrate.py --dry-run
        """
    )

    parser.add_argument(
        '--config', help='Configuration YAML file'
    )
    parser.add_argument(
        '--step', type=int, nargs='+', choices=[1, 2, 3, 4, 5],
        help='Run specific steps (default: all)'
    )
    parser.add_argument(
        '--output-dir', help='Override output directory'
    )
    parser.add_argument(
        '-v', '--verbose', action='store_true', help='Verbose output'
    )
    parser.add_argument(
        '--dry-run', action='store_true', help='Print commands without executing'
    )

    args = parser.parse_args()

    # Create orchestrator
    orchestrator = WorkflowOrchestrator(
        config_file=args.config,
        output_dir=args.output_dir,
        verbose=args.verbose
    )

    # Run workflow
    if args.dry_run:
        orchestrator.logger.info("DRY RUN: Would execute workflow steps")
        return 0

    if args.step:
        # Run specific steps
        step_functions = {
            1: orchestrator.run_genomad,
            2: orchestrator.run_quality_filtering,
            3: orchestrator.run_clustering,
            4: orchestrator.run_ani_validation,
            5: orchestrator.run_refinement
        }

        for step in args.step:
            return_code = step_functions[step]()
            if return_code != 0:
                return 1
        return 0
    else:
        # Run all steps
        return orchestrator.run_all()


if __name__ == "__main__":
    sys.exit(main())
