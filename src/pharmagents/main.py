from __future__ import annotations

import argparse
import json

from pharmagents.core.config import PipelineConfig
from pharmagents.workflows.pipeline import VirtualPharmaPipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description='PharmAgents production-ready working project')
    sub = parser.add_subparsers(dest='command', required=True)

    run = sub.add_parser('run', help='Run full pipeline')
    run.add_argument('--disease', required=True, type=str)
    run.add_argument('--top-k-targets', type=int, default=3)
    run.add_argument('--n-leads', type=int, default=12)
    run.add_argument('--optimization-rounds', type=int, default=4)
    run.add_argument('--real-fetch', action='store_true', help='Enable live UniProt/RCSB fetching')
    run.add_argument('--docking-backend', choices=['auto', 'heuristic', 'vina', 'gnina'], default='auto')
    run.add_argument('--receptor-path', type=str, default=None, help='Optional prepared receptor file path')
    run.add_argument('--pdb-path', type=str, default=None, help='Optional raw PDB receptor file path for auto-preparation')
    run.add_argument('--center', nargs=3, type=float, default=None, metavar=('X', 'Y', 'Z'))
    run.add_argument('--box-size', nargs=3, type=float, default=None, metavar=('SX', 'SY', 'SZ'))
    run.add_argument('--json', action='store_true')

    history = sub.add_parser('history', help='Show recent experiment runs')
    history.add_argument('--limit', type=int, default=10)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.command == 'run':
        config = PipelineConfig(use_real_fetch=args.real_fetch, docking_backend=args.docking_backend)
        pipe = VirtualPharmaPipeline(config=config)
        overrides = {
            'receptor_path': args.receptor_path,
            'pdb_path': args.pdb_path,
            'pocket_center': list(args.center) if args.center else None,
            'box_size': list(args.box_size) if args.box_size else None,
        }
        result = pipe.run(
            disease=args.disease,
            top_k_targets=args.top_k_targets,
            n_leads=args.n_leads,
            optimization_rounds=args.optimization_rounds,
            target_overrides=overrides,
        )
        if args.json:
            print(json.dumps(result, indent=2))
        else:
            print(result['final_report'])
            print(f"\nRun ID: {result.get('run_id')}")
            print(f"Artifacts: {result.get('artifacts', {})}")
    elif args.command == 'history':
        pipe = VirtualPharmaPipeline()
        print(json.dumps(pipe.experiment_logger.list_runs(limit=args.limit), indent=2))


if __name__ == '__main__':
    main()
