#!/usr/bin/env python3

"""
This is the main entry point for the workflow.
"""

from argparse import ArgumentParser

from cpg_flow import workflow

from cpg_exomiser.stages import CombineExomiserGeneTsvs, CombineExomiserVariantTsvs


def cli_main() -> None:
    """
    CLI entrypoint - starts up the workflow
    """
    parser = ArgumentParser()
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    stages = [CombineExomiserGeneTsvs, CombineExomiserVariantTsvs]

    workflow.run_workflow(stages=stages, dry_run=args.dry_run)


if __name__ == '__main__':
    cli_main()
