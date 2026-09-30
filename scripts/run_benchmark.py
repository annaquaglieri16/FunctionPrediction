#!/usr/bin/env python3
"""CLI entrypoint for the PLM-on-incomplete-proteins benchmark."""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import yaml

from plm_benchmark.evaluate import run_benchmark


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--organism-config",
        default="configs/organism_ecoli_k12.yaml",
        help="Path to an organism config YAML.",
    )
    parser.add_argument(
        "--enzyme-config",
        default="configs/enzymes.yaml",
        help="Path to an enzyme/digestion config YAML.",
    )
    parser.add_argument("--out-dir", default="results")
    args = parser.parse_args()

    with open(args.organism_config) as fh:
        organism_cfg = yaml.safe_load(fh)
    with open(args.enzyme_config) as fh:
        enzyme_cfg = yaml.safe_load(fh)

    run_benchmark({**organism_cfg, **enzyme_cfg}, out_dir=args.out_dir)


if __name__ == "__main__":
    main()
