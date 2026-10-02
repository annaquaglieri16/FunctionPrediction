#!/usr/bin/env python3
"""CLI entrypoint for fetching a model organism FASTA file from Uniprot."""

import argparse
from pathlib import Path
import yaml

from plm_benchmark.data_io import fetch_uniprot_proteome


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--organism-config",
        default="configs/organism_ecoli_k12.yaml",
        help="Path to an organism config YAML.",
    )
    parser.add_argument(
        "--out-path",
        required=True,
        help="Output path where to save the downloaded FASTA file.",
    )

    # This is a argument that does not require a value after it - it becomes a flag
    parser.add_argument("--force", action="store_true", help="Force download.")

    args = parser.parse_args()

    with open(args.organism_config) as fh:
        organism_cfg = yaml.safe_load(fh)
        
    reviewed_only = organism_cfg.get("reviewed_only", True)
    proteome_id = organism_cfg.get("uniprot_proteome_id")
    if not proteome_id:
        raise ValueError(f"'uniprot_proteome_id' missing from {args.organism_config}.")

    force = args.force
    fasta_path = Path(args.out_path)

    if force or not fasta_path.exists():
        fasta_path = fetch_uniprot_proteome(
            proteome_id=proteome_id,
            out_path=fasta_path,
            reviewed_only=reviewed_only,
        )
        print("downloaded ->", fasta_path)
    else:
        print("already present ->", fasta_path)


if __name__ == "__main__":
    main()
