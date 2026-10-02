#!/usr/bin/env python3
"""CLI entrypoint for digesting a proteome."""

import argparse
from pathlib import Path
import yaml

from plm_benchmark.data_io import load_fasta
from plm_benchmark.digest import digest_proteome

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fasta-path",
        required=True,
        help="Path to organism proteome FASTA file.",
    )

    parser.add_argument(
        "--enzyme-config",
        default="configs/enzymes.yaml",
        help="Path to an enzyme settings YAML.",
    )

    parser.add_argument(
        "--out-path",
        default = "peptide.parquet",
        help = "Output parquet file with the digested proteome."
    )

    args = parser.parse_args()

    out_path = Path(args.out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(args.enzyme_config) as fh:
        enzyme_cfg = yaml.safe_load(fh)
        
    enzyme = enzyme_cfg.get("default_enzyme", "trypsin")
    missed_cleavages = enzyme_cfg.get("missed_cleavages", [2])
    min_length = enzyme_cfg.get("min_length", 6)
    max_length = enzyme_cfg.get("max_length", 40)


    fasta_path = Path(args.fasta_path)

    sequences = load_fasta(fasta_path)
    peptide_df = digest_proteome(
        sequences,
        enzyme=enzyme,
        missed_cleavages=max(missed_cleavages),
        min_length=min_length,
        max_length=max_length,
    )

    peptide_df.to_parquet(out_path, index=False)

    print(f"{len(peptide_df):,} peptides from {len(sequences):,} proteins -> {out_path}")


if __name__ == "__main__":
    main()
