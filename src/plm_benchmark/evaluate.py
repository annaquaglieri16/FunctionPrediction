"""Orchestrates the end-to-end benchmark:

  load proteome + labels -> full-length embeddings -> train probe
  -> for each fragmentation strategy, embed fragments -> evaluate probe
  -> collect a table of performance vs. fragment length/coverage.

Label: EC top-level class (`plm_benchmark.labels`), chosen for Phase 1 --
covers ~40% of the E. coli proteome (enzymes only) but needs no free-text
parsing and is a reasonably balanced 7-way classification problem.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from . import digest, fragments
from .labels import accession_from_header


def build_labeled_dataset(
    proteins: dict[str, str], ec_labels: dict[str, str]
) -> pd.DataFrame:
    """Join protein sequences with EC top-level class labels.

    `proteins` is keyed by FASTA header (as `data_io.load_fasta` returns),
    `ec_labels` by bare UniProt accession (as `labels.load_ec_labels`
    returns). Proteins with no EC annotation (non-enzymes) are dropped.
    """
    rows = []
    for protein_id, seq in proteins.items():
        accession = accession_from_header(protein_id)
        ec_class = ec_labels.get(accession)
        if ec_class is None:
            continue
        rows.append(
            {
                "protein_id": protein_id,
                "accession": accession,
                "sequence": seq,
                "ec_class": ec_class,
            }
        )
    return pd.DataFrame(rows)


def build_fragment_dataset(
    proteins: dict[str, str],
    enzyme: str = "trypsin",
    missed_cleavages: int = 2,
    min_length: int = 6,
    max_length: int = 40,
) -> pd.DataFrame:
    """Digest every protein and return a long table: one row per peptide,
    with a back-reference to its parent protein and coverage fraction."""
    rows = []
    for protein_id, seq in proteins.items():
        peptides = digest.digest(
            seq,
            enzyme=enzyme,
            missed_cleavages=missed_cleavages,
            min_length=min_length,
            max_length=max_length,
        )
        for pep in peptides:
            rows.append(
                {
                    "protein_id": protein_id,
                    "peptide": pep.sequence,
                    "start": pep.start,
                    "end": pep.end,
                    "missed_cleavages": pep.missed_cleavages,
                    "coverage": len(pep) / len(seq),
                    "enzyme": enzyme,
                }
            )
    return pd.DataFrame(rows)


def build_truncation_dataset(
    proteins: dict[str, str],
    fractions=(0.1, 0.25, 0.5, 0.75, 1.0),
    termini=("N", "C", "center"),
) -> pd.DataFrame:
    rows = []
    for protein_id, seq in proteins.items():
        for terminus in termini:
            for fraction in fractions:
                frag = fragments.truncate(seq, fraction=fraction, terminus=terminus)
                rows.append(
                    {
                        "protein_id": protein_id,
                        "fragment": frag.sequence,
                        "start": frag.start,
                        "end": frag.end,
                        "coverage": fraction,
                        "terminus": terminus,
                    }
                )
    return pd.DataFrame(rows)


def run_benchmark(config: dict, out_dir: str | Path) -> None:
    """Entry point wired up by scripts/run_benchmark.py.

    TODO:
      1. load_fasta(proteome_path) -> proteins
      2. labels.load_ec_labels(ec_labels_path) -> ec_labels; build_labeled_dataset(proteins, ec_labels)
      3. embed_sequences(full-length proteins) -> X_full; train_probe(X_full, y)
      4. build_fragment_dataset / build_truncation_dataset
      5. embed_sequences(fragments) -> X_frag; evaluate_probe(..., X_frag, y)
      6. aggregate accuracy/F1 vs. coverage bucket, save to out_dir
    """
    raise NotImplementedError(
        "Not wired up yet -- see notebooks/04_labels.ipynb for the pieces "
        "prototyped so far, and the Phase 1 checklist in README.md."
    )
