"""Fetch and parse EC (Enzyme Commission) number labels from UniProt.

EC numbers only annotate enzymes, so this label covers a subset of the
proteome (~40% for E. coli K-12) -- non-enzymatic proteins have no EC
number and are simply dropped by `load_ec_labels`. The label used for the
probe is the top-level EC class (the first of the four dot-separated
digits, e.g. "1.1.1.44" -> "1"), which is a reasonably balanced 7-way
classification problem for this proteome.
"""

from __future__ import annotations

import csv
from pathlib import Path

EC_CLASS_NAMES: dict[str, str] = {
    "1": "oxidoreductase",
    "2": "transferase",
    "3": "hydrolase",
    "4": "lyase",
    "5": "isomerase",
    "6": "ligase",
    "7": "translocase",
}


def fetch_ec_labels(
    proteome_id: str,
    out_path: str | Path,
    reviewed_only: bool = True,
) -> Path:
    """Download accession -> EC number annotations from the UniProt REST API.

    Requires network access -- run this from a normal terminal / notebook,
    not necessarily through a sandboxed bridge that may have no egress.
    """
    import requests

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    query = f"proteome:{proteome_id}"
    if reviewed_only:
        query += "+AND+reviewed:true"
    url = (
        "https://rest.uniprot.org/uniprotkb/stream"
        f"?query={query}&fields=accession,ec&format=tsv"
    )
    response = requests.get(url, timeout=120)
    response.raise_for_status()
    out_path.write_text(response.text)
    return out_path


def load_ec_labels(path: str | Path, level: str = "top") -> dict[str, str]:
    """Parse a UniProt accession/EC-number TSV (as written by
    `fetch_ec_labels`) into {accession: EC label}.

    Proteins with no EC number (non-enzymes) are omitted. Proteins with
    several EC numbers (~10% of annotated ones -- multifunctional enzymes)
    keep only the first-listed one, since UniProt lists the primary
    activity first.

    level: "top" (default) returns just the top-level class digit (e.g.
    "1.1.1.44" -> "1"), a coarse 7-way category -- used by the probe
    benchmark (notebooks 04-05). "full" returns the complete primary EC
    number unchanged (e.g. "1.1.1.44") -- a much more specific "same
    enzymatic function" label, used by the peptide-pair similarity pilot
    (notebook 06), where "same broad category" would be too loose a
    definition of "same function".
    """
    if level not in ("top", "full"):
        raise ValueError("level must be 'top' or 'full'")
    # "full_set" is deliberately not handled here -- see load_ec_label_sets,
    # which returns frozensets and so cannot share this dict[str, str] type.

    labels: dict[str, str] = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            ec_field = (row.get("EC number") or "").strip()
            if not ec_field:
                continue
            accession = row["Entry"]
            primary_ec = ec_field.split(";")[0].strip()
            labels[accession] = primary_ec if level == "full" else primary_ec.split(".")[0]
    return labels


def load_ec_label_sets(path: str | Path) -> dict[str, frozenset[str]]:
    """Parse the same TSV into {accession: frozenset of ALL its EC numbers}.

    `load_ec_labels` keeps only the primary EC number, so two multifunctional
    enzymes sharing a *secondary* activity are labelled different-function
    (174 entries in the E. coli K-12 file carry more than one EC). With
    same-function prevalence at ~0.3% that barely moves a pooled metric, but
    it is a genuine label error and it matters once positives are restricted
    to non-homologous pairs, where every positive counts.

    Pair up with `pair_eval.same_function_from_label_sets`, which defines
    same-function as a non-empty intersection.
    """
    label_sets: dict[str, frozenset[str]] = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            ec_field = (row.get("EC number") or "").strip()
            if not ec_field:
                continue
            numbers = {part.strip() for part in ec_field.split(";") if part.strip()}
            if numbers:
                label_sets[row["Entry"]] = frozenset(numbers)
    return label_sets


def accession_from_header(header: str) -> str:
    """Extract the UniProt accession from a FASTA header of the form
    "sp|P0AD86|LPT_ECOLI" (or "tr|...|..."), i.e. the same header format
    `data_io.load_fasta` uses as its dict keys. Returns the header
    unchanged if it isn't pipe-delimited."""
    parts = header.split("|")
    return parts[1] if len(parts) >= 2 else header
