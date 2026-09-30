"""Fetch and parse reference proteomes from UniProt."""

from __future__ import annotations

from pathlib import Path


def fetch_uniprot_proteome(
    proteome_id: str,
    out_path: str | Path,
    reviewed_only: bool = True,
) -> Path:
    """Download a reference proteome FASTA from the UniProt REST API.

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
        f"?query={query}&format=fasta&compressed=false"
    )
    response = requests.get(url, timeout=120)
    response.raise_for_status()
    out_path.write_text(response.text)
    return out_path


def load_fasta(path: str | Path) -> dict[str, str]:
    """Minimal FASTA parser: header (up to the first whitespace) -> sequence."""
    sequences: dict[str, str] = {}
    header: str | None = None
    chunks: list[str] = []

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    sequences[header] = "".join(chunks)
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            sequences[header] = "".join(chunks)

    return sequences
