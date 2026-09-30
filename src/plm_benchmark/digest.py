"""In-silico enzymatic digestion of protein sequences.

Cleavage rules come from pyteomics' ExPASy table (`expasy_rules`) rather than
being maintained here, so the biochemistry is curated upstream. Each pattern
matches the *residue* whose C-terminal peptide bond breaks -- "trypsin" is
"([KR](?=[^P]))|..." meaning cleave after a K or R unless followed by P. The
cleavage point is therefore the match *end*, not its start.

Enzyme names are the ExPASy spellings (hyphens and spaces, e.g. "arg-c",
"glutamyl endopeptidase", "chymotrypsin high specificity"); print
`sorted(CLEAVAGE_RULES)` for the full list.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from pyteomics.parser import expasy_rules

# Cleavage-site rules for common proteases used in bottom-up proteomics,
# sourced entirely from pyteomics (ExPASy). Kept as a module-level alias so
# callers have a stable name to import.
CLEAVAGE_RULES: dict[str, str] = expasy_rules


@dataclass(frozen=True)
class Peptide:
    sequence: str
    start: int  # 0-based, inclusive, position in the parent protein
    end: int  # 0-based, exclusive
    missed_cleavages: int

    def __len__(self) -> int:
        return len(self.sequence)


def _cleavage_sites(sequence: str, pattern: str) -> list[int]:
    """Sorted cleavage positions (0-based, peptide boundary *at* this index),
    always including the sequence end.

    ExPASy patterns match the residue whose C-terminal bond breaks, so the
    boundary is `m.end()`. This also holds for rules that match the residue
    *before* the cut (e.g. asp-n is r"\\w(?=D)"), and for zero-width
    lookaround patterns, where start and end coincide.
    """
    sites = {m.end() for m in re.finditer(pattern, sequence)}
    sites.add(len(sequence))
    return sorted(sites)


def digest(
    sequence: str,
    enzyme: str = "trypsin",
    missed_cleavages: int = 2,
    min_length: int = 1,
    max_length: int | None = None,
) -> list[Peptide]:
    """Digest a protein sequence in silico with the given enzyme rule.

    Returns every peptide obtainable with 0..missed_cleavages missed
    cleavage sites, deduplicated, filtered by length.
    """
    if enzyme not in CLEAVAGE_RULES:
        raise ValueError(
            f"Unknown enzyme {enzyme!r}. Known enzymes: {sorted(CLEAVAGE_RULES)}"
        )
    if missed_cleavages < 0:
        raise ValueError("missed_cleavages must be >= 0")

    sequence = sequence.strip().upper()
    sites = _cleavage_sites(sequence, CLEAVAGE_RULES[enzyme])
    boundaries = [0] + sites

    peptides: dict[tuple[int, int], Peptide] = {}
    n_sites = len(boundaries) - 1
    for i in range(n_sites):
        for mc in range(0, missed_cleavages + 1):
            j = i + 1 + mc
            if j >= len(boundaries):
                break
            start, end = boundaries[i], boundaries[j]
            if start == end:
                continue
            pep_seq = sequence[start:end]
            if len(pep_seq) < min_length:
                continue
            if max_length is not None and len(pep_seq) > max_length:
                continue
            peptides[(start, end)] = Peptide(pep_seq, start, end, mc)

    return sorted(peptides.values(), key=lambda p: (p.start, p.end))


def digest_multi_enzyme(
    sequence: str,
    enzymes: list[str],
    missed_cleavages: int = 2,
    min_length: int = 1,
    max_length: int | None = None,
) -> list[Peptide]:
    """Union of peptides from digesting with several enzymes independently
    (models a combined/parallel multi-enzyme proteomics workflow)."""
    seen: dict[tuple[int, int], Peptide] = {}
    for enzyme in enzymes:
        for pep in digest(sequence, enzyme, missed_cleavages, min_length, max_length):
            seen[(pep.start, pep.end)] = pep
    return sorted(seen.values(), key=lambda p: (p.start, p.end))
