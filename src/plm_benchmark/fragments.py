"""Non-enzymatic ways to generate incomplete protein fragments.

These complement digest.py's enzymatic digestion with a controlled
ablation: fixed-length windows / terminal truncations, independent of any
protease's cleavage bias. Useful for a clean "how much sequence length is
needed" dose-response curve.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Fragment:
    sequence: str
    start: int  # 0-based, inclusive
    end: int  # 0-based, exclusive

    def __len__(self) -> int:
        return len(self.sequence)


def sliding_windows(sequence: str, size: int, stride: int = 1) -> list[Fragment]:
    """All fixed-length windows of `size` residues, stepped by `stride`."""
    if size <= 0:
        raise ValueError("size must be positive")
    if stride <= 0:
        raise ValueError("stride must be positive")
    sequence = sequence.strip().upper()
    if size > len(sequence):
        return []
    frags = []
    for start in range(0, len(sequence) - size + 1, stride):
        end = start + size
        frags.append(Fragment(sequence[start:end], start, end))
    return frags


def truncate(sequence: str, fraction: float, terminus: str = "N") -> Fragment:
    """Keep `fraction` of the sequence, truncating from one terminus.

    terminus: "N" keeps the N-terminal portion (removes the C-terminal
    tail), "C" keeps the C-terminal portion, "center" keeps a centered
    window of that fraction.
    """
    if not 0 < fraction <= 1:
        raise ValueError("fraction must be in (0, 1]")
    sequence = sequence.strip().upper()
    n = len(sequence)
    keep = max(1, round(n * fraction))

    if terminus == "N":
        return Fragment(sequence[:keep], 0, keep)
    if terminus == "C":
        return Fragment(sequence[n - keep:], n - keep, n)
    if terminus == "center":
        start = (n - keep) // 2
        end = start + keep
        return Fragment(sequence[start:end], start, end)
    raise ValueError("terminus must be 'N', 'C', or 'center'")


def random_window(sequence: str, size: int, rng) -> Fragment:
    """A single random internal window of `size` residues.

    `rng` is a `random.Random` instance (pass one explicitly for
    reproducibility rather than relying on the global random module).
    """
    sequence = sequence.strip().upper()
    if size > len(sequence):
        raise ValueError("size larger than sequence length")
    start = rng.randint(0, len(sequence) - size)
    end = start + size
    return Fragment(sequence[start:end], start, end)
