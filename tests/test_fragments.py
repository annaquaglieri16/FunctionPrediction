import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pytest

from plm_benchmark.fragments import random_window, sliding_windows, truncate


def test_sliding_windows_basic():
    seq = "ABCDEFGHIJ"  # length 10
    windows = sliding_windows(seq, size=4, stride=1)
    assert len(windows) == 7  # 10 - 4 + 1
    assert windows[0].sequence == "ABCD"
    assert windows[-1].sequence == "GHIJ"
    for w in windows:
        assert seq[w.start:w.end] == w.sequence


def test_sliding_windows_stride():
    seq = "ABCDEFGHIJ"
    windows = sliding_windows(seq, size=4, stride=3)
    assert [w.sequence for w in windows] == ["ABCD", "DEFG", "GHIJ"]


def test_sliding_windows_size_larger_than_sequence():
    assert sliding_windows("ABC", size=10) == []


def test_truncate_n_terminus():
    seq = "ABCDEFGHIJ"
    frag = truncate(seq, fraction=0.5, terminus="N")
    assert frag.sequence == "ABCDE"
    assert (frag.start, frag.end) == (0, 5)


def test_truncate_c_terminus():
    seq = "ABCDEFGHIJ"
    frag = truncate(seq, fraction=0.5, terminus="C")
    assert frag.sequence == "FGHIJ"
    assert (frag.start, frag.end) == (5, 10)


def test_truncate_center():
    seq = "ABCDEFGHIJ"
    frag = truncate(seq, fraction=0.4, terminus="center")
    assert len(frag.sequence) == 4
    assert seq[frag.start:frag.end] == frag.sequence


def test_truncate_invalid_fraction():
    with pytest.raises(ValueError):
        truncate("ABCDEFG", fraction=0)
    with pytest.raises(ValueError):
        truncate("ABCDEFG", fraction=1.5)


def test_random_window_reproducible():
    seq = "ABCDEFGHIJ" * 3
    rng1 = random.Random(42)
    rng2 = random.Random(42)
    f1 = random_window(seq, size=5, rng=rng1)
    f2 = random_window(seq, size=5, rng=rng2)
    assert f1 == f2
    assert len(f1.sequence) == 5
