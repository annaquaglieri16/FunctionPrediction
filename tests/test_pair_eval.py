import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd
import pytest

from plm_benchmark.pair_eval import (
    build_cross_protein_pairs,
    evaluate_pairs,
    pairwise_cosine_similarity,
    precision_at_k,
    sample_peptides_per_protein,
)


def test_sample_peptides_per_protein_caps_each_group():
    df = pd.DataFrame(
        {
            "protein_id": ["A", "A", "A", "A", "B", "B", "C"],
            "peptide": ["p1", "p2", "p3", "p4", "p5", "p6", "p7"],
        }
    )
    out = sample_peptides_per_protein(df, max_per_protein=2, rng=np.random.default_rng(0))

    counts = out["protein_id"].value_counts()
    assert counts["A"] == 2  # capped from 4
    assert counts["B"] == 2  # already exactly 2, kept
    assert counts["C"] == 1  # already fewer than cap, kept
    assert set(out["protein_id"]) == {"A", "B", "C"}


def test_pairwise_cosine_similarity_identical_vectors():
    embeddings = np.array([[1.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    sim = pairwise_cosine_similarity(embeddings)

    assert sim.shape == (3, 3)
    assert sim[0, 1] == pytest.approx(1.0)  # identical direction
    assert sim[0, 2] == pytest.approx(0.0)  # orthogonal
    assert np.allclose(np.diag(sim), 1.0)


def test_pairwise_cosine_similarity_scale_invariant():
    embeddings = np.array([[1.0, 0.0], [3.0, 0.0]])
    sim = pairwise_cosine_similarity(embeddings)
    assert sim[0, 1] == pytest.approx(1.0)


def test_build_cross_protein_pairs_excludes_same_protein():
    protein_ids = np.array(["A", "A", "B", "B"])
    function_labels = np.array(["f1", "f1", "f1", "f2"])

    i_idx, j_idx, same_function = build_cross_protein_pairs(protein_ids, function_labels)

    # Only cross-protein pairs: (0,2) (0,3) (1,2) (1,3) -- (0,1) and (2,3)
    # are same-protein and excluded.
    assert len(i_idx) == 4
    for i, j in zip(i_idx, j_idx):
        assert protein_ids[i] != protein_ids[j]


def test_build_cross_protein_pairs_labels_match_function():
    protein_ids = np.array(["A", "B", "C"])
    function_labels = np.array(["f1", "f1", "f2"])

    i_idx, j_idx, same_function = build_cross_protein_pairs(protein_ids, function_labels)
    pairs = dict(zip(zip(i_idx, j_idx), same_function))

    assert pairs[(0, 1)] == 1  # A,B share f1
    assert pairs[(0, 2)] == 0  # A,C differ
    assert pairs[(1, 2)] == 0  # B,C differ


def test_precision_at_k_perfect_ranking():
    # Highest-similarity pairs are all same-function.
    similarity = np.array([0.9, 0.8, 0.7, 0.2, 0.1])
    same_function = np.array([1, 1, 1, 0, 0])
    assert precision_at_k(similarity, same_function, k=3) == pytest.approx(1.0)
    assert precision_at_k(similarity, same_function, k=5) == pytest.approx(0.6)


def test_precision_at_k_raises_when_k_exceeds_pairs():
    with pytest.raises(ValueError):
        precision_at_k(np.array([0.1, 0.2]), np.array([0, 1]), k=5)


def test_evaluate_pairs_perfect_separation_gives_auroc_one():
    # same-function pairs all score strictly higher than different-function pairs.
    similarity = np.array([0.9, 0.85, 0.8, 0.3, 0.2, 0.1])
    same_function = np.array([1, 1, 1, 0, 0, 0])

    metrics = evaluate_pairs(similarity, same_function, k_values=(2, 6))

    assert metrics["auroc"] == pytest.approx(1.0)
    assert metrics["average_precision"] == pytest.approx(1.0)
    assert metrics["n_pairs"] == 6
    assert metrics["n_same_function"] == 3
    assert metrics["prevalence"] == pytest.approx(0.5)
    assert metrics["precision_at_k"][2] == pytest.approx(1.0)


def test_evaluate_pairs_drops_k_values_larger_than_pairs():
    similarity = np.array([0.9, 0.1])
    same_function = np.array([1, 0])
    metrics = evaluate_pairs(similarity, same_function, k_values=(1, 100))
    assert set(metrics["precision_at_k"]) == {1}
