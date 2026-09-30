"""Tests for the post-review additions to pair_eval (baselines, blocked
similarity, group stratification, group-level bootstrap, set-valued labels)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pytest

from plm_benchmark.pair_eval import (
    bootstrap_auroc_by_group,
    center_embeddings,
    composition_features,
    cosine_similarity_pairs,
    evaluate_pairs,
    length_match_score,
    pairwise_cosine_similarity,
    precision_at_k,
    same_function_from_label_sets,
    sample_cross_protein_pairs,
    stratify_pairs,
)


def test_cosine_similarity_pairs_matches_full_matrix():
    rng = np.random.default_rng(0)
    embeddings = rng.normal(size=(40, 8)).astype(np.float32)
    i_idx, j_idx = np.triu_indices(40, k=1)
    full = pairwise_cosine_similarity(embeddings)[i_idx, j_idx]
    blocked = cosine_similarity_pairs(embeddings, i_idx, j_idx, chunk=7)
    assert np.allclose(full, blocked, atol=1e-6)


def test_cosine_similarity_pairs_respects_chunk_boundaries():
    embeddings = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 0.0]], dtype=np.float32)
    i_idx = np.array([0, 0, 1])
    j_idx = np.array([1, 2, 2])
    out = cosine_similarity_pairs(embeddings, i_idx, j_idx, chunk=1)
    assert out == pytest.approx([0.0, 1.0, 0.0], abs=1e-6)


def test_center_embeddings_zero_mean():
    embeddings = np.array([[1.0, 5.0], [3.0, 7.0]])
    centered = center_embeddings(embeddings)
    assert centered.mean(axis=0) == pytest.approx([0.0, 0.0])
    # centring changes cosine geometry -- that is the point
    assert pairwise_cosine_similarity(centered)[0, 1] == pytest.approx(-1.0)


def test_composition_features_are_fractions():
    features = composition_features(["AAAC", "GGGG"])
    assert features.shape == (2, 20)
    assert features.sum(axis=1) == pytest.approx([1.0, 1.0])
    alanine = "ACDEFGHIKLMNPQRSTVWY".index("A")
    assert features[0, alanine] == pytest.approx(0.75)


def test_composition_features_ignore_nonstandard_residues():
    # X is not in the 20-letter alphabet, so the fractions no longer sum to 1
    features = composition_features(["AAXX"])
    assert features.sum() == pytest.approx(0.5)


def test_length_match_score_is_negative_absolute_difference():
    peptides = ["AAAA", "AA", "AAAAAA"]
    scores = length_match_score(peptides, np.array([0, 0]), np.array([1, 2]))
    assert scores == pytest.approx([-2.0, -2.0])


def test_sample_cross_protein_pairs_excludes_same_protein():
    protein_ids = np.array(["A", "A", "B", "C"])
    function_labels = np.array(["f1", "f1", "f1", "f2"])
    i_idx, j_idx, same_function = sample_cross_protein_pairs(
        protein_ids, function_labels, n_pairs=500, rng=np.random.default_rng(0)
    )
    assert len(i_idx) > 0
    assert np.all(protein_ids[i_idx] != protein_ids[j_idx])
    assert set(np.unique(same_function)) <= {0, 1}


def test_sample_cross_protein_pairs_prevalence_approximates_enumeration():
    rng = np.random.default_rng(1)
    protein_ids = np.array([f"p{i}" for i in range(60)])
    function_labels = np.array(["f1" if i % 6 == 0 else f"f{i}" for i in range(60)])
    _, _, sampled = sample_cross_protein_pairs(
        protein_ids, function_labels, n_pairs=200_000, rng=rng
    )
    i_idx, j_idx = np.triu_indices(60, k=1)
    exact = (function_labels[i_idx] == function_labels[j_idx]).mean()
    assert sampled.mean() == pytest.approx(exact, abs=0.01)


def test_stratify_pairs_splits_on_shared_group():
    group_ids = np.array([0, 0, 1, 2])
    i_idx = np.array([0, 0, 2])
    j_idx = np.array([1, 2, 3])
    masks = stratify_pairs(group_ids, i_idx, j_idx)
    assert masks["same_group"].tolist() == [True, False, False]
    assert masks["different_group"].tolist() == [False, True, True]
    assert (masks["same_group"] | masks["different_group"]).all()


def test_same_function_from_label_sets_uses_intersection():
    label_sets = [
        frozenset({"1.3.5.1", "1.3.5.4"}),
        frozenset({"1.3.5.4"}),  # shares only the secondary EC
        frozenset({"2.7.1.1"}),
    ]
    i_idx = np.array([0, 0, 1])
    j_idx = np.array([1, 2, 2])
    same = same_function_from_label_sets(label_sets, i_idx, j_idx)
    assert same.tolist() == [1, 0, 0]


def test_same_function_from_label_sets_disagrees_with_primary_only():
    # the case load_ec_labels gets wrong: shared secondary activity
    label_sets = [frozenset({"1.1.1.1", "1.2.3.4"}), frozenset({"9.9.9.9", "1.2.3.4"})]
    same = same_function_from_label_sets(label_sets, np.array([0]), np.array([1]))
    primary_only = int(sorted(label_sets[0])[0] == sorted(label_sets[1])[0])
    assert same.tolist() == [1]
    assert primary_only == 0


def test_evaluate_pairs_reports_enrichment_over_prevalence():
    similarity = np.array([0.9, 0.8, 0.2, 0.1])
    same_function = np.array([1, 0, 0, 0])
    metrics = evaluate_pairs(similarity, same_function, k_values=(1, 4))
    assert metrics["prevalence"] == pytest.approx(0.25)
    assert metrics["enrichment_at_k"][1] == pytest.approx(4.0)
    assert metrics["enrichment_at_k"][4] == pytest.approx(1.0)


def test_precision_at_k_handles_k_equal_to_n():
    similarity = np.array([0.9, 0.1, 0.5])
    same_function = np.array([1, 0, 1])
    assert precision_at_k(similarity, same_function, k=3) == pytest.approx(2 / 3)


def test_bootstrap_auroc_by_group_brackets_perfect_separation():
    rng = np.random.default_rng(0)
    n_groups = 30
    group_i = rng.integers(0, n_groups, 4000)
    group_j = rng.integers(0, n_groups, 4000)
    same_function = rng.integers(0, 2, 4000)
    # perfectly separable scores
    similarity = same_function + rng.normal(0, 0.01, 4000)
    result = bootstrap_auroc_by_group(
        similarity, same_function, group_i, group_j, n_boot=25, rng=rng
    )
    assert result["n_boot"] > 0
    assert result["auroc_mean"] > 0.99
    assert result["ci95"][0] <= result["auroc_mean"] <= result["ci95"][1]


def test_bootstrap_auroc_by_group_returns_nan_when_no_usable_resample():
    group_i = np.array([0, 0])
    group_j = np.array([1, 1])
    same_function = np.array([1, 0])
    result = bootstrap_auroc_by_group(
        similarity=np.array([0.9, 0.1]),
        same_function=same_function,
        group_i=group_i,
        group_j=group_j,
        n_boot=5,
        min_positives=50,
    )
    assert result["n_boot"] == 0
    assert np.isnan(result["auroc_mean"])
