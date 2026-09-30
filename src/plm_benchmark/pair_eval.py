"""Pairwise same-function vs. different-function peptide similarity eval.

Tests a narrower claim than the probe benchmark (`evaluate.py`/`probe.py`):
not "can a classifier recover function from an embedding" but "do two
peptides from proteins with the same function land closer together in
embedding space than two peptides from proteins with different function" --
i.e. is function signal present in the raw geometry, no classifier fitted.

Pairs are built across *different* parent proteins only (same-protein pairs
are excluded -- two peptides from one protein are often overlapping
substrings, which would trivially look "similar" regardless of any real
functional signal). Scored by cosine similarity; labelled 1 if the two
parent proteins share the same function annotation (e.g. full EC number),
0 otherwise. Same-function pairs are typically a small minority, so AUROC
alone is optimistic -- report average precision / precision@k alongside it.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def sample_peptides_per_protein(
    peptide_df: pd.DataFrame,
    max_per_protein: int,
    protein_col: str = "protein_id",
    rng: np.random.Generator | None = None,
) -> pd.DataFrame:
    """Cap the number of peptides kept per protein.

    Pairwise comparisons grow quadratically in the number of peptides, and a
    protein digested with several missed cleavages contributes many
    overlapping peptides; capping keeps compute tractable and stops any one
    protein from dominating the pair pool.
    """
    if rng is None:
        rng = np.random.default_rng(0)

    kept_indices = []
    for _, group in peptide_df.groupby(protein_col, sort=False):
        if len(group) <= max_per_protein:
            kept_indices.append(group.index.to_numpy())
        else:
            kept_indices.append(
                rng.choice(group.index.to_numpy(), size=max_per_protein, replace=False)
            )
    return peptide_df.loc[np.concatenate(kept_indices)].reset_index(drop=True)


def pairwise_cosine_similarity(embeddings: np.ndarray) -> np.ndarray:
    """Full pairwise cosine similarity matrix, shape (n, n).

    Materialises n^2 floats (466 MB at n=10.8k, 4.6 GB at n=34k), so it caps
    the peptide budget by memory. Prefer `cosine_similarity_pairs`, which
    scores only the pairs you asked for, in blocks.
    """
    norms = np.linalg.norm(embeddings, axis=1, keepdims=True)
    normalized = embeddings / np.clip(norms, 1e-12, None)
    return normalized @ normalized.T


def center_embeddings(embeddings: np.ndarray) -> np.ndarray:
    """Subtract the dataset mean embedding.

    Contextual-embedding spaces are anisotropic: every sequence sits in a
    narrow cone, so raw cosines are uniformly high (pilot 06: 1st percentile
    0.74, median 0.95) and the differences that carry signal are compressed.
    Centring is the cheapest correction and it moved the pilot's
    same-function AUROC from 0.468 (below chance) to 0.501. Removing the top
    one or two principal components ("all-but-the-top") is the next lever.
    """
    return embeddings - embeddings.mean(axis=0, keepdims=True)


def cosine_similarity_pairs(
    embeddings: np.ndarray,
    i_idx: np.ndarray,
    j_idx: np.ndarray,
    chunk: int = 250_000,
) -> np.ndarray:
    """Cosine similarity for the given index pairs only, computed in blocks.

    Memory is O(chunk * dim) instead of O(n^2), so the number of peptides is
    no longer bounded by the similarity matrix.
    """
    norms = np.linalg.norm(embeddings, axis=1, keepdims=True)
    normalized = (embeddings / np.clip(norms, 1e-12, None)).astype(np.float32)
    out = np.empty(len(i_idx), dtype=np.float32)
    for start in range(0, len(i_idx), chunk):
        stop = min(start + chunk, len(i_idx))
        out[start:stop] = np.einsum(
            "ij,ij->i", normalized[i_idx[start:stop]], normalized[j_idx[start:stop]]
        )
    return out


AMINO_ACIDS = tuple("ACDEFGHIKLMNPQRSTVWY")


def composition_features(peptides) -> np.ndarray:
    """Amino-acid composition (fractions, 20-dim) per peptide.

    The baseline floor for any embedding claim: if a PLM embedding does not
    beat composition, it is adding nothing over counting letters. Report
    Delta(embedding - composition), not the embedding metric alone.
    """
    return np.array(
        [[peptide.count(aa) / max(len(peptide), 1) for aa in AMINO_ACIDS] for peptide in peptides],
        dtype=np.float32,
    )


def length_match_score(peptides, i_idx: np.ndarray, j_idx: np.ndarray) -> np.ndarray:
    """-|length difference| as a pair score.

    A second baseline, and a diagnostic: mean-pooled embeddings that include
    special tokens encode length, so a length score that ranks pairs similarly
    to the embedding score is evidence of that artefact rather than of signal.
    """
    lengths = np.array([len(peptide) for peptide in peptides], dtype=np.float32)
    return -np.abs(lengths[i_idx] - lengths[j_idx])


def build_cross_protein_pairs(
    protein_ids: np.ndarray, function_labels: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """All i<j index pairs from *different* proteins, plus the same-function
    label for each pair.

    Returns (i_idx, j_idx, same_function) where same_function[k] is 1 if
    peptide i_idx[k] and peptide j_idx[k] come from proteins sharing the
    same function label, 0 otherwise.

    Enumerates all n(n-1)/2 pairs; use `sample_cross_protein_pairs` once n
    makes that impractical.
    """
    n = len(protein_ids)
    i_idx, j_idx = np.triu_indices(n, k=1)
    cross_protein = protein_ids[i_idx] != protein_ids[j_idx]
    i_idx, j_idx = i_idx[cross_protein], j_idx[cross_protein]
    same_function = (function_labels[i_idx] == function_labels[j_idx]).astype(int)
    return i_idx, j_idx, same_function


def sample_cross_protein_pairs(
    protein_ids: np.ndarray,
    function_labels: np.ndarray,
    n_pairs: int,
    rng: np.random.Generator | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """A uniform random sample of cross-protein pairs.

    Same distribution as `build_cross_protein_pairs`, so metrics are
    unbiased estimates of the full-enumeration values -- at 3M sampled pairs
    the pilot's prevalence estimate (0.329%) matched the full 58.4M-pair
    value (0.324%). Use this to scale past the point where enumerating every
    pair is affordable.
    """
    if rng is None:
        rng = np.random.default_rng(0)
    n = len(protein_ids)
    i_idx = rng.integers(0, n, n_pairs, dtype=np.int64)
    j_idx = rng.integers(0, n, n_pairs, dtype=np.int64)
    keep = protein_ids[i_idx] != protein_ids[j_idx]
    i_idx, j_idx = i_idx[keep], j_idx[keep]
    same_function = (function_labels[i_idx] == function_labels[j_idx]).astype(int)
    return i_idx, j_idx, same_function


def same_function_from_label_sets(
    label_sets, i_idx: np.ndarray, j_idx: np.ndarray
) -> np.ndarray:
    """Same-function label for set-valued annotations (any shared member).

    Multifunctional enzymes carry several EC numbers. Comparing only the
    primary one (as `load_ec_labels` does) mislabels two proteins that share
    a *secondary* activity as different-function. This compares full sets via
    a sparse membership matrix, so it stays vectorised over tens of millions
    of pairs.
    """
    from scipy import sparse

    vocabulary = {label: index for index, label in enumerate(sorted({l for s in label_sets for l in s}))}
    rows, cols = [], []
    for row, labels in enumerate(label_sets):
        for label in labels:
            rows.append(row)
            cols.append(vocabulary[label])
    membership = sparse.csr_matrix(
        (np.ones(len(rows), dtype=np.int8), (rows, cols)),
        shape=(len(label_sets), len(vocabulary)),
    )
    shared = membership[i_idx].multiply(membership[j_idx]).sum(axis=1)
    return (np.asarray(shared).ravel() > 0).astype(int)


def stratify_pairs(
    group_ids: np.ndarray, i_idx: np.ndarray, j_idx: np.ndarray
) -> dict[str, np.ndarray]:
    """Boolean masks splitting pairs by whether the two members share a group.

    Pass sequence-identity cluster ids (e.g. from `mmseqs easy-cluster
    --min-seq-id 0.5 -c 0.8`) to separate the two questions the pooled metric
    confounds:

      "different_group" -- can the embedding tell that two *non-homologous*
          peptides come from same-function proteins? This is the question
          worth answering; homologous positives are why the pooled pilot
          showed 91x enrichment at k=10 with ~all hits cross-species.
      "same_group" -- expected easy; a positive control that the pipeline can
          detect anything at all, and a job alignment already does better.

    Also useful with organism ids, to separate "same genome, shared
    amino-acid usage" from "different genome, shared function".
    """
    same = group_ids[i_idx] == group_ids[j_idx]
    return {"same_group": same, "different_group": ~same}


def bootstrap_auroc_by_group(
    similarity: np.ndarray,
    same_function: np.ndarray,
    group_i: np.ndarray,
    group_j: np.ndarray,
    n_boot: int = 200,
    min_positives: int = 50,
    rng: np.random.Generator | None = None,
) -> dict:
    """Bootstrap CI for the pair AUROC, resampling *groups* not pairs.

    Pairs are strongly dependent: the pilot's 58.4M pairs come from 5,403
    proteins, and with a 2-peptides-per-protein cap one protein pair
    contributes up to 4 peptide pairs. A pair-level interval is therefore far
    too narrow. Resample protein (or cluster) ids with replacement and keep
    pairs whose *both* members are in the resample. On the pilot this gave
    0.501 [0.494, 0.508], i.e. indistinguishable from chance.
    """
    from sklearn.metrics import roc_auc_score

    if rng is None:
        rng = np.random.default_rng(0)
    groups = np.unique(np.concatenate([group_i, group_j]))
    lookup = np.zeros(groups.max() + 1, dtype=bool)
    aucs = []
    for _ in range(n_boot):
        lookup[:] = False
        lookup[rng.choice(groups, size=len(groups), replace=True)] = True
        mask = lookup[group_i] & lookup[group_j]
        if same_function[mask].sum() >= min_positives and same_function[mask].mean() < 1:
            aucs.append(roc_auc_score(same_function[mask], similarity[mask]))
    aucs = np.asarray(aucs)
    if aucs.size == 0:
        return {"auroc_mean": float("nan"), "ci95": (float("nan"), float("nan")), "n_boot": 0}
    return {
        "auroc_mean": float(aucs.mean()),
        "ci95": tuple(float(x) for x in np.quantile(aucs, [0.025, 0.975])),
        "n_boot": int(aucs.size),
    }


def precision_at_k(similarity: np.ndarray, same_function: np.ndarray, k: int) -> float:
    """Precision among the top-k highest-similarity pairs.

    Uses `argpartition` (O(n)) rather than a full sort of every pair score.
    """
    if k > len(similarity):
        raise ValueError(f"k={k} exceeds the number of pairs ({len(similarity)})")
    if k == len(similarity):
        return float(same_function.mean())
    top_k = np.argpartition(-similarity, k)[:k]
    return float(same_function[top_k].mean())


def evaluate_pairs(
    similarity: np.ndarray, same_function: np.ndarray, k_values: tuple[int, ...] = (10, 50, 100, 500)
) -> dict:
    """AUROC + average precision + precision@k for same-function pairs
    ranked by cosine similarity."""
    from sklearn.metrics import average_precision_score, roc_auc_score

    n_pos = int(same_function.sum())
    n_total = len(same_function)

    prevalence = n_pos / n_total
    at_k = {k: precision_at_k(similarity, same_function, k) for k in k_values if k <= n_total}

    return {
        "n_pairs": n_total,
        "n_same_function": n_pos,
        "prevalence": prevalence,
        "auroc": roc_auc_score(same_function, similarity),
        "average_precision": average_precision_score(same_function, similarity),
        "precision_at_k": at_k,
        # precision@k is hard to read against a 0.3% base rate; enrichment is
        # the same number divided by prevalence. Pilot 06 read 1.0x AUROC-wise
        # but 91x at k=10 -- the top of the ranking carried all the (homology)
        # signal, which precision@k alone made easy to miss.
        "enrichment_at_k": {k: (p / prevalence if prevalence else float("nan")) for k, p in at_k.items()},
    }
