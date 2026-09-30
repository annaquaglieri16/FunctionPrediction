"""Re-run the notebook-06 pair evaluation with the post-review fixes.

What this adds over the original notebook run:

  1. special-token-masked mean pooling (`embed.py`), re-embedding the same
     peptide table;
  2. amino-acid composition and peptide-length baselines, so the embedding
     metric is interpretable as a difference rather than an absolute;
  3. mean-centred embeddings alongside raw, since anisotropy alone pushed the
     original AUROC below chance;
  4. uniform random pair sampling + blocked cosine, so the peptide budget is
     not capped by an n-by-n matrix;
  5. protein-level bootstrap CIs, because pairs are not independent.

Sequence-identity clustering is deliberately NOT included yet -- when cluster
ids are available, pass them to `pair_eval.stratify_pairs` and report the
different-cluster stratum as the headline.

Usage (from the repo root, with the repo venv):
    .venv/bin/python scripts/rerun_pilot_eval.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from plm_benchmark.embed import DEFAULT_MODEL, embed_sequences  # noqa: E402
from plm_benchmark.pair_eval import (  # noqa: E402
    bootstrap_auroc_by_group,
    center_embeddings,
    composition_features,
    cosine_similarity_pairs,
    evaluate_pairs,
    length_match_score,
    sample_cross_protein_pairs,
    stratify_pairs,
)

PROCESSED = REPO_ROOT / "data" / "processed"
RESULTS = REPO_ROOT / "results"
N_PAIRS = 3_000_000
K_VALUES = (10, 100, 1_000, 10_000)


def main() -> None:
    digest = pd.read_csv(PROCESSED / "pilot_peptide_digest.csv")
    peptides = digest["peptide"].tolist()
    protein_ids = digest["protein_id"].to_numpy()
    ec_full = digest["ec_full"].to_numpy()
    organism_codes = pd.factorize(digest["organism"])[0]
    protein_codes = pd.factorize(digest["protein_id"])[0]

    fixed_path = PROCESSED / "pilot_peptide_embeddings_fixed.npy"
    if fixed_path.exists():
        fixed = np.load(fixed_path)
    else:
        fixed = embed_sequences(peptides, model_name=DEFAULT_MODEL, pooling="mean")
        np.save(fixed_path, fixed)
    old = np.load(PROCESSED / "pilot_peptide_embeddings.npy")

    rng = np.random.default_rng(0)
    i_idx, j_idx, same_function = sample_cross_protein_pairs(
        protein_ids, ec_full, n_pairs=N_PAIRS, rng=rng
    )

    representations = {
        "esm2_special_tokens_in_pool_raw": lambda: cosine_similarity_pairs(old, i_idx, j_idx),
        "esm2_special_tokens_in_pool_centred": lambda: cosine_similarity_pairs(
            center_embeddings(old), i_idx, j_idx
        ),
        "esm2_masked_pool_raw": lambda: cosine_similarity_pairs(fixed, i_idx, j_idx),
        "esm2_masked_pool_centred": lambda: cosine_similarity_pairs(
            center_embeddings(fixed), i_idx, j_idx
        ),
        "aa_composition": lambda: cosine_similarity_pairs(
            composition_features(peptides), i_idx, j_idx
        ),
        "peptide_length_match": lambda: length_match_score(peptides, i_idx, j_idx),
    }

    organism_masks = stratify_pairs(organism_codes, i_idx, j_idx)
    rows = []
    for name, score_fn in representations.items():
        score = score_fn()
        metrics = evaluate_pairs(score, same_function, k_values=K_VALUES)
        boot = bootstrap_auroc_by_group(
            score, same_function, protein_codes[i_idx], protein_codes[j_idx],
            n_boot=40, rng=np.random.default_rng(7),
        )
        row = {
            "representation": name,
            "auroc": metrics["auroc"],
            "auroc_ci_low": boot["ci95"][0],
            "auroc_ci_high": boot["ci95"][1],
            "average_precision": metrics["average_precision"],
            "prevalence": metrics["prevalence"],
            "length_corr": float(
                np.corrcoef(score, -length_match_score(peptides, i_idx, j_idx))[0, 1]
            ),
        }
        row.update({f"enrichment_at_{k}": v for k, v in metrics["enrichment_at_k"].items()})
        for stratum, mask in organism_masks.items():
            label = {"same_group": "within_organism", "different_group": "cross_organism"}[stratum]
            if same_function[mask].sum() > 50:
                row[f"auroc_{label}"] = evaluate_pairs(
                    score[mask], same_function[mask], k_values=()
                )["auroc"]
        rows.append(row)
        print(f"{name:36s} AUROC {row['auroc']:.4f} "
              f"[{row['auroc_ci_low']:.4f}, {row['auroc_ci_high']:.4f}]  "
              f"AP {row['average_precision']:.5f}  "
              f"enrich@10 {row.get('enrichment_at_10', float('nan')):.1f}x")

    table = pd.DataFrame(rows)
    RESULTS.mkdir(exist_ok=True)
    table.to_csv(RESULTS / "pilot06_rerun_metrics.csv", index=False)
    print(f"\npairs {len(same_function):,}  prevalence {same_function.mean():.4%}  "
          f"proteins {len(np.unique(protein_codes)):,}")
    print("wrote", RESULTS / "pilot06_rerun_metrics.csv")


if __name__ == "__main__":
    main()
