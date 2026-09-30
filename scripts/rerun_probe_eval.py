"""Re-run the notebook-05 probe benchmark with special-token-masked pooling.

Reproduces the accuracy-vs-coverage table from saved state, adds the two
reference points the original table lacked (full-length test baseline and the
majority-class rate), and reports both pooling variants side by side.
Writes results/accuracy_vs_coverage_v2.csv.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from plm_benchmark.embed import DEFAULT_MODEL, embed_sequences
from plm_benchmark.probe import evaluate_probe, train_probe

P = ROOT / "data" / "processed"
prot = pd.read_csv(P / "labeled_proteins_split.csv")
pep = pd.read_csv(P / "labeled_peptide_digest.csv")
X_prot_old = np.load(P / "labeled_full_embeddings.npy")
X_pep_old = np.load(P / "labeled_peptide_embeddings.npy")

for name, seqs, cache in [("prot", prot["sequence"].tolist(), "labeled_full_embeddings_fixed.npy"),
                          ("pep", pep["peptide"].tolist(), "labeled_peptide_embeddings_fixed.npy")]:
    if not (P / cache).exists():
        print(f"embedding {name}: {len(seqs)} sequences", flush=True)
        np.save(P / cache, embed_sequences(seqs, model_name=DEFAULT_MODEL, pooling="mean"))
X_prot_new = np.load(P / "labeled_full_embeddings_fixed.npy")
X_pep_new = np.load(P / "labeled_peptide_embeddings_fixed.npy")

buckets = pd.cut(pep["coverage"], bins=[0, 0.1, 0.25, 0.5, 0.75, 1.0], include_lowest=True)
pep = pep.assign(coverage_bucket=buckets)
is_train_prot = (prot["split"] == "train").to_numpy()
pep_test = (pep["split"] == "test").to_numpy()

rows = []
for tag, Xp, Xq in [("special_tokens_in_pool", X_prot_old, X_pep_old), ("masked_pool", X_prot_new, X_pep_new)]:
    scaler, clf = train_probe(Xp[is_train_prot], prot.loc[is_train_prot, "ec_class"])
    full = evaluate_probe(scaler, clf, Xp[~is_train_prot], prot.loc[~is_train_prot, "ec_class"])
    rows.append({"pooling": tag, "stratum": "full-length proteins (test)", "accuracy": full["accuracy"],
                 "f1_macro": full["f1_macro"], "n": full["n"], "median_length": int(prot.loc[~is_train_prot, "sequence"].str.len().median())})
    sub = pep[pep_test]
    for bucket, grp in sub.groupby("coverage_bucket", observed=True):
        m = evaluate_probe(scaler, clf, Xq[pep_test][sub.index.get_indexer(grp.index)], grp["ec_class"])
        rows.append({"pooling": tag, "stratum": f"peptides, coverage {bucket}", "accuracy": m["accuracy"],
                     "f1_macro": m["f1_macro"], "n": m["n"], "median_length": int(grp["peptide"].str.len().median())})
    print(tag, "done", flush=True)

majority = prot.loc[~is_train_prot, "ec_class"].value_counts(normalize=True).iloc[0]
pep_major = pep[pep_test]["ec_class"].value_counts(normalize=True).iloc[0]
rows.append({"pooling": "reference", "stratum": "majority-class rate (test proteins)", "accuracy": majority,
             "f1_macro": np.nan, "n": int((~is_train_prot).sum()), "median_length": np.nan})
rows.append({"pooling": "reference", "stratum": "majority-class rate (test peptides)", "accuracy": pep_major,
             "f1_macro": np.nan, "n": int(pep_test.sum()), "median_length": np.nan})

out = pd.DataFrame(rows)
out.to_csv(ROOT / "results" / "accuracy_vs_coverage_v2.csv", index=False)
print(out.to_string(index=False))
