# Review of notebooks 05 and 06 — what to change before the rerun

Reviewed 2026-09-24 against `notebooks/05_probe_benchmark.ipynb`,
`notebooks/06_function_similarity_pilot.ipynb`, `src/plm_benchmark/{pair_eval,embed,labels}.py`
and `results/pilot_summary_metrics.csv`.

All numbers below were **recomputed from `data/processed/pilot_peptide_embeddings.npy` and
`pilot_peptide_digest.csv`** on a uniform random sample of 3.0M cross-protein pairs
(prevalence 0.329%, matching the full 58.4M-pair prevalence of 0.324%). Figure:
`results/pilot06_diagnosis.png`.

## 1. The reported AUROC of 0.474 is an artefact, not a result

| representation | AUROC | note |
| --- | --- | --- |
| ESM-2 mean-pooled, raw (as run) | **0.468** | below chance |
| ESM-2 mean-pooled, mean-centred | **0.501** | protein-level bootstrap 95% CI [0.494, 0.508] |
| amino-acid composition (20-dim) | 0.490 | baseline that was missing |
| peptide-length match (−\|Δlength\|) | 0.501 | length alone says nothing about EC |

Subtracting the dataset mean embedding before computing cosine moves the AUROC from 0.468 to
0.501, and a bootstrap **resampled over proteins** (not pairs) puts the CI across 0.5. So the
correct statement is "no detectable global signal", not "anti-correlated". The distinction
matters: a sub-chance number invites a story about the embedding actively misleading, when the
real cause is a geometric offset.

**Candidate mechanism, and what the rerun actually showed.** `embed.py` mean-pooled over
`attention_mask`, which **includes the `<cls>` and `<eos>` tokens**. For a 6-mer the pooled vector
is 2/8 = 25% two constant vectors; for a 40-mer, 2/42 = 5%. That constant share is a function of
peptide length, so the pooled embedding partly encodes length, and cosine similarity becomes partly
a length-similarity detector: measured `corr(cosine, |Δlength|) = −0.330`.

Re-embedding with the special tokens masked out confirms the first half of that chain and refutes
the second:

- the embeddings do change, length-dependently: 4.9% relative change for a 27mer rising to
  **15.3% for a 7mer**;
- length coupling drops but does not go away: `corr(cosine, |Δlength|)` −0.330 → **−0.272**;
- **the AUROC barely moves: 0.470 → 0.467 raw, 0.5025 → 0.5019 centred.**

So the special tokens were *not* the cause of the below-chance AUROC. **Centring is what fixes it**
(0.467 → 0.502, CI crossing 0.5), which means the driver is the large shared mean component of the
embedding space generally, not the two constant token vectors specifically. The pooling fix is still
correct and worth keeping — it removes a length artefact that would confound any length-stratified
analysis, and it does measurably sharpen the top of the ranking (enrichment at k=100 rises from
21.5× to 30.7×) — but it should not be advertised as the explanation for the headline number. The
narrow cosine cone (1st percentile 0.743, median 0.945) is ordinary anisotropy.

## 2. There was no baseline, so the negative result was not falsifiable

Composition (0.490) and length (0.501) are both at chance too. That means the current design
cannot distinguish "peptide embeddings lost the function signal" from "this pair construction has
no signal for anything to find". Every rerun must carry the composition baseline through the
identical pipeline; the quantity of interest is Δ(embedding − composition), not the AUROC itself.

## 3. The only signal present is homology, and it is at the extreme top of the ranking

Same-EC enrichment over prevalence, mean-centred embeddings:

| k | precision@k | enrichment | fraction of the same-EC hits that are cross-species |
| --- | --- | --- | --- |
| 10 | 0.300 | **91×** | 100% |
| 100 | 0.070 | 21× | 86% |
| 1,000 | 0.014 | 4.3× | 79% |
| 10,000 | 0.005 | 1.6× | 61% |

Strong enrichment confined to the top ~100 pairs, almost all of it cross-organism. That is the
signature of **orthologues**: the same enzyme in *E. coli* and *P. aeruginosa* yields near-identical
tryptic peptides, so the pair is "same EC" *and* "same sequence". The notebook's write-up reports
the histogram overlap and concludes anisotropy, but does not mention that precision@10 is 91× — a
reader would take "no signal" away from a result that in fact contains a large, uninteresting one.

(Only 4 peptides in the set have an identical sequence elsewhere, and no sampled pair shared a
peptide sequence, so exact-duplicate peptides are *not* the cause. That suspicion was wrong.)

## 4. "Form pairs only across sequence-identity clusters" — what it means and why

The pair filter currently excludes only pairs from the *same protein*. It does not exclude pairs of
*different but homologous* proteins, which is where the whole §3 signal comes from. Cosine
similarity between peptides of two orthologues is high because the sequences are similar; the fact
that they also share an EC number is a consequence, not something the embedding discovered. With
homologues in the positive class, a high AUROC would mean "PLM embeddings detect sequence
similarity", which is already known and which alignment does better and cheaper.

Concretely:

1. Pool the parent proteins of all four proteomes into one FASTA.
2. `mmseqs easy-cluster proteins.fasta clu tmp --min-seq-id 0.5 -c 0.8 --cov-mode 1`
   (50% identity, 80% coverage) and read `clu_cluster.tsv` as `protein_id → cluster_id`.
3. Attach `cluster_id` to every peptide via its parent protein.
4. Build pairs as now, then **split them into two strata** rather than pooling:
   - **different-cluster pairs** — the question worth answering: can the embedding tell that two
     *non-homologous* peptides come from proteins with the same enzymatic function?
   - **same-cluster pairs** — expected to be easy; report separately as a positive control that the
     pipeline can detect anything at all.
5. Report AUROC/AP/precision@k for each stratum, with the protein- or cluster-level bootstrap CI.

If the different-cluster stratum sits at chance with a tight CI, that is a clean, defensible
negative result: mean-pooled short-peptide embeddings carry no non-homologous function signal. The
current run cannot make that claim because the two strata are mixed.

An identity threshold of 0.5 with 0.8 coverage is a reasonable default; 0.3 is the stricter
"remote homology" line if the different-cluster positives are still numerous enough.

## 5. Should the rerun use two organisms instead of four?

No — this does not address the confound, and it costs positives. All of the current signal is
cross-organism (§3), so dropping organisms removes signal without removing homology: paralogues
*within* one genome are homologous too, and the enzyme families that dominate EC classes are
exactly the ones with many in-genome relatives. Keep four proteomes, add the cluster stratification,
and additionally report within-organism and cross-organism strata separately (the notebook's own
"Next" list already proposes this, and it is worth doing — it separates "same genome, shared
amino-acid usage" from "different genome, shared function").

## 6. Label construction: `;`-separated EC numbers are truncated to the first

`load_ec_labels` keeps only the first-listed EC number for multifunctional enzymes (174 entries in
the *E. coli* file alone). Two proteins sharing a *secondary* activity are then labelled
different-function. With prevalence at 0.3% this barely moves the metrics, but it is a label error
that will matter once positives are restricted to different-cluster pairs. Parse the field into a
**set** and define same-function as non-empty intersection.

## 7. Uncertainty must be estimated over proteins, not pairs

58.4M pairs come from 5,403 proteins, and the peptide cap means one protein pair contributes up to
4 peptide pairs. Pairs are strongly dependent, so a pair-level CI is far too narrow. The
protein-level bootstrap used above ([0.494, 0.508] on 40 resamples) is the honest interval; with
cluster stratification, resample clusters.

## 8. Notebook 05: the split is protein-level but not identity-clustered

The train/test split is `train_test_split(..., stratify=ec_class)` over proteins, with peptides
inheriting their parent's split. That correctly prevents same-protein leakage — good, and better
than most published probes — but homologous proteins still land on both sides, so probe accuracy is
optimistic relative to genuinely novel sequences. Use `GroupShuffleSplit` (or `StratifiedGroupKFold`)
with `groups=cluster_id` from the same MMseqs2 run.

Also add the **matched full-length control**: score the probe on held-out *full-length* test
proteins as well as on their peptides. The drop between those two numbers is the cost of
fragmentation specifically; without it, the peptide result confounds "fragment carries less signal"
with "probe generalises imperfectly to held-out proteins". Notebook 05's framing already draws this
distinction in prose — the control makes it a number.

Both notebooks share the `embed.py` pooling bug in §1, so 05's fragment-length curve should be
regenerated after the fix: the constant-token share is largest for the shortest fragments, which is
exactly where the curve collapses.

## 9. Smaller things

- `precision_at_k` calls `np.argsort(-similarity)` (full sort of 58M floats) — `np.argpartition`
  is O(n) and is what the k-largest needs.
- `pairwise_cosine_similarity` materialises the full n×n matrix: 466 MB at n=10.8k, 4.6 GB at
  n=34k. Compute pair similarities in blocks over the pair index instead, so the peptide budget
  isn't capped by memory.
- The "Next" note suggests CLS pooling as an anisotropy fix. Careful: ESM-2 has no
  sentence-level pretraining objective, so its `<cls>` representation is not trained to summarise
  the sequence and is usually *worse* than mean pooling for retrieval. Fix the special-token
  masking and centring first; whitening (removing the top 1-2 principal components) is the more
  reliable next lever.
- Consider a larger checkpoint only after the above: `esm2_t12_35M` is a weak model, but the
  current result is dominated by artefacts, so scaling up now would not be interpretable.

## Code changes made (2026-09-24)

Implemented in `src/`, with tests; sequence-identity clustering deliberately left out for now,
but the hooks are in place so it drops in without further code changes.

**`embed.py`**
- `pooling="mean"` (the default) now pools over **residue positions only**, excluding padding and
  the `<cls>`/`<eos>` special tokens (`return_special_tokens_mask=True`, then
  `attention_mask * (1 - special_tokens_mask)`).
- `pooling="mean_with_special"` preserves the old behaviour so pre-fix results stay reproducible.
- New `layer` argument to pool an intermediate hidden layer instead of the final one.
- Verified on four peptides: the embedding change versus the old pooling scales as 1/length —
  4.9% relative change for a 27mer, 6.3% for 19aa, 8.2% for 14aa, **15.3% for a 7mer** — which is
  the length-coupling mechanism from §1, measured directly.

**`pair_eval.py`**
- `center_embeddings` — the anisotropy correction.
- `cosine_similarity_pairs` — blocked scoring of just the requested pairs; memory is
  O(chunk × dim) instead of O(n²), so the peptide budget is no longer memory-capped.
  `pairwise_cosine_similarity` is kept (and documented as the memory-heavy path).
- `composition_features` (20-dim AA fractions) and `length_match_score` — the two baselines.
- `sample_cross_protein_pairs` — uniform random pair sampling for scaling past full enumeration.
- `stratify_pairs(group_ids, ...)` — returns `same_group` / `different_group` masks; takes MMseqs2
  cluster ids when they exist, and organism ids in the meantime.
- `same_function_from_label_sets` — set-valued same-function via a sparse membership matrix,
  vectorised over tens of millions of pairs.
- `bootstrap_auroc_by_group` — CI by resampling proteins/clusters, not pairs.
- `precision_at_k` uses `argpartition`; `evaluate_pairs` additionally reports `enrichment_at_k`
  (precision@k ÷ prevalence), the framing that makes the 91× top-of-ranking effect visible.

**`labels.py`**
- `load_ec_label_sets` returns `{accession: frozenset(all EC numbers)}`. `load_ec_labels` is
  unchanged, so notebooks 04-05 keep their current behaviour.

**`scripts/rerun_pilot_eval.py`** — re-embeds the saved peptide table with the fixed pooling and
scores six representations (old pooling raw/centred, fixed pooling raw/centred, composition,
length) with bootstrap CIs and within/cross-organism strata, writing
`results/pilot06_rerun_metrics.csv`.

**`tests/`** — 15 new tests in `test_pair_eval_v2.py` (blocked cosine equals the full matrix,
centring, composition fractions including the non-standard-residue case, sampled-pair prevalence
matching enumeration, group stratification, set-label intersection versus primary-only,
enrichment reporting, bootstrap behaviour including the degenerate case) plus label-set coverage.
Suite: 51 passed.

## Post-fix results (`results/pilot06_rerun_metrics.csv`)

3.0M sampled cross-protein pairs, prevalence 0.326%, 5,403 proteins; AUROC CIs are
protein-level bootstraps (40 resamples).

| representation | AUROC [95% CI] | AP | corr(cos, \|Δlen\|) | enrich@10 | enrich@100 | AUROC within-org | AUROC cross-org |
| --- | --- | --- | --- | --- | --- | --- | --- |
| special tokens in pool, raw | 0.470 [0.457, 0.483] | 0.00302 | −0.330 | 0× | 0× | 0.468 | 0.468 |
| special tokens in pool, centred | 0.502 [0.496, 0.508] | 0.00375 | −0.491 | 123× | 21.5× | 0.514 | 0.493 |
| **masked pool, raw** | 0.467 [0.454, 0.479] | 0.00299 | −0.272 | 0× | 3.1× | 0.466 | 0.466 |
| **masked pool, centred** | 0.502 [0.495, 0.507] | 0.00382 | −0.426 | 123× | 30.7× | 0.515 | 0.491 |
| amino-acid composition | 0.495 [0.484, 0.504] | 0.00323 | −0.025 | 0× | 0× | 0.509 | 0.485 |
| peptide-length match | 0.504 [0.495, 0.513] | 0.00330 | −1.000 | 0× | 0× | 0.507 | 0.502 |

Reading it:

1. **No global function signal, and none beyond composition.** Best embedding AUROC 0.502
   [0.495, 0.507] against composition's 0.495 [0.484, 0.504] — overlapping, both at chance.
   Δ(embedding − composition) ≈ 0.007 with CIs that overlap, so the pilot's headline question
   ("do peptide embeddings carry function signal beyond counting letters?") answers **no** at this
   model size and pooling.
2. **The raw-cosine CIs exclude 0.5 from below** (0.454-0.483) in both pooling variants, so the
   anti-predictive behaviour is real and reproducible, and centring is the correction.
3. **The pooling fix improves retrieval, not ranking.** Global AUROC unchanged; enrichment at
   k=100 rises 21.5× → 30.7× and at k=1,000 5.5× → 6.7×. If the eventual application is
   "retrieve candidate same-function peptides", that is the regime that matters.
4. **Within-organism > cross-organism**, consistently (0.515 vs 0.491 for the fixed centred
   embeddings; composition shows the same pattern, 0.509 vs 0.485). Shared amino-acid usage within
   a genome behaves like weak function signal, which is exactly why the organism stratification
   was worth adding — and why the sequence-identity stratification is still the missing control.
5. **Length is not the confound it looked like.** Length-match alone gives 0.504 [0.495, 0.513]:
   chance. So length coupling degrades the embedding score without itself predicting function.

Net: the negative result now has its baseline and its interval, and it is an honest negative —
with the caveat that homologous and non-homologous positives are still pooled, so the claim it
supports is "no signal beyond composition", not yet "no *non-homologous* function signal".

## How much of a short-peptide embedding *is* composition (added 2026-09-24)

A direct measurement of the information ceiling, needing no model runs — ridge regression from
amino-acid composition + log length (21 features) onto the 480-dim embedding, scored on a
protein-grouped held-out 30%:

| input | median length | embedding variance explained by composition + length |
| --- | --- | --- |
| tryptic peptides (n=10,804) | 18 aa | **R² = 0.536** |
| full-length proteins (n=1,740) | 338 aa | **R² = 0.111** |

A *linear* map recovers over half the variance of a peptide embedding, and only a ninth of a
full-protein embedding — and a linear map is a lower bound on what composition determines, so the
true share is higher still. This is the mechanism behind the pilot's null: at tryptic-peptide
length the embedding is largely a re-encoding of composition plus length, so a composition
baseline is hard to beat by construction. At protein length the same features explain almost
nothing, which is why PLM embeddings are useful there.

(The full-protein embeddings come from the pre-fix pooling, but at 338 aa the two special tokens
are 0.6% of the pool, so the comparison is unaffected.)

Practical consequence for the benchmark: the peptide-level ceiling is low for reasons that more
model capacity will not change. If the project continues, the interesting variants are the ones
that add *context* (peptide embedded in the parent-protein frame, or a trained peptide→protein
projection) rather than a bigger checkpoint.

## Notebook 05 rerun: accuracy vs coverage, with the two missing reference points

`results/accuracy_vs_coverage_v2.csv`, from `scripts/rerun_probe_eval.py`. The probe is logistic
regression trained on full-length **train**-split embeddings; peptides scored are **test**-split
only. 7 EC top-level classes.

| stratum | median length | n | accuracy (masked pool) | accuracy (old pooling) | macro F1 (masked) |
| --- | --- | --- | --- | --- | --- |
| full-length proteins (test) | 342 aa | 348 | **0.698** | 0.698 | 0.646 |
| peptides, coverage ≤0.1 | 17 aa | 25,622 | **0.205** | 0.183 | 0.096 |
| peptides, coverage 0.1-0.25 | 31 aa | 2,449 | **0.302** | 0.275 | 0.185 |
| peptides, coverage 0.25-0.5 | 34 aa | 81 | **0.346** | 0.309 | 0.250 |
| *majority-class rate, test proteins* | — | 348 | *0.319* | — | — |
| *majority-class rate, test peptides* | — | 28,152 | *0.292* | — | — |

Three things this now shows that the original three-row table could not:

1. **The trivial baseline beats the probe on 91% of the peptides.** Always predicting the most
   common EC class scores 0.292 on the test peptides. The largest bucket (coverage ≤0.1, 25,622
   peptides, median 17 aa) scores **0.205 — 0.087 *below* that baseline**. The 0.1-0.25 bucket is
   +0.010 and the 81-peptide bucket +0.054. So the rising 0.18 → 0.27 → 0.31 curve in the original
   figure is not a weak-but-real dose-response; on the bulk of the data the probe is worse than
   guessing the majority class. Macro F1 of 0.096 says the same thing differently: it is
   effectively predicting one or two classes.
2. **The cost of fragmentation, quantified.** Same probe, same proteins: 0.698 on full-length
   sequences (2.2× the majority rate) against ~0.21 on their own tryptic peptides. That gap is
   the information lost by digestion, and it is the ceiling any peptide-level functional method
   works under.
3. **The pooling fix helps exactly where predicted and nowhere else.** Peptides gain +0.023,
   +0.027, +0.037 accuracy; full-length proteins are unchanged to three decimals (0.698, macro F1
   0.6466 → 0.6459). The special-token share of the pooled mean scales as 1/length, so a 342-aa
   protein cannot care and a 17-aa peptide does.

Caveat retained from §6: the train/test split is by protein but not by sequence-identity cluster,
so homologues sit on both sides and 0.698 is optimistic for genuinely novel proteins.

## Suggested order for the rerun

1. Fix special-token pooling in `embed.py`; re-embed (both notebooks).
2. Add composition and length baselines to `pair_eval`; report Δ against them.
3. MMseqs2 cluster the pooled parent proteins; attach `cluster_id`.
4. Recompute 06 stratified by cluster (different-cluster = headline, same-cluster = positive
   control) × organism (within/cross), with cluster-level bootstrap CIs.
5. Re-split 05 by cluster group, add the full-length test control, regenerate the
   accuracy-vs-coverage curve.
6. Only then vary the checkpoint size or pooling scheme.
