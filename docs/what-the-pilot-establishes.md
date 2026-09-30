# What the pilot establishes, and what it leaves open

Written 2026-09-24, after the notebook 05/06 review and rerun. Numbers from
`results/pilot06_rerun_metrics.csv` and `docs/pilot-06-review.md`.

## The motivating questions

1. *Can this fragment help infer the function of its parent protein?*
2. *Can this fragment help infer the functions present in the microbiome under study?*

These differ in their estimand. (1) asks for a per-item label; (2) asks for a distribution over a
population. A method can fail at (1) and still serve (2), because aggregation tolerates per-item
error that an individual assignment cannot.

## Claim 1 — a scoped negative

For in-silico tryptic peptides (median 18 aa) from four well-annotated bacterial proteomes,
**frozen `esm2_t12_35M` embeddings, mean-pooled and compared by raw cosine, do not rank same-EC
peptide pairs above different-EC pairs any better than amino-acid composition does.**

AUROC 0.502 [0.495, 0.507] versus composition 0.495 [0.484, 0.504], protein-level bootstrap;
average precision 0.0038 against a prevalence of 0.0033.

The three qualifiers are load-bearing and should always be said out loud: *frozen* (nothing was
trained), *mean-pooled* (no other pooling or layer tested), *raw cosine* (the weakest possible
probe of an embedding space). What is refuted is the specific, cheap version of the idea — the one
a reviewer would ask about first, and the one most likely to be tried by someone who assumes it
works.

One thing the result is **not** fragile to: the pair set includes the easy homologous positives,
and still shows nothing. Those are ~100 of ~26,000 positives, so excluding them cannot rescue the
metric. The sequence-identity control (still outstanding) will sharpen the statement, not reverse
it.

## Claim 2 — the mechanism, which is the transferable part

Composition + log length (21 features, linear ridge, protein-grouped held-out set) explains
**R² = 0.536** of the variance of a tryptic-peptide embedding, but only **R² = 0.111** of a
full-length protein embedding (median 338 aa).

A masked language model earns its keep by conditioning each residue on its context. A tryptic
peptide is ~1-5% of its protein and carries almost none, so at that length the representation
collapses toward a re-encoding of composition and length — and a composition baseline becomes hard
to beat by construction. This predicts that **more model capacity will not fix it**, which is a
stronger and more useful statement than the AUROC alone. It also predicts where the headroom is:
anything that restores context.

## Claim 3 — the methodological findings

Four traps that make a naive version of this experiment look like it works, or like it fails for
the wrong reason. Each was measured here, not assumed:

- **Anisotropy produces a below-chance AUROC.** Raw cosine gave 0.467-0.470 with CIs excluding 0.5
  from below. Mean-centring moves it to 0.502. Reporting "0.47" without centring invites a story
  about the embedding being actively misleading, which is false.
- **Homologous positives create spectacular top-of-ranking enrichment that is not function.**
  123× at k=10, 31× at k=100, of which 79-100% are cross-species pairs — the same enzyme in two
  genomes, i.e. sequence identity, which alignment detects better and cheaper. A study that
  reported only precision@k would claim success here.
- **Special tokens in the pooled mean encode length.** `<cls>`/`<eos>` are 25% of a 6-mer's pooled
  vector and 5% of a 40-mer's; masking them changes a 7mer's embedding by 15.3% and a 27mer's by
  4.9%, and reduces `corr(cosine, |Δlength|)` from −0.330 to −0.272. It did *not* change the AUROC —
  worth knowing, because the fix is correct but is not the explanation.
- **Pair-level confidence intervals are roughly an order of magnitude too narrow.** 58.4M pairs come
  from 5,403 proteins; the interval has to be bootstrapped over proteins.

For an interview, Claim 3 is the part that demonstrates method judgement: the pipeline was tested
against its own artefacts before its result was believed.

## Mapping back to the two questions

**Question 1 (parent-protein function).** The pilot tested the hardest possible version: one
18-mer, alone, no parent context, no training. It says that a single tryptic peptide does not
carry standalone functional signal at this representation. It does **not** say the question is
unanswerable, because the realistic unit is different in two ways:

- In a real metaproteomics run you usually observe **several peptides from the same protein**. The
  natural unit is the peptide *set*, and joint evidence from k peptides is a different statistical
  object from one peptide. (Directly testable with data already on disk: notebook 05's digest has
  every peptide per protein, not the 2-peptide cap used in 06.)
- The incumbent route does not ask the peptide to carry function at all: **align the peptide to a
  reference protein, then inherit that protein's annotation** (NovoLign-style DIAMOND alignment +
  LCA). That works, and the pilot is not evidence against it. Where it fails is the population the
  de novo peptides were sequenced to reach — novel organisms with no reference match — and that is
  the real open problem.

**Question 2 (microbiome-level function).** Untouched by this pilot, and better posed. The
estimand is a *distribution over functions in the sample*, not a label per peptide, which changes
what "good enough" means: unbiased or correctable per-peptide error can still support a population
estimate, coverage matters more than per-item accuracy, and the deliverable is a functional profile
with its denominator (what fraction of signal is annotatable at all). This is the same logic
phylopeptidomics uses for taxonomy — do not force an assignment on each ambiguous peptide, model
the aggregate — transplanted to function.

## What to do next, in order of expected value

1. **Peptide sets, not single peptides.** Does same-EC AUROC rise with the number of peptides
   aggregated per protein? Cheap, uses existing embeddings, and maps exactly onto Question 1. If it
   rises steeply, the negative result is about single fragments rather than about peptides.
2. **Restore context.** Embed the peptide inside its candidate parent frame, or train the
   asymmetric map (dual encoder on peptide/parent-protein pairs with in-batch negatives) so the
   peptide representation is learned to point at a protein rather than assumed comparable to one.
3. **Reformulate at the community level** (Question 2): given a set of de novo peptides, estimate
   the sample's functional composition with coverage and confidence, validated against a paired
   metagenome-derived database search. The comparator is "what did the database search miss", which
   is also NovoLign's most interesting use.
4. **Use the PLM where it is additive rather than redundant.** Sequence plausibility
   (pseudo-perplexity) needs no database and no function labels, and improves the quality of the
   peptide set feeding 1-3. Alignment remains better for containment and identity; a PLM is not
   competing there.
5. **Finish the controls** on anything that does show signal: sequence-identity clusters
   (`mmseqs easy-cluster --min-seq-id 0.5 -c 0.8`), the composition baseline, the non-homologous
   stratum as the headline.

## "Why bother with a PLM when BLAST exists?"

Claim 1 found that embedding proximity tracks sequence homology. That is the regime where
alignment is *already* strong, so using embeddings as a homology detector is largely redundant:
you pay a forward pass per sequence to approximate, without a calibrated significance measure,
something DIAMOND computes exactly and cheaply. Two commonly-offered advantages do not survive
scrutiny:

- **"No database search needed."** Not for retrieval. To find similar sequences by embedding you
  must still embed the reference set and search it — you have swapped a k-mer index for a vector
  index, made it approximate (ANN), and made it more expensive to build. The database is still
  there. What genuinely needs no database is **single-sequence scoring**, which is not homology
  detection.
- **"It's pre-trained, so I just feed it a sequence."** True, but BLAST needs no training *and* no
  weights, no GPU, no framework; it is deterministic and ships an E-value. On ease of use,
  alignment wins. Convenience is not the differentiator.

What a PLM does add, in decreasing order of how defensible it is:

1. **Scoring a sequence with no reference at all.** Pseudo-perplexity under masked-LM scoring
   answers "does this read like real protein", and variant scoring answers "is this substitution
   tolerated". Alignment structurally cannot do this: no hit means no output. For de novo peptides
   this is the one clearly additive use, and it needs no function labels and no database.
2. **Remote homology, below the alignment twilight zone.** Under roughly 25-30% identity, pairwise
   E-values lose discriminative power. Note the fair comparator here is not BLAST but **profile**
   methods — PSI-BLAST, jackhmmer, HHblits — which already recover much of that ground by building
   an MSA. The PLM's advantage is reaching comparable sensitivity *from a single sequence*, with no
   MSA to construct, which matters when the family is too sparse to profile.
   (Specific benchmark numbers should be checked before being quoted; the verified datapoint in
   this vault's reading is Empathi, doi:10.1038/s41467-025-64177-5, which reports better
   sensitivity than homology-based tools for phage protein function — at protein level.)
3. **Information alignment does not represent.** Embeddings encode residue covariation, hence
   contacts and fold; identity encodes none of that. So embeddings can group sequences whose
   function is conserved past the point where their letters are recognisably related.
4. **Fixed-length features for supervised models.** A 480-dim vector drops into a classifier. Worth
   noting the honest baseline: k-NN over alignment bit-scores is also a usable feature and is
   frequently competitive.

What alignment keeps: **containment and coordinates** (where in this protein does this peptide
sit — a similarity scalar cannot answer that), **calibrated significance** (an E-value has a null
model; cosine similarity has none, which is exactly why centring changed our AUROC by 0.035),
**short queries**, and cost.

**The consequence for this project.** The remote-homology argument is a protein-length argument.
You cannot establish remote homology from 18 residues with *any* method — there is not enough
information in the query — so the main reason to prefer a PLM over alignment does not transfer to
tryptic peptides. For peptide→parent mapping, alignment is the right tool; the PLM's additive role
in a de novo pipeline is sequence plausibility, upstream of the mapping.

**The experiment that would settle it** turns the outstanding control into a positive question:
stratify pairs by parent-protein sequence identity, then ask whether embedding proximity recovers
same-function pairs *that alignment misses* — i.e. positives in the <30% identity stratum. If it
does, the method has a defensible niche; if it does not, embeddings are a more expensive route to
what DIAMOND already returns. This is the same `mmseqs easy-cluster` output needed for the control,
used to define the regime rather than to remove it.

## One-sentence version

The pilot establishes that a single tryptic peptide's frozen PLM embedding carries no functional
signal beyond its amino-acid composition, and explains why — at that length the embedding largely
*is* composition — which redirects the question from "label each peptide" to "aggregate evidence
across peptides and across the community".
