# PLM Robustness to Incomplete Proteins — Benchmark

Question: how well do protein language model (PLM) embeddings retain
predictive signal when computed from an incomplete protein (a peptide
fragment) instead of the full sequence? This is the feasibility check
behind applying PLMs to metaproteomics-derived peptides to infer
properties of their parent proteins.

Full design notes live in the "Protein Modeling" Claude project
(`plm-benchmark-brief.md`).

## Phase 1 scope
- Organism: *E. coli* K-12 MG1655 (UniProt reference proteome UP000000625)
- Fragmentation: in-silico trypsin digestion (+ missed cleavages) and
  controlled sliding-window / terminal truncation
- Embeddings: frozen ESM2 embeddings (mean-pooled), no fine-tuning
- Evaluation: simple probe (logistic regression) trained on full-length
  embeddings, evaluated on fragment embeddings; performance vs. fragment
  length/coverage is the headline result

## Project layout
```
configs/        organism + enzyme/fragmentation configs
data/raw        downloaded proteome FASTA + label sources (not committed)
data/processed  digested/truncated peptide sets, embeddings (not committed)
src/plm_benchmark/
  digest.py     in-silico enzymatic digestion (trypsin, Lys-C, Glu-C, ...)
  fragments.py  sliding-window / terminal truncation utilities
  data_io.py    fetch + parse UniProt reference proteomes
  labels.py     fetch + parse EC (Enzyme Commission) number labels
                (top-level class or full number)
  embed.py      frozen PLM embedding extraction (ESM2 via HF transformers)
  probe.py      simple sklearn probe: train/evaluate on embeddings
  evaluate.py   orchestrates the end-to-end benchmark run
  pair_eval.py  pairwise same-function-vs-different-function similarity
                eval (AUROC / average precision / precision@k)
scripts/run_benchmark.py   CLI entrypoint
tests/          pytest tests for the deterministic logic (digestion, fragments, labels, pair_eval)
docker/         containerized benchmark environment
docs/           agent/testing workflow notes
results/        benchmark outputs (metrics, plots) (not committed)
```

## Setup (local development)
```bash
cd ~/Projects/ProteinPrediction
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements-core.txt   # digestion, tests, data IO
# pip install -r requirements-ml.txt   # only once ready to embed (torch + transformers, large)
PYTHONPATH=src pytest tests/ -v
```

## Setup (Docker)
```bash
docker compose -f docker/docker-compose.yml build
docker compose -f docker/docker-compose.yml run --rm benchmark pytest tests/ -v
```

## Status
- [x] Repo scaffold
- [x] Trypsin + alternative-enzyme in-silico digestion, tested
- [x] Sliding-window / truncation fragment generation, tested
- [x] Fetch E. coli K-12 reference proteome
- [x] Pick a label: EC top-level class (`plm_benchmark.labels`, notebook 04) —
      covers ~1,740/4,403 proteins (enzymes only)
- [x] Frozen ESM2 embedding extraction (notebook 03, subsampled so far)
- [x] Probe training + fragment-vs-full-length evaluation (notebook 05):
      full-length accuracy 0.70 vs. 0.18-0.31 on tryptic peptides —
      accuracy degrades sharply with coverage
- [x] Results plot: accuracy vs. fragment length/coverage
      (`results/accuracy_vs_coverage.png`)
- [x] Pilot: cross-organism peptide function-similarity via pairwise AUROC
      (notebook 06, `plm_benchmark.pair_eval`) — 4 bacterial proteomes,
      same-vs-different EC number peptide pairs. AUROC ~0.47 (chance),
      driven by embedding anisotropy (mean-pooled cosine similarities all
      cluster in 0.85-1.0 regardless of function) — corroborates the
      probe result from a model-free angle
      (`results/pilot_function_similarity.png`)
- [ ] Dockerize end to end
- [ ] Set up `run_benchmark` as a Nextflow pipeline (digest -> embed -> probe ->
      evaluate -> report), one process per stage, containerised, with a
      `test` profile
- [ ] Extend to more organisms (Phase 2/3), then PTMs (Phase 3+)
- [ ] Follow-ups suggested by the pilot: CLS-token pooling / embedding
      whitening (anisotropy fix), a larger ESM2 checkpoint, within- vs.
      cross-organism AUROC breakdown

See `docs/agent-workflow.md` for how to use Claude subagents (research /
build / verify) alongside pytest while developing this.
