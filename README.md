# PLM Robustness to Incomplete Proteins — Benchmark

Question: how well do protein language model (PLM) embeddings retain
predictive signal when computed from an incomplete protein (a peptide
fragment) instead of the full sequence? This is the feasibility check
behind **applying PLMs to metaproteomics-derived peptides** to infer
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
main.nf         Nextflow workflow (WIP): fetch -> digest
nextflow.config Nextflow params + profiles (docker)
pyproject.toml  package definition; registers the plm-* CLI commands
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
  cli/          thin command-line wrappers, one per pipeline stage
                (plm_fetch.py, plm_digest.py; embed/probe/eval still to come)
scripts/run_benchmark.py   CLI entrypoint (stub; superseded by src/plm_benchmark/cli/)
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
python -m pip install -e .             # installs the package + plm-fetch / plm-digest commands
pytest tests/ -v
```
`make venv` builds the same environment and registers the "Python (plm-benchmark)"
Jupyter kernel.

## Setup (Docker)
The image contains only the installed package (`pip install .` from
`pyproject.toml`); configs and data are passed in at run time.
```bash
docker build -t plm-benchmark -f docker/Dockerfile .
docker run --rm plm-benchmark plm-digest --help
```
The old `docker compose ... pytest` flow no longer works (tests are not in the
image); run tests from the local venv.

## Running the stages by hand (CLI)
Each stage reads files and writes a file, so it can be run on its own:
```bash
plm-fetch  --organism-config configs/organism_ecoli_k12.yaml --out-path data/raw/ecoli.fasta
plm-digest --fasta-path data/raw/ecoli.fasta --enzyme-config configs/enzymes.yaml \
           --out-path data/processed/peptides.parquet
```
The science code stays in `src/plm_benchmark/`; the CLIs only parse arguments,
read the input, call the library function and write the output.

## Nextflow pipeline (work in progress)
`main.nf` wires the CLIs together, one process per stage. Currently: `FETCH` ->
`DIGEST`. Embed, probe and evaluation stages are not yet added.

Requires Java 17+ and Nextflow.
```bash
nextflow run main.nf                   # uses plm-* commands on your PATH (activate .venv first)
nextflow run main.nf -profile docker   # runs each process in the plm-benchmark image
nextflow run main.nf -resume           # reuse cached steps
```
- Parameters (`params.organism_config`, `params.enzyme_config`, `params.outdir`)
  are set in `nextflow.config` and can be overridden, e.g. `--outdir my_run`.
- Published outputs go to `results/` (`raw/`, `processed/`); `work/` is
  Nextflow's cache and is not committed.
- `-resume` does not detect changes to the Python package. After editing code,
  rerun without `-resume` (or use a new image tag), and rebuild the image if
  using `-profile docker`.
- Adding a step: write the function and a CLI in `src/plm_benchmark/cli/`,
  register it under `[project.scripts]` in `pyproject.toml`, reinstall
  (`python -m pip install -e .`), add a process in `main.nf`, and rebuild the image.

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
- [ ] Dockerize end to end (image covers fetch + digest so far)
- [ ] Set up `run_benchmark` as a Nextflow pipeline (digest -> embed -> probe ->
      evaluate -> report), one process per stage, containerised, with a
      `test` profile. **In progress:** `fetch` and `digest` stages done
      (`main.nf`, CLIs in `src/plm_benchmark/cli/`); embed, probe, labels and
      pair-eval stages still to do
- [ ] Extend to more organisms (Phase 2/3), then PTMs (Phase 3+)
- [ ] Follow-ups suggested by the pilot: CLS-token pooling / embedding
      whitening (anisotropy fix), a larger ESM2 checkpoint, within- vs.
      cross-organism AUROC breakdown

See `docs/agent-workflow.md` for how to use Claude subagents (research /
build / verify) alongside pytest while developing this.
