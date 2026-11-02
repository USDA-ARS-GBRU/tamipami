# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

TamiPami identifies PAM/TAM recognition sequences of Cas/TnpB endonucleases from paired-end FASTQ of a control (uncut) and an experimental (cut) plasmid library. It ships as a CLI (`tamipami`) and a Streamlit web app. The version comes from git tags via `setuptools_scm` (`tamipami/_version.py` is generated; don't edit it). `CHANGELOG.md` tracks unreleased v0.4.0 changes that alter results.

## Setup and commands

Dependencies are split deliberately: conda (`conda-requirements.txt`, includes **BBMap** which provides the external `bbmerge.sh` binary) then pip (`pip-requirements.txt`: streamlit, ortools, ckmeans, logomaker), then `pip install .`. Install order matters on Apple Silicon (Abseil conflicts with ortools). Python 3.13.

```bash
python -m pytest -q                                   # full suite (CI runs this before building Docker)
python -m pytest tests/test_release_regressions.py -q # one file
python -m pytest tests/fastq/test_pam_stream.py::test_name -q   # one test
tamipami process --cont1 c_R1.fq.gz --cont2 c_R2.fq.gz --exp1 e_R1.fq.gz --exp2 e_R2.fq.gz --library RTW574 --outfile out.h5
tamipami predict --input out.h5 --predict_out outdir [--cutoff '{"3":0.7}'] [--max_fraction 0.5]
tamipami serve                                        # Streamlit app on :8501
```

`process` and `predict` can be piped (HDF5 bytes over STDOUT/STDIN). There is no linter configured. Benchmarks live in `benchmark/` (`run_benchmark.py` for CLI timing, `benchmark_app_chrome.py` via Selenium, `degenerate_solver/` for solver oracle comparisons) and need data from `data/fastq` (untracked, large).

## Architecture

Pipeline: FASTQ → k-mer counts → z-scores/cutoff → selected k-mers → minimal degenerate IUPAC set.

- **`fastq.py`** – Shells out to `bbmerge.sh` (args from `assets/config.yaml`, loaded by `config.py`) to merge read pairs, streams the output, and counts PAM/TAM k-mers next to the spacer (orientation `5prime` = PAM before target, `3prime` = target then PAM like SpCas9). Addgene library presets (RTW554/555/572/574) supply spacer/orientation. Input paths are symlinked into a temp dir because `bbmerge.sh` evals its args (see `_safe_link`).
- **`pam.py`** – `pamSeqExp` holds control/experimental counts per k-mer length (`multikmerdict`, lengths 3–8), builds the multi-length tables, computes CLR-based z-scores with Dirichlet sampling (fixed seed 0 for determinism) and BH-adjusted p-values. `find_breakpoint` picks the automatic cutoff (ckmeans / KDE); `limit_cutoff` / `MAX_SELECTED_FRACTION` (0.5) cap how many top-z k-mers can be selected. `entropy.py` supports the entropy-based breakpoint.
- **`degenerate.py` + `candidates.py`** – Convert the selected k-mer set into the fewest IUPAC "rectangles" that exactly cover it. `candidates.py` enumerates exact rectangle candidates (dendrogram-induced + positional subsets, with size caps raising `CandidateLimitExceeded`); `degenerate.py` solves a three-stage lexicographic OR-tools CP-SAT problem (min rectangles → min degenerate positions → min allowed bases) with a deterministic time budget. Results carry `proven_optimal`; a timeout returns the best exact cover, never a fallback. Design notes are in `docs/degenerate_methods.md`.
- **`tpio.py`** – HDF5 import/export (keys `k_<length>`), STDIN-to-tempfile handling, Altair histogram.
- **`cli.py`** – argparse entry (`process`, `predict`, `serve`); `cli_utils.py` holds argument validators (cutoff JSON for lengths 3–5, max_fraction). `serve` launches Streamlit on `app.py`.
- **`app.py`** – Streamlit UI; uploads are written to a temp dir and processed through the same `fastq`/`pam` code as the CLI. Does not cache FASTQ processing. `GA_MEASUREMENT_ID` env var enables optional Google Analytics.

CLI and web app must stay behaviorally compatible (`tests/test_cli_app_compatibility.py`); results-affecting changes should be recorded in `CHANGELOG.md`.
