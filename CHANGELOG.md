# Changelog

## v0.4.0 (unreleased)

Changes since v0.3.4. Several of them change results, so z-scores, cutoffs and
degenerate sequences from v0.3.4 runs will not match exactly.

### Changes that affect results

- **Depth-aware statistics.** The CLR model now draws Dirichlet samples to estimate
  error from sequencing depth. P-values and Benjamini–Hochberg adjusted p-values are
  calculated more accurately. The sampler uses a fixed seed by default (0), so
  identical inputs give identical results.
- **Low-depth data.** At very low depth with strong depletion, z-scores were NaN.
  They are now finite.
- **Read counting.** The total read count was off by one. 5′ PAMs truncated by the
  end of a read are now skipped instead of counted.
- **Degenerate sequences.**
  - Candidate generation now enumerates every exact degenerate pattern.
  - The three-stage CP-SAT solve uses a deterministic time budget, so the same
    selection gives the same sequences on any machine.
  - A solve that reaches its budget returns the best exact cover found. It no
    longer fails or falls back to one sequence per k-mer.
  - Selecting all k-mers of a length returns `NNN…` instead of an error.
- **Top-50% selection limit.** By default, at most the top 50% of k-mers by z-score
  can be selected; a lower cutoff is raised to that limit. Even near-PAMless nucleases
  such as SpRY prefer NRN, and smaller selections give degenerate sequences that are
  found faster and proven minimal. The automatic cutoff never selects more than the
  top 40%, so it is unaffected.

### New

- **Optimality reporting.** When the degenerate sequences are an exact cover but not
  proven to be the smallest set:
  - the web app shows a note under the sequences;
  - `tamipami predict` logs a warning and sets `proven_optimal` to `False` in
    `degenerate_pam_tam.csv` (a new column).
- **Web app.** A "Limit selection to the top 50% of k-mers" checkbox, on by default.
- **CLI.**
  - New `predict --max_fraction` option (default 0.5; use `1` to disable the limit).
  - `predict` reads from STDIN when `--input` is omitted.
  - `--cutoff` is optional; a cutoff is computed automatically for each length when
    it is omitted.
- **Python API.** New `degenerate.solve_degenerates()`, which returns the sequences and
  a `proven_optimal` flag. `seqs_to_degenerates()` is unchanged. New
  `pamSeqExp.limit_cutoff()` and `pamSeqExp.top_fraction_cutoff()`.
- **Benchmarks.** `benchmark/degenerate_solver/` compares the degenerate solver with an
  independent global-optimum oracle on 2,365 test sets, and times selections of 5–100%
  of k-mers.

### Changed

- **CLI options aligned with the web app.** `--length` takes 3–5 (default 5). Logos use
  k-mers with z-score ≥ cutoff.
- **HDF5 export overwrites the output file instead of appending to it.**
- **Logs.**
  - The CLI log file is overwritten on each run, at INFO level by default.
  - The web app logs to `TAMIPAMI_LOGFILE` (default: `tamipami_streamlit.log` in the
    system temp directory).
- **Logo files.** Logos can also be saved as SVG, EPS, PS, PNG, JPEG or TIFF.
- **Dependencies.** Minimum versions are declared in `pyproject.toml`, including
  `ortools>=9.15,<10` and `streamlit>=1.50,<2`. `numpy` is now listed explicitly.

### Fixed

- **Security.**
  - Paths containing shell metacharacters could be interpreted by `bbmerge.sh`.
    Input files are now passed through safe symlinks.
  - The spacer search used a regular expression open to catastrophic backtracking.
  - Spacer sequences are validated (A, C, G and T only).
- **Error handling.**
  - A failed BBMerge run now raises an error. Before, it silently produced empty
    counts.
  - A missing input file raises `FileNotFoundError`.
- **Web app.**
  - Uploads are written to a private temporary directory that is deleted after
    processing, instead of the package directory.
  - Plain FASTQ uploads are no longer re-compressed.
  - The HDF5 export is built once per run.
  - Figures are closed after display, which fixes a memory leak.
  - A failure in the degenerate-sequence calculation shows a warning instead of
    stopping the app.
- **`tamipami serve`** works from any working directory.

### Docker

- **Hardening.** The image runs as a non-root user. `compose.yaml` runs the container
  with a read-only root filesystem, a `tmpfs` at `/tmp`, all capabilities dropped,
  `no-new-privileges` and a 4 GB memory limit.
- **Image.** The image is smaller (3.3 GB, down from 3.9 GB) and has a health check.
  Upload size, error detail and solver threads (4) are set via environment variables.
- **CI.** The test suite must pass before an image is published to GHCR.
